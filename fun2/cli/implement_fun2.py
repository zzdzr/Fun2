#!/usr/bin/env python3
"""
Master Pipeline for integrating O/E matrix generation, cooler transformation,
SoN signal generation, summits combination, and planning.

This file defines:
    - PipelineConfig: A dataclass that encapsulates all configuration parameters.
    - GenerateOEMatrix: Generates O/E matrices per chromosome.
    - TransformToCooler: Transforms the generated O/E matrices into a .cooler file.
    - GenerateSoNSignal: Generates SoN signals based on the .cooler file and given parameters.
    - CombineSummits: Combines multiple SoN output files (summits) into a merged result.
    - Planning: Uses MCTS for fountains identification.
    - MasterPipeline: Orchestrates the complete workflow sequentially.
    
Each class and function is documented using Google style docstrings.
"""

import os
import logging
import argparse
import yaml
import warnings
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Any

from fun2.planning.planning import Fun2Planning
from fun2.logging import log_step, skip_if_exists
from fun2.perform_command import run_single_command
from fun2.preprocessing.find_summits import SummitsFinder

# -----------------------------------------------------------------------------
# Global Logging Setup
# -----------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logging.captureWarnings(True)
logging.getLogger("py.warnings").setLevel(logging.ERROR)

warnings.filterwarnings("ignore", module="numpy")
warnings.filterwarnings("ignore", module="pandas")
# -----------------------------------------------------------------------------
# Pipeline Configuration (using dataclass)
# -----------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    """Configuration for the pipeline execution.

    Attributes:
        working_dir (str): Working directory for pipeline execution.
        juicer_jar (str): Path to the Juicer tools jar file.
        chrom_size (str): Path to the chromosome size file.
        full_name (str): Identifier for the current run.
        input_hic (str): Path to the input hic file.

        generate_cooler_script (str): Script for converting OE matrix to cooler.
        oe_matrix_dir (str): Directory to store OE matrix files.
        oe_output_template (str): Template for per-chromosome OE matrix output.

        cooler_output_dir (str): Directory for storing cooler file.
        assembly (str): Genome assembly identifier (e.g., hg19).

        summits_threshold (float): Threshold for summits combination.
        padding (int): Padding to use when combining summits.
        summits_output_dir (str): Directory for storing summits output.

        calSoN_script (str): Script for generating SoN signals.
        SoN_params (List[Dict[str, Any]]): List of dictionaries with SoN parameters.
        resolution (int): Resolution for OE and SoN calculations.
        oe_normalization (str): Normalization method for OE matrix.

        ini_height (int): Initial height parameter for planning.
        ini_width (int): Initial width parameter for planning.
        ini_angle (int): Initial angle parameter for planning.
        layer_height (int): Layer height for planning.
        max_iter (int): Maximum iterations for planning.
        exploration_constant (float): Exploration constant.
        alpha (float): Alpha parameter for planning.
        es_popsize (int): Population size for evolutionary strategy.
        es_sigma (float): Sigma for evolutionary strategy.
        momentum (float): Momentum parameter.
        max_rollout_depth (int): Maximum rollout depth.
        eta (float): Eta parameter.
        initial_lr (float): Initial learning rate.
        min_lr (float): Minimum learning rate.
        decay_factor (float): Decay factor for learning rate.
        gamma (float): Gamma parameter.
        action_limits (List[Tuple[float, float]]): Limits for actions.
        fountain_path (str): Path to save fountain results.
        n_process (int): Number of parallel processes for planning.
    """

    # Basic settings
    working_dir: str = "."
    juicer_jar: str = ""
    chrom_size: str = "" # Path to a chromosome size file (e.g., a text file with chrom names)
    full_name: str = ""
    input_hic: str = ""

    # O/E matrix generation settings
    generate_cooler_script: str = ""
    oe_matrix_dir: str = ""
    oe_output_template: str = "{oe_matrix_dir}/{chrom}.oe_mat.txt"  # Template for output per chromosome

    # Cooler conversion settings
    cooler_output_dir: str = ""
    assembly: str = "hg19"

    # Summits combination settings
    summits_threshold: float = 0.0
    padding: int = 0
    summits_output_dir: str = ""

    # SoN calculation settings
    calSoN_script: str = ""
    SoN_params: List[Dict[str, Any]] = field(default_factory=lambda: [
        {"height": 100, "width": 25, "angle": 45, "output": None},
        {"height": 200, "width": 25, "angle": 45, "output": None},
        {"height": 300, "width": 25, "angle": 45, "output": None}
    ])
    resolution: int = 5000
    oe_normalization: str = "VC_SQRT"
    chrom: str = None

    # Derived fields (populated in __post_init__)
    cooler_output: str = ""
    combined_summits_output: str = ""

    # For planning
    boundary_margin: int = 300
    expansion_margin: int = 300
    ini_height: int = 180
    ini_width: int = 25
    ini_angle: int = 45
    layer_height: int = 1
    edge_width: int = 4
    max_iter: int = 250
    exploration_constant: float = 5
    alpha: float = 0.1
    es_popsize: int = 8
    es_sigma: float = 5.0
    es_lr: float = 2.0
    momentum: float = 0.5
    max_rollout_depth: int = 3
    eta: float = 0.1

    seed: int = 1
    gamma: float = 0.9
    mode: str = "default"
    action_limits: List[Tuple[float, float]] = field(default_factory = lambda: [(-2, 2), (-50, 50), (-2,2), (-2,2)])
    angle_boundary: List[float] = field(default_factory = lambda: [40, 50])
    fountain_path: str = ""
    n_process: int = 3

    merged_threshold: int = 2
    width_threshold: float = 5

    def __post_init__(self):
        """Post-initialization processing to compute derived fields.

        This method computes the full paths for the cooler output file and the
        combined summits file. It also sets default output paths for each SoN parameter
        if they are not provided in the configuration.
        """
        # Compute derived paths using os.path.join for portability
        self.cooler_output = os.path.join(self.cooler_output_dir, f"{self.full_name}_oe.cooler")
        self.combined_summits_output = os.path.join(self.summits_output_dir, f"{self.full_name}_summits_merged.bed")
        # Set default output paths for each SoN parameter if not provided
        for param in self.SoN_params:
             if param.get("output") is None:
                param["output"] = os.path.join(
                    self.summits_output_dir,
                    f"{self.full_name}_height_{param['height']}kb_width_{param['width']}kb_angle_{param['angle']}"
                )          

# -----------------------------------------------------------------------------
# Pipeline Modules
# -----------------------------------------------------------------------------
class GenerateOEMatrix:
    """Module to generate the Observed/Expected (O/E) matrix for each chromosome.

    This module reads the chromosome names from the provided chromosome size file,
    then iterates over each chromosome to call an external Java tool (Juicer) that
    computes the O/E matrix. The output for each chromosome is saved according to
    a template defined in the configuration.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the GenerateOEMatrix module.

        Args:
            config (PipelineConfig): The configuration object containing paths and parameters.
        """
        self.config = config

    @log_step
    # @skip_if_exists(lambda self: self.config.oe_matrix_dir)
    def run(self) -> None:
        """Generates OE matrices for all chromosomes.

        Returns:
            str: The directory where OE matrices are stored.
        """
        with open(self.config.chrom_size, "r") as f:
            chroms = [line.split()[0] for line in f if line.strip()]
        for chrom in chroms:
            output_path = self.config.oe_output_template.format(
                oe_matrix_dir=self.config.oe_matrix_dir,
                chrom=chrom
            )
            command = (
                f"java -jar {self.config.juicer_jar} dump oe {self.config.oe_normalization} "
                f"{self.config.input_hic} {chrom} {chrom} BP {self.config.resolution} {output_path}"
            )
            logging.info("Running dump OE for chromosome %s: %s", chrom, command)
            run_single_command(command, working_dir=self.config.working_dir)
        return self.config.oe_matrix_dir

class TransformToCooler:
    """Module to convert the generated O/E matrices into a .cooler file.

    This module calls an external Python script (generate_cooler_script) to transform
    the OE matrices (stored in a directory) into a single cooler file that represents
    the contact matrix at the specified resolution.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the TransformToCooler module.

        Args:
            config (PipelineConfig): The configuration object.
        """
        self.config = config
    
    @log_step
    # @skip_if_exists(lambda self: self.config.cooler_output_dir)
    def run(self) -> str:
        """Transforms OE matrices to a cooler file.

        Returns:
            str: The path to the generated cooler file.
        """
        command = (
            f"python {self.config.generate_cooler_script} "
            f"--chrom_size_path {self.config.chrom_size} "
            f"--output {self.config.cooler_output} "
            f"--pixels_dir {self.config.oe_matrix_dir} "
            f"--resolution {self.config.resolution} "
            f"--assembly {self.config.assembly}"
        )
        logging.info("Transforming O/E matrix to .cooler with command: %s", command)
        run_single_command(command, working_dir=self.config.working_dir)
        return self.config.cooler_output

class GenerateSoNSignal:
    """Module to generate SoN (Signal-over-Noise) signals from a cooler file.

    This module iterates over a set of SoN parameter configurations and for each one,
    calls an external Python script (calSoN_script) to compute the SoN signal. The output
    file for each parameter set is stored as specified in the configuration.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the GenerateSoNSignal module.

        Args:
            config (PipelineConfig): The configuration object.
        """
        self.config = config
    
    @log_step
    # @skip_if_exists(lambda self: self.config.summits_output_dir)
    def run(self) -> List[str]:
        """Generates SoN signals for each provided parameter set.

        Returns:
            List[str]: A list of output paths for the SoN signals.
        """
        for param in self.config.SoN_params:
            # Transform gunine length to pixel length
            
            command = (
                f"python {self.config.calSoN_script} {self.config.cooler_output} "
                f"--chrom {self.config.chrom} "
                f"--chrom_size_path {self.config.chrom_size} --output {param['output']} "
                f"--height {param['height']} --width {param['width']} "
                f"--angle {param['angle']} --merged_threshold {self.config.merged_threshold} --width_threshold {self.config.width_threshold} "
                f"--resolution {self.config.resolution} "
            )
            logging.info("Generating SoN signal with command: %s", command)
            run_single_command(command, working_dir=self.config.working_dir)

class CombineSummits:
    """Module to combine summits from multiple SoN signal outputs.

    This module reads the SoN signal output files (assumed to be in BED format)
    for different parameter sets, then calls a function from an external module
    (utils.find_summits.SummitsFinder.combine_summits) to merge overlapping or adjacent
    summits based on specified thresholds. The merged result is saved to a combined output file.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the CombineSummits module.

        Args:
            config (PipelineConfig): The configuration object.
        """
        
        import pandas as pd

        self.combine_summits = SummitsFinder.combine_summits
        self.SoN_summits_path: List[str] = []
        self.config = config

        for param in self.config.SoN_params:
            self.SoN_summits_path.append(param['output']+'.bed')

        self.summits_dfs = [pd.read_csv(path, sep='\t') for path in self.SoN_summits_path]

    @log_step
    def run(self) -> None:
        """Combines multiple summits and saves the merged result.

        Returns:
            str: The path to the combined summits output file.
        """
        summits = self.combine_summits(
            dfs = self.summits_dfs, threshold = self.config.summits_threshold, padding = self.config.padding
        )

        combined_output = self.config.combined_summits_output
        summits.to_csv(combined_output, sep='\t', index=False, header=False)
        logging.info("Summits combined and saved to %s", combined_output)
    
class Planning:
    """Module for the planning step (e.g., fountains identification).

    This module instantiates a planning algorithm (here assumed to be Fun2Planning
    from an external module) with various parameters from the configuration. It then
    runs the planning process in parallel and saves the results.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the Planning module.

        Args:
            config (PipelineConfig): The configuration object.
        """
        self.config = config
        self.mcts_planner = Fun2Planning(
            weight_name = False,
            clr_path = self.config.cooler_output,
            boundary_margin = self.config.boundary_margin,
            expansion_margin = self.config.expansion_margin,
            chrom_size_path = self.config.chrom_size,
            summits_path = self.config.combined_summits_output,
            resolution = self.config.resolution,
            ini_height = self.config.ini_height,
            ini_width = self.config.ini_width,
            ini_angle = self.config.ini_angle,
            layer_height = self.config.layer_height,
            edge_width = self.config.edge_width,
            max_iter = self.config.max_iter,
            exploration_constant = self.config.exploration_constant,
            alpha = self.config.alpha,
            action_limits = self.config.action_limits,
            es_popsize = self.config.es_popsize,
            es_sigma = self.config.es_sigma,
            es_lr = self.config.es_lr,
            gamma = self.config.gamma,
            momentum = self.config.momentum,
            mode = self.config.mode,
            angle_boundary=self.config.angle_boundary,
            seed = self.config.seed,
            eta = self.config.eta
        )
    
    @log_step
    def run(self) -> None:
        """Runs the planning algorithm in parallel and saves the results."""
        self.mcts_planner.run_parallel(
            save_path = self.config.fountain_path,
            processes = self.config.n_process
        )
        

class MasterPipeline:
    """Master pipeline that orchestrates the complete workflow.

    This pipeline executes the following steps sequentially:
        1. Generate O/E matrix.
        2. Transform O/E matrix to a cooler file.
        3. Generate SoN signals.
        4. Combine summits.
        5. (Optional) Run planning.

    Each step checks if its output already exists (via decorators) to enable
    "breakpoint resume" behavior.
    """
    def __init__(self, config: PipelineConfig):
        """Initializes the MasterPipeline with the given configuration.

        Args:
            config (PipelineConfig): The configuration object.
        """
        self.config = config
        # Ensure necessary directories exist
        for path in [
            self.config.working_dir,
            self.config.oe_matrix_dir,
            self.config.cooler_output_dir,
            self.config.summits_output_dir
        ]:
            os.makedirs(path, exist_ok=True)
            logging.info("Ensured directory exists: %s", path)

    @log_step        
    def run(self) -> None:
        """Executes the full pipeline workflow.

        Returns:
            str: The path to the combined summits output file.
        """
        logging.info("Starting Master Pipeline execution...")

        # Step 1: Generate O/E matrix.
        oe_generator = GenerateOEMatrix(self.config)
        oe_generator.run()

        # Step 2: Transform O/E matrix to .cooler
        transformer = TransformToCooler(self.config)
        transformer.run()

        # Step 3: Generate SoN signal.
        son_generator = GenerateSoNSignal(self.config)
        son_generator.run()

        # Step 4: Combine summits.
        combiner = CombineSummits(self.config)
        combiner.run()

        # Step 5: Planning
        planner = Planning(self.config)
        planner.run()
        
        logging.info("Master Pipeline execution finished successfully.")

# -----------------------------------------------------------------------------
# Configuration File Loading
# -----------------------------------------------------------------------------
def load_config(config_file: str) -> PipelineConfig:
    """Loads pipeline configuration from a YAML file.

    Args:
        config_file (str): Path to the YAML configuration file.

    Returns:
        PipelineConfig: The populated configuration object.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file '{config_file}' not found.")
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    return PipelineConfig(**config)

# -----------------------------------------------------------------------------
# Main: Parse Arguments and Run Pipeline
# -----------------------------------------------------------------------------
def main():
    """Main function to parse arguments, load configuration, and run the master pipeline."""
    parser = argparse.ArgumentParser(
        description="Run the integrated pipeline for O/E matrix generation, "
                    "cooler transformation, SoN signal generation, summit combination, "
                    "and fountains identification."
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Path to the YAML configuration file (default: config.yaml)"
    )
    args = parser.parse_args()
    
    # Load configuration from the YAML file.
    config = load_config(args.config)
    logging.info("Loaded configuration from %s", args.config)

    # Instantiate and run the master pipeline.
    master_pipeline = MasterPipeline(config)
    final_result = master_pipeline.run()
    logging.info("Final result stored at: %s", final_result)

if __name__ == '__main__':
    main()