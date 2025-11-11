#!/usr/bin/env python3
"""
Module for running Monte Carlo Tree Search (MCTS) on the CPU for fountain detection,
and for sampling box operations (e.g. computing intensities, rotations, etc.).

This module includes:
    - MCTSRunnerCPU: Class to load genomic matrices, run MCTS on chunks of data, and
      process results.
    - A set of numba-optimized functions for geometric computations.
    - SamplingBoxAgent: Class for representing and manipulating sampling boxes.

All functions and classes are documented using Google style docstrings.
"""

import logging
import numpy as np
import pandas as pd
import cooler
from tqdm import tqdm
from numba import njit
from typing import List, Tuple, Any

# Import external modules
from fun2.planning.mcts_continuous_action_space import MCTS
from fun2.planning.sampling_box import SamplingBoxAgent
from fun2.utils import adjust_p_values
from fun2.transform import mapping_coverage, pixel_to_genuine

# Setup global logging configuration
logging.basicConfig(
    level = logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# -----------------------------------------------------------------------------
# MCTSRunnerCPU: Main class for running MCTS on the CPU for fountain detection.
# -----------------------------------------------------------------------------
class MCTSRunnerCPU:
    """Runs Monte Carlo Tree Search (MCTS) on the CPU for fountain detection.

    This class loads a genomic contact matrix for a given chromosome from a Cooler file,
    processes summit data, divides the matrix into chunks, and for each chunk applies MCTS
    to detect fountain signals. The results are compiled into a DataFrame with adjusted p-values.

    Attributes:
        summits_path (str): Path to the summit data file.
        clr_path (str): Path to the cooler file.
        weight_name (Union[str, bool]): Weight parameter name or flag for balancing.
        chrom (str): Chromosome identifier.
        resolution (int): Resolution in base pairs.
        width (int): Initial sampling box width.
        height (int): Initial sampling box height.
        angle (int): Initial sampling box orientation angle.
        layer_height (int): Height of each layer in the sampling box.
        max_iter (int): Maximum iterations for MCTS.
        exploration_constant (float): Exploration constant used in UCB1.
        boundary_margin (int): Margin for boundary adjustment.
        expansion_margin (int): Margin for expanding the sampling region.
        alpha (float): Progressive widening exponent.
        action_limits (List[Tuple[float, float]]): Allowed ranges for action updates.
        es_popsize (int): Evolutionary strategy sample size.
        es_sigma (float): Standard deviation of noise in ES.
        es_lr (float): Initial learning rate for ES.
        momentum (float): Momentum factor for ES gradient update.
        max_rollout_depth (int): Maximum rollout depth.
        eta (float): Parameter controlling penalty for fountain extension.
        initial_lr (float): Initial dynamic learning rate.
        min_lr (float): Minimum dynamic learning rate.
        decay_factor (float): Decay factor for dynamic learning rate.
        gamma (float): Discount rate for future rewards.
    """

    def __init__(
        self,
        summits_path: str,
        clr_path: str,
        chrom: str,
        resolution: int,
        width: int,
        height: int,
        angle: int,
        layer_height: int,
        edge_width: int,
        max_iter: int,
        exploration_constant: float,
        weight_name: str = False,
        boundary_margin: int = 50,
        expansion_margin: int = 400,
        alpha: float = 0.1,
        action_limits: List[Tuple[float, float]] = [(-2, 2), (-50, 50), (-2, 2), (-1, 1)],
        angle_boundary: List[float] = [40, 50],
        es_popsize: int = 8,
        es_sigma: float = 5.0,
        es_lr: float = 2.0,
        momentum: float = 0.5,
        max_rollout_depth: int = 3,
        eta: float = 0.5,
        gamma: float = 0.9,
        seed: int = 1,
        mode: str = "default"
    ) -> None:
        
        self.summits_path: str = summits_path
        self.clr_path: str = clr_path
        self.weight_name: str = weight_name
        self.chrom: str = chrom
        self.resolution: int = resolution
        self.width: int = width
        self.height: int = height
        self.angle: int = angle
        self.layer_height: int = layer_height
        self.edge_width: int = edge_width
        self.max_iter: int = max_iter
        self.exploration_constant: float = exploration_constant
        self.boundary_margin: int = boundary_margin
        self.expansion_margin: int = expansion_margin
        self.alpha: float = alpha
        self.action_limits: list = action_limits
        self.es_popsize: int = es_popsize
        self.es_sigma: float = es_sigma
        self.es_lr: float = es_lr
        self.momentum: float = momentum
        self.max_rollout_depth: int = max_rollout_depth
        self.eta: float = eta
        self.gamma: float = gamma
        self.mode: str = mode
        self.seed:int = seed
        self.angle_boundary: List[float] = angle_boundary

        self.matrix: np.ndarray = self._load_matrices()
    

    def get_submatrix_bounds(self, bin_start: int, bin_end: int, summits: pd.DataFrame) -> Tuple[int, int]:
        """Adjusts submatrix boundaries based on proximity to summit centers.

        Args:
            bin_start (int): Starting bin index.
            bin_end (int): Ending bin index.
            summits (pd.DataFrame): DataFrame of summit data.

        Returns:
            Tuple[int, int]: Adjusted (bin_start, bin_end) values.
        """
        if not summits.empty:
            min_summit = summits['center'].min()
            max_summit = summits['center'].max()

            # if self.mode == 'default':
            if min_summit - bin_start < self.boundary_margin:
                bin_start = max(0, bin_start - self.expansion_margin)

            if bin_end - max_summit < self.boundary_margin:
                bin_end = min(bin_end + self.expansion_margin, self.matrix.shape[0])
        
        return bin_start, bin_end
    
    def _load_matrices(self) -> np.ndarray:
        """Loads the contact matrix for the specified chromosome.

        Returns:
            np.ndarray: The upper-triangular contact matrix.
        """

        clr = cooler.Cooler(self.clr_path)
        if self.weight_name == "False":
            self.weight_name = False

        matrix = clr.matrix(balance=self.weight_name).fetch(self.chrom)
        matrix = np.triu(matrix)

        return matrix

    def load_file(self) -> pd.DataFrame:
        """
        Loads summit data from file and calculates the center positions of each summit.

        Returns:
            pd.DataFrame: Dataframe containing chrom, start, end, width, and center of each summit.
        """
        summits = pd.read_table(self.summits_path, sep='\t', 
        names = ['chrom', 'start', 'end', 'width' , 'SoN'])
        summits['center'] = (summits['start'] + summits['end']) // 2 // self.resolution
        summits = summits[summits['chrom'] == self.chrom].reset_index(drop=True)

        return summits[['chrom', 'start', 'end', 'width', 'center']]

    def run_on_cpu(self, chunk_size=500) -> pd.DataFrame:
        """Processes the contact matrix in chunks using MCTS and compiles results.

        Args:
            chunk_size (int, optional): The size of each chunk in bins. Defaults to 500.

        Returns:
            pd.DataFrame: A DataFrame containing the results for each chunk, with adjusted p-values.
        """
        summits = self.load_file()
        results: List = []
        bin_start: int = 0

        total_n = len(summits)
        pbar = tqdm(total = total_n, desc="Processing summits", leave=False)
        # Process matrix in chunks until the end.
        while bin_start < self.matrix.shape[0]:
            # Adjust bin boundaries based on summit proximity.
            bin_start_adjusted, bin_end_adjusted = self.get_submatrix_bounds(
                bin_start, bin_start + chunk_size, summits = summits
            )

            # Extract submatrix for the current chunk
            submatrix = self.matrix[
                bin_start_adjusted:bin_end_adjusted, 
                bin_start_adjusted:bin_end_adjusted
            ]

            # Filter summits that fall within the current chunk.
            chunk_summits = summits[(summits['center'] >= bin_start) & (summits['center'] < bin_start + chunk_size)]
            chunk_summits = chunk_summits.copy()
            chunk_summits['center'] = chunk_summits['center'].values - bin_start_adjusted

            if not chunk_summits.empty:
                # Perform MCTS for each summit center in the chunk
                for row in chunk_summits.itertuples(index=False):
                    chrom, start, end, _, center = row
                    sampling_box = SamplingBoxAgent(
                        image = submatrix, 
                        center = (center, center),
                        width = self.width,
                        height = self.height,
                        angle = self.angle,
                        layer_height = self.layer_height,
                        edge_width = self.edge_width
                    )

                    # Execute MCTS
                    mcts_solver = MCTS(
                        sampling_box, 
                        max_iter=self.max_iter, 
                        exploration_c=self.exploration_constant,
                        alpha_pw = self.alpha,
                        action_limits = self.action_limits,
                        angle_boundary = self.angle_boundary,
                        es_pop = self.es_popsize,
                        es_sigma = self.es_sigma,
                        lr0 = self.es_lr,
                        momentum = self.momentum,
                        rollout_depth = self.max_rollout_depth,
                        gamma = self.gamma,
                        eta = self.eta,
                        seed = self.seed,
                        mode = self.mode   
                    )
                    node, reward, quality, intensity, sum_interactions, p_up, p_dn = mcts_solver.search()

                    if node is None:
                        logging.info("Skipping current summit due to invalid MCTS result.")
                        # Append results to the list
                        results.append(
                            (chrom, start, end, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
                            intensity, quality, sum_interactions, p_up, p_dn)
                        )
                    else:
                        # Extract results from the best node
                        height = pixel_to_genuine(
                            self.resolution, 
                            node.agent.height, 
                            type_ = 'height'
                        ) # notice the orientation and angle transformation

                        width = pixel_to_genuine(
                            self.resolution, 
                            node.agent.width, 
                            type_ = 'width'
                        )

                        center = (node.agent.center[0] + bin_start_adjusted) * self.resolution
                        elongation_up, elongation_down, width_up, width_down = mapping_coverage(node.agent, bin_start_adjusted, self.resolution)

                        # Append results to the list
                        results.append(
                            (chrom, start, end, center, height, width, node.agent.angle, elongation_up, elongation_down, width_up, width_down,
                            reward, intensity, quality, sum_interactions, p_up, p_dn)
                        )

                    # Release memory
                    del mcts_solver

                    # update pbar
                    pbar.update(1)

            # Move to the next chunk
            bin_start += chunk_size  

        # Compile results into a DataFrame and adjust p-values
        df = pd.DataFrame(results, columns=['chrom', 'summit_start', 'summit_end', 'identified_center', 'height', 'width', 'angle',
         'elongation_up', 'elongation_down', 'width_up', 'width_down', 'reward', 'intensity', 'quality', 'sum_interactions', 'pval_upstream', 'pval_downstream'])
        df = adjust_p_values(df)

        return df
