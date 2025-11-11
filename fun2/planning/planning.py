import os
import argparse
import pandas as pd
import multiprocessing as mp
from typing import Dict, List, Tuple
from functools import partial

from fun2.planning.mcts import MCTSRunnerCPU
from fun2.transform import transform

class Fun2Planning:

    def __init__(
        self,
        clr_path: str,
        weight_name: str,
        chrom_size_path: str,
        summits_path: str,
        resolution: int,
        ini_height: float,
        ini_width: float,
        ini_angle: float,
        layer_height: int,
        edge_width: int,
        max_iter: int,
        exploration_constant: float,
        alpha: float,
        action_limits: List[Tuple[float, float]],
        es_popsize: int,
        es_sigma: float,
        es_lr: float,
        momentum: float,
        angle_boundary: List[float] = [40, 50],
        boundary_margin: int = 50,
        expansion_margin: int = 400,
        max_rollout_depth: int = 3,
        eta: float = 0.5,
        gamma: float = 0.9,
        mode: str = "default",
        seed: int = 1,
    ) -> None:

        self.clr_path: str = clr_path
        self.weight_name: str = weight_name
        self.chrom_size_path: str = chrom_size_path
        self.summits_path: str = summits_path
        self.resolution: int = resolution
        self.ini_angle: float = ini_angle
        self.ini_height: float = transform(self.resolution, ini_height, type_ = 'height')
        self.ini_width: float = transform(self.resolution, ini_width, type_ = 'width')
        self.layer_height: int = layer_height
        self.edge_width: int = edge_width
        self.max_iter: int = max_iter
        self.exploration_constant: float = exploration_constant
        self.alpha: float = alpha
        self.action_limits: List[Tuple[float, float]] = action_limits
        self.es_popsize: int = es_popsize
        self.es_sigma: float = es_sigma
        self.es_lr: float = es_lr
        self.momentum: float = momentum

        self.boundary_margin: int = boundary_margin
        self.expansion_margin: int = expansion_margin
        self.max_rollout_depth: int = max_rollout_depth
        self.eta: float = eta
        self.gamma: float = gamma
        self.mode: str = mode
        self.seed: int = seed
        self.angle_boundary: List[float] = angle_boundary
        self.chrom_sizes: Dict[str, int] = self._load_chrom_sizes()


    def _load_chrom_sizes(self) -> Dict[str, int]:
        """Loads chromosome sizes from the provided file.

        Returns:
            Dict[str, int]: Mapping from chromosome names to sizes.
        """
        chrom_sizes: Dict[str, int] = {}
        with open(self.chrom_size_path, 'r') as f:
            for line in f:
                chrom, size = line.strip().split()
                chrom_sizes[chrom] = int(size)
        return chrom_sizes
    
    def _process_single_chrom(self, chrom: str, save_path: str) -> pd.DataFrame:
        """Processes a single chromosome using MCTSRunnerCPU.

        Args:
            chrom (str): Chromosome identifier.
            save_path (str): Path to save chromosome-specific results.

        Returns:
            pd.DataFrame: Fountain detection results for the given chromosome.
        """
        
        print(f"Processing chromosome: {chrom}...")
        mcts_cpu = MCTSRunnerCPU(
            summits_path=self.summits_path,
            clr_path=self.clr_path,
            weight_name=self.weight_name,
            chrom=chrom,
            resolution=self.resolution,
            width=self.ini_width,
            height=self.ini_height,
            angle=self.ini_angle,
            layer_height=self.layer_height,
            edge_width = self.edge_width,
            max_iter=self.max_iter,
            exploration_constant=self.exploration_constant,
            boundary_margin=self.boundary_margin,
            expansion_margin=self.expansion_margin,
            alpha=self.alpha,
            action_limits=self.action_limits,
            es_popsize=self.es_popsize,
            es_sigma=self.es_sigma,
            es_lr = self.es_lr,
            momentum=self.momentum,
            max_rollout_depth=self.max_rollout_depth,
            eta=self.eta,
            gamma=self.gamma,
            mode = self.mode,
            angle_boundary = self.angle_boundary,
            seed = self.seed,

        )
        fountains = mcts_cpu.run_on_cpu()
        tmp_save_name = f"{os.path.splitext(save_path)[0]}_{chrom}.txt"
        fountains.to_csv(tmp_save_name, sep='\t', header=True, index=False)
        print(f"Results for chromosome {chrom} saved to {tmp_save_name}")
        return fountains

    def run_parallel(self, save_path: str, processes: int = None) -> None:
        """Runs fountain detection in parallel for all chromosomes and saves the final results.

        Args:
            save_path (str): Path to save the final concatenated results.
            processes (int, optional): Number of parallel processes. Defaults to CPU count.
        """
        chrom_list: List[str] = list(self.chrom_sizes.keys())
        with mp.Pool(processes=processes or mp.cpu_count()) as pool:
            func = partial(self._process_single_chrom, save_path=save_path)
            results: List[pd.DataFrame] = pool.map(func, chrom_list)

        fountains_total = pd.concat(results, ignore_index=True)
        fountains_total.to_csv(save_path, sep = '\t', header=True, index=False)
        print(f"Final fountain detection results saved to {save_path}")

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate SoN values and detect summits for genomic matrices."
    )
    parser.add_argument("--clr_path", type=str, required=True, help="Path to the main .cool matrix file.")
    parser.add_argument("--weight_name", type=str, default=None, help="Name of the balancing weight to apply.")
    parser.add_argument("--chrom_size_path", type=str, required=True, help="Path to the chromosome size file.")
    parser.add_argument("--summits_path", type=str, required=True, help="Base path for output files.")
    parser.add_argument("--resolution", type=int, required=True, help="Matrix resolution in base pairs.")
    parser.add_argument("--ini_height", type=float, required=True, help="Initial height of the sampling box.")
    parser.add_argument("--ini_width", type=float, required=True, help="Initial width of the sampling box.")
    parser.add_argument("--ini_angle", type=float, required=True, help="Initial angle for sampling box orientation.")
    parser.add_argument("--layer_height", type=int, default=5, help="Height of each layer (default: 5).")
    parser.add_argument("--max_iter", type=int, default=500, help="Maximum MCTS iterations (default: 500).")
    parser.add_argument("--exploration_constant", type=float, default=1.5, help="Exploration constant (default: 1.5).")
    parser.add_argument("--alpha", type=float, required=True, help="Progressive widening exponent.")
    parser.add_argument("--es_popsize", type=int, required=True, help="ES sample size.")
    parser.add_argument("--es_sigma", type=float, required=True, help="ES noise standard deviation.")
    parser.add_argument("--es_lr", type=float, default=2.0, help="Learning Ratio for ES.")
    parser.add_argument("--momentum", type=float, default=0.5, help="Momentum factor for ES (default: 0.5).")
    parser.add_argument("--boundary_margin", type=int, default=50, help="Boundary margin for submatrix adjustment (default: 50).")
    parser.add_argument("--expansion_margin", type=int, default=400, help="Expansion margin for submatrix adjustment (default: 400).")
    parser.add_argument("--max_rollout_depth", type=int, default=3, help="Maximum rollout depth (default: 3).")
    parser.add_argument("--eta", type=float, default=0.5, help="Penalty parameter for fountain extension (default: 0.5).")
    parser.add_argument("--gamma", type=float, default=0.9, help="Discount rate for future rewards (default: 0.9).")
    parser.add_argument("--save_path", type=str, required=True, help="Path to save the final concatenated fountain detection results.")
    parser.add_argument("--processes", type=int, default=None, help="Number of parallel processes (default: CPU count).")
    parser.add_argument("--mode", type=str, default=None, help="Mode to used: default for fountains and fixed_angle for stripes.")
    parser.add_argument("--edge_width", type=int, default=None, help="Edge_width used to split edge of sampling box.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for MCTS planning.")
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    planner = Fun2Planning(
        clr_path=args.clr_path,
        weight_name=args.weight_name,
        chrom_size_path=args.chrom_size_path,
        summits_path=args.summits_path,
        resolution=args.resolution,
        ini_height=args.ini_height,
        ini_width=args.ini_width,
        ini_angle=args.ini_angle,
        layer_height=args.layer_height,
        edge_width=args.edge_width,
        max_iter=args.max_iter,
        exploration_constant=args.exploration_constant,
        alpha=args.alpha,
        action_limits=[(-2, 2), (-50, 50), (-2, 2), (-2, 2)],
        es_popsize=args.es_popsize,
        es_sigma=args.es_sigma,
        es_lr=args.es_lr,
        momentum=args.momentum,
        boundary_margin=args.boundary_margin,
        expansion_margin=args.expansion_margin,
        max_rollout_depth=args.max_rollout_depth,
        eta=args.eta,
        gamma=args.gamma,
        seed=args.seed,
        mode=args.mode
    )
    
    # Choose parallel run
    planner.run_parallel(save_path=args.save_path, processes=args.processes)