from typing import Tuple, List
from fun2.preprocessing.find_summits import SummitsFinder

class SummitsCombine:
    def __init__(self, dfs_path: List[str], threshold: float, padding: int):
        self.combine_summits = SummitsFinder.combine_summits
        self.dfs_path = dfs_path
        self.threshold = threshold
        self.padding = padding

    def combine(self):
        summits = self.combine_summits(
            dfs = dfs_path, threshold = threshold, padding = padding
        ) 

        return summits

def parse_args():
    """Parses command-line arguments for SoNCalculator.
    
    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Calculate SoN values for a given matrix and find summits.")
    parser.add_argument("matrix_path", type=str, help="Path to the main matrix file.")
    parser.add_argument("--control_matrix_path", type=str, default=None, help="Path to the control matrix file.")
    parser.add_argument("--chrom_size_path", type=str, help="Path to the chromosome size file.")
    parser.add_argument("--output", type=str, required=True, help="Base path to save output `.bedgraph` and `.bed` files.")
    parser.add_argument("--height", type=int, default=200, help="Height of the sampling box.")
    parser.add_argument("--width", type=int, default=15, help="Width of the sampling box.")
    parser.add_argument("--angle", type=float, default=45, help="Rotation angle for the sampling box.")
    parser.add_argument("--resolution", type=int, default=5000, help="Matrix resolution in base pairs.")
    parser.add_argument("--merge_threshold", type=int, default=2, help="Threshold for merging nearby summit bins.")
    parser.add_argument("--width_threshold", type=int, default=5, help="Minimum width for summit regions.")

    return parser.parse_args()

def main():
    """Main function to initialize SoNCalculator and calculate SoN matrix."""
    args = parse_args()
    
    # Initialize SoNCalculator with command-line arguments
    calculator = SoNCalculator(
        matrix_path=args.matrix_path,
        control_matrix_path=args.control_matrix_path,
        chrom_size_path=args.chrom_size_path,
        output=args.output,
        height=args.height,
        width=args.width,
        angle=args.angle,
        resolution=args.resolution
    )

    # Calculate SoN values for whole chromosome
    SoN_df, summits_df, regions_df = calculator.calculate_SoN_all_chromosomes()
    # Save results to output files
    calculator.save_results(SoN_df, summits_df, regions_df)

if __name__ == '__main__':
    main()

