import argparse
import cooler
import numpy as np
import pandas as pd

from fun2.utils import rotate_point, compute_mean_intensity, generate_box
from fun2.transform import transform
from fun2.preprocessing.find_summits import SummitsFinder

class SoNCalculator:
    """
    Calculator for Signal-over-Noise (SoN) values in a matrix and outputs results in `.bedgraph` and `.bed` files.
    """

    def __init__(self, matrix_path: str, control_matrix_path: str = None, weight_name: str = None, chrom_size_path: str = None, output: str = None,
                 height: int = 200, width: int = 5, angle: float = 45, resolution: int = 5000, merged_threshold: int = 2, width_threshold: float = 5):
        """
        Initialize SoNCalculator with paths to input files and calculation parameters.
        
        Args:
            matrix_path (str): Path to the main matrix file.
            control_matrix_path (str, optional): Path to the control matrix file.
            chrom_size_path (str, optional): Path to the chromosome size file.
            output (str): Base path to save output files.
            height (int): Height of the sampling box.
            width (int): Width of the sampling box.
            angle (float): Rotation angle for sampling box.
            resolution (int): Matrix resolution.
        """
        self.matrix_path = matrix_path
        self.control_matrix_path = control_matrix_path
        self.weight_name = weight_name
        # self.chrom = chrom
        self.chrom_size_path = chrom_size_path
        self.output = output
        self.resolution = resolution
        self.angle = angle
        self.height = transform(self.resolution, height, type_ = 'height')
        self.width = transform(self.resolution, width, type_ = 'width')
        self.rad_angle = np.deg2rad(self.angle)
        self.merged_threshold = merged_threshold
        self.width_threshold = width_threshold

        # Load chromosome sizes
        self.chrom_sizes = self._load_chrom_sizes() if chrom_size_path else {}
        
    def _load_chrom_sizes(self) -> dict:
        """Load chromosome sizes from a chrom-size file.

        Returns:
            dict: Chromosome names and their sizes.
        """
        chrom_sizes = {}
        with open(self.chrom_size_path, 'r') as f:
            for line in f:
                chrom, size = line.strip().split()
                chrom_sizes[chrom] = int(size)
        return chrom_sizes

    def _load_matrices(self, chrom: str) -> np.ndarray:
        """Load main and control matrices for a specific chromosome.

        Args:
            chrom (str): Chromosome name to fetch matrix data for.

        Returns:
            np.ndarray: Matrix with control subtracted if provided, otherwise main matrix.
        """ 
        clr = cooler.Cooler(self.matrix_path)
        clr_control = cooler.Cooler(self.control_matrix_path) if self.control_matrix_path else None

        matrix = clr.matrix(balance=self.weight_name).fetch(chrom)
        matrix = np.triu(matrix) # only triu matrix

        if clr_control:
            control_matrix = clr.matrix(balance=self.weight_name).fetch(chrom)
            subtracted_matrix = matrix - control_matrix
            subtracted_matrix[subtracted_matrix < 0] = 0
            
            return subtracted_matrix

        return matrix

    def _calculate_single_SoN(self, matrix: np.ndarray, idx: int) -> float:
        """Calculate the SoN value for a single index position.

        Args:
            matrix (np.ndarray): The matrix to calculate SoN values on.
            idx (int): Index position to calculate SoN for.

        Returns:
            float: SoN value at the given index position.
        """

        sampling_box = generate_box(
            (idx, idx), self.width, self.height, self.angle
        )

        upstream_center = rotate_point(
            idx - self.width, idx, idx, idx, self.rad_angle
        )
        downstream_center = rotate_point(
            idx + self.width, idx, idx, idx, self.rad_angle
        )

        upstream_background_box = generate_box(
            upstream_center, self.width, self.height, self.angle
        )
        downstream_background_box = generate_box(
            downstream_center, self.width, self.height, self.angle
        )

        #Compute mean intensities for signal and background
        sampling_mean = compute_mean_intensity(
            matrix, sampling_box
        )
        upstream_mean = compute_mean_intensity(
            matrix, upstream_background_box
        )
        downstream_mean = compute_mean_intensity(
            matrix, downstream_background_box
        )

        # Compute median for signal and background
        # sampling_median = compute_median_intensity(
        #     matrix, sampling_box
        # )
        # upstream_median = compute_median_intensity(
        #     matrix, upstream_background_box
        # )
        # downstream_median = compute_median_intensity(
        #     matrix, downstream_background_box
        # )

        #print(f"sampling median is {sampling_median}, upstream median is {upstream_median}, downstream median is {downstream_median}")
        
        avg_background_mean = (upstream_mean + downstream_mean) / 2 if (upstream_mean + downstream_mean) != 0 else 1e-10
        # use median for robust estimate
        # maximum_background = max(upstream_median, downstream_median)
        # if maximum_background == 0:
        #     maximum_background = 1e-10

        return sampling_mean * np.log(sampling_mean / avg_background_mean) if sampling_mean > 0 else 0

    def cal_SoN(self, chrom: str) -> np.ndarray:
        """Calculate SoN values for a specific chromosome.

        Args:
            chrom (str): Chromosome name to calculate SoN values for.

        Returns:
            np.ndarray: Array of SoN values for each index.
        """
        matrix = self._load_matrices(chrom)
        
        SoN_vals = np.zeros(matrix.shape[0])
        for idx in range(matrix.shape[0]):
            SoN_vals[idx] = self._calculate_single_SoN(matrix, idx)

        return SoN_vals
    
    def calculate_SoN_all_chromosomes(self) -> tuple:
        """Calculate SoN values and summits for all chromosomes in the chromosome size file.

        Returns:
            tuple: DataFrames for SoN values and summits across all chromosomes.
        """
        SoN_all_df = pd.DataFrame()
        summits_all_df = pd.DataFrame()
        regions_all_df = pd.DataFrame()
        
        for chrom in self.chrom_sizes.keys():
            print(f"Processing chromosome {chrom}...")
            SoN_vals = self.cal_SoN(chrom)
            SoN_df = self._format_SoN_df(SoN_vals, chrom)
            summits_df, regions_df = self._format_summits_df(SoN_vals, chrom, merged_threshold=self.merged_threshold, width_threshold=self.width_threshold)

            SoN_all_df = pd.concat([SoN_all_df, SoN_df], ignore_index=True)
            summits_all_df = pd.concat([summits_all_df, summits_df], ignore_index=True)
            regions_all_df = pd.concat([regions_all_df, regions_df], ignore_index=True)
        
        return SoN_all_df, summits_all_df, regions_all_df

    def _format_SoN_df(self, SoN_vals: np.ndarray, chrom: str) -> pd.DataFrame:
        """Format SoN values into a DataFrame for `.bedgraph` output.

        Args:
            SoN_vals (np.ndarray): Array of calculated SoN values.
            chrom (str): Chromosome name.

        Returns:
            pd.DataFrame: Formatted DataFrame for SoN values.
        """
        # Generate the .bedgraph DataFrame with chromosome, start, end, and SoN values
        SoN_vals[np.isnan(SoN_vals)] = 0

        n_bins = len(SoN_vals)
        chrom_size = self.chrom_sizes.get(
            chrom, n_bins * self.resolution
        )
        end_positions = np.clip(
            np.arange(self.resolution, (n_bins+1) * self.resolution, self.resolution),
            0, chrom_size
        )

        return pd.DataFrame({
            'chrom': [chrom] * n_bins,
            'start': np.arange(0, n_bins * self.resolution, self.resolution),
            'end': end_positions,
            'SoN': SoN_vals
        })
    
    def _format_summits_df(self, SoN_vals: np.ndarray, chrom: str, merged_threshold: int=2, width_threshold: int=5) -> pd.DataFrame:
        """Format summit regions into a DataFrame for `.bed` output.

        Args:
            SoN_vals (np.ndarray): Array of SoN values.
            chrom (str): Chromosome name.
            merged_threshold (int): Distance threshold for merging adjacent bins.
            width_threshold (int): Minimum width for summit regions.

        Returns:
            pd.DataFrame: Formatted DataFrame for summit regions.
        """
        summit_finder = SummitsFinder(SoN_vals, self.resolution, merged_threshold, width_threshold)
        merged_summits, widths, starts, ends, max_SoNs = summit_finder.find_summits()

        return pd.DataFrame({
            'chrom': [chrom] * len(merged_summits), 
            'start': [s[0] for s in merged_summits], 
            'end': [s[1] for s in merged_summits], 
            'width': widths,
            'SoN': max_SoNs
        }), pd.DataFrame({
            'chrom': [chrom] * len(merged_summits),
            'start': starts,
            'end': ends,
            'SoN': max_SoNs
        })

    def save_results(self, SoN_df: pd.DataFrame, summits_df: pd.DataFrame, regions_df: pd.DataFrame) -> None:
        """Save the SoN and summits data to `.bedgraph` and `.bed` files.

        Args:
            SoN_df (pd.DataFrame): DataFrame with SoN values.
            summits_df (pd.DataFrame): DataFrame with summit regions.
        """
        SoN_df.to_csv(self.output + '.bedgraph', sep = '\t', header = False, index = False)
        summits_df.to_csv(self.output + '.bed', sep = '\t', header = True, index = False)

        print(f"SoN values saved to {self.output}.bedgraph and summits saved to {self.output}.bed")

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
    parser.add_argument("--merged_threshold", type=int, default=2, help="Threshold for merging nearby summit bins.")
    parser.add_argument("--width_threshold", type=float, default=5, help="Minimum width for summit regions.")
    # parser.add_argument("--chrom", "-c", type=str, default=None, help="If set, only process this chromosome (must be present in your --chrom_size_path).")

    return parser.parse_args()

def main():
    """Main function to initialize SoNCalculator and calculate SoN matrix."""
    args = parse_args()
    
    # Initialize SoNCalculator with command-line arguments
    calculator = SoNCalculator(
        matrix_path=args.matrix_path,
        control_matrix_path=args.control_matrix_path,
        chrom_size_path=args.chrom_size_path,
        # chrom = args.chrom,
        output=args.output,
        height=args.height,
        width=args.width,
        angle=args.angle,
        resolution=args.resolution,
        width_threshold=args.width_threshold,
        merged_threshold=args.merged_threshold
    )

    # if args.chrom:
    #     chrom = args.chrom
    #     SoN_vals = calculator.cal_SoN(chrom)
    #     SoN_df = calculator._format_SoN_df(SoN_vals, chrom)
    #     summits_df, regions_df = calculator._format_summits_df(SoN_vals, chrom)

    # else:
        # Calculate SoN values for whole chromosome
    SoN_df, summits_df, regions_df = calculator.calculate_SoN_all_chromosomes()
        # Save results to output files
    
    calculator.save_results(SoN_df, summits_df, regions_df)

if __name__ == "__main__":
    main()








