import numpy as np
import argparse
import pandas as pd
from typing import List, Tuple

from fun2.preprocessing.calculate_son import SoNCalculator

def calculate_overlap(summits_df: pd.DataFrame, target_zones: pd.DataFrame) -> Tuple[int, float]:
    """
    Calculate the overlap between summit regions and target zones.

    Args:
        summits_df (pd.DataFrame): DataFrame containing summit regions.
        target_zones (pd.DataFrame): DataFrame containing target zones.

    Returns:
        Tuple[int, float]: Overlap count and overlap ratio.
    """

    overlap_count = 0

    for chrom in summits_df['chrom'].unique():
        summits_chrom = summits_df[summits_df['chrom'] == chrom]
        target_chrom = target_zones[target_zones['chrom'] == chrom]

        # Iterate through summit intervals and check overlap
        for _, zone in target_chrom.iterrows():
            overlap = summits_chrom.apply(
                lambda summit: not (summit['end'] < zone['start'] or summit['start'] > zone['end']),
                axis = 1
            )

            overlap_count += overlap.sum()

    overlap_ratio = overlap_count / len(summits_df) if len(summits_df) > 0 else 0

    return overlap_count, overlap_ratio

def grid_search(
    matrix_path: str,
    target_zones: str,
    chrom_size_path: str,
    output: str,
    heights: List[float],
    widths: List[float],
    angle: float = 45,
    resolution: int = 5000,
    merge_threshold: int = 0,
    width_threshold: int = 0
) -> pd.DataFrame:
    """
    Perform grid search over heights and widths for SoN calculations.

    Args:
        matrix_path (str): Path to the main matrix file.
        target_zones (pd.DataFrame): DataFrame of target zones.
        chrom_size_path (str): Path to the chromosome size file.
        output (str): Base path to save output files.
        heights (List[float]): List of heights to iterate over.
        widths (List[float]): List of widths to iterate over.
        angle (float): Rotation angle for the sampling box. Default is 45 degrees.
        resolution (int): Matrix resolution in base pairs. Default is 5000.
        merge_threshold (int): Threshold for merging nearby summit bins. Default is 0.
        width_threshold (int): Minimum width for summit regions. Default is 0.

    Returns:
        pd.DataFrame: DataFrame containing SoN statistics for all height-width combinations.
    """
    
    results = []
    for height in heights:
        for width in widths:
            print(f"This is height:{height}, width:{width}.")
            # Initialize SoNCalculator with command-line arguments
            calculator = SoNCalculator(
                matrix_path=matrix_path,
                control_matrix_path=None,
                chrom_size_path=chrom_size_path,
                output=output,
                height=height,
                width=width,
                angle=angle,
                resolution=resolution
            )

            # Calculate SoN values for whole chromosome
            SoN_df, summits_df, regions_df = calculator.calculate_SoN_all_chromosomes()
            # Save results to output files
            #calculator.save_results(SoN_df, summits_df, regions_df)
            SoN_df.to_csv(output + f'_height_{height}_width_{width}.bedgraph', sep = '\t', header = False, index = False)
            summits_df.to_csv(output + f'_height_{height}_width_{width}.bed', sep = '\t', header = False, index = False)
            regions_df.to_csv(output + f'_height_{height}_width_{width}_merged_regions.bed', sep = '\t', header = False, index = False)

            # calculate number summits
            num_summits = len(summits_df)
            
            # number of overlapping
            n, alpha = calculate_overlap(summits_df, target_zones)

            results.append(
                {
                    'height': height,
                    'width': width,
                    'num_summits': num_summits,
                    'n_overlapped': n,
                    'perc_overlapped': alpha,
                }
            )

            

    # Convert results to a DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv(
        output + '_grid_search_summary.csv', index=False
    )

    return results_df

def parse_args():
    parser = argparse.ArgumentParser(description="Perform SoN grid search with adjustable parameters.")
    
    parser.add_argument("--matrix_path", type=str, required=True, help="Path to the main matrix file.")
    parser.add_argument("--chrom_size_path", type=str, required=True, help="Path to the chromosome size file.")
    parser.add_argument("--output", type=str, required=True, help="Base path to save output files.")
    parser.add_argument("--target_zones", type=str, required=True, help="Path to the target zones file (CSV format).")
    parser.add_argument("--heights", type=str, required=True, help="Comma-separated list of heights for grid search.")
    parser.add_argument("--widths", type=str, required=True, help="Comma-separated list of widths for grid search.")
    parser.add_argument("--angle", type=float, default=45, help="Rotation angle for the sampling box (default: 45).")
    parser.add_argument("--resolution", type=int, default=5000, help="Matrix resolution in base pairs (default: 5000).")

    return parser.parse_args()


def main():
    args = parse_args()

    matrix_path = args.matrix_path
    chrom_size_path = args.chrom_size_path
    output = args.output
    target_zones_path = args.target_zones
    heights = [float(h) for h in args.heights.split(',')]
    widths = [float(w) for w in args.widths.split(',')]
    angle = args.angle
    resolution = args.resolution

    target_zones = pd.read_table(target_zones_path, names=['chrom','start','end'])
    grid_search(
        matrix_path=matrix_path,
        target_zones=target_zones,
        chrom_size_path=chrom_size_path,
        output=output,
        heights=heights,
        widths=widths,
        angle=angle,
        resolution=resolution
    )


if __name__ == "__main__":
    main()