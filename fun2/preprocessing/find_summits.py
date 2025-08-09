import numba
import numpy as np
import pandas as pd
from intervaltree import IntervalTree
from typing import Tuple, List
from numba.typed import List

from fun2.utils import generate_box
from fun2.validate import validate_summits_df

@numba.njit
def merge_adjacent_bins(pos_arr, resolution: int, threshold: int = 1) -> list:
    """Merges adjacent bins based on a threshold distance.

    Args:
        pos_arr (np.ndarray): Array of positions to be merged.
        resolution (int): Resolution scale for merging positions.
        threshold (int, optional): Maximum allowable distance between bins for merging. Defaults to 1.

    Returns:
        list: List of tuples, each containing the start and end of a merged bin.
    """
    idx, coordinates = 0, []

    while idx < len(pos_arr) - 1:
        start = pos_arr[idx]
        merge_idx = idx + 1
        # merge adjacent bin
        while merge_idx < len(pos_arr) and pos_arr[merge_idx] - pos_arr[merge_idx-1] <= threshold:
            merge_idx += 1
        # record end position
        end = pos_arr[merge_idx - 1]
        coordinates.append((start * resolution, (end+1) * resolution)) # consider two edges of each bin
        idx = merge_idx
    
    return list(coordinates)

class SummitsFinder:
    """Identifies summits and merges adjacent bins based on provided SoN values.

    Attributes:
        SoN_vals (np.ndarray): Array of SoN values for each position.
        resolution (int): Resolution of the matrix.
        merge_threshold (int): Threshold for merging adjacent bins.
        width_threshold (int): Minimum width for a summit.
    """

    def __init__(self, SoN_vals: np.ndarray, resolution: int, merge_threshold: int, width_threshold: int):
        """
        Args:
            SoN_vals (np.ndarray): Array of Signal-over-Noise (SoN) values.
            resolution (int): Resolution of the matrix.
            merge_threshold (int): Maximum distance between bins to consider them adjacent.
            width_threshold (int): Minimum width (in bins) for a summit to be valid.
        """
        self.SoN_vals = SoN_vals
        self.resolution = resolution
        self.merge_threshold = merge_threshold
        self.width_threshold = width_threshold

    def _get_zones(self) -> list:
        """Obtains zones with valid SoN values and merges adjacent bins.

        Returns:
            list: Merged coordinates as a list of (start, end) tuples.
        """
        valid_positions = np.where(self.SoN_vals > 0)[0]

        return merge_adjacent_bins(valid_positions, self.resolution, self.merge_threshold)
    
    def estimate_width(self) -> list:
        """Estimates widths of merged zones.

        Returns:
            list: List of widths for each merged zone.
        """
        merged_coordinates = self.get_zones()
        estimated_widths = [i[1] - i[0] for i in merged_coordinates]

        return estimated_widths

    def find_summits(self) -> tuple:
        """Finds summits within SoN values, merging adjacent bins as needed.

        Returns:
            tuple: Contains:
                - peaks (list of tuple): Coordinates of peaks (start, end) for each summit.
                - widths (list): Widths of each merged summit.
        """    
        merged_coordinates = self._get_zones()
        peaks, widths = [], []
        starts, ends = [], []
        max_SoNs = []
        
        for start, end in merged_coordinates:
            # Convert to bin indices
            s_bin, e_bin = start // self.resolution, end // self.resolution
            
            # Check if the merged region width meets the threshold
            if e_bin - s_bin > self.width_threshold:
                SoN_segments = self.SoN_vals[s_bin: e_bin+1]
                relative_peak_idx, max_SoN = SummitsFinder.find_local_maximum(SoN_segments)
                peak_idx = s_bin + relative_peak_idx
                peaks.append(
                    (peak_idx * self.resolution, (peak_idx + 1) * self.resolution)
                )
                widths.append(end - start)
                starts.append(start)
                ends.append(end)
                max_SoNs.append(max_SoN)
        
        return peaks, widths, starts, ends, max_SoNs

    @staticmethod
    def find_local_maximum(signal, window_size: int = 3) -> Tuple[int, float]:
        """Finds the index of the local maximum in a smoothed signal.

        Args:
            signal (np.ndarray): The signal array in which to find the maximum.
            window_size (int, optional): Size of the smoothing window. Defaults to 3.

        Returns:
            int: Index of the local maximum in the smoothed signal.
        """
  
        smoothed_signal =  np.convolve(
            signal, np.ones(window_size) / window_size, mode = 'same'
        )

        return np.argmax(smoothed_signal), np.max(signal)
     
    @staticmethod
    def combine_summits(dfs: List[pd.DataFrame], threshold: float, padding: int) -> pd.DataFrame:
        def build_interval_tree(df: pd.DataFrame, padding: int) -> IntervalTree:
            """Build IntervalTree for fast overlap

            Args:
                df (pd.DataFrame): dataframe of detected summits. Dataframe contains 
                columns (chrom, start, end)
                padding (int): width of padding (kb) for each of summit
                threshold (float): threshold filter for summit's SoN

            Returns:
                intervaltree: IntervalTrees for summits

            Raises:
                ValueError: If the input dataframe does not contain the required columns.
            """
            tree = IntervalTree()
            for idx, row in df.iterrows():
                chrom, start, end = row['chrom'], row['start'] - padding, row['end'] + padding
                tree[start: end] = (idx, chrom)
            
            return tree

        if not dfs or len(dfs) < 2:
            raise ValueError("At least two dataframes.")
        
        # check validation of dataframe
        valid_dfs = [validate_summits_df(df, required_columns = ['chrom', 'start', 'end']) for df in dfs]
        valid_dfs = [df[df.iloc[:, -1]>threshold].reset_index(drop=True) for df in valid_dfs]
        df_main = valid_dfs[0]
        other_dfs = valid_dfs[1:]

        interval_trees = [build_interval_tree(df, padding) for df in other_dfs]
        valid_indices = []

        for idx_main, row in df_main.iterrows():
            chrom_main, start_main, end_main = row['chrom'], row['start'], row['end']
            range_start, range_end = start_main, end_main

            valid = True
            for df, tree in zip(other_dfs, interval_trees):
                overlaps = tree[range_start:range_end]

                current_values = [
                    df.iloc[overlap.data[0], -1]
                    for overlap in overlaps if overlap.data[1] == chrom_main
                ]

                if not current_values or max(current_values) <= threshold:
                    valid = False
            
            if valid:
                valid_indices.append(idx_main)
        
        return df_main.loc[valid_indices].reset_index(drop=True)