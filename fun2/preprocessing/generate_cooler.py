import re
import os
import cooler
import argparse
import numpy as np
import pandas as pd

class CoolerGenerator:

    def __init__(self, chrom_size_path: str, pixels_dir: str, resolution: int) -> None:
        self.chrom_size_path = chrom_size_path
        self.pixels_dir = pixels_dir
        self.resolution = resolution

    def _load_sparse_matrices(self) -> pd.DataFrame:
        """load sparse matrix and generate pixel-pairs 

        Returns:
            pd.DataFrame: dataframe containing pixel-pairs for all chroms
        """
        sparse_matrix_li = [file for file in os.listdir(self.pixels_dir) if file.endswith('.txt')]

        pixels = []
        for file in sparse_matrix_li:
            # get chromsome information
            chrom_match = re.search(r"\bchr[0-9XYM]+\b", file)
            if chrom_match:
                chrom = chrom_match.group()
                file_path = os.path.join(self.pixels_dir, file)

                with open(file_path, 'r') as f:
                    for line in f:
                        start, end, count = line.strip().split()
                        pixels.append((chrom, int(start), int(end), float(count)))
        
        pixels = pd.DataFrame(pixels, columns = ['chrom', 'start', 'end', 'count'])
        
        return pixels

    def _load_chrom_sizes(self) -> dict:
        """Loads chromosome sizes from a chrom-size file.

        Returns:
            dict: Chromosome names and their sizes.
        """
        chrom_sizes = {}
        with open(self.chrom_size_path, 'r') as f:
            for line in f:
                chrom, size = line.strip().split()
                chrom_sizes[chrom] = int(size)
        return chrom_sizes

    def _generate_bins(self, chrom_sizes: dict) -> pd.DataFrame:
        """Generate dataframe containing bins

        Args:
            chrom_sizes (dict): dict containing chromosome names and their sizes

        Returns:
            pd.DataFrame: pd.DataFrame: dataframe containing bins
        """
        bins = []
        for chrom, size in chrom_sizes.items():
            bins.append(pd.DataFrame({
                'chrom': [chrom] * ((size // self.resolution)+1),
                'start': range(0, size, self.resolution),
                'end': range(self.resolution, size + self.resolution, self.resolution)
            }))

        bins = pd.concat(bins, ignore_index=True)

        return bins
 
    def _generate_refined_pixels(self, bins: pd.DataFrame, pixels: pd.DataFrame) -> pd.DataFrame:
        """transform coordinates(start/end) to index of bin

        Args:
            bins (pd.DataFrame): dataframe containing bins
            pixels (pd.DataFrame): dataframe containing pixel-value pair

        Returns:
            pd.DataFrame: dataframe containing transformed pixel-value pair
        """
        bins = bins.reset_index().rename(
            columns = {'index': 'bin_id'}
        )

        pixels = pixels.merge(
            bins[['chrom', 'start', 'bin_id']],
            left_on = ['chrom', 'start'],
            right_on = ['chrom', 'start'],
            how = 'left'
        ).rename(columns = {'bin_id': 'bin1_id'})

        pixels = pixels.merge(
            bins[['chrom', 'start', 'bin_id']],
            left_on = ['chrom', 'end'],
            right_on = ['chrom', 'start'],
            how = 'left'
        ).rename(columns = {'bin_id':'bin2_id'})

        
        if len(pixels)>0:
            pixels = pixels.loc[:, ['bin1_id', 'bin2_id', 'count']]
            
        return pixels

    def generate_cooler_obj(self, bins: pd.DataFrame, pixels: pd.DataFrame, output: str, assembly: str) -> cooler:
        """generate cooler object from pixels dataframe

        Args:
            bins (pd.DataFrame): dataframe containing bins
            pixels (pd.DataFrame): dataframe containing pixel-value pairs
            output (string): path for output file
            assembly (string): genome version

        Returns:
            cooler: cooler object containing matrix for all chromosomes
        """
        
        cooler.create_cooler(
            cool_uri = output,
            bins = bins,
            pixels = pixels,
            assembly = assembly,
            dtypes={'count': 'float64'}
        )

def main(args) -> None:
    """generate cooler object

    Args:
        args (parser.parse_args): parameters help to generate cooler object
    """
    generator = CoolerGenerator(
        args.chrom_size_path, args.pixels_dir, args.resolution
    )
    
    chrom_sizes = generator._load_chrom_sizes()
    pixels = generator._load_sparse_matrices()
    bins = generator._generate_bins(chrom_sizes)

    refined_pixels = generator._generate_refined_pixels(bins, pixels)

    generator.generate_cooler_obj(
        bins = bins, pixels = refined_pixels, output = args.output, assembly = args.assembly 
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate .cool file from sparse matrix data.")
    parser.add_argument("--chrom_size_path", type=str, required=True, 
                        help="Path to the chromosome sizes file.")
    parser.add_argument("--pixels_dir", type=str, required=True, 
                        help="Directory containing sparse matrix files.")
    parser.add_argument("--resolution", type=int, required=True, 
                        help="Resolution (bin size) for the Hi-C data.")
    parser.add_argument("--output", type=str, required=True, 
                        help="Path for the output .cool file.")
    parser.add_argument("--assembly", type=str, default="hg19", 
                        help="Genome assembly version (default: hg19).")

    args = parser.parse_args()
    main(args)