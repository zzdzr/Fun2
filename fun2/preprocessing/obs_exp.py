import os
import numba
import argparse
import numpy as np
from numba import njit

def parse_args():
    """
    Parse command-line arguments using argparse.
    
    Returns:
        argparse.Namespace: The parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description = "Process genomic data and generate matrices."
    )

    parser.add_argument(
        '--resolution', type = int, default = 5000,
        help='Resolution for the matrix (default: 5000)'
    )

    parser.add_argument(
        '--input_directory', type=str, required=True,
        help='Path to the directory containing input .txt files'
    )

    parser.add_argument(
        '--chrom_size_file', type = str, required=True,
        help = 'Path to chrom-size file'
    )
    parser.add_argument(
        '--output_directory', type = str, required=True,
        help = 'Path to the directory where the output .npy files will be saved'
    )

    return parser.parse_args()


class MatrixProcesser:

    def __init__(self, resolution, input_directory, chrom_size_file, output_directory):
        """
        Initializes the matrix processor with the given configuration.
        
        Args:
            resolution (int): The resolution for the matrix.
            input_directory (str): The directory containing input .txt files.
            chrom_size_file (str): Path to the chrom-size file.
            output_directory (str): The directory where output .npy files will be saved.
        """
        self.resolution = resolution
        self.input_directory = input_directory
        self.chrom_size_file = chrom_size_file
        self.output_directory = output_directory

    def extract_chrom_from_filename(self, filename):
        """
        Extracts the chromosome name from a file name. Assumes the file name format is
        'MiC-WT_merge-rp2.chrX.oe_mat.txt'.

        Args:
            filename (str): The input file name.
        
        Returns:
            str: The chromosome name extracted from the file name.
        """
        chrom = filename.split('.')[-3]
        return chrom

    def load_chrom_sizes(self):
        """
        Loads chromosome sizes from a chrom-size file and stores them in a dictionary.
        The file should be in two-column format: chromosome name and size.

        Args:
            chrom_size_file (str): The path to the chrom-size file.
        """
        self.chrom_sizes = {}
        with open(self.chrom_size_file, 'r') as f:
            for line in f:
                chrom, size = line.strip().split()
                size = int(size)
                self.chrom_sizes[chrom] = size

    def load_data(self, file_path):
        """
        Loads triplet data (i, j, value) from a text file.

        Args:
            file_path (str): Path to the text file containing the data.

        Returns:
            list: A list of tuples representing (i, j, value) for the matrix data.
        """
        data = []
        with open(file_path, 'r') as f:
            for line in f:
                i, j, value = map(float, line.strip().split())
                data.append((int(i), int(j), value))
            return data 
    
    def initialize_matrix(self, chrom):
        """
        Initializes an empty matrix for a specific chromosome based on its size.
        
        Args:
            chrom (str): The chromosome name.
        """
        matrix_size = self.chrom_sizes[chrom] // self.resolution + 1
        self.matrix = np.zeros((matrix_size, matrix_size))

    @staticmethod
    @njit
    def fill_upper_triangle(matrix, data, resolution):
        """
        Fills the upper triangular part of the matrix with the provided triplet data.

        Args:
            matrix (ndarray): The matrix to be filled.
            data (list): A list of tuples (i, j, value) representing the matrix data.
            resolution (int): The resolution used to map genomic coordinates to matrix indices.

        Returns:
            ndarray: The matrix after filling the upper triangle with the data.
        """
        for i, j, value in data:
            idx_i = i // resolution
            idx_j = j // resolution
            if idx_i <= idx_j:
                matrix[idx_i, idx_j] = value

        return matrix

    def save_npy(self, filename):
        """
        Saves the matrix to a .npy file.

        Args:
            filename (str): The path to the file where the matrix will be saved.
        """
        np.save(filename, self.matrix)

    def process(self):
        """
        Processes all .txt files in the directory, fills the matrices based on the data,
        and saves the resulting matrices as .npy files.

        """

        # load chromsome size information
        self.load_chrom_sizes()

        for filename in os.listdir(self.input_directory):
            if filename.endswith('.txt'):
                file_path = os.path.join(self.input_directory, filename)

                # extract chromosome
                chrom = self.extract_chrom_from_filename(filename)
                # load data
                data = self.load_data(file_path)

                # Initialize the matrix for chromosome
                self.initialize_matrix(chrom)

                # fill matrix 
                self.matrix = self.fill_upper_triangle(self.matrix, data, self.resolution)
                # output
                output_path = os.path.join(
                    self.output_directory, f"matrix_{chrom}.npy"
                )
                self.save_npy(output_path)

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()
    
    # Initialize the processor with command-line arguments
    processor = MatrixProcesser(
        resolution=args.resolution,
        input_directory=args.input_directory,
        chrom_size_file=args.chrom_size_file,
        output_directory=args.output_directory
    )
    
    # Run the processing
    processor.process()