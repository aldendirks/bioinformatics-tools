#!/usr/bin/env python3

"""
Print a sequence from a multifasta file by a specified position. 

Author: Alden Dirks
Date: November 14, 2025
Version: 0.1

"""

import argparse

def parse_positions(position_args):
    """
    Convert a list of position strings or ranges into a sorted list of unique positions.
    E.g., ['1', '3', '5-7'] -> [1, 3, 5, 6, 7]
    """
    positions = set()
    for arg in position_args:
        if '-' in arg:
            start, end = arg.split('-')
            positions.update(range(int(start), int(end) + 1))
        else:
            positions.add(int(arg))
    return sorted(positions)

def read_fasta(fasta_file):
    """Read a FASTA file and return headers and sequences lists."""
    headers = []
    sequences = []
    current_seq = []

    with open(fasta_file) as file:
        for line in file:
            # If it's a header, add header to headers list, 
            # join sequences in current_seq list together
            # (if any sequences in current_seq), 
            # append that concatenated sequence to the 
            # sequences list, and reset current_seq list
            if line.startswith('>'):
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
                headers.append(line.strip())
            # If it's a sequence line, add to current_seq list
            else:
                current_seq.append(line.strip())

        # Append the last sequence if exists
        if current_seq:
            sequences.append(''.join(current_seq))

    return headers, sequences

def print_sequences(fasta_file, positions):
    headers, sequences = read_fasta(fasta_file)
    for pos in positions:
        if 1 <= pos <= len(sequences):
            print(headers[pos - 1])
            print(sequences[pos - 1])
        else:
            print(f"Position {pos} is out of range. There are only {len(sequences)} sequences.")

def main():
    parser = argparse.ArgumentParser(
        description="Print sequences from a multifasta file based on its position in the file."
    )
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("positions", nargs="+", help="Sequence positions or ranges to print (e.g., 1 3 5-7)")

    args = parser.parse_args()
    positions = parse_positions(args.positions)

    print_sequences(args.fasta_file, positions)

if __name__ == "__main__":
    main()