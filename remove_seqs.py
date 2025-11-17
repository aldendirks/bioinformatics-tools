#!/usr/bin/env python3

"""
Remove sequences from a multifasta file smaller than a specified length. 

Author: Alden Dirks
Date: November 14, 2025
Version: 0.1

"""

import argparse
import os

def remove_short_sequences(fasta_file, length):
    with open(fasta_file) as file:
        headers = []
        sequences = []
        current_seq = []
        
        for line in file:
            if line.startswith('>'): 
                headers.append(line.strip())
                if current_seq:
                    sequences.append(''.join(current_seq))
                    current_seq = []
            else: 
                current_seq.append(line.strip())

        if current_seq:
            sequences.append(''.join(current_seq))
    
    total_sequences = len(sequences)
    kept = 0
    removed = 0

    base = os.path.splitext(fasta_file)[0]  # removes file extension
    output_path = base + "_length-filtered.fasta"

    with open(output_path, "w") as outfile:
        for header, sequence in zip(headers, sequences):
            if len(sequence) >= length:
                outfile.write(f"{header}\n{sequence}\n")
                kept += 1
            else:
                removed += 1
    
    print(f"Input file: {fasta_file}")
    print(f"Output file: {output_path}")
    print(f"Minimum length: {length}")
    print(f"Total sequences: {total_sequences}")
    print(f"Kept: {kept}")
    print(f"Removed: {removed}")

def main():
    parser = argparse.ArgumentParser(
        description="Remove sequences smaller than a specified length from a multifasta file."
    )
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("length", type=int, help="Minimum length of sequences to keep")

    args = parser.parse_args()

    remove_short_sequences(args.fasta_file, args.length)

if __name__ == "__main__":
    main()