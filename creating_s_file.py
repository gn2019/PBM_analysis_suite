import os
import pandas as pd
import argparse

def find_header(file_path, keywords=("Block", "Column", "Row", "Name", "ID", "RefNumber")):
    try:
        with open(file_path, 'r') as file:
            for row_index, line in enumerate(file, start=0):
                if all(keyword in line for keyword in keywords):  # If a line contains all keywords
                    return row_index
        print("No header found in this file!")
        return None  # Explicitly return None if no header is found
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def main(gal_path, seqFile_path, output_name):
    header = find_header(gal_path)
    
    if header is None:
        print("Exiting: No valid header found in the provided GAL file.")
        return

    try:
        gal = pd.read_table(gal_path, header=header)
        gal = gal[gal.Block ==1]
        seqFile = pd.read_table(seqFile_path)
        gal_with_seqs = gal.merge(seqFile, how='left', left_on='ID', right_on='ProbeID').loc[:, ["Column", "Row", "Name", "ID", "Sequence"]]
        output_path = os.path.join(os.path.dirname(gal_path), output_name)
        gal_with_seqs.to_csv(output_path, index=None, sep='\t')
        print(f"Output saved to {output_path}")
    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GAL and sequence files and generate output.")
    parser.add_argument("gal_path", help="Path to the GAL file")
    parser.add_argument("seqFile_path", help="Path to the SequenceFile")
    parser.add_argument("output_name", help="Name of the output file")

    args = parser.parse_args()
    main(args.gal_path, args.seqFile_path, args.output_name)
