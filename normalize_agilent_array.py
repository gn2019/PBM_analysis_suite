#!/usr/bin/env python3
"""
(c)2006-2008 The Brigham and Women's Hospital, Inc and President and Fellows of Harvard College
"""

import sys
import argparse
from statistics import median

import pandas as pd

from analyze_agilent import cy3_regression, alexa488_normalize, median_neighborhood

###########################################################
###
### normalize_agilent_array.pl
### Normalize and spatially detrend single subgrid of 4x44K
###   or 8x15K array
###
### M. Berger, 12/12/2006
###   modified 2/10/2007
###   modified 7/5/2007
###   modified 4/3/2008
###
### Takes as input (1) an Alexa488 PBM GPR file (raw or
###   Masliner-adjusted), (2) optionally the Cy3 DNA
###   GPR file from the corresponding microarray subgrid, and
###   (3) the probe ID's and sequences, and outputs the following:
###
### (1) "Raw data" file
### 	Column 1: Column
### 	Column 2: Row
### 	Column 3: Spot Name
### 	Column 4: Spot ID
### 	Column 5: Sequence
### 	Column 6: Cy3 ADJBSI (if Cy3 not included, 1)
### 	Column 7: Cy3 Flags (if Cy3 not included, 0)
### 	Column 8: Alexa488 ADJBSI
### 	Column 9: Alexa488 Flags
### (Row 1 is a header. Data starts in row 2.)
###
### (2) "All data" file [same as above but with calculated values for
###     expected Cy3; obs/exp Cy3; normalized Alexa488; spatially-
###     adjusted Alexa488]
###
### (3) "Combinatorial" file [ranked list of normalized intensities and
###	sequences of combinatorial 'all k-mer' probes, with control spots
###     removed. This can be used directly as input for k-mer scoring
###     and motif finding software.]
###
### (4) "Regression" file [results of linear regression over Cy3
###     probe intensities and sequences, including R^2 value indicating
###     quality of fit]
###
### Performs the following steps to normalize and
###    spatially detrend the signal on the array...
### Step 1: Calculate (Observed/Expected) Cy3
### Step 2: Normalize Alexa488 by (Obs/Exp) Cy3
### Step 3: Spatially adjust normalized PBM
###
### Array format (4x44K or 8x15K) is inferred from the sequence file
###
############################################################


############################################################
## 'Wavelength' numerical values from GenePix 6.0 software
##    Should match column in GPR file labeled 'F Median - B'
############################################################

DIM_TO_ARRAY_TYPE = {
    (170, 266): "4x44k",
    (164, 96): "8x15k",
    (328, 192): "8x60k",
    (340, 532): "4x180k",
}
ARRAY_TYPE_TO_DIM = {v: k for k, v in DIM_TO_ARRAY_TYPE.items()}

COLS = [
    "Column",
    "Row",
    "Name",
    "ID",
    "Sequence",
    "Cy3",
    "Cy3Flags",
    "Alexa488",
    "Alexa488Flags",
    "Cy3Exp",
    "Obs/Exp",
    "Alexa488Norm",
    "Alexa488Median",
    "Alexa488Adjusted"
]


#################################################
### Read in PBM, Cy3 information for each protein
#################################################

def read_sequence_file(sequence_file):
    data_matrix = {}
    with open(sequence_file, 'r') as sequence_file_handle:
        next(sequence_file_handle)  # Skip header
        for line in sequence_file_handle:
            entry = line.strip('\r\n ').split('\t')
            col, row = map(int, entry[:2])  # Column, Row
            if col not in data_matrix:
                data_matrix[col] = {}
            data_matrix[col][row] = [entry[2], entry[3], entry[4]]  # Probe Name, Probe ID, Sequence

    last_column, last_row = col, row

    array_type = DIM_TO_ARRAY_TYPE.get((last_row, last_column), None)
    if array_type:
        print(f"Array format recognized as '{array_type}'.")
    else:
        sys.exit("Cannot recognize array format by dimensions")

    return data_matrix, array_type


def read_input_file(input_file, data_matrix, array_type, gp, is_cy3=False):
    print(f"Reading {'cy3' if is_cy3 else 'Input'} File {input_file}.")
    fluorescence_col = f'F{gp} Median - B{gp}' if args.background else f'F{gp} Median'

    spot_check = False
    cols_dict = {
        'flag': 0,
        'adjbsi': 0,
        'bsi': 0,
    }

    # cy3_spots = []

    with open(input_file, 'r') as pbm_file_handle:
        for line in pbm_file_handle:
            line = line.strip().split('\t')
            if spot_check and line[1]:
                col, row = map(int, line[1:3])
                cell = data_matrix[col][row]
                if not is_cy3:
                    cell.extend([1, 0])  # Cy3 BSI and Flags default (overwrite later)
                    if cols_dict['adjbsi'] > 0:
                        cell.append(int(line[cols_dict['adjbsi']]))  # Alexa488 ADJBSI
                    else:
                        cell.append(int(line[cols_dict['bsi']]))  # Cy3 BSI
                    cell.append(int(line[cols_dict['flag']]))  # Alexa488 Flags
                else:
                    if cols_dict['adjbsi'] > 0:
                        cell[3] = int(line[cols_dict['adjbsi']])  # Cy3 ADJBSI
                    else:
                        cell[3] = int(line[cols_dict['bsi']])  # Cy3 BSI
                    cell[4] = int(line[cols_dict['flag']])  # Cy3 Flags
                    # in source: cy3 file only
                    # if data_matrix[col][row] > -100 and \
                    #         ("dBr" in data_matrix[col][row][1] or "Ctrl" in data_matrix[col][row][1]):
                    #     cy3_spots.append(data_matrix[col][row][3])

            if "Block" in line[0]:
                spot_check = True
                for tally in range(len(line)):
                    if line[tally] == "Flags" or line[tally] == '"Flags"':
                        cols_dict['flag'] = tally
                    if line[tally] == "ADJBSI" or line[tally] == '"ADJBSI"':
                        cols_dict['adjbsi'] = tally
                    if line[tally] == fluorescence_col or line[tally] == f'"{fluorescence_col}"':
                        cols_dict['bsi'] = tally
                if cols_dict['flag'] < 1:
                    raise ValueError(f"Unable to find 'Flags' column in {input_file}")
                if cols_dict['adjbsi'] < 1 and cols_dict['bsi'] < 1:
                    raise ValueError(
                        f"Unable to find column labeled 'ADJBSI' or '{fluorescence_col}' in {input_file}")
    if (row, col) != ARRAY_TYPE_TO_DIM[array_type]:
        raise ValueError(f"Columns and rows are uneven between sequence file and {input_file}")


def perform_cy3_regression(output_prefix, data_matrix, array_type):
    """
    ###################################################
    ### STEP 1: COMPUTE REGRESSION COEFFICIENTS FOR CY3
    ###################################################

    Use combinatorial (unflagged) spots to calculate regression coeffecients:
    intercept, AAAA, AAAC, AAAG, AAAT, ..., ATTT
    """
    final_index = 6
    regression_order = 3
    regression_output = output_prefix + "_regression.txt"

    if args.cy3_file:
        print(f"Computing Expected Cy3 Signals. 4^{regression_order} Parameters.")
        final_index = cy3_regression(data_matrix, regression_order, array_type, regression_output, final_index)
    else:
        last_row, last_col = ARRAY_TYPE_TO_DIM[array_type]
        for row in range(1, last_row + 1):
            for col in range(1, last_col + 1):
                data_matrix[col][row].extend([1, 1])
        final_index = final_index + 2

    return final_index


def perform_alexa488_normalization(data_matrix, array_type, final_index):
    """
    Appends expected Cy3 intensity to data matrix.
    Calculates (Observed/Expected) Cy3 signal for all non-Agilent sequences.
    If no Cy3 normalization, appends (Observed/Expected) = 1.

    ###################################################
    ### STEP 2: NORMALIZE ALEXA488 BY CY3
    ###################################################
    """
    if args.cy3_file:
        print("Normalizing Alexa488 Intensities by Cy3.")
    norm_alexa488 = []
    final_index = alexa488_normalize(data_matrix, norm_alexa488, array_type, final_index)
    alexa488_median = median(norm_alexa488)
    return alexa488_median, final_index


def perform_median_neighborhood_analysis(data_matrix, array_type, alexa488_median, final_index):
    """
    ###################################################
    ### STEP 3: ADJUST ALEXA488 BY LOCATION ON ARRAY
    ###################################################

    Calculate median neighborhood (block radius R) for each spot.
    """
    radius = 7
    print(f"Calculating Neighborhood Median for Alexa488. Radius = {radius}.")
    final_index = median_neighborhood(data_matrix, radius, "Alexa488", array_type, alexa488_median, final_index)
    return final_index


### Appends block median to data matrix.
### For Alexa488, divides by block median and multiplies by array median.

###################################################
### PRINT RESULTS TO OUTPUT FILES
###################################################

def write_output_files(output_prefix, data_matrix, final_index, array_type):
    print(f"Printing to output files for {output_prefix}.\n")

    intensity = {}
    all_output = f"{output_prefix}_alldata.txt"
    combinatorial_output = f"{output_prefix}_combinatorial.txt"
    raw_output = f"{output_prefix}_rawdata.txt"

    with open(all_output, "w") as all_output_handler, open(raw_output, "w") as raw_output_handler:
        all_output_handler.write('\t'.join(COLS) + '\n')
        raw_output_handler.write('\t'.join(COLS[:9]) + '\n')

        last_row, last_col = ARRAY_TYPE_TO_DIM[array_type]

        for j in range(1, last_row + 1):
            for k in range(1, last_col + 1):
                all_output_handler.write(f"{k}\t{j}")
                raw_output_handler.write(f"{k}\t{j}")

                for item in range(final_index + 1):
                    all_output_handler.write(f"\t{data_matrix[k][j][item]}")
                    if item <= 6:
                        raw_output_handler.write(f"\t{data_matrix[k][j][item]}")

                all_output_handler.write("\n")
                raw_output_handler.write("\n")

                if "dBr" in data_matrix[k][j][1] and data_matrix[k][j][final_index] != "NA":
                    intensity[data_matrix[k][j][2]] = data_matrix[k][j][final_index]

    sorted_spots = sorted(intensity, key=lambda x: intensity[x], reverse=True)

    with open(combinatorial_output, "w") as combinatorial_output_handler:
        for spot_sequence in sorted_spots:
            combinatorial_output_handler.write(f"{intensity[spot_sequence]}\t{spot_sequence}\n")


def create_sequence_file(tdt_file, seqlist_file, output_prefix):
    tdt = pd.read_csv(tdt_file, sep='\t', header=0)
    seqlist = pd.read_csv(seqlist_file, sep='\t', header=0)
    seqlist.rename(columns={'ProbeID': 'ID'}, inplace=True)
    # seqlist.loc[seqlist['Sequence'] == '#N/A', 'Sequence'] = ''

    seqs = pd.merge(tdt[['Column', 'Row', 'Name', 'ID']], seqlist[['ID', 'Sequence']], on='ID', how='left')
    fn = f'{output_prefix}_sequence_file.txt'
    seqs.to_csv(fn, sep='\t', index=False)
    return fn


def parse_args():
    parser = argparse.ArgumentParser(description='''Normalize and spatially detrend Agilent microarray data.

    EXAMPLE: perl normalize_agilent_array.pl -i pbm.gpr -c cy3.gpr -s 4x44k_v1_sequences.txt -o ./prot_name
    Will output the following files: (1) ./prot_name_rawdata.txt
                                     (2) ./prot_name_alldata.txt
                                     (3) ./prot_name_combinatorial.txt
                                     (4) ./prot_name_regression.txt''')
    parser.add_argument('-i', required=True, dest='pbm_file', help='''.GPR file from GenePix (either raw or adjusted by Masliner)
    If generated by Masliner, will use the column labeled 'ADJBSI'.
    Otherwise, should have a column labeled 'F1 Median - B1'.''')
    parser.add_argument('-c', dest='cy3_file', help='''.GPR file from GenePix for Cy3 channel (optional, either exact output or adjusted by Masliner)
    If generated by Masliner, will use the column labeled 'ADJBSI'.
    Otherwise, should have a column labeled 'F1 Median - B1'.
    **If this file is not included, only spatial detrending (not Cy3 normalization) will be performed.''')

    seq_file_subparser = parser.add_subparsers(dest='seq_file_subparser', required=True)

    # Subparser for requiring one argument
    seq_file_parser = seq_file_subparser.add_parser('s')
    seq_file_parser.add_argument('sequence_file', type=str, help='''Text file containing probe sequences and coordinates for the corresponding array design.
    Contains five columns (with header): Column / Row / Probe Name / Probe ID / Sequence''')

    # Subparser for requiring two arguments
    tdt_seqlist_parser = seq_file_subparser.add_parser('ts')
    tdt_seqlist_parser.add_argument('tdt_file', type=str, help='TDT file downloaded from Agilent (DNAFront_BCBottom)')
    tdt_seqlist_parser.add_argument('seqlist_file', type=str, help='SequenceList file downloaded from Agilent')

    # parser.add_argument('-s', required=True, dest='sequence_file', help='''Text file containing probe sequences and coordinates for the corresponding array design.
    # Contains five columns (with header): Column / Row / Probe Name / Probe ID / Sequence''')

    parser.add_argument('-o', required=True, dest='output_prefix', help='''Prefix for output filenames.
    If priting files to a subdirectory, the subdirectory must already exist''')
    parser.add_argument('--cy3_gp', type=int, default=1, help='Cy3 wavelength, relevant to column names (default = 1).')
    parser.add_argument('--pbm_gp', type=int, default=1, help='PBM wavelength, relevant to column names (default = 1).')
    parser.add_argument('-b', '--background', action='store_true', help='Consider "F? Median - B?" column instead of "F? Median"')
    args = parser.parse_args()
    return args


def main():
    if args.seq_file_subparser == 'ts':
        args.sequence_file = create_sequence_file(args.tdt_file, args.seqlist_file, args.output_prefix)
    data_matrix, array_type = read_sequence_file(args.sequence_file)

    read_input_file(args.pbm_file, data_matrix, array_type, args.pbm_gp, is_cy3=False)
    if args.cy3_file:
        read_input_file(args.cy3_file, data_matrix, array_type, args.cy3_gp, is_cy3=True)
    final_index = perform_cy3_regression(args.output_prefix, data_matrix, array_type)
    alexa488_median, final_index = perform_alexa488_normalization(data_matrix, array_type, final_index)
    final_index = perform_median_neighborhood_analysis(data_matrix, array_type, alexa488_median, final_index)
    write_output_files(args.output_prefix, data_matrix, final_index, array_type)


if __name__ == '__main__':
    args = parse_args()
    main()
