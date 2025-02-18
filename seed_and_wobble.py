##################################################################
### seed_and_wobble.pl
###
### Takes as input a two-column list of all "combinatorial" probe
###   intensities and sequences.  MUST BE PRE-ORDERED FROM BRIGHTEST
###   TO DIMMEST.
### Also takes file of gapped patterns (e.g., 11111.111) to
###   consider as candidate seeds.
### Also takes file of all gapped patterns (e.g., 111..11.1.11) that
###   are covered evenly on the universal array design.  This may be
###   the same as or longer than the above file of candidate seeds.
###
### Outputs list of all possible k-mers and their corresponding
###   enrichment score, median signal, and Z score (separate file for each seed).
### Also outputs single integrated file for all seeds above a cutoff.
### Also outputs file of top "N" position weight matrices.
###
### M. Berger 05/02/07 (adapted from A. Philippakis)
###   Modified 04/04/08
##################################################################
import os
import argparse

from seed_and_wobble_modules import *

# Parameters
SPOT_LENGTH = 36  # total number of nucleotides in probe after common primer
START_POSITION = 2  # position from the end of probe to consider, starting at 1
ESCORE_CUTOFF = 0.25  # to store for the integrated list of top k-mers from all seeds
PRINT_SEEDS = True  # "yes" if the output file for each seed is desired
PRINT_TOP_SEEDS = True  # "yes" if the output file for top seeds (above cutoff) is desired
TOP_N = 3  # number of top k-mer seeds to use for seed-and-wobble PWM construction


def parse_args():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Seed and Wobble Script")

    # Add the command-line arguments
    parser.add_argument("intensity_file", type=str, help="PBM data file (sorted by intensities) | file that ends with '.._Normalized_deBruijn.txt'")
    parser.add_argument("order", type=int, help="Width of k-mers to inspect | '7' for 8x15k, '8' for 4x44k or 8x60k")
    parser.add_argument("seed_pattern_file", type=str, help="File List of candidate seed patterns (e.g., 'patterns_7of9.txt') or number of gaps to generate automaticaly, default=2")
    parser.add_argument("total_pattern_file", type=str, help="File List of all covered gapped patterns (e.g., 'patterns_8x15k_all_7mer.txt') or number of gaps to generate automaticaly, default=8 (for 7-mers put 6)")
    parser.add_argument("output_prefix", type=str, help="Output file prefix")

    # Parse the command-line arguments
    args = parser.parse_args()
    return args


def generate_patterns(number_of_gaps, order):
    assert number_of_gaps.isnumeric(), "number of gaps must be a positive number"
    number_of_gaps = int(number_of_gaps)
    patterns = []
    for i in range(int('1' * order, 2), 2 ** (order + number_of_gaps), 2):
        pattern = bin(i)[2:]
        if pattern.count("1") == order:
            patterns.append(pattern.replace("0", "."))
    return patterns


def read_patterns(patterns_file, order):
    """
    Read in list of gapped patterns to consider for seeds
    """
    if os.path.exists(patterns_file):
        seed_patterns = list(map(str.strip, open(patterns_file, "r")))
        if any(lambda x: len(x.replace('.', '')) != order, seed_patterns):
            raise ValueError(
                    f"Number of positions in seed in {patterns_file} does not agree with order = {order}")
        return seed_patterns
    else:
        print(f"{patterns_file} doesn't exists, try to read as number of gaps")
        return generate_patterns(patterns_file, order)


def read_intensity(intensity_file):
    """
    Read in intensities and create array of sequences
    """
    data_matrix = []
    with open(intensity_file, 'r') as f:
        for line in f:
            line_data = line.strip().split("\t")
            data_matrix.append([float(line_data[0]), line_data[1]])  # spot intensity, spot sequence
            if len(data_matrix) > 1 and data_matrix[-2][0] < data_matrix[-1][0]:
                raise ValueError("Probes in the input file are not sorted from brightest to dimmest.")
    return data_matrix


def remove_duplicates(seed_patterns):
    """
    Remove duplicate seeds from list of seeds, reversed words are duplicates
    """
    unique_patterns = []
    for pattern in seed_patterns:
        if pattern not in unique_patterns and pattern[::-1] not in unique_patterns:
            unique_patterns.append(pattern)
    return unique_patterns


def seed_and_wobble(data_matrix, order, seed_patterns, output_prefix):
    """
    Calculate median intensity and enrichment score and Z score for each 8-mer
    """
    seed_patterns = remove_duplicates(seed_patterns)

    number_spots_array = len(data_matrix)  # number of spots on array
    tally = 0
    top_kmer_areas = {}
    top_kmer_medians = {}
    top_kmer_zscores = {}
    kmer_ranks = {}
    output_file_1 = ""
    keep_fraction = 0.5

    while tally <= len(seed_patterns) - 1:

        if PRINT_SEEDS or tally == 0:
            if PRINT_SEEDS:
                output_file_1 = f"{output_prefix}_{order}mers_{seed_patterns[tally]}.txt"
            else:
                output_file_1 = f"{output_prefix}_{order}mers.txt"

            with open(output_file_1, "w") as output1:
                top_label = f"{order}-mer"
                output1.write(f"{top_label}\t{top_label}\tE-score\tMedian\tZ-score\n")

        kmer_ranks.clear()  # data structure that stores ranks of every k-mer
        kmer_intensities = {}  # data structure that stores intensities of every k-mer
        kmer_areas = {}  # stores area for each gapped k-mer

        # get the ranks of every k-mer of the spaced seed in seed_patterns
        observed_patterns = {}
        spaced_seed = seed_patterns[tally]
        rev_spaced_seed = spaced_seed[::-1]
        spaced_seed_width = len(spaced_seed)

        print(f"Currently on spaced seed: {spaced_seed}")

        for spot_number, key in enumerate(data_matrix):
            spot_sequence = key[1]
            for i in range(START_POSITION - 1, SPOT_LENGTH - spaced_seed_width + 1):
                current_string = spot_sequence[i:i + spaced_seed_width]
                fwd_element = ""
                rev_element = ""

                for counter in range(spaced_seed_width):
                    if spaced_seed[counter] == "1":
                        if len(current_string) > counter:
                            fwd_element += current_string[counter]
                    else:
                        fwd_element += "."
                rc_fwd_element = rc(fwd_element)

                zero_element = fwd_element if fwd_element < rc_fwd_element else rc_fwd_element
                observed_patterns[zero_element] = 0

                if spaced_seed != rev_spaced_seed:
                    for counter in range(spaced_seed_width):
                        if rev_spaced_seed[counter] == "1":
                            if len(current_string) > counter:
                                rev_element += current_string[counter]
                        else:
                            rev_element += "."
                    rc_rev_element = rc(rev_element)

                    zero_element = rev_element if rev_element < rc_rev_element else rc_rev_element
                    observed_patterns[zero_element] = 0

            for key in observed_patterns:
                if key not in kmer_ranks:
                    kmer_ranks[key] = []
                    kmer_intensities[key] = []
                kmer_ranks[key].append(spot_number + 1)
                kmer_intensities[key].append(data_matrix[spot_number][0])

            observed_patterns.clear()

            if not spot_number % 1000:
                print(spot_number)

        get_kmer_truncated_area(kmer_ranks, kmer_areas, number_spots_array, keep_fraction)

        ###########################################################
        ### Calculate median of all log (median), median absolute deviation
        ###########################################################

        log_median = []
        bit_seed = spaced_seed.replace(".", "0")[::-1]
        decimal_bit_seed = int(bit_seed, 2)

        for k in range(4 ** order):
            word = gapped_convert_to_letters(k, order, decimal_bit_seed, spaced_seed_width)
            rev_comp = rc(word)

            if word <= rev_comp or spaced_seed != rev_spaced_seed:
                if word > rev_comp:
                    word = rev_comp
                if kmer_intensities.get(word):
                    median_intensity = median(kmer_intensities[word])
                    log = math.log(median_intensity)
                    log_median.append(log)

        median_log_median = median(log_median)
        deviations = []

        for log in log_median:
            dev = abs(log - median_log_median)
            deviations.append(dev)

        median_abs_dev = median(deviations)

        ##############################################################
        ### Print: Word / Enrichment / Median Intensity / Z-Score
        ##############################################################

        with open(output_file_1, "a") as output1:
            for k in range(4 ** order):
                word = gapped_convert_to_letters(k, order, decimal_bit_seed, spaced_seed_width)
                rev_comp = rc(word)

                if word <= rev_comp or spaced_seed != rev_spaced_seed:
                    if PRINT_SEEDS or tally == 0:
                        output1.write(f"{word}\t{rev_comp}\t")
                    if word > rev_comp:
                        word = rev_comp
                    if word in kmer_intensities:
                        median_intensity = median(kmer_intensities[word])
                        log = math.log(median_intensity)
                        zscore = (log - median_log_median) / (1.4826 * median_abs_dev)

                        if PRINT_SEEDS or tally == 0:
                            output1.write(f"{kmer_areas[word]:.5f}\t{median_intensity:.2f}\t{zscore:.4f}\n")
                        if kmer_areas[word] > ESCORE_CUTOFF:
                            top_kmer_areas[word] = kmer_areas[word]
                            top_kmer_medians[word] = median_intensity
                            top_kmer_zscores[word] = zscore
                    else:
                        if PRINT_SEEDS or tally == 0:
                            output1.write("NA\tNA\tNA\n")

        tally += 1

    return top_kmer_areas, top_kmer_medians, top_kmer_zscores


def find_top_n_seeds(top_kmer_areas, top_kmer_medians, top_kmer_zscores, output_prefix, order, top_n=TOP_N):
    """
    Find top N seeds (adapted from A. Philippakis)
    """
    print(f"Finding top {top_n} seeds.\n")

    top_n_elements = []
    element_values = []
    top_n_counter = 0

    for key, value in top_kmer_areas.items():
        element_values.append({
            'element': key,
            'value': value,
            'median': top_kmer_medians[key],
            'zscore': top_kmer_zscores[key]
        })
        top_n_counter += 1

    if top_n_counter < top_n:
        top_n = top_n_counter

    element_values = sorted(element_values, key=lambda x: x['value'], reverse=True)

    for N in range(top_n):
        top_n_elements.append({
            'element': element_values[N]['element'],
            'value': element_values[N]['value']
        })

    if PRINT_TOP_SEEDS:
        outputfile2 = f"{output_prefix}_{order}mers_top_enrichment.txt"
        with open(outputfile2, "w") as OUTPUT2:
            toplabel = f"{order}-mer"
            OUTPUT2.write(f"{toplabel}\t{toplabel}\tE-score\tMedian\tZ-score\n")
            N = 0
            while N < len(element_values):
                rc_element = rc(element_values[N]['element'])
                OUTPUT2.write(
                    f"{element_values[N]['element']}\t{rc_element}\t{element_values[N]['value']:.5f}\t{element_values[N]['median']:.2f}\t{element_values[N]['zscore']:.4f}\n")
                N += 1

    return top_n_elements


def construct_pwm(intensity_file, top_n_elements, top_kmer_areas, total_patterns, output_prefix, order):
    """
    Seed-and-Wobble PWM construction
    """
    area_pwm = {}
    kmer_ranks = {}
    array = intensity_file

    outputfile3 = f"{output_prefix}_{order}mers_pwm.txt"
    with open(outputfile3, "w") as output3:
        for z, top_elem in enumerate(top_n_elements):
            ranking = z + 1
            print(f"Currently on element ranked: {ranking}")
            seed = top_elem['element']
            print(f"{seed}\t{top_elem['value']}")
            output3.write(f"{ranking}\t{seed}\t{top_kmer_areas[seed]}\n\n")
            seed = "......." + seed + "......."
            area_pwm[seed] = {}
            wobble_seed_rerank(kmer_ranks, seed, area_pwm[seed], array, SPOT_LENGTH, START_POSITION, 1)

            minimum_info_pos = find_minimum_info_pos(area_pwm[seed], seed, math.log(10000))

            extend_seed_allpatterns_rerank(kmer_ranks, seed, minimum_info_pos, area_pwm[seed], array, SPOT_LENGTH,
                                           START_POSITION, 1, total_patterns)

            while len(area_pwm[seed]['A']) > 0 and area_pwm[seed]['A'][0] == 0 and \
                    len(area_pwm[seed]['C']) > 0 and area_pwm[seed]['C'][0] == 0 and \
                    len(area_pwm[seed]['G']) > 0 and area_pwm[seed]['G'][0] == 0 and \
                    len(area_pwm[seed]['T']) > 0 and area_pwm[seed]['T'][0] == 0:
                area_pwm[seed]['A'].pop(0)
                area_pwm[seed]['C'].pop(0)
                area_pwm[seed]['G'].pop(0)
                area_pwm[seed]['T'].pop(0)

            # Convert Perl while loop to Python while loop for removing trailing zeros
            while len(area_pwm[seed]['A']) > 0 and area_pwm[seed]['A'][-1] == 0 and \
                    len(area_pwm[seed]['C']) > 0 and area_pwm[seed]['C'][-1] == 0 and \
                    len(area_pwm[seed]['G']) > 0 and area_pwm[seed]['G'][-1] == 0 and \
                    len(area_pwm[seed]['T']) > 0 and area_pwm[seed]['T'][-1] == 0:
                area_pwm[seed]['A'].pop()
                area_pwm[seed]['C'].pop()
                area_pwm[seed]['G'].pop()
                area_pwm[seed]['T'].pop()

            output3.write("Enrichment score matrix\n\n")
            for key in sorted(area_pwm[seed]):
                output3.write(f"{key}:")
                for y in range(len(area_pwm[seed][key])):
                    output3.write(f"\t{area_pwm[seed][key][y]}")
                output3.write("\n")

            output3.write("\nEnergy matrix for enoLOGOS\n\n")
            output3.write("PO")
            for counter in range(1, len(area_pwm[seed]['A']) + 1):
                output3.write(f"\t{counter}")
            output3.write("\n")
            for key in sorted(area_pwm[seed]):
                output3.write(f"{key}:")
                for y in range(len(area_pwm[seed][key])):
                    logscaled = area_pwm[seed][key][y] * (-math.log(10000))
                    output3.write(f"\t{logscaled}")
                output3.write("\n")

            #	print output3 "\nReverse complement matrix for enoLOGOS\n\n";
            #	print output3 "PO";
            #	for (my $counter=1; $counter<=($#{$area_pwm{$seed}{A}}+1); $counter++) {
            #		print output3 "\t$counter";
            #	}
            #	print output3 "\n";
            #	foreach $key (sort keys %{$area_pwm{$seed}}) {
            #		my $compkey;
            #		if ($key eq "A") {$compkey="T";}
            #		if ($key eq "C") {$compkey="G";}
            #		if ($key eq "G") {$compkey="C";}
            #		if ($key eq "T") {$compkey="A";}
            #		print output3 "$compkey:";
            #		for (my $y=$#{$area_pwm{$seed}{$key}}; $y>=0; $y--) {
            #		    my $logscaled = $area_pwm{$seed}{$key}[$y]*(-log(10000));
            #		    print output3 "\t$logscaled";
            #		}
            #		print output3 "\n";
            #	}

            output3.write("\nProbability matrix\n\n")
            for key in sorted(area_pwm[seed]):
                output3.write(f"{key}:")
                for y in range(len(area_pwm[seed][key])):
                    numerator = math.exp(math.log(10000) * area_pwm[seed][key][y])
                    denominator = math.exp(math.log(10000) * area_pwm[seed]['A'][y]) + math.exp(
                        math.log(10000) * area_pwm[seed]['C'][y]) + math.exp(
                        math.log(10000) * area_pwm[seed]['G'][y]) + math.exp(math.log(10000) * area_pwm[seed]['T'][y])
                    probability = numerator / denominator
                    output3.write(f"\t{probability}")
                output3.write("\n")
            output3.write("\n\n\n")


def main():
    patterns = read_patterns(args.seed_pattern_file, args.order)
    data_matrix = read_intensity(args.intensity_file)
    top_kmer_areas, top_kmer_medians, top_kmer_zscores = \
        seed_and_wobble(data_matrix, args.order, patterns, args.output_prefix)
    top_n_elements = find_top_n_seeds(top_kmer_areas, top_kmer_medians, top_kmer_zscores, args.output_prefix,
                                      args.order)
    # Read in list of all gapped patterns covered in universal PBM design (for wobble step)
    total_patterns = read_patterns(args.total_pattern_file, args.order)
    construct_pwm(args.intensity_file, top_n_elements, top_kmer_areas, total_patterns,
                  args.output_prefix, args.order)


if __name__ == '__main__':
    args = parse_args()
    main()
