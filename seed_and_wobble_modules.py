import math
import re

import numpy as np


def element_builder(word, current_order, start_order, stop_order, last_dot_pos, word_length, spaced_words_Aref):
    """
    Recursive function to build all gapped patterns of a certain order for a given input word.
    Here, currentorder is the number of non-dot positions at the current function call, desiredorder
    is the ultimate number of non-dot positions that are desired.
    written by A. Philippakis, 2006

    word: word used in current function call
    current_order: number of non-. positions in current word
    start_order: minimum number of non-. positions in order for it to be stored
    stop_order: maximum number of non-. positions in order for it to be stored
    last_dot_pos: location of where last dot is
    word_length: length of starting word
    spaced_words_Aref: array that holds all of the spaced words
    """

    if current_order == stop_order:
        spaced_words_Aref.append(word)
    else:
        if current_order <= start_order:
            spaced_words_Aref.append(word)

        for i in range(last_dot_pos + 1, word_length - 1):
            current_word = "".join([word[j] if i != j else "." for j in range(word_length)])
            element_builder(current_word, current_order - 1, start_order, stop_order, i, word_length, spaced_words_Aref)


def get_kmer_truncated_area(kmer_ranks_Href, kmer_trunc_areas_Href, number_spots, keep_fraction):
    """
    Program that takes all k-mers and their ranks and, for each k-mer,
    calculates the "area" statistic after dropping a certain fraction from
    the foreground and background.
    Recall that area = (1 / (B + F)) * (average background rank - average foreground rank).
    written by A. Philippakis, 2006

    kmer_ranks_Href: ranks of every kmer
    kmer_trunc_areas_Href: will store truncated areas for every kmer
    number_spots: number of spots on the array
    keep_fraction: fraction of spots to keep.  For example 0.75 drops bottom quarter, 0.5 drops bottom half
    """
    for key in sorted(kmer_ranks_Href):
        kmer_trunc_areas_Href[key] = truncated_enrichmentU_one_array(kmer_ranks_Href[key], number_spots, keep_fraction)


def truncated_enrichmentU_one_array(foreground_Aref, n_points, keep_fraction):
    """
    Function to compute the area statistic in the case of one array, but with the
    proviso that a fraction of the bottom of the foreground and background will be dropped.
    This is based on the idea that there will be outliers in the bottom that
    one wants to avoid.

    Written by A. Philippakis, 2006.

    foreground_Aref: array holding foreground elements
    n_points: number of data points
    keep_fraction: fraction of data points to drop. For example, if we want to drop bottom half do 0.5, bottom quarter do 0.75.
    """
    foreground_Aref = sorted(foreground_Aref)
    original_fgsize = len(foreground_Aref)
    fgsize = int(original_fgsize * keep_fraction)  # Foreground size after truncation
    bgsize = int((n_points - original_fgsize) * keep_fraction)  # Background size after truncation
    rank_sum = 0
    result = 0.0

    last_rank = bgsize

    for i in range(fgsize):
        rank_sum += last_rank + 1 if foreground_Aref[i] - i > bgsize else foreground_Aref[i]
        last_rank += 1

    if fgsize > 0:
        result = ((fgsize * (fgsize + 1) / 2) + (fgsize * bgsize / 2) - rank_sum) / (fgsize * bgsize)

    return result


def convert_to_letters(number, k):
    """
    Function takes a number and a "k" and returns the ACGT-equivalent of that number
    for length k.

    written by A. Philippakis, 2006
    """
    letters = "ACGT"
    string = "".join([letters[(number >> i) & 0x3] for i in range((k - 1) * 2, -1, -2)])
    return string


def median(a):
    return np.median(a)


def get_top_X_elements(element_dict, top_X):
    """
    Function to take the hash that stores the area for every k-mer and find the
    top X values (where "X" is the passed variable $topX).  These are stored in an
    array of hashes that is sorted by these areas.

    written by A. Philippakis, 2006

    elementHref: hash storing area of every element
    topX: top X elements to return
    topelementsAref: array
    """
    topelements_list = []
    element_values = []

    for key, value in element_dict.items():
        element_values.append({
            'element': key,
            'value': value
        })

    element_values.sort(key=lambda x: x['value'], reverse=True)

    for i in range(top_X):
        topelements_list.append({
            'element': element_values[i]['element'],
            'value': element_values[i]['value']
        })

    return topelements_list


def extend_element(element, max_seed_order):
    """
    Simple function to take an element of a spaced seed and extend it by adding .'s.
    The number of positions added is such that it is possible to wobble all the
    elements up to a given max-width
    written by A. Philippakis, 2006

    element: word to extend
    max_seed_order: max seed order collected in this run
    """
    iterations = max_seed_order - len(element)
    element = "." * iterations + element + "." * iterations
    return element


def wobble_seed_rerank(kmer_ranks, center, pwm_matrix, array_file, spot_length, start_position, keep_fraction):
    """
    Function to build a pwm based on area medians
    written by A. Philippakis, 2006

    center: element to wobble
    pwm_matrix: ref to HofA storing PWM entries
    """
    pwm_matrix["A"], pwm_matrix["C"], pwm_matrix["G"], pwm_matrix["T"] = [[0 for _ in center] for _ in range(4)]

    for i, char in enumerate(center):
        if char == ".":
            continue

        A_variant, C_variant, G_variant, T_variant = get_variants(center, i)

        for var in (A_variant, C_variant, G_variant, T_variant):
            if var not in kmer_ranks:
                add_element(var, kmer_ranks, array_file, spot_length, start_position)

        collapsed_ranks = {}
        number_collapsed_obs = collapse_rankings(A_variant, C_variant, G_variant, T_variant, kmer_ranks,
                                                 collapsed_ranks)

        vals = [get_reranked_truncated_enrichmentU(
            kmer_ranks[var], collapsed_ranks, number_collapsed_obs, keep_fraction)
            for var in (A_variant, C_variant, G_variant, T_variant)]

        pwm_matrix["A"][i], pwm_matrix["C"][i], pwm_matrix["G"][i], pwm_matrix["T"][i] = vals


def get_variants(center, wobble_position):
    """
    Quick function to get all of the variants at a given position for an
    input word;  NOTE THAT THIS RETURNS VARIANTS THAT ARE LOWER IN LEXICOGRAPHIC
    ORDER WITH RESPECT TO REVERSE COMPLEMENTATION.  NOTE THAT THIS ALSO TRIMS
    OFF .'S THAT BEGIN AND END THE STRING

    written by A. Philippakis, 2006

    wobble_position: position to wobble
    """
    start_pos, stop_pos = wobble_position, wobble_position
    for j, char in enumerate(center):
        if char != ".":
            start_pos, stop_pos = min(start_pos, j), max(stop_pos, j)

    bef = center[start_pos:wobble_position]
    aft = center[wobble_position + 1:stop_pos + 1]
    A_variant = bef + "A" + aft
    C_variant = bef + "C" + aft
    G_variant = bef + "G" + aft
    T_variant = bef + "T" + aft

    A_variant = min(A_variant, rc(A_variant))
    C_variant = min(C_variant, rc(C_variant))
    G_variant = min(G_variant, rc(G_variant))
    T_variant = min(T_variant, rc(T_variant))

    return A_variant, C_variant, G_variant, T_variant


def add_element(element, kmer_ranks, array_file, spot_length, start_pos):
    """
    Function takes a k-mer that may have been ignored in the first pass, scans through the
    array, and adds it as needed.

    written by A. Philippakis, 2006
    modified by M. Berger (4/24/07) to discard positions closest to end of spot

    element: element which it is desired to get the rank orderings of
    kmerranksHref: place where rank orderings are stored
    arrayfile: name of arrayfile
    spotlength: length of spot to consider
    startpos: position from end of strand to consider, starting at 1
    """
    rank = 1
    rc_element = rc(element)
    min_element = min(element, rc_element)

    with open(array_file, "r") as array_file:
        for line in array_file:
            sequence = line.strip().split("\t")[-1][start_pos - 1: spot_length]
            if re.search(element, sequence, re.IGNORECASE) or \
                    re.search(rc_element, sequence, re.IGNORECASE):
                kmer_ranks.setdefault(min_element, []).append(rank)
            rank += 1


def collapse_rankings(Avariant, Cvariant, Gvariant, Tvariant, kmer_ranks, collapsed_ranks):
    """
    Function to take ranks of A variant, C variant G variant and T variant
    and collapse them into a hash that stores their rankings with respect to each other

    written by A. Philippakis, 2006

    kmerranks: ref to ranks
    """
    collapsed_ranks.clear()
    observed_ranks = set()

    def check_and_add_rank(variant):
        nonlocal observed_ranks
        if variant not in kmer_ranks:
            return
        for rank in kmer_ranks[variant]:
            if rank not in observed_ranks:
                observed_ranks.add(rank)

    for var in (Avariant, Cvariant, Gvariant, Tvariant):
        check_and_add_rank(var)

    observed_ranks = sorted(observed_ranks)
    for i, rank in enumerate(observed_ranks):
        collapsed_ranks[rank] = i + 1

    return len(observed_ranks)


def get_reranked_truncated_enrichmentU(variant_ranks, collapsed_ranks, number_obs, keep_fraction):
    """
    Function to get the truncated area for a given variant AFTER DOING THE RERANKING PROCEDURE

    written by A. Philippakis, 2006

    variantranksAref: array that stores the ranks in original array of a given word variant
    collapsedranksHref: stores the ranks for the wobbles at that position
    numberobs: total number of data points in the foreground and background
    """
    collapsed_median = [collapsed_ranks[rank] for rank in variant_ranks]
    return_val = truncated_enrichmentU_one_array(collapsed_median, number_obs, keep_fraction)
    return return_val


def find_minimum_info_pos(area_pwm, seed, scaling_factor):
    """
    Function to take a nascent pwm (i.e., the entries are in areas), and find the
    position of minimum information content. Here it converts areas to a boltzman
    using a scaling factor that is passed to the function.  It then computes entropy at
    each position using a uniform background, and returns the position that is the least info.
    Note that this also takes as input a seed that tells it which positions to ignore.

    written by A. Philippakis, 2006

    areapwm: nascent pwm that has areas as entries
    """
    seed_length = len(seed)
    for l in 'ACGT':
        if l not in area_pwm:
            area_pwm[l] = []
    tmp_pwm = {key: [0] * seed_length for key in area_pwm}

    if not len(seed) == len(area_pwm['A']) == len(area_pwm['C']) == len(area_pwm['G']) == len(area_pwm['T']):
        raise ValueError("pwm rows do not have the same number of entries in find_minimum_info_pos")

    min_entropy = 5  # value greater than 2 to start it at
    min_entropy_pos = None

    for i in range(seed_length):
        if seed[i] in "ACGTacgt":
            denominator = 0
            entropy = 0
            for key in area_pwm:
                tmp_pwm[key][i] = math.exp(area_pwm[key][i] * scaling_factor)
                denominator += tmp_pwm[key][i]
            for key in area_pwm:
                tmp_pwm[key][i] = tmp_pwm[key][i] / denominator
                entropy += tmp_pwm[key][i] * math.log(tmp_pwm[key][i], 2)
            entropy = 2 + entropy
            if entropy < min_entropy:
                min_entropy = entropy
                min_entropy_pos = i

    return min_entropy_pos


def extend_seed_rerank(kmer_ranks, center, minimum_info_pos, pwm, array_file, spot_length, start_pos, keep_fraction):
    """
    Variation on wobble_seed_rerank, but this takes a pwm that has already
    been built and extends it.

    written by A. Philippakis, 2006

    center: element to wobble
    minimuminfopos: position to now wobble in center when looking at the others
    pwm: ref to HofA storing PWM entries
    arrayfile: ref to arrayfile
    keepfraction: fraction of foreground and background to keep in getting area
    """
    new_center = "".join(center[i] if i != minimum_info_pos else "." for i in range(len(center)))
    for l in 'ACGT':
        if l not in pwm:
            pwm[l] = []

    for i in range(len(center)):
        if center[i] == "." and i != minimum_info_pos:
            Avariant, Cvariant, Gvariant, Tvariant = get_variants(new_center, i)

            for var in (Avariant, Cvariant, Gvariant, Tvariant):
                if var not in kmer_ranks:
                    add_element(var, kmer_ranks, array_file, spot_length, start_pos)

            collapsed_ranks = {}
            number_collapsed_obs = collapse_rankings(Avariant, Cvariant, Gvariant, Tvariant, kmer_ranks, collapsed_ranks)

            vals = [get_reranked_truncated_enrichmentU(
                kmer_ranks[var], collapsed_ranks, number_collapsed_obs, keep_fraction)
                for var in (Avariant, Cvariant, Gvariant, Tvariant)]

            pwm['A'][i], pwm['C'][i], pwm['G'][i], pwm['T'][i] = vals


def extend_seed_allpatterns_rerank(kmer_ranks, center, minimum_info_pos, pwm, array_file, spot_length, start_pos,
                                   keep_fraction, patterns):
    """
    Variation on extend_seed_rerank, but this first discards the minimum
     information position and examines all other positions corresponding to
     spaced seeds covered in the array design.
    Requires a separate file of gapped patterns specific to that array.

    adapted from extend_seed_rerank by M. Berger (4/24/07)

    center: element to wobble
    minimuminfopos: position to now wobble in center when looking at the others
    pwm: ref to HofA storing PWM entries
    arrayfile: ref to arrayfile
    keepfraction: fraction of foreground and background to keep in getting area
    patterns: reference to array of patterns covered [e.g. 11111111, 1.1111111, 11.111111, etc.]
    """
    for l in 'ACGT':
        if l not in pwm:
            pwm[l] = []

    new_center = "".join(center[i] if i != minimum_info_pos else "." for i in range(len(center)))
    center_seed = "".join('1' if l != '.' else '.' for l in new_center)

    for i, c in enumerate(center):
        if c == "." and i != minimum_info_pos:
            query_seed = trim_seed(center_seed[:i] + "1" + center_seed[i+1:])

            seed_present = False
            for pattern in patterns:
                reverse_pattern = pattern[::-1]
                if query_seed == pattern or query_seed == reverse_pattern:
                    seed_present = True
                    break

            if seed_present:
                Avariant, Cvariant, Gvariant, Tvariant = get_variants(new_center, i)

                for var in (Avariant, Cvariant, Gvariant, Tvariant):
                    if var not in kmer_ranks:
                        add_element(var, kmer_ranks, array_file, spot_length, start_pos)

                collapsed_ranks = {}
                number_collapsed_obs = collapse_rankings(Avariant, Cvariant, Gvariant, Tvariant, kmer_ranks,
                                                       collapsed_ranks)

                vals = [get_reranked_truncated_enrichmentU(
                    kmer_ranks[var], collapsed_ranks, number_collapsed_obs, keep_fraction)
                    for var in (Avariant, Cvariant, Gvariant, Tvariant)]
                pwm['A'][i], pwm['C'][i], pwm['G'][i], pwm['T'][i] = vals


def gapped_convert_to_letters(number, k, seed, seedwidth):
    """
    Very similar to convert to letters, but also takes a gapped pattern
    to generate the corresponding word

    seed: gapped pattern represented as a bit string
    seedwidth: width of gapped pattern
    """
    letters = "ACGT"
    string = ""
    seed_position_holder = 0

    for i in range((k - 1) * 2, -1, -2):
        while not (1 << seed_position_holder) & seed:
            string += "."
            seed_position_holder += 1

        string += letters[(number >> i) & 0x3]
        seed_position_holder += 1

    return string


def rc(string):
    """
    Takes a string and returns its reverse complement
    """
    string = string.translate(str.maketrans('acgtACGT', 'tgcaTGCA'))
    string = string[::-1]
    return string


def trim_seed(fullseed):
    """
    Takes a spaced seed element (e.g. ..111.1..111..) and
    trims uninformative positions from end.

    written by M. Berger (04/24/07)
    """
    return fullseed.strip('.')


def rerank_spots_unexpectbed_signals(ranked_seq_Aref, reranked_seq_Aref, pwm_Href, spot_length):
    """
    Takes a given PWM as well as a ranked list of sequences and finds
    the max-scoring match to the motif on each spot of the array.

    Let i be the rank (1=highest) of a given spot by signal
    intensity, and let j be the rank by sequence socre.
    We re-sort the spots in ascending order by (i/j).

    written by A. Philippakis, 2006; adapted by M. Berger (5/14/07)

    rankedseqAref : array of sequences and signals, ranked in descending order
    rerankedseqAref : array where output is sent after this re-ranking procedure
    pwmHref: pwm being used on this iteration
    spotlength: amount of sequence excluding primer
    """
    array_data = []
    number_spots = 0

    for seq in ranked_seq_Aref:
        number_spots += 1
        sequence = seq[1]
        sequence_score = score_sequence_against_pwm(pwm_Href, sequence, spot_length)
        array_data.append({
            'sequence': sequence,
            'observedrank': number_spots,
            'sequence_score': sequence_score,
            'rerankvalue': 0
        })

        if (number_spots % 1000) == 0:
            print(f"reviewing spot {number_spots}")

    array_data = sorted(array_data, key=lambda x: x['sequence_score'], reverse=True)
    for i, d in enumerate(array_data):
        d['rerankvalue'] = d['observedrank'] / (i + 1)

    array_data = sorted(array_data, key=lambda x: x['rerankvalue'])
    for seq, d in zip(reranked_seq_Aref, array_data):
        seq[:2] = d['rerankvalue'], d['sequence']


def score_sequence_against_pwm(pwm_Href, sequence, spot_length):
    """
    Takes input PWM and sequence, and finds the score
    of the maximum scoring sub-sequence.

    spotlength:  amount of spot to consider after removing primer
    """
    for l in 'ACGT':
        if l not in pwm_Href:
            pwm_Href[l] = []
    pwm_len = len(pwm_Href['A'])

    rc_pwm = {
        'A': pwm_Href['T'][::-1],
        'C': pwm_Href['G'][::-1],
        'G': pwm_Href['C'][::-1],
        'T': pwm_Href['A'][::-1]
    }

    max_score = rcmax_score = 0

    for i in range(pwm_len - 1, -1, -1):
        sum_score = rcsum_score = 0
        for j in range(i, pwm_len):
            if (j - i) < spot_length:
                sum_score += pwm_Href[sequence[j - i]][j]
                rcsum_score += rc_pwm[sequence[j - i]][j]
        if i == (pwm_len - 1):
            max_score = sum_score
            rcmax_score = rcsum_score
        else:
            if sum_score > max_score:
                max_score = sum_score
            if rcsum_score > rcmax_score:
                rcmax_score = rcsum_score

    for i in range(spot_length):
        sum_score = rcsum_score = 0
        for j in range(pwm_len):
            if (i + j) < spot_length:
                sum_score += pwm_Href[sequence[i + j]][j]
                rcsum_score += rc_pwm[sequence[i + j]][j]
        if sum_score > max_score:
            max_score = sum_score
        if rcsum_score > rcmax_score:
            rcmax_score = rcsum_score

    return max(max_score, rcmax_score)


def rerank_spots_unexpected_signals_gomer(ranked_seq_Aref, reranked_seq_Aref, pwm_Href, spot_length):
    """
    Takes a given *Probability* PWM as well as a ranked list of
    sequences and finds the max-scoring match to the motif on each spot
    of the array *Using the GOMER framework* adapted by Chen and Morris
    in Bioinformatics 2007 (RankMotif).

    Let i be the rank (1=highest) of a given spot by signal
    intensity, and let j be the rank by sequence socre.
    We re-sort the spots in ascending order by (i/j).

    written by A. Philippakis, 2006; adapted by M. Berger (9/5/07)

    rankedseqAref : array of sequences and signals, ranked in descending order
    rerankedseqAref : array where output is sent after this re-ranking procedure
    pwmHref: pwm being used on this iteration
    spotlength: amount of sequence to consider
    """
    array_data = []
    number_spots = 0

    for i in range(len(ranked_seq_Aref)):
        number_spots += 1
        sequence = ranked_seq_Aref[i][1]
        sequence_score = score_sequence_against_pwm_gomer(pwm_Href, sequence, spot_length)
        array_data.append({
            'sequence': sequence,
            'observedrank': number_spots,
            'sequence_score': sequence_score,
            'rerankvalue': 0
        })

        if number_spots % 1000 == 0:
            print("reviewing spot", number_spots)

    array_data = sorted(array_data, key=lambda x: x['sequence_score'], reverse=True)

    for i, d in enumerate(array_data):
        d['rerankvalue'] = d['observedrank'] / (i + 1)

    array_data = sorted(array_data, key=lambda x: x['rerankvalue'])

    for seq, d in zip(reranked_seq_Aref, array_data):
        seq[:2] = d['rerankvalue'], d['sequence']


def score_sequence_against_pwm_gomer(pwm_Href, sequence, spot_length):
    """
    Takes input *Probability* PWM and sequence, and
    finds the score of the maximum scoring sub-sequence
    using the GOMER framework adapted by Chen and Morris
    Bioinformatics 2007 (RankMotif).

    Adapted by M. Berger (9/5/07) from S. Jaeger and A. Philippakis

    spotlength:  amount of spot to consider
    """
    value_gomer = 1
    for l in 'ACGT':
        if l not in pwm_Href:
            pwm_Href[l] = []

    seq_array = list(sequence)

    pwm_len = len(pwm_Href['A'])
    rc_pwm = {
        'A': [pwm_Href['T'][pwm_len - i - 1] for i in range(pwm_len)],
        'C': [pwm_Href['G'][pwm_len - i - 1] for i in range(pwm_len)],
        'G': [pwm_Href['C'][pwm_len - i - 1] for i in range(pwm_len)],
        'T': [pwm_Href['A'][pwm_len - i - 1] for i in range(pwm_len)]
    }

    for i in range(pwm_len - 1, 0, -1):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwm_len):
            if j < i:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            elif j - i > len(seq_array) - 1:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= pwm_Href[seq_array[j - i]][j]
                prod_gomer_rc *= rc_pwm[seq_array[j - i]][j]
        value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)

    for i in range(spot_length):
        prod_gomer = 1
        prod_gomer_rc = 1
        for j in range(pwm_len):
            if j + i >= spot_length:
                prod_gomer *= 0.25
                prod_gomer_rc *= 0.25
            else:
                prod_gomer *= pwm_Href[seq_array[j + i]][j]
                prod_gomer_rc *= rc_pwm[seq_array[j + i]][j]
        value_gomer *= (1 - prod_gomer) * (1 - prod_gomer_rc)

    value_gomer = 1 - value_gomer
    return value_gomer


def collapse_dinucleotide_rankings(AAvariant, ACvariant, AGvariant, ATvariant,
                                   CAvariant, CCvariant, CGvariant, CTvariant,
                                   GAvariant, GCvariant, GGvariant, GTvariant,
                                   TAvariant, TCvariant, TGvariant, TTvariant,
                                   kmer_ranks_Href, collapsed_ranks_Href):
    """
    Function to take ranks of 16 dinucleotide variants and collapse
    them into a hash that stores their rankings with respect to each other

    kmer_ranks_Href: ref to ranks
    """
    collapsed_ranks_Href.clear()
    observed_ranks = []

    for kmerranks in [kmer_ranks_Href[AAvariant], kmer_ranks_Href[ACvariant], kmer_ranks_Href[AGvariant],
                      kmer_ranks_Href[ATvariant],
                      kmer_ranks_Href[CAvariant], kmer_ranks_Href[CCvariant], kmer_ranks_Href[CGvariant],
                      kmer_ranks_Href[CTvariant],
                      kmer_ranks_Href[GAvariant], kmer_ranks_Href[GCvariant], kmer_ranks_Href[GGvariant],
                      kmer_ranks_Href[GTvariant],
                      kmer_ranks_Href[TAvariant], kmer_ranks_Href[TCvariant], kmer_ranks_Href[TGvariant],
                      kmer_ranks_Href[TTvariant]]:
        for rank in kmerranks:
            if rank not in observed_ranks:
                observed_ranks.append(rank)

    observed_ranks.sort()

    for i, rank in enumerate(observed_ranks):
        collapsed_ranks_Href[rank] = i + 1

    return len(observed_ranks)
