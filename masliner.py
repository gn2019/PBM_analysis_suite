#!/usr/bin/env python3
# masliner: microarray scaling through linear regression
# takes as parameters:
#  two genepix files (-g1 -g2)
#  linear range (-ll -lh) (i.e., linear low and linear high)
#  fluorophores (-f1 -f2)
#  mode (-m) (type of linear regression estimate)
#  for some modes: -inc (increment from -ll and -hl to be processed)

import re
import sys
import math
import argparse

VERSION = 1.02  # v1.01 = 11/29/2001, v1.02 = 03/25/2002

# mode identifiers
MODE_STRAIGHT = 0
MODE_DYNOPT = 1
MODE_LABEL = ["straight", "dynopt"]
MODE_DESC = [
    f"      {MODE_LABEL[MODE_STRAIGHT]}: -lh and -ll considered to define linear range directly",
    f"      {MODE_LABEL[MODE_DYNOPT]}: -lh and -ll considered to be within linear range but are\n" +
    f"         extended dynamically to a larger range using -inc and -q parameters."
]

# status identifiers
STATUS_NEWFILE = 0
STATUS_ADJFILE = 1

# parameter defaults
Dlinlow = 2000  # default linear range low
Dlinhigh = 60000  # default linear range high
Dfluor1 = 1  # default fluorophore file 1
Dfluor2 = 1  # default fluorophore file 2
Dmode = MODE_STRAIGHT
Dinc = 1000
Dq = 0.05

# DEBUG mode
DEBUG = 0


class Genepix:
    def __init__(self, id_, scanbsi, bkg, adjbsi, status, adjcount, regerr, satflag):
        self.id = id_
        self.scanbsi = scanbsi
        self.bkg = bkg
        self.adjbsi = adjbsi
        self.status = status
        self.adjcount = adjcount
        self.regerr = regerr
        self.satflag = satflag
        self.bkgstdev = None


class RegInfo:
    def __init__(self, x, y):
        if len(x) < 2:
            return

        self.N = len(x)
        self.sumx = sum(x)
        self.sumy = sum(y)
        self.avgx = self.sumx / self.N
        self.avgy = self.sumy / self.N
        self.sumx2 = sum([i * i for i in x])
        self.sumy2 = sum([i * i for i in y])
        self.sumxy = sum([i * j for i, j in zip(x, y)])
        self.sumxm2 = self.sumx2 - self.N * self.avgx * self.avgx
        self.sumym2 = self.sumy2 - self.N * self.avgy * self.avgy
        self.sumxym = self.sumxy - self.N * self.avgx * self.avgy
        self.stdevx = math.sqrt(self.sumxm2 / (self.N - 1))
        self.stdevy = math.sqrt(self.sumym2 / (self.N - 1))
        self.df = self.N - 2
        self.r = self.sumxym / (self.stdevx * self.stdevy * (self.N - 1))
        self.r2 = self.r * self.r
        self.b = self.sumxym / self.sumxm2
        self.a = self.avgy - self.b * self.avgx
        self.ssnexp = (
                self.sumy2
                + self.N * self.a * self.a
                + self.sumx2 * self.b * self.b
                - 2 * self.b * self.sumxy
                - 2 * self.a * self.sumy
                + 2 * self.a * self.b * self.sumx
        )
        self.ssexp = self.N * (self.a - self.avgy) * (self.a - self.avgy) + \
            self.b * self.b * self.sumx2 + \
            2 * (self.a - self.avgy) * self.b * self.sumx
        self.mse = math.sqrt(self.ssnexp / self.df)
        self.mse_over_avgy = self.mse / self.avgy
        self.F = self.ssexp / (self.ssnexp / self.df)
        self.stdevb = math.sqrt(self.ssnexp / (self.df * self.sumxm2))
        self.stdeva = math.sqrt(self.ssnexp * ((1 / self.N) + (self.avgx * self.avgx / self.sumxm2)) / self.df)

    def __str__(self):
        return f"""
Regression information:
  N={self.N}
  a={self.a} (regression y-intercept)
  b= {self.b} (regression slope)
  stdevb={self.stdevb}
  stdeva={self.stdeva}
  r={self.r}
  mse={self.mse} (mean unexplained error)
  mse/avgy={self.mse_over_avgy} (quality indicator [cf. -q parameter])
  R^2={self.r2} (coefficient of determination)
  ssexp={self.ssexp} (explained sum of squares)
  ssnexp={self.ssnexp} (unexplained sum of squares)
  df={self.df}
  F={self.F}
  avgx={self.avgx}
  stdevx={self.stdevx}
  sumx2={self.sumx2} (sum x^2)
  sumxm2={self.sumxm2} (sum (x-avgx)^2)
  avgy={self.avgy}
  stdevy={self.stdevy}
  sumy2={self.sumy2} (sum y^2)
  sumym2={self.sumym2} (sum (y-avgy)^2)
  sumxy={self.sumxy} (sum xy)
  sumxym={self.sumxym} (sum (x-avgx)*(y-avgy))
"""


def compute_stdev(data):
    n = 0
    total = 0
    total_squared = 0

    for value in data:
        if value != "":
            n += 1
            total += value
            total_squared += value ** 2

    if n <= 1:
        stdev = 0
    else:
        avg = total / n
        stdev = math.sqrt((total_squared - n * avg ** 2) / (n - 1))

    return stdev


def perform_adjustment(gp1, gp2):
    reginfo = None
    new_linlow = None
    new_linhigh = None
    lbkg1, lbkg2 = [], []

    # do regression analysis
    if args.mode == MODE_STRAIGHT:
        reginfo = perform_straight_regression(gp1, gp2)
        new_linlow = args.linlow
        new_linhigh = args.linhigh
    elif args.mode == MODE_DYNOPT:
        reginfo, new_linlow, new_linhigh = perform_dynopt_regression(gp1, gp2)
    else:
        print("\nPerformAdjustment internal error: Unrecognized regression mode", args.mode)
        print("Program terminating.")
        exit()

    # Now compute stdev background in linear range for both scans 1 and 2.
    # These values are used in computing initial regerr values. 11/29/2001
    for blkrowcol in gp1.id.keys():
        if new_linlow <= gp1.scanbsi[blkrowcol] <= new_linhigh:
            lbkg1.append(gp1.bkg[blkrowcol])
        if new_linlow <= gp2.scanbsi[blkrowcol] <= new_linhigh:
            lbkg2.append(gp2.bkg[blkrowcol])

    gp1.bkgstdev = compute_stdev(lbkg1)
    gp2.bkgstdev = compute_stdev(lbkg2)  # used in WriteAdjustedFile()

    # extrapolate using regression and generate output, saving adjusted values in %gp2.adjbsi.
    # %gp2.regerr, %gp2.adjcount, %gp2.satflag modified for any generated %adjvals entry.
    # 11/29/2001 - logic reformulated.  Now trigger for adjustment is %gp2.scanbsi value
    #   and %gp2.adjbsi gets (a) all adjusted values, (b) copies of unadjusted values with
    #  %adjcount set to 0.

    if DEBUG:
        print(
            "adjusted?\tblkcolrow\tid\tlow-scan-scanbsi\tlow-scan-adjcount\tlow-scan-adjbsi\t"
            "low-scan-regerr\thigh-scan-scanbsi\thigh-scan-adjcount\thigh-scan-adjbsi\t"
            "high-scan-regerr"
        )

    numadjcount = 0

    for blkcolrow, scanbsi in gp2.scanbsi.items():
        if scanbsi > new_linhigh or gp1.adjcount[blkcolrow] > 0:
            if gp1.adjcount[blkcolrow] > 0:
                gp2.adjbsi[blkcolrow], gp2.regerr[blkcolrow] = \
                    regression_prediction(reginfo, gp1.adjbsi[blkcolrow], gp1.regerr[blkcolrow])
                gp2.adjcount[blkcolrow] = gp1.adjcount[blkcolrow] + 1
                if gp1.satflag[blkcolrow] == 1:
                    gp2.satflag[blkcolrow] = 1
            else:
                gp2.adjbsi[blkcolrow], gp2.regerr[blkcolrow] = \
                    regression_prediction(reginfo, gp1.scanbsi[blkcolrow], gp1.bkgstdev)
                gp2.adjcount[blkcolrow] = 1
                if gp1.scanbsi[blkcolrow] > new_linhigh:
                    gp2.satflag[blkcolrow] = 1
            gp2.adjbsi[blkcolrow] = int(gp2.adjbsi[blkcolrow] + 0.5)
            numadjcount += 1
            adjust_ind = 1
        else:
            gp2.adjbsi[blkcolrow] = gp2.scanbsi[blkcolrow]
            gp2.regerr[blkcolrow] = gp2.bkgstdev
            adjust_ind = 0

        if DEBUG:
            print(
                adjust_ind,
                blkcolrow,
                gp2.id[blkcolrow],
                gp1.scanbsi[blkcolrow],
                gp1.adjcount[blkcolrow],
                gp1.adjbsi[blkcolrow],
                gp1.regerr[blkcolrow],
                gp2.scanbsi[blkcolrow],
                gp2.adjcount[blkcolrow],
                gp2.adjbsi[blkcolrow],
                gp2.regerr[blkcolrow],
            )

    print(f"\n{numadjcount} high intensity scan BSI values adjusted.\n")
    print("\nRewriting high intensity scan array file with adjusted values.")

    write_adjusted_file(gp2)


def perform_straight_regression(gp1, gp2):
    x = []
    y = []

    print(f"\nPerforming linear regression of data points with intensities between {args.linlow} and {args.linhigh}.")

    # collect genes within low and high range in both slides
    num_adjcount_in_range = 0
    for blkcolrow in sorted(gp1.id.keys()):
        if args.linlow <= gp2.scanbsi[blkcolrow] <= args.linhigh:
            if gp1.adjcount[blkcolrow] == 0 and args.linlow <= gp1.scanbsi[blkcolrow] <= args.linhigh:
                x.append(gp1.scanbsi[blkcolrow])
                y.append(gp2.scanbsi[blkcolrow])
            elif gp1.adjcount[blkcolrow] > 0 and args.linlow <= gp1.adjbsi[blkcolrow] <= args.linhigh:
                x.append(gp1.adjbsi[blkcolrow])
                y.append(gp2.adjbsi[blkcolrow])
                num_adjcount_in_range += 1

    n = len(x)
    if n <= 3:
        print("Too few spots (", n, ") in linear range to perform regression (n<=3).")
        print("Program terminating.")
        exit()

    print(f"\n{n} spots will be used to perform interpolation. {num_adjcount_in_range} were previously adjusted.")
    if n <= 10:
        print("Warning: Number of spots (", n, ") in linear range is small (n<=10).")

    # perform regression
    reginfo = calc_regression(x, y)
    print_regression_info(reginfo)

    return reginfo


def perform_dynopt_regression(gp1, gp2):
    regstorage = []
    regrangedata = []
    num_regressions_n_too_small = 0
    num_regressions_failing_q_thresh = 0
    num_regressions = 0
    lowest_q_over_thresh = -1

    print("\nAnalyzing data for largest BSI range of regression containing", args.linlow, "to", args.linhigh)
    print("  meeting regression quality threshold", args.q)

    # find maximum and minimum intensities >= 0
    maxint = ""
    minint = ""
    for blkcolrow in gp1.id.keys():
        if gp1.adjcount[blkcolrow] > 0:
            val1 = gp1.adjbsi[blkcolrow]
        else:
            val1 = gp1.scanbsi[blkcolrow]
        if val1 <= 0 or gp2.scanbsi[blkcolrow] <= 0:
            continue
        maxint = max(maxint, max(val1, gp2.scanbsi[blkcolrow]))
        minint = min(minint, min(val1, gp2.scanbsi[blkcolrow]))
    minint = max(minint, 0)
    maxint = max(maxint, 0)
    if args.linlow < minint:
        start_linlow = args.linlow
    else:
        if (args.linlow - minint) % args.inc == 0:
            start_linlow = minint
        else:
            start_linlow = args.linlow - ((args.linlow - minint) // args.inc + 1) * args.inc
    if maxint < args.linhigh:
        end_linhigh = args.linhigh
    else:
        if (maxint - args.linhigh) % args.inc == 0:
            end_linhigh = maxint
        else:
            end_linhigh = args.linhigh + ((maxint - args.linhigh) // args.inc + 1) * args.inc

    # get statistics
    for curr_linlow in range(start_linlow, args.linlow + 1, args.inc):
        for curr_linhigh in range(args.linhigh, end_linhigh + 1, args.inc):
            x = []
            y = []
            num_adjcount_in_range = 0
            for blkcolrow in gp1.id.keys():
                if gp1.adjcount[blkcolrow] > 0:
                    val1 = gp1.adjbsi[blkcolrow]
                else:
                    val1 = gp1.scanbsi[blkcolrow]
                if curr_linlow <= val1 <= curr_linhigh \
                        and curr_linlow <= gp2.scanbsi[blkcolrow] <= curr_linhigh:
                    x.append(val1)
                    y.append(gp2.scanbsi[blkcolrow])
                    if gp1.adjcount[blkcolrow] > 0:
                        num_adjcount_in_range += 1
            num_regressions += 1
            if len(x) <= 3:
                num_regressions_n_too_small += 1
            else:
                reginfo = calc_regression(x, y)
                if reginfo.mse_over_avgy > args.q:
                    num_regressions_failing_q_thresh += 1
                    if lowest_q_over_thresh == -1:
                        lowest_q_over_thresh = reginfo.mse_over_avgy
                    elif lowest_q_over_thresh > reginfo.mse_over_avgy:
                        lowest_q_over_thresh = reginfo.mse_over_avgy
                else:
                    regstorage.append(reginfo)
                    intrange = curr_linhigh - curr_linlow
                    regqual = reginfo.mse_over_avgy
                    regrangedata.append(
                        f"{intrange}\t{regqual}\t{num_adjcount_in_range}\t{curr_linlow}"
                        f"\t{curr_linhigh}\t{len(regstorage) - 1}"
                    )

    print("\nOptimal regression range analysis complete:")
    print("    ", num_regressions, "ranges considered.")
    print("    ", num_regressions_n_too_small, "could not be computed (N<=3).")
    print("    ", num_regressions_failing_q_thresh, "failed regression quality test.")
    print("    ", len(regstorage), "regressions passed regression quality test.")

    if len(regstorage) == 0:
        print("\nNo suitable linear regressions found.")
        if num_regressions_failing_q_thresh > 0 and lowest_q_over_thresh > -1:
            print("Best quality measure for any regression considered was", lowest_q_over_thresh)
        print("\nProgram terminating.")
        exit()

    from functools import cmp_to_key
    regrangedata.sort(key=cmp_to_key(by_range_and_qual))
    (
        intrange,
        regqual,
        num_adjcount_in_range,
        curr_linlow,
        curr_linhigh,
        i,
    ) = regrangedata[0].split("\t")
    reginfo = regstorage[int(i)]

    print("\nRegression range adjcount to", curr_linlow, "to", curr_linhigh)
    print("\n", reginfo.N, "spots used to perform interpolation. ", num_adjcount_in_range,
          "were previously adjusted.")
    if reginfo.N <= 10:
        print("Warning: Number of spots (", reginfo.N, ") in linear range is small (N<=10).")

    print_regression_info(reginfo)

    return reginfo, curr_linlow, curr_linhigh


def by_range_and_qual(a, b):
    adata = a.split("\t")
    bdata = b.split("\t")
    if adata[0] == bdata[0]:
        return int(adata[1]) - int(bdata[1])
    else:
        return int(bdata[0]) - int(adata[0])


# def max_val(a, b):
#     if a == "" and b == "":
#         return ""
#     elif a == "" and b != "":
#         return b
#     elif a != "" and b == "":
#         return a
#     return max(a, b)


# def min_val(a, b):
#     if a == "" and b == "":
#         return ""
#     elif a == "" and b != "":
#         return b
#     elif a != "" and b == "":
#         return a
#     return min(a, b)


def read_genepix_file(gfile, fluor):
    id_hash = {}
    scanbsi_hash = {}
    bkg_hash = {}
    adjbsi_hash = {}
    adjcount_hash = {}
    regerr_hash = {}
    satflag_hash = {}

    labels_found = False
    for line in gfile:
        if re.match(r'^"?Block"?\t"?Column"?\t"?Row"?\t"?Name"?', line):  # 03/25/2002 (make "" optional)
            labels_found = True
            break

    if not labels_found:
        print("Error: Genepix file", gfile, "has invalid format (no label line).")
        print("Program terminating.")
        exit()

    datalabel = "F" + str(fluor) + " Median"  # 03/25/2002
    bkglabel = "B" + str(fluor) + " Median"  # 03/25/2002

    line = line.rstrip()
    cols = re.split(r'\s*\t\s*', line)
    datacol = next((i for i, col in enumerate(cols) if re.match(r'^"?' + re.escape(datalabel) + r'"?$', col)),
                   None)  # 03/25/2002
    bkgcol = next((i for i, col in enumerate(cols) if re.match(r'^"?' + re.escape(bkglabel) + r'"?$', col)),
                  None)  # 03/25/2002

    for col, label in ((datacol, datalabel), (bkgcol, bkglabel)):
        if col is None:
            print("Error: Genepix file", gfile, "has invalid format (no", label, "column).")
            print("Program terminating.")
            exit()

    adjcountlabel = "ADJCOUNT"  # label for adjustment count column added during prior masliner run
    adjbsilabel = "ADJBSI"  # label for adjusted BSI column added during prior masliner run
    regerrlabel = "REGERR"  # label for regression error column added during prior run
    satflaglabel = "SATURATION-FLAG"  # label for saturation flag column added prior run

    if cols[-4:] == [adjcountlabel, adjbsilabel, regerrlabel, satflaglabel]:
        adjcountcol = len(cols) - 4
        adjbsicol = len(cols) - 3  # 11/20/2001
        regerrcol = len(cols) - 2
        satflagcol = len(cols) - 1
        extracols_present = True
        status = STATUS_ADJFILE
    else:
        extracols_present = False
        status = STATUS_NEWFILE

    for line in gfile:
        line = line.strip()
        if line == "":
            continue
        cols = re.split(r'\t', line)
        blkcolrow = ",".join(cols[:3])
        id_hash[blkcolrow] = cols[4]
        id_hash[blkcolrow] = id_hash[blkcolrow].lstrip('"').rstrip('"')
        scanbsi_hash[blkcolrow] = float(cols[datacol]) - float(cols[bkgcol])  # 11/29/2001
        bkg_hash[blkcolrow] = float(cols[bkgcol])  # 11/29/2001
        if extracols_present:
            adjcount_hash[blkcolrow] = float(cols[adjcountcol])
            adjbsi_hash[blkcolrow] = float(cols[adjbsicol])  # 11/29/2001
            regerr_hash[blkcolrow] = float(cols[regerrcol])
            satflag_hash[blkcolrow] = float(cols[satflagcol])
        else:
            adjcount_hash[blkcolrow] = 0
            adjbsi_hash[blkcolrow] = 0  # 11/29/2001
            regerr_hash[blkcolrow] = 0  # 11/29/2001
            satflag_hash[blkcolrow] = 0

    return Genepix(id_hash, scanbsi_hash, bkg_hash, adjbsi_hash, status, adjcount_hash, regerr_hash,
                   satflag_hash)  # 11/29/2001


def write_adjusted_file(gp):
    # writes out a version of $gfile2 that has adusted values.
    # Revised: 11/29/2001
    # Four extra columns at the end maintain the adjusted values and indicators
    args.gfile2.seek(0)
    for line in args.gfile2:
        if line.startswith('"Block"\t"Column"\t"Row"\t"Name"'):  # 03/25/2002 (make "" optional)
            break
        else:
            args.outfile.write(line)

    newlabelline = line
    endws = re.search(r'(\s*)$', newlabelline).group(1)
    newlabelline = re.sub(r'\s*$', '', newlabelline)

    cols = re.split(r' *\t *', line.rstrip())

    if gp.status == STATUS_ADJFILE:
        adjcountcol = len(cols) - 4
        adjbsicol = len(cols) - 3
        regerrcol = len(cols) - 2
        satflagcol = len(cols) - 1
    elif gp.status == STATUS_NEWFILE:
        adjcountcol = len(cols)
        adjbsicol = len(cols) + 1
        regerrcol = len(cols) + 2
        satflagcol = len(cols) + 3
    else:
        print("\nWriteAdjustedFile internal error: Unrecognized high intensity scan "
              "column status", gp.status)
        print("Program terminating.")
        exit()

    if gp.status == STATUS_NEWFILE:
        newlabelline += "\tADJCOUNT\tADJBSI\tREGERR\tSATURATION-FLAG" + endws
        args.outfile.write(newlabelline)
    else:
        args.outfile.write(newlabelline + endws)

    for line in args.gfile2:
        if line.strip() == "":
            args.outfile.write(line)
            continue

        cols = line.rstrip('\r\n').split('\t')
        blkcolrow = ",".join(cols[:3])
        cols.extend([None] * 4)
        cols[adjcountcol] = gp.adjcount[blkcolrow]
        cols[regerrcol] = f'{gp.regerr[blkcolrow]:.12f}'
        cols[adjbsicol] = int(gp.adjbsi[blkcolrow])
        cols[satflagcol] = gp.satflag[blkcolrow]
        newline = "\t".join(map(str, cols))
        args.outfile.write(newline + endws)


def calc_regression(x, y):
    # regression calculations are based on Sokal. R.R. and Rohlf, 1995, F.J.
    #  Biometry, 3rd ed., (Freeman and Company, NY), chapter 14.
    #
    return RegInfo(x, y)


def print_regression_info(reginfo):
    print(reginfo)


def regression_prediction(reginfo, x, stdevx):
    # stdevy is adapted from Sokal & Rohlf (1995), Biometry, 3rd ed.,
    #  box 14.2 (p. 468) for standard error of prediction from linear
    #  regression.  This formula differs by incorporating a term for
    #  the variance of the x term on which the prediction is
    #  calculated.  This is the $stdevx*$stdevx term.  This formula
    #  may be derived by considering the expected value of the variance
    #  of the prediction given by box 14.2.
    y = reginfo.b * x + reginfo.a
    stdevy = reginfo.mse * math.sqrt(
        (1 + 1 / reginfo.N + (stdevx * stdevx + (x - reginfo.avgx) * (x - reginfo.avgx))
        / reginfo.sumxm2
    ))
    return y, stdevy


# Function to analyze parameters
def analyze_parameters(sys_args):
    if len(sys_args) == 1:
        usage()
        exit()

    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-g1', '--gfile1', required=True, type=argparse.FileType('r'), help='Genepix file 1')
    parser.add_argument('-g2', '--gfile2', required=True, type=argparse.FileType('r'), help='Genepix file 2')
    parser.add_argument('-f1', '--fluor1', required=True, default="", type=str, help='Fluor used in Genepix file 1')
    parser.add_argument('-f2', '--fluor2', required=True, default="", type=str, help='Fluor used in Genepix file 2')
    parser.add_argument('-ll', '--linlow', default=Dlinlow, type=int, help='Low limit of linear adjbsi range')
    parser.add_argument('-lh', '--linhigh', default=Dlinhigh, type=int, help='High limit of linear adjbsi range')
    parser.add_argument('-inc', default=Dinc, type=int, help='Dynamic regression range increment')
    parser.add_argument('-m', '--modestr', default=MODE_LABEL[Dmode], type=str, help='Regression calculation mode')
    parser.add_argument('-q', default=Dq, type=float, help='Dynamic regression range quality')
    parser.add_argument('-o', '--outfile', required=True, type=argparse.FileType('w'),
                        help='Output adjusted high intensity scan file')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')

    args_ = parser.parse_args(sys_args[1:])
    assert args_.linlow <= args_.linhigh, \
        f"Error: linear range low parameter {args_.linlow} value exceeds linear range high value {args_.linhigh}"
    args_.mode = MODE_LABEL.index(args_.modestr)
    assert args_.inc > 0, f"Error: -lh/-ll increment value (-inc) {args_.inc} not > 0."
    assert args_.q > 0, f"Error: regression quality value (-q) {args_.q} not > 0."

    print(f"""Parameters for current masliner execution:
  -g1: {args_.gfile1.name} (genepix file 1)
  -g2: {args_.gfile2.name} (genepix file 2)
  -f1: {args_.fluor1} (fluor used in genepix file 1)
  -f2: {args_.fluor2} (fluor used in genepix file 2)
  -m:  {MODE_LABEL[args_.mode]} (regression calculation mode)
  -ll: {args_.linlow} (low limit of linear adjbsi range)
  -lh: {args_.linhigh} (high limit of linear adjbsi range)
  -inc: {args_.inc} (dynamic regression range increment)
  -q:  {args_.q} (dynamic regression range quality)
  -o: {args_.outfile.name} (output adjusted high intensity scan file)
""")

    return args_


def usage():
    modelabels = ", ".join(MODE_LABEL)
    modedescs = "\n".join(MODE_DESC)

    print(f'''masliner: MicroArray Spot LINEar Regression.  Takes two microarray genepix
  output files and uses linear regression to compute adjusted background-
  subtracted intensity (BSI) values for spots above the scanner's linear
  range on the higher intensity scan on the basis of the BSI values for 
  these spots on the lower intensity scan array, thereby correcting for 
  saturation and constructing a common linear scale for a broader intensity
  range.  Output is an adjusted genepix-format file for the higher 
  intensity scan  that can be used as input to masliner to adjust scans 
  that are even more intense.  Adjusted BSI values and related information
  are placed in four extra columns of each row of the output file:

    ADJCOUNT: the number of masliner adjustments performed for this spot
    ADJBSI:   adjusted BSI value computed by masliner
    REGERR:   estimated error for computed ADJBSI value (or standard 
              deviation of background for unadjusted values)
    SATURATION-FLAG: 1 if ADJBSI was computed based on a lower intensity
              scan value that was itself saturated; 0 otherwise.

Version: {VERSION}

Syntax:
  masliner -g1 genepix-file-1 -g2 genepix-file-2 -o output-file
     [-m calibration-mode] 
     [-f1 fluor-genepix-file-1] [-f2 fluor-genepix-file-2]
     [-ll linear-range-low-limit] [-lh linear-range-high-limit]
     [-inc rangelimit-increment-factor] [-q regression-quality]

Notes:
  output-file: adjusted genepix-value for higher BSI array with values 
    above -lh reomputed recomputed based on linear regression. 
  fluor-genepix-file-1 (2): Indicates whether BSI values in 
    genepix-file-1 (2) are from the "F1" or "F2" columns.  Defaults are 
    {Dfluor1} and {Dfluor2}, respectively (used by ScanArray [GSI Lumonics]).
  linear-range-low-limit, linear-range-high-limit: BSI values 
    considered to be within linear range of scanner. Defaults are {Dlinlow}
    and {Dlinhigh}, respectively.
  calibration-mode: {modelabels}
{modedescs}
    Default = {MODE_LABEL[Dmode]}.
  rangelimit-increment-factor: increments of -lh and decrements of -ll
    that are considered when extending linear range for regression for
    "{MODE_LABEL[MODE_DYNOPT]}" processing.  Largest range that meets
    the quality thresholds in -q is selected. Default = {Dinc}.
  regression_quality: threshold for RMS error of regression over 
    average BSI of higher BSI values in -ll to -ll range (an
    estimate of coefficient of variation of unexplained regression error).  
    Default = {Dq}
  Regression calculations are based on Sokal, R.R., and Rohlf, F.S., 1995
    Biometry, 3rd ed., (W.H. Freeman and Company, New York), chapter 14.
''')


def main():
    # Read Genepix files
    gp1 = read_genepix_file(args.gfile1, args.fluor1)
    gp2 = read_genepix_file(args.gfile2, args.fluor2)

    # Validate files
    for blkcolrow1, name1 in gp1.id.items():
        if gp2.id[blkcolrow1] != name1:
            print(f"Error: Genepix files disagree on gene in spot ({blkcolrow1}): {name1} vs. {gp2.id[blkcolrow1]}.")
            print("Program terminating.")
            exit()
    for blkcolrow2, name2 in gp2.id.items():
        if gp1.id[blkcolrow2] != name2:
            print(f"Error: Genepix files disagree on gene in spot ({blkcolrow2}): {gp1.id[blkcolrow2]} vs. {name2}.")
            print("Program terminating.")
            exit()

    # Get total scan intensity and interchange slides if necessary
    sum1 = sum2 = 0
    for blkcolrow, scanbsi in gp1.scanbsi.items():
        sum1 += scanbsi
        sum2 += gp2.scanbsi[blkcolrow]
    if sum1 > sum2:
        args.gfile1, args.gfile2 = args.gfile2, args.gfile1
        args.fluor1, args.fluor2 = args.fluor2, args.fluor1
        gp1.scanbsi, gp2.scanbsi = gp2.scanbsi, gp1.scanbsi
        gp1.bkg, gp2.bkg = gp2.bkg, gp1.bkg
        gp1.adjbsi, gp2.adjbsi = gp2.adjbsi, gp1.adjbsi
        gp1.adjcount, gp2.adjcount = gp2.adjcount, gp1.adjcount
        gp1.regerr, gp2.regerr = gp2.regerr, gp1.regerr
        gp1.satflag, gp2.satflag = gp2.satflag, gp1.satflag
        gp1.status, gp2.status = gp2.status, gp1.status

    # Check file status
    if gp2.status == STATUS_ADJFILE:
        print("\nWARNING: Higher adjbsi scan file", args.gfile2, "has been previously adjusted.")
        print("         Results of previous adjustment will be overwritten.")

    if DEBUG == 1:
        print("\nLow intensity scan file background standard deviation =", gp1.bkgstdev)
        print("\nHigh intensity scan file background standard deviation =", gp2.bkgstdev)

    print("\nCommencing regression analysis and adjustment.")
    perform_adjustment(gp1, gp2)


if __name__ == '__main__':
    # Analyze parameters from command-line arguments
    args = analyze_parameters(sys.argv)
    main()
