#! /usr/bin/env python3

from pathlib import Path
import argparse
import re

# parser
parser = argparse.ArgumentParser(description="""""")
subparsers = parser.add_subparsers()

parser_data  = subparsers.add_parser("histo_mk", help="")
parser_data.add_argument('input',
                    metavar='input(s)',
                    type=str,
                    action='store',
                    nargs='+',
                    help="mcfm .dat file(s) to be processed")
parser_data.add_argument('-o', '--output',
                    dest='output',
                    type=str,
                    action='store',
                    default='histogram',
                    help="""The prefix that is attached at the start of the
                    output filename""")
parser_data.add_argument('--average',
                    action='store_true',
                    help="""Average over the input files, i.e. report
                    histograms as sum over inputs divide the number of inputs.
                    One should use this if runs were performed in parallel, but
                    with the same setup""")

parser_scale = subparsers.add_parser("scale_var", help="")

parser_match = subparsers.add_parser("match", help="")

parser_delta = subparsers.add_parser("delta", help="")


class mcfmhisto(object):
    def __init__(self, obs, nbins, xmin, xmax, bins, xsecs):
        self.obs = obs
        self.nbins = nbins
        self.xmin = xmin
        self.xmax = xmax
        self.bins = bins
        self.xsecs = xsecs

    def __str__(self):
        print(self.obs)
        print(self.nbins)
        print(self.xmin)
        print(self.xmax)
        print(self.bins)
        print(self.xsecs)
        return ""

####  mcfm .dat output  ####

## Header information
## Histogram Block
   # newline
   # HIST = number
   # observable name
   # data
   # newline
   # AVG RMS INTEGRAL
   # newline
   # ENTRIES O`FLOW U`FLOW
   # newline
## Histogram Block
   # ...

def parse_histo(histo_in):
    try:
        obs_name, *data, integ, shots = histo_in

        bins = []
        xsecs = []

        for line in data:
            bin, xsec, xsecerr = line.split()
            bins.append(float(bin))
            xsecs.append(float(xsec))

        nbins = len(bins)
        xmin = bins[0]
        xmax = bins[-1]

        return mcfmhisto(obs_name, nbins, xmin, xmax, bins, xsecs)

    except:
        return mcfmhisto(None, None, None, None, None, None)

# helper functions
def nonempty_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def per_section(it, is_delimiter=lambda line: line.startswith("HIST")):
    ret = []
    for line in it:
        # remove leading and trailing spaces
        line = line.rstrip().lstrip()
        if is_delimiter(line):
            yield ret
            ret = []
        else:
            ret.append(line)
    yield ret

if __name__ == "__main__":

    # args = parser.parse_args()

    # parse input into header and separate histograms
    with open("mcfmtest.dat", "r") as f:
        header, *histograms = list(per_section(nonempty_lines(f)))

    # parse header

    # parse histograms
    print(parse_histo(histograms[0]))

    # with open("mcfmtest.dat", "r") as f:
    #     for line in nonempty_lines(f):
    #         # print(type(line), line)
    #         # s = re.findall(r'\s\(.*?\)', line)    # this is what it should be if parens appeared properly
    #         s = re.findall(r'\s\(.*?$', line)
    #         if s:
    #             print(s.pop())

    #         if "HIST" in line:
    #             hist = line
    #             obs_name = f.readline()
    #             print(hist)
    #             print(obs_name)

    #             while True:
    #                 try:
    #                     bin, xsec, xsecerr = f.readline().split()
    #                     print(bin, xsec, xsecerr)
    #                 except:
    #                     break

    #             avg, rms, integ = re.findall(r"[+-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)", f.readline())
    #             print(avg, rms, integ)
    #             f.readline()
    #             entries, uflow, oflow = re.findall(r"\d+", f.readline())
    #             print(entries, uflow, oflow)
