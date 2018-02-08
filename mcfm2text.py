#! /usr/bin/env python3

from pathlib import Path
import argparse
import numbers
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

class MixedObs(Exception):
    pass

class mcfmhisto(object):
    def __init__(self, obs, nbins, xmin, xmax, bins, xsecs):
        self.obs = obs
        self.nbins = nbins
        self.xmin = xmin
        self.xmax = xmax
        self.bins = bins
        self.xsecs = xsecs

    def __repr__(self):
        pass

    def __str__(self):
        print(self.obs)
        print(self.nbins)
        print(self.xmin)
        print(self.xmax)
        print(self.bins)
        print(self.xsecs)
        return ""

    def __len__(self):
        assert len(self.bins) == len(self.xsecs)
        return len(self.bins)

    def __add__(self, other):
        if isinstance(other, numbers.Real):
            xsecs = [other + x for x in self.xsecs]
            return mcfmhisto(self.obs, self.bins, self.xmin, self.xmax, self.bins, xsecs)
        elif isinstance(other, mcfmhisto):
            try:
                if self.obs != other.obs or None:
                    raise MixedObs
                xsecs = [x + y for x, y in zip(self.xsecs, other.xsecs)]
                return mcfmhisto(self.obs, self.bins, self.xmin, self.xmax, self.bins, xsecs)
            except TypeError:
                if self.xsecs is None and other.xsecs is None:
                    return self
                elif self.xsecs is None:
                    return other
                elif other.xsecs is None:
                    return self
                else:
                    return Exception
            except MixedObs:
                return "You can't mix observables, this doesn't make sense"
        else:
            return NotImplemented

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

def parse_header(header_in):
    pass

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

        # currently unused
        # avg, rms, integral = integ.split()
        # entries, uflow, oflow = shots.split()

        return mcfmhisto(obs_name, nbins, xmin, xmax, bins, xsecs)

    except ValueError:
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
    print(parse_histo(histograms[0]) + parse_histo(histograms[0]))
