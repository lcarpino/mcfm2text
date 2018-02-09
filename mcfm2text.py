#! /usr/bin/env python3

from pathlib import Path
from functools import reduce
from operator import add
import argparse
import numbers
import re




############################################################
####                 Main Methods                       ####
############################################################

def histo_mk(args):
    inputs = [Path(path) for path in args.input if Path(path).exists()]
    ninputs = len(inputs)

    # get mcfm output in list
    mcfm_histos = [read_hists(mcfm_in) for mcfm_in in inputs]

    try:
        # get nested list of all histograms
        all_histos = map(mcfm2hist, mcfm_histos)
        # transpose for reduce call
        all_histos_tp = list(map(lambda *sl : list(sl), *all_histos))
        # Add together histograms with like observables
        histos = [reduce(add, hist) for hist in all_histos_tp]
    except TypeError:
        print("No valid inputs found")
        exit()

    return histos

def scale_var(args):
    pass

def match(args):
    pass

def delta(args):
    pass




############################################################
####                       Parser                       ####
############################################################

parser = argparse.ArgumentParser(description="""""")
subparsers = parser.add_subparsers(title="",
                                   description="",
                                   help="")

parser_data  = subparsers.add_parser("histo_mk",
                                     description="",
                                     help="")
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
parser_data.set_defaults(func=histo_mk)

parser_scale = subparsers.add_parser("scale_var",
                                     description="",
                                     help="")
parser_scale.add_argument('input',
                    metavar='input(s)',
                    type=str,
                    action='store',
                    nargs='+',
                    help="mcfm histogram file(s) to be processed")
parser_scale.add_argument('-o', '--output',
                    dest='output',
                    type=str,
                    action='store',
                    default='histogram',
                    help="""The prefix that is attached at the start of the
                    output filename""")
parser_scale.add_argument('-c', '--central',
                    dest='central',
                    type=str,
                    action='store',
                    help="the histograms of the central scale")
parser_scale.set_defaults(func=scale_var)

parser_match = subparsers.add_parser("match",
                                     description="",
                                     help="")
parser_match.set_defaults(func=match)
match_subparsers = parser_match.add_subparsers(title="subcommands",
                                   description="des",
                                   help="help me")

match_nlo = match_subparsers.add_parser("nlo",
                                        description="",
                                        help="")
match_nlo.add_argument('-lo', '--leading_order',
                    dest='lo',
                    type=str,
                    action='store',
                    nargs='+',
                    help='')
match_nlo.add_argument('-v', '--virt',
                    required=True,
                    dest='virt',
                    type=str,
                    action='store',
                    nargs='+',
                    help='')
match_nlo.add_argument('-r', '--real',
                    required=True,
                    dest='real',
                    type=str,
                    action='store',
                    nargs='+',
                    help='')

match_nlo_nll = match_subparsers.add_parser("nlo+nll",
                                            description="",
                                            help="")
match_nlo_nll.add_argument('-lo', '--leading_order',
                    dest='lo',
                    type=str,
                    action='store',
                    help='')
match_nlo_nll.add_argument('-v', '--virt',
                    dest='virt',
                    type=str,
                    action='store',
                    help='')
match_nlo_nll.add_argument('-r', '--real',
                    dest='real',
                    type=str,
                    action='store',
                    help='')
match_nlo_nll.add_argument('--nll',
                    dest='nll',
                    type=str,
                    action='store',
                    help='')
match_nlo_nll.add_argument('--nllexpd',
                    dest='nllexpd',
                    type=str,
                    action='store',
                    help='')

match_nlo_nnll = match_subparsers.add_parser("nlo+nnll",
                                            description="",
                                            help="")
match_nlo_nnll.add_argument('-lo', '--leading_order',
                    dest='lo',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('-v', '--virt',
                    dest='virt',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('-r', '--real',
                    dest='real',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('--nnll',
                    dest='nnll',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('--nnllexpd',
                    dest='nnllexpd',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('--lumi0',
                    dest='lumi0',
                    type=str,
                    action='store',
                    help='')
match_nlo_nnll.add_argument('--lumi1',
                    dest='lumi1',
                    type=str,
                    action='store',
                    help='')


parser_delta = subparsers.add_parser("delta", help="")




############################################################
####                 Histogram Class                    ####
############################################################

####  mcfm .dat output 

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
        pretty_histo = """## BEGIN HEADER
# observable name:
# {obs}
# nbins\txmin\txmax
# {nbins}        {xmin}        {xmax}
## END HEADER
## BEGIN HISTOGRAM
""".format(obs=self.obs, nbins=self.nbins, xmin=self.xmin, xmax=self.xmax)
        for bin, xsec in zip(self.bins, self.xsecs):
            pretty_histo += "{:08.3f}\t{:014.10f}\n".format(bin, xsec)
        pretty_histo += "## END HISTOGRAM"

        return pretty_histo

    def __format__(self):
        pass

    def __len__(self):
        assert len(self.bins) == len(self.xsecs)
        return len(self.bins)

    def __add__(self, other):
        if isinstance(other, numbers.Real):
            xsecs = [other + x for x in self.xsecs]
            return mcfmhisto(self.obs, self.bins, self.xmin, self.xmax, self.bins, xsecs)
        elif isinstance(other, mcfmhisto):
            try:
                if self.obs == None:
                    pass
                elif other.obs == None:
                    pass
                elif self.obs != other.obs:
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




############################################################
####                 Histogram Parsing                  ####
############################################################

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




############################################################
####                 Helper Functions                   ####
############################################################

def mcfm2hist(mcfm_hist):
    histlist = []
    for enum, hist_elem in enumerate(mcfm_hist):
        histlist.append(parse_histo(hist_elem))
    return histlist

def read_hists(mcfm_in):
    with open(str(mcfm_in), "r") as f:
        header, *histograms = per_section(nonempty_lines(f))
    return histograms

def read_mcfmhisto(mcfmhisto):
    with open(str(mcfmhisto), "r") as f:
        header, histogram = per_section2(nonempty_lines(f))
    print(header)
    print(histogram)

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

def per_section2(it, is_head=lambda line: line.startswith("##"), is_tail=lambda line: line.startswith("##")):
    ret = []
    record_mode = False
    for line in it:
        # remove leading and trailing spaces
        line = line.rstrip().lstrip()
        if not record_mode:
            if is_head(line):
                record_mode = True
        elif is_tail(line):
            yield ret
            ret = []
            record_mode = False
        else:
            ret.append(line)


if __name__ == "__main__":

    read_mcfmhisto("histogram-ptll_full")

    # args = parser.parse_args()

    # histograms = args.func(args)

    # for hist in histograms:
    #     with open(args.output + "-" + hist.obs, "w") as f:
    #         f.write(str(hist))
