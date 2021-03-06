#! /usr/bin/env python3

from pathlib import Path
from functools import reduce
from operator import add
from math import sqrt
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
    central   = [Path(path) for path in args.central_scale
                 if Path(path).exists()]
    variation = [Path(path) for path in args.variation_scales
                 if Path(path).exists()]

    try:
        cen_hist = [read_mcfmhisto(hist_in) for hist_in in central]
        var_hist = [read_mcfmhisto(hist_in) for hist_in in variation]

        # match up like observables
        all_hist = [ [cen] + [var for var in var_hist if var.obs == cen.obs]
                     for cen in cen_hist ]
    except:
        # pass
        print("Something nasty happened processing histograms for the scale variations")
        exit()

    # central scale
    cen_scale = [ cen.xsecs for cen in cen_hist]

    # extract min and max values
    min_scale = [ [min(xsec_tup) for xsec_tup in
                   zip(*[hist.xsecs for hist in obs])]
                  for obs in all_hist ]

    max_scale = [ [max(xsec_tup) for xsec_tup in
                   zip(*[hist.xsecs for hist in obs])]
                  for obs in all_hist ]

    # if we want symmetric error bands need to compute the maximum that the
    # envelope differs from the central value.
    if args.symmetric:
        # calculate the absolute error of the symmetric bands
        sym_scale = [ [max(abs(mxb - cb), abs(mnb - cb)) for cb, mnb, mxb
                       in zip(cen, mini, maxi)]
                       for cen, mini, maxi in zip(cen_scale,
                                                  min_scale, max_scale) ]

        # reconstruct min_scale with the symmetric error band
        min_scale = [ [cb - mnb for cb, mnb in zip(cen, sym)]
                      for cen, sym in zip(cen_scale, sym_scale)]

        max_scale = [ [cb + mxb for cb, mxb in zip(cen, sym)]
                      for cen, sym in zip(cen_scale, sym_scale)]

        result = [scalesymhisto(cen.obs, cen.nbins, cen.xmin, cen.xmax,
                          cen.bins, cen.xsecs, err)
            for cen, err in zip(cen_hist, sym_scale)]

    else:
        result = [scalehisto(cen.obs, cen.nbins, cen.xmin, cen.xmax,
                             cen.bins, cen.xsecs, mini, maxi)
                  for cen, mini, maxi in zip(cen_hist, min_scale, max_scale)]

    return result


def match_nlo(args):
    try:
        # inputs
        lo_in   = [Path(path) for path in args.lo if Path(path).exists()]
        virt_in = [Path(path) for path in args.virt if Path(path).exists()]
        real_in = [Path(path) for path in args.real if Path(path).exists()]

        # number of inputs for consistency checks
        num_lo_in   = len(lo_in)
        num_virt_in = len(virt_in)
        num_real_in = len(real_in)
    except TypeError:
        print("there need to be inputs for lo, virt and real")
        exit()
    except:
        print("something weird and unexpected happened")
        exit()

    if num_lo_in == num_virt_in == num_real_in:
        # histograms from inputs
        lo_histos   = [read_mcfmhisto(hist_in) for hist_in in lo_in]
        virt_histos = [read_mcfmhisto(hist_in) for hist_in in virt_in]
        real_histos = [read_mcfmhisto(hist_in) for hist_in in real_in]

        # derived histograms
        nlo_histos  = [lo + virt + real for lo, virt, real
                       in zip(lo_histos, virt_histos, real_histos)]
    else:
        print("Unequal numbers of inputs")
        exit()

    return nlo_histos

def match_nlo_nll(args):
    try:
        # inputs
        lo_in      = [Path(path) for path in args.lo if Path(path).exists()]
        virt_in    = [Path(path) for path in args.virt if Path(path).exists()]
        real_in    = [Path(path) for path in args.real if Path(path).exists()]
        nll_in     = [Path(path) for path in args.nll if Path(path).exists()]
        nllexpd_in = [Path(path) for path in args.nllexpd
                      if Path(path).exists()]

        # number of inputs for consistency checks
        num_lo_in      = len(lo_in)
        num_virt_in    = len(virt_in)
        num_real_in    = len(real_in)
        num_nll_in     = len(nll_in)
        num_nllexpd_in = len(nllexpd_in)
    except TypeError:
        print("there need to be inputs for lo, virt and real")
        exit()
    except:
        print("something weird and unexpected happened")
        exit()

    if num_lo_in == num_virt_in == num_real_in == num_nll_in == num_nllexpd_in:
        # histograms from inputs
        lo_histos   = [read_mcfmhisto(hist_in) for hist_in in lo_in]
        virt_histos = [read_mcfmhisto(hist_in) for hist_in in virt_in]
        real_histos = [read_mcfmhisto(hist_in) for hist_in in real_in]
        nll_histos  = [read_mcfmhisto(hist_in) for hist_in in nll_in]
        nll1_histos = [read_mcfmhisto(hist_in) for hist_in in nllexpd_in]

        # derived histograms
        nlo1_histos = [virt + real for virt, real
                       in zip(virt_histos, real_histos)]

        matched_histos = [nll*(1 + (nlo1-nll1)/lo) for lo, nlo1, nll, nll1 in
                          zip(lo_histos, nlo1_histos, nll_histos, nll1_histos)]
    else:
        print("Unequal numbers of inputs")
        exit()

    return matched_histos

def match_nlo_nnll(args):
    # assume the default of scheme a throughout
    try:
        # inputs
        lo_in       = [Path(path) for path in args.lo
                       if Path(path).exists()]
        virt_in     = [Path(path) for path in args.virt
                       if Path(path).exists()]
        real_in     = [Path(path) for path in args.real
                       if Path(path).exists()]
        nnll_in     = [Path(path) for path in args.nnll
                       if Path(path).exists()]
        nnllexpd_in = [Path(path) for path in args.nnllexpd
                       if Path(path).exists()]
        lumi0_in    = [Path(path) for path in args.lumi0
                       if Path(path).exists()]
        lumi1_in    = [Path(path) for path in args.lumi1
                       if Path(path).exists()]
        if args.scheme == "b":
            real_incl_in = [Path(path) for path in args.real_inclusive
                            if Path(path).exists]

        # number of inputs for consistency checks
        num_lo_in       = len(lo_in)
        num_virt_in     = len(virt_in)
        num_real_in     = len(real_in)
        num_nnll_in     = len(nnll_in)
        num_nnllexpd_in = len(nnllexpd_in)
        num_lumi0_in    = len(lumi0_in)
        num_lumi1_in    = len(lumi1_in)
        if args.scheme == "b":
            num_real_incl_in = len(real_incl_in)
    except TypeError:
        print("there need to be inputs for lo, virt and real")
        exit()
    except:
        print("something weird and unexpected happened")
        exit()

    if (num_lo_in == num_virt_in == num_real_in == num_nnll_in == 
        num_nnllexpd_in == num_lumi0_in == num_lumi1_in):
        # histograms from inputs
        lo_histos    = [read_mcfmhisto(hist_in) for hist_in in lo_in]
        virt_histos  = [read_mcfmhisto(hist_in) for hist_in in virt_in]
        real_histos  = [read_mcfmhisto(hist_in) for hist_in in real_in]
        nnll_histos  = [read_mcfmhisto(hist_in) for hist_in in nnll_in]
        nnll1_histos = [read_mcfmhisto(hist_in) for hist_in in nnllexpd_in]
        lumi0_histos = [read_mcfmhisto(hist_in) for hist_in in lumi0_in]
        lumi1_histos = [read_mcfmhisto(hist_in) for hist_in in lumi1_in]

        # derived histograms
        dlumi_histos   = [lumi1/lumi0 for lumi0, lumi1
                          in zip(lumi0_histos, lumi1_histos)]
        Sigma1_histos  = [virt + real for virt, real
                          in zip(virt_histos, real_histos)]

        matched_histos = [nnll*(1 + (Sigma1-nnll1)/(lo*(1+dlumi)))
                          for lo, Sigma1, nnll, nnll1, dlumi
                          in zip(lo_histos, Sigma1_histos, nnll_histos,
                                 nnll1_histos, dlumi_histos)]
    else:
        print("Unequal numbers of inputs")
        exit()

    if args.scheme == "b":
        if num_real_incl_in == num_lo_in:
            real_incl_histos = [read_mcfmhisto(hist_in) for hist_in
                                in real_incl_in]
            sigma1_histos    = [virt + real for virt, real
                                in zip(virt_histos, real_incl_histos)]
            matched_histos   = [ma - sigma1*nnll/(lo*(1+dlumi))
                                for ma, lo, sigma1, nnll, dlumi in
                                zip(matched_histos, lo_histos, sigma1_histos,
                                    nnll_histos, dlumi_histos)]

    return matched_histos

def match_nnll_pure(args):
    # assume the default of scheme a throughout
    try:
        # inputs
        nnll_in     = [Path(path) for path in args.nnll
                       if Path(path).exists()]
        if args.scheme == "b":
            lo_in       = [Path(path) for path in args.lo
                           if Path(path).exists()]
            virt_in     = [Path(path) for path in args.virt
                           if Path(path).exists()]
            real_incl_in = [Path(path) for path in args.real_inclusive
                            if Path(path).exists]

        # number of inputs for consistency checks
        num_nnll_in     = len(nnll_in)
        if args.scheme == "b":
            num_lo_in        = len(lo_in)
            num_virt_in      = len(virt_in)
            num_real_incl_in = len(real_incl_in)
    except TypeError:
        print("there need to be inputs for nnll")
        exit()
    except:
        print("something weird and unexpected happened")
        exit()

    # histograms from inputs
    nnll_histos  = [read_mcfmhisto(hist_in) for hist_in in nnll_in]

    # derived histograms
    matched_histos = [nnll for nnll in nnll_histos]

    if args.scheme == "b":
        if num_nnll_in == num_lo_in == num_virt_in == num_real_incl_in:
            lo_histos    = [read_mcfmhisto(hist_in) for hist_in in lo_in]
            virt_histos  = [read_mcfmhisto(hist_in) for hist_in in virt_in]
            real_incl_histos = [read_mcfmhisto(hist_in) for hist_in
                                in real_incl_in]
            sigma1_histos    = [virt + real for virt, real
                                in zip(virt_histos, real_incl_histos)]
            matched_histos   = [ma - sigma1*nnll/lo
                                for ma, lo, sigma1, nnll in
                                zip(matched_histos, lo_histos, sigma1_histos,
                                    nnll_histos)]

    return matched_histos


def efficiency(args):
    vetoed    = [Path(path) for path in args.vetoed
                 if Path(path).exists()]
    inclusive = [Path(path) for path in args.inclusive
                 if Path(path).exists()]

    try:
        veto_hist = [read_mcfmhisto(hist_in) for hist_in in vetoed]
        incl_hist = [read_mcfmhisto(hist_in) for hist_in in inclusive]

        # have a check to make sure observables align and there are
        # equal numbers of inputs?
    except:
        pass

    # efficiency
    eff = [veto/incl for veto, incl in zip(veto_hist, incl_hist)]

    return eff

def jet_veto_efficiency(args):
    efficiency = [Path(path) for path in args.efficiency
                  if Path(path).exists()]
    inclusive  = [Path(path) for path in args.inclusive_xsec
                  if Path(path).exists()]

    try:
        eff_hist  = [read_scalesymhisto(hist_in) for hist_in in efficiency]
        incl_hist = [read_scalesymhisto(hist_in) for hist_in in inclusive]

        # have a check to make sure observables align and there are
        # equal numbers of inputs?
    except:
        pass

    # JVE
    jve = [incl*eff for incl, eff in zip(incl_hist, eff_hist)]

    return jve

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
parser_data.add_argument(
    'input',
    metavar='input(s)',
    type=str,
    action='store',
    nargs='+',
    help="mcfm .dat file(s) to be processed")
parser_data.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_data.add_argument(
    '--average',
    action='store_true',
    help="""Average over the input files, i.e. report histograms as sum over
    inputs divide the number of inputs. One should use this if runs were
    performed in parallel, but with the same setup""")
parser_data.set_defaults(func=histo_mk)

parser_scale = subparsers.add_parser(
    "scale_var",
    description="",
    help="")
parser_scale.add_argument(
    '-cs', '--central_scale',
    type=str,
    action='store',
    nargs='+',
    help="""mcfm histogram file(s) for a given set of
    observables that consistitute the central scale""")
parser_scale.add_argument(
    '-vs', '--variation_scales',
    type=str,
    action='store',
    nargs='*',
    help="""mcfm histogram file(s) that constitute the scale
    variations to be performed""")
parser_scale.add_argument(
    '-s', '--symmetric',
    dest='symmetric',
    action='store_true',
    help="""Symmetrise the uncertainty bands by taking the maximum of the upper
    and lower uncertainties""")
parser_scale.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_scale.set_defaults(func=scale_var)

parser_match = subparsers.add_parser(
    "match",
    description="",
    help="")
match_subparsers = parser_match.add_subparsers(
    title="subcommands",
    description="des",
    help="help me")
parser_match_nlo = match_subparsers.add_parser(
    "nlo",
    description="",
    help="")
parser_match_nlo.add_argument(
    '-l', '--lo',
    dest='lo',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo.add_argument(
    '-v', '--virt',
    required=True,
    dest='virt',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo.add_argument(
    '-r', '--real',
    required=True,
    dest='real',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_match_nlo.set_defaults(func=match_nlo)

parser_match_nlo_nll = match_subparsers.add_parser(
    "nlo+nll",
    description="",
    help="")
parser_match_nlo_nll.add_argument(
    '-l', '--lo',
    dest='lo',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nll.add_argument(
    '-v', '--virt',
    dest='virt',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nll.add_argument(
    '-r', '--real',
    dest='real',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nll.add_argument(
    '--nll',
    dest='nll',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nll.add_argument(
    '--nllexpd',
    dest='nllexpd',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nll.add_argument(
    '-ri', '--real_inclusive',
    dest='real_inclusive',
    type=str,
    action='store',
    nargs='+',
    help="""The inclusive real cross section with no constraints on QCD
    radiation applied, only necessary for scheme 'b'""")
parser_match_nlo_nll.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_match_nlo_nll.set_defaults(func=match_nlo_nll)

parser_match_nlo_nnll = match_subparsers.add_parser(
    "nlo+nnll",
    description="",
    help="")
parser_match_nlo_nnll.add_argument(
    '-s', '--scheme',
    dest='scheme',
    type=str,
    action='store',
    default='a',
    choices=['a', 'b'],
    help="""""")
parser_match_nlo_nnll.add_argument(
    '-l', '--lo',
    dest='lo',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '-v', '--virt',
    dest='virt',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '-r', '--real',
    dest='real',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '--nnll',
    dest='nnll',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '--nnllexpd',
    dest='nnllexpd',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '--lumi0',
    dest='lumi0',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '--lumi1',
    dest='lumi1',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nlo_nnll.add_argument(
    '-ri', '--real_inclusive',
    dest='real_inclusive',
    type=str,
    action='store',
    nargs='+',
    help="""The inclusive real cross section with no constraints on QCD
    radiation applied, only necessary for scheme 'b'""")
parser_match_nlo_nnll.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_match_nlo_nnll.set_defaults(func=match_nlo_nnll)

parser_match_nnll_pure = match_subparsers.add_parser(
    "nnll_pure",
    description="",
    help="")
parser_match_nnll_pure.add_argument(
    '-s', '--scheme',
    dest='scheme',
    type=str,
    action='store',
    default='a',
    choices=['a', 'b'],
    help="""""")
parser_match_nnll_pure.add_argument(
    '--nnll',
    dest='nnll',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nnll_pure.add_argument(
    '-l', '--lo',
    dest='lo',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nnll_pure.add_argument(
    '-v', '--virt',
    dest='virt',
    type=str,
    action='store',
    nargs='+',
    help='')
parser_match_nnll_pure.add_argument(
    '-ri', '--real_inclusive',
    dest='real_inclusive',
    type=str,
    action='store',
    nargs='+',
    help="""The inclusive real cross section with no constraints on QCD
    radiation applied, only necessary for scheme 'b'""")
parser_match_nnll_pure.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_match_nnll_pure.set_defaults(func=match_nnll_pure)

parser_efficiency = subparsers.add_parser(
    "efficiency",
    description="",
    help="")
parser_efficiency.add_argument(
    '-c', '--with_cuts',
    dest='vetoed',
    type=str,
    action='store',
    nargs='+',
    help="""""")
parser_efficiency.add_argument(
    '-i', '--inclusive',
    dest='inclusive',
    type=str,
    action='store',
    nargs='+',
    help="""""")
parser_efficiency.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_efficiency.set_defaults(func=efficiency)

parser_jve = subparsers.add_parser(
    "JVE",
    description="",
    help="")
parser_jve.add_argument(
    '-e', '--efficiency',
    dest='efficiency',
    type=str,
    action='store',
    nargs='+',
    help="""""")
parser_jve.add_argument(
    '-x', '--inclusive_xsec',
    dest='inclusive_xsec',
    type=str,
    action='store',
    nargs='+',
    help="""""")
parser_jve.add_argument(
    '-o', '--output',
    dest='output',
    type=str,
    action='store',
    default='histogram',
    help="""The prefix that is attached at the start of the output filename""")
parser_jve.set_defaults(func=jet_veto_efficiency)

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
        return "hist blah"
        # pass

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
            pretty_histo += "{:08.3f}\t{:014.10E}\n".format(bin, xsec)
        pretty_histo += "## END HISTOGRAM"

        return pretty_histo

    def __format__(self):
        pass

    def __len__(self):
        assert len(self.bins) == len(self.xsecs)
        return len(self.bins)

    def __neg__(self):
        xsecs = [-x for x in self.xsecs]
        return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                         self.bins, xsecs)

    def __pos__(self):
        xsecs = [+x for x in self.xsecs]
        return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                         self.bins, xsecs)

    def __pow__(self, other):
        if isinstance(other, numbers.Real):
            # catch negative powers of zero that would give divide by zero
            # exceptions we tactically assume that this errors always comes
            # from the matching where we have expressions like 0*(1 + 0/0)
            # which is well defined
            xsecs = [0.0 if x == 0 and other <= 0
                     else x**other for x in self.xsecs]
            return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                             self.bins, xsecs)
        else:
            return NotImplemented

    def __add__(self, other):
        if isinstance(other, numbers.Real):
            xsecs = [other + x for x in self.xsecs]
            return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                             self.bins, xsecs)
        elif isinstance(other, mcfmhisto):
            try:
                if self.obs == None:
                    pass
                elif other.obs == None:
                    pass
                elif self.obs != other.obs:
                    raise MixedObs
                xsecs = [x + y for x, y in zip(self.xsecs, other.xsecs)]
                return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                                 self.bins, xsecs)
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

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        if isinstance(other, numbers.Real):
            xsecs = [other * x for x in self.xsecs]
            return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                             self.bins, xsecs)
        elif isinstance(other, mcfmhisto):
            try:
                if self.obs == None:
                    pass
                elif other.obs == None:
                    pass
                elif self.obs != other.obs:
                    raise MixedObs
                xsecs = [x * y for x, y in zip(self.xsecs, other.xsecs)]
                return mcfmhisto(self.obs, self.nbins, self.xmin, self.xmax,
                                 self.bins, xsecs)
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

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self * other**-1

    def __rtruediv__(self, other):
        return self**-1 * other

class scalehisto(object):
    def __init__(self, obs, nbins, xmin, xmax, bins, central_scale,
                 min_scale, max_scale):
        self.obs = obs
        self.nbins = nbins
        self.xmin = xmin
        self.xmax = xmax
        self.bins = bins
        self.central_scale = central_scale
        self.min_scale = min_scale
        self.max_scale = max_scale

    def __str__(self):
        pretty_histo = """## BEGIN HEADER
# observable name:
# {obs}
# nbins\txmin\txmax
# {nbins}        {xmin}        {xmax}
## END HEADER
## BEGIN HISTOGRAM
""".format(obs=self.obs, nbins=self.nbins, xmin=self.xmin, xmax=self.xmax)
        for bin, cen, mini, maxi in zip(self.bins,
                                        self.central_scale,
                                        self.min_scale,
                                        self.max_scale):
            pretty_histo += "{:08.3f}\t{:014.10E}\t{:014.10E}\t{:014.10E}\n".format(bin, cen, mini, maxi)
        pretty_histo += "## END HISTOGRAM"

        return pretty_histo

class scalesymhisto(object):
    def __init__(self, obs, nbins, xmin, xmax, bins, central_scale, sym_scale):
        self.obs = obs
        self.nbins = nbins
        self.xmin = xmin
        self.xmax = xmax
        self.bins = bins
        self.central_scale = central_scale
        self.sym_scale = sym_scale

    def __str__(self):
        pretty_histo = """## BEGIN HEADER
# observable name:
# {obs}
# nbins\txmin\txmax
# {nbins}        {xmin}        {xmax}
## END HEADER
## BEGIN HISTOGRAM
""".format(obs=self.obs, nbins=self.nbins, xmin=self.xmin, xmax=self.xmax)
        for bin, cen, err in zip(self.bins,
                                 self.central_scale,
                                 self.sym_scale):
            pretty_histo += "{:08.3f}\t{:014.10E}\t{:014.10E}\n".format(bin, cen, err)
        pretty_histo += "## END HISTOGRAM"

        return pretty_histo

    def __mul__(self, other):
        if isinstance(other, numbers.Real):
            central_scales = [other * x for x in self.central_scale]
            sym_scales     = [other * x for x in self.sym_scale]
            return scalesymhisto(self.obs, self.nbins, self.xmin, self.xmax,
                                 self.bins, central_scales, sym_scales)
        elif isinstance(other, scalesymhisto):
            try:
                if self.obs == None:
                    pass
                elif other.obs == None:
                    pass
                elif self.obs != other.obs:
                    raise MixedObs
                central_scales = [x * y for x, y in zip(self.central_scale, other.central_scale)]
                sym_scales = [err_quad(sc, sdc, oc, odc) for sc, sdc, oc, odc
                              in zip(self.central_scale,
                                     self.sym_scale,
                                     other.central_scale,
                                     other.sym_scale)]

                # sym_scales = [abs(sc*oc)*sqrt( (sdc/sc)**2 + (odc/oc)**2 )
                #               for sc, sdc, oc, odc
                #               in zip(self.central_scale,
                #                      self.sym_scale,
                #                      other.central_scale,
                #                      other.sym_scale)]

                return scalesymhisto(self.obs, self.nbins, self.xmin, self.xmax,
                                     self.bins, central_scales, sym_scales)
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

    def __rmul__(self, other):
        return self.__mul__(other)



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

def err_quad(x, dx, y, dy):

    if x == 0:
        dxonxsq = 0.0
    else:
        dxonxsq = (dx/x)**2

    if y == 0:
        dyonysq = 0.0
    else:
        dyonysq = (dy/y)**2

    return abs(x*y)*sqrt(dxonxsq + dyonysq)

def mcfm2hist(mcfm_hist):
    histlist = []
    for enum, hist_elem in enumerate(mcfm_hist):
        histlist.append(parse_histo(hist_elem))
    return histlist

def read_hists(mcfm_in):
    with open(str(mcfm_in), "r") as f:
        header, *histograms = per_section(nonempty_lines(f))
    return histograms

def read_mcfmhisto(mcfmhisto_in):
    with open(str(mcfmhisto_in), "r") as f:
        header, histogram = per_section2(nonempty_lines(f))

    ## Header
    # removes "# " from strings, maybe I can do this without
    # magic numbers
    obs_name = header[1][2:]
    nbins, xmin, xmax = header[3][2:].split()
    nbins = int(nbins)
    xmin = float(xmin)
    xmax = float(xmax)

    ## Data
    bins = []
    xsecs = []
    for line in histogram:
        bin, xsec = line.split()
        bins.append(float(bin))
        xsecs.append(float(xsec))

    ## Footer

    return mcfmhisto(obs_name, nbins, xmin, xmax, bins, xsecs)

def read_scalesymhisto(mcfmhisto_in):
    with open(str(mcfmhisto_in), "r") as f:
        header, histogram = per_section2(nonempty_lines(f))

    ## Header
    # removes "# " from strings, maybe I can do this without
    # magic numbers
    obs_name = header[1][2:]
    nbins, xmin, xmax = header[3][2:].split()
    nbins = int(nbins)
    xmin = float(xmin)
    xmax = float(xmax)

    ## Data
    bins = []
    central_scales = []
    sym_scales = []
    for line in histogram:
        bin, central_scale, sym_scale = line.split()
        bins.append(float(bin))
        central_scales.append(float(central_scale))
        sym_scales.append(float(sym_scale))

    ## Footer

    return scalesymhisto(obs_name, nbins, xmin, xmax, bins, central_scales, sym_scales)

def nonempty_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def per_section(it,
                is_delimiter=lambda line: line.startswith("HIST")):
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

def per_section2(it,
                 is_head=lambda line: line.startswith("##"),
                 is_tail=lambda line: line.startswith("##")):
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

    args = parser.parse_args()

    histograms = args.func(args)

    for histnum, hist in enumerate(histograms):
        try:
            with open(args.output + "-" + hist.obs, "w") as f:
                f.write(str(hist))
        except:
            # print("Histogram {} failed".format(histnum))
            pass
