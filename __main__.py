import os
import ROOT
import pstats
import cProfile
from glob import glob
from argparse import ArgumentParser

from eventdisplay import EventDisplay
from constants import *

def parseargs():

    parser = ArgumentParser()
    parser.add_argument('-input-dir', '-i')
    parser.add_argument('--diffuser', choices=['barefibre', 'collimator', 'diffuser', 'monitor'], default='barefibre')
    parser.add_argument('--injector', choices=['B1', 'B2', 'B3', 'B4', 'B5'], default='B1')
    parser.add_argument('-b', '--batch', action='store_true', default=False)
    parser.add_argument('-p', '--profile', action='store_true', default=False)
    args = parser.parse_args()

    return args

def run(args):

    if args.batch:
        ROOT.gROOT.SetBatch(True)

    chain, run = load_data(args)

    EventDisplay(chain, run, wvar='occ', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, walltime_cut=False, correct=True)
    EventDisplay(chain, run, wvar='charge', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, walltime_cut=False, correct=True)
    #EventDisplay(chain, run, wvar='occ', tof_cut_override=RunInfo.Runs[run].time_sel, fit=True, walltime_cut=False, logz=True)
    #EventDisplay(chain, run, wvar='charge', tof_cut_override=RunInfo.Runs[run].time_sel, fit=True, walltime_cut=False)
    #EventDisplay(chain, run, wvar='charge', tof_cut_override=RunInfo.Runs[run].time_sel, fit=True, walltime_cut=False, logz=True)

    return

def main():

    args = parseargs()

    if args.profile:
        profile_func(run, args=[args])
    else:
        run(args)

    return

def load_data(args):

    print '\nLoading data...'

    flist, sel_run = expand_file_list(args.input_dir, args.diffuser, args.injector)
    chain = load_chain(flist)

    return chain, sel_run

def expand_file_list(input_dir, diffuser, injector, sel_run=None):

    for run in RunInfo.Runs:

        if injector == Injector.tostr(RunInfo.Runs[run].injector) and diffuser == Source.tostr(RunInfo.Runs[run].source):
            sel_run = run
            break

    if sel_run is not None:
        runs = '*%s.root' % str(sel_run)
    else:
        runs = '*.root'

    flist = {sel_run: f for f in glob(os.path.join(input_dir, diffuser, runs))}

    print '\tLoaded files:'

    for f in flist:
        print '\t\t', flist[f]

    return flist, sel_run

def load_chain(flist):
    
    chain = ROOT.TChain("tqtree")
    
    for run in flist:
        chain.AddFile(flist[run])
        
    return chain

def run_profiler(func, *args, **kwargs):
    profile = cProfile.Profile()
    profile.enable()
    func(*args, **kwargs)
    profile.disable()
    return profile

def profile_func(func, args=[], kwargs={}, num=20):
    prof = run_profiler(func, *args, **kwargs)
    ps = pstats.Stats(prof)
    for key in ["time", "cumulative"]:
        print "--- top {} sorted by {}".format(num, key)
        ps.sort_stats(key).print_stats(num)
    return

if __name__ == "__main__":
    main()