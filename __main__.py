import os
import ROOT
import pstats
import cProfile
from glob import glob
from argparse import ArgumentParser

from eventdisplay import EventDisplay
from constants import *

def parseargs():
    """Does command line parsing.

    Args:
        None

    Returns:
        argparse.Namespace: The parsed commandline arguments.

    """
    parser = ArgumentParser()
    parser.add_argument('-input-dir', '-i')
    parser.add_argument('--diffuser', choices=['barefibre', 'collimator', 'diffuser', 'monitor'], default='barefibre')
    parser.add_argument('--injector', choices=['B1', 'B2', 'B3', 'B4', 'B5'], default='B1')
    parser.add_argument('-b', '--batch', action='store_true', default=False)
    parser.add_argument('-p', '--profile', action='store_true', default=False)
    args = parser.parse_args()

    return args

def run(args):
    """Runs the event display generation.

    Args:
        args (argparse.Namespace): The parsed commandline arguments.

    Returns:
        None

    """
    if args.batch:
        ROOT.gROOT.SetBatch(True)

    chain, run = load_data(args)

    EventDisplay(chain, run, wvar='occ', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, walltime_cut=False, correct=False)
    EventDisplay(chain, run, wvar='charge', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, walltime_cut=False, correct=True)

def main():
    """Parse args and run event display generation.

    Args:
        None

    Returns:
        None

    """

    args = parseargs()

    if args.profile:
        profile_func(run, args=[args])
    else:
        run(args)

def load_data(args):
    """Gets the TChain for the run specified in the commandline arguments. 

    Args:
        args (argparse.Namespace): The parsed commandline arguments.

    Returns:
        ROOT:TChain: The TChain requested. 
        int: The ID of the run selected.
    """
    print '\nLoading data...'

    flist, sel_run = expand_file_list(args.input_dir, args.diffuser, args.injector)
    chain = load_chain(flist)

    return chain, sel_run

def expand_file_list(input_dir, diffuser, injector, sel_run=None):
    """Gets the correct files or all files in the specified directory.

    Args:
        input_dir (str): The root directory of the ROOT files.
        diffuser (str): The name of the diffuser type requested.
        injector (str): The name of the injector position requested.
        sel_run (int): The ID of the run requested.

    Returns:
        list: The file list.
        int: The ID of the run selected.
    """

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
    """Loads the ROOT files into a ROOT chain.

    Args:
        flist (list): The list of ROOT files.

    Returns:
        ROOT:TChain: The TChain of data requested.

    """

    chain = ROOT.TChain("tqtree")
    
    for run in flist:
        chain.AddFile(flist[run])
        
    return chain

def run_profiler(func, *args, **kwargs):
    """Loads the ROOT files into a ROOT chain.

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        cProfile.Profile

    """

    profile = cProfile.Profile()
    profile.enable()
    func(*args, **kwargs)
    profile.disable()

    return profile

def profile_func(func, args=[], kwargs={}, num=20):
    """Loads the ROOT files into a ROOT chain.

    Args:
        args (list): Variable length argument list.
        kwargs (dict): Arbitrary keyword arguments.
        num (int): Number of process stats to show.

    Returns:
        None

    """
    prof = run_profiler(func, *args, **kwargs)
    ps = pstats.Stats(prof)
    for key in ["time", "cumulative"]:
        print "--- top {} sorted by {}".format(num, key)
        ps.sort_stats(key).print_stats(num)

if __name__ == "__main__":
    main()