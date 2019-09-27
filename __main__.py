import os
import ROOT
import pstats
import cProfile
from glob import glob
from argparse import ArgumentParser
from collections import OrderedDict

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
    parser.add_argument('--flist', '-i', nargs="+")
    parser.add_argument('--source', choices=Source._enums.keys(), default=Source._names[Source.COLLIMATOR])
    parser.add_argument('--injector', choices=Injector._enums.keys(), default=Injector._names[Injector.B1])
    parser.add_argument('--run-period', choices=RunPeriod._enums.keys(), default=RunPeriod._names[RunPeriod.LIVERPOOL_LASER])
    parser.add_argument('--run', default=None)
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

    sel_chains = load_data(args)

    for run, chain in sel_chains.items():
        EventDisplay(chain, run, wvar='occ', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, event_id_cut=False, correct=False)
        EventDisplay(chain, run, wvar='charge', tof_cut_override=RunInfo.Runs[run].time_sel, fit=False, event_id_cut=False, correct=True)

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

    if args.run is None:
        run = find_requested_runs(run_period=args.run_period, injector=args.injector, source=args.source)
    else:
        run = args.run

    sel_files = get_runs(args.flist, *run)
    sel_chains = load_chain(sel_files)

    return sel_chains

def get_runs(flist, *runs):
    """Select runs from flist by run number.

    Args:
        flist (list): The list of all files passed at command line.
        *runs: The run IDs.

    Returns:
        OrderedDict: A dict of run numbers as key and file path as value.
    """

    run_file_map = OrderedDict()

    for run in runs:
        for file in flist:
            if str(run) in str(file):
                run_file_map[run] = file

    return run_file_map

def find_requested_runs(run_period, injector, source):
    """Finds runs matching the run period, injector and source queries.

    Args:
        run_period (str): The run period.
        injector (str): The injector.
        source (str): The source.

    Returns:
        list: The list of runs.
    """

    runs = []

    print 'Finding runs during run period \'%s\' with %s %s...' % (run_period, injector, source)

    for run_id, info in RunInfo.Runs.items():
        if RunPeriod.tostr(info.runperiod) == run_period and Injector.tostr(info.injector) == injector and Source.tostr(info.source) == source:
            runs.append(info.runnum)

    print '\tFound run(s) %s.' % ','.join(map(str, runs))

    return runs

def load_chain(sel_files):
    """Loads the ROOT files into a ROOT chain.

    Args:
        sel_files (OrderedDict): A dict of run numbers as key and paths as value.

    Returns:
        OrderedDict: A dict of run numbers as key and ROOT:TChain as value.
    """

    sel_chains = OrderedDict()

    for run, file in sel_files.items():
        chain = ROOT.TChain("tqtree")
        chain.AddFile(file)

        sel_chains[run] = chain
        
    return sel_chains

def run_profiler(func, *args, **kwargs):
    """Runs the cProfiler.

    Args:
        func (callable): The function to profile.
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        cProfile.Profile: The profiler results.
    """

    profile = cProfile.Profile()
    profile.enable()
    func(*args, **kwargs)
    profile.disable()

    return profile

def profile_func(func, args=[], kwargs={}, num=20):
    """Runs the profiler.

    Args:
        func (callable): The function to profile.
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