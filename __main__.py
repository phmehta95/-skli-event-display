import os
import ROOT
from glob import glob
from argparse import ArgumentParser

from eventdisplay import EventDisplay
from constants import *

data_dir = '/Users/billyvinning/Work/hk/skli_analysis/data/korean_laser_feb19'

def parseargs():

    parser = ArgumentParser()
    parser.add_argument('--diffuser', choices=['barefibre', 'collimator', 'diffuser'], default='barefibre')
    parser.add_argument('--injector', choices=['B1', 'B2', 'B3', 'B4', 'B5'], default='B1')
    parser.add_argument('-b', '--batch', action='store_true', default=False)
    args = parser.parse_args()

    return args

def main():

    args = parseargs()

    if args.batch:
        ROOT.gROOT.SetBatch(True)

    chain, run = load_data(args)

    EventDisplay(chain, run, rot=0.)

    return

def load_data(args):
    print '\nLoading data...'
    flist, sel_run = expand_file_list(args.diffuser, args.injector)
    chain = load_chain(flist)

    return chain, sel_run

def expand_file_list(diffuser, injector, sel_run=None):

    for run in RunInfo.Runs:

        if injector == Injector.tostr(RunInfo.Runs[run].injector) and diffuser == Source.tostr(RunInfo.Runs[run].source):
            sel_run = run
            break

    if sel_run is not None:
        runs = '*%s.root' % str(sel_run)
    else:
        runs = '*.root'

    flist = {sel_run: f for f in glob(os.path.join(data_dir, diffuser, runs))}

    print '\tLoaded files:'

    for f in flist:
        print '\t\t', flist[f]

    return flist, sel_run

def load_chain(flist):
    
    chain = ROOT.TChain("tqtree")
    
    for run in flist:
        chain.AddFile(flist[run])
        
    return chain


if __name__ == "__main__":
    main()
