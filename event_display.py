import os
import ROOT
import pandas as pd
from array import array
import numpy as np
from glob import glob
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import root_numpy
from collections import namedtuple
from scipy.ndimage.filters import median_filter

from scipy import ndimage
ROOT.gROOT.SetBatch(True)
Pos = namedtuple("Pos", ["X", "Y", "Z"])
Inj = namedtuple("Inj", ["Pos", "Tar"])

class Source:
    BARE = 1
    COLLIMATOR = 2
    DIFFUSER = 3

    _names = {BARE:"barefibre",
              COLLIMATOR:"collimator",
              DIFFUSER:"diffuser",
    }

    ALL = [BARE, COLLIMATOR, DIFFUSER]

    @classmethod
    def tostr(cls, code):
        return cls._names[code]

class Injector:
    #OLD_TOP = Pos(-35.3, 777.7, 1802.7) # Location of Jan. SK deployment (AS ON WIKI)
    OLD_TOP = Inj(Pos(-38., 700., 1610.), Pos(-38., 700., -1810.)) # Location of Jan. SK deployment (AS IN THIS ANALYSIS CODE (?))
    NEW_TOP = Inj(Pos(-70.7, -777.7, 1802.7), Pos(-25.00, -694.5, -1810.0)) # Default for vertical laser analysis
    B1 = Inj(Pos(1490.73, 768.14, 1232.25 + 70.7), Pos(-1474.44, -825.362, 1243.0 + 70.7)) # UK injectors, add or subtract one PMT spacing
    B2 = Inj(Pos(1490.73, 768.14, 595.95 + 70.7), Pos(-1453.88, -860.984, 600.0 + 70.7)) # as we are either below or above the existing Korean
    B3 = Inj(Pos(1490.73, 768.14, -40.35 - 70.7), Pos(-1494.65, -788.019, -99.00 - 70.7)) # injectors depending on depth.
    B4 = Inj(Pos(1490.73, 768.14, -605.95 - 70.7), Pos(-1459.59, -851.269, -565.00 - 70.7))
    B5 = Inj(Pos(1490.73, 768.14, -1242.25 - 70.7), Pos(-1427.93, -903.300, -1232.00 - 70.7))
    BOTTOM = Inj(Pos(-70.7, 777.7, -1802.7), Pos(-70.7, 777.7, 1802.7)) # For completeness, do not use.

    _names = {OLD_TOP: "oldtop",
                NEW_TOP: "newtop",
                B1: "B1",
                B2: "B2",
                B3: "B3",
                B4: "B4",
                B5: "B5",
                BOTTOM: "BOTTOM"}

    ALL_BARREL = [B1, B2, B3, B4, B5]
    ALL = [OLD_TOP, B1, B2, B3, B4, B5]

    @classmethod
    def tostr(cls, code):
        return cls._names[code]

class RunInfo:
    RunInfo = namedtuple("RunInfo", ["runnum", "injector", "source", "intensity", "quality"])
    RI = RunInfo
    Runs = {r.runnum:r for r in [
        # Tuesday 23rd January 2018
        RI(77480, Injector.OLD_TOP, Source.BARE, 5, False),
        RI(77481, Injector.OLD_TOP, Source.BARE, 0, False),
        RI(77483, Injector.OLD_TOP, Source.BARE, 6, True),
        RI(77484, Injector.OLD_TOP, Source.BARE, 7, True),
        RI(77485, Injector.OLD_TOP, Source.BARE, 5, True),
        RI(77486, Injector.OLD_TOP, Source.BARE, 8, True),
        RI(77488, Injector.OLD_TOP, Source.COLLIMATOR, 5, True),
        RI(77489, Injector.OLD_TOP, Source.COLLIMATOR, 6, True),
        RI(77490, Injector.OLD_TOP, Source.COLLIMATOR, 7, True),
        # Wednesday 24th January 2018
        RI(77496, Injector.OLD_TOP, Source.COLLIMATOR, 8, True),
        RI(77497, Injector.OLD_TOP, Source.DIFFUSER, 10, True),
        RI(77498, Injector.OLD_TOP, Source.DIFFUSER, 15, True),
        RI(77499, Injector.OLD_TOP, Source.DIFFUSER, 20, True),
        RI(77500, Injector.OLD_TOP, Source.DIFFUSER, 25, True),
        # Tuesday 5th February 2019
        RI(80174, Injector.B1, Source.COLLIMATOR, None, True),
        RI(80175, Injector.B2, Source.COLLIMATOR, None, True),
        RI(80176, Injector.B3, Source.COLLIMATOR, None, True),
        RI(80177, Injector.B4, Source.COLLIMATOR, None, True),
        RI(80178, Injector.B5, Source.COLLIMATOR, None, True),
        RI(80180, Injector.B1, Source.DIFFUSER, None, True),
        RI(80181, Injector.B2, Source.DIFFUSER, None, True),
        RI(80182, Injector.B3, Source.DIFFUSER, None, True),
        RI(80183, Injector.B4, Source.DIFFUSER, None, True),
        RI(80184, Injector.B5, Source.DIFFUSER, None, True),
        RI(80186, Injector.B1, Source.BARE, None, True),
        RI(80187, Injector.B2, Source.BARE, None, True),
        RI(80188, Injector.B3, Source.BARE, None, True),
        RI(80189, Injector.B4, Source.BARE, None, True),
        RI(80190, Injector.B5, Source.BARE, None, True),
    ]}

data_dir = '/Users/billyvinning/Work/hk/skli_analysis/data/korean_laser_feb19'

m = 1.
cm = 1e-2
mm = 1e-3

class sk_constants:
    #-1689.998779296875 1689.998779296875
    #-1689.6700439453125 1689.62939453125
    #-1810.0 1810.0

    IDPMTRadius = .254*m
    WCIDDiameter          = 33.6815*m #16.900*2*cos(2*pi*rad/75)*m; //inner detector diameter
    WCIDHeight            = 36.200*m #"" "" height
    WCBarrelPMTOffset     = 0.0715*m #offset from vertical
    WCBarrelNumPMTHorizontal  = 150 
    WCBarrelNRings        = 17.
    WCPMTperCellHorizontal= 4
    WCPMTperCellVertical  = 3 
    WCCapPMTSpacing       = 0.707*m # distance between centers of top and bottom pmts
    WCCapEdgeLimit        = 16.9*m
    WCBlackSheetThickness = 2.0*cm
    WCIDCircumference = np.pi*WCIDDiameter

def parseargs():

    parser = ArgumentParser()
    parser.add_argument('--diffuser', choices=['barefibre', 'collimator', 'diffuser'], default='barefibre')
    parser.add_argument('--injector', choices=['B1', 'B2', 'B3', 'B4', 'B5'], default='B1')
    args = parser.parse_args()

    return args

def main():

    args = parseargs()
    chain, run = load_data(args)

    EventDisplay(chain, run, rot=0.)

    return

def load_data(args):

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

    print os.path.join(data_dir, diffuser, runs)
    flist = {sel_run: f for f in glob(os.path.join(data_dir, diffuser, runs))}

    print 'Loaded files:'

    for f in flist:
        print '\t', flist[f]

    return flist, sel_run

def load_chain(flist):
    
    chain = ROOT.TChain("tqtree")
    
    for run in flist:
        chain.AddFile(flist[run])
        
    return chain

class EventDisplay():

    def __init__(self, tree, run, wvar='', fit=True, norm='area', rot=0.):
        self._run = run
        self._diffuser = Source.tostr(RunInfo.Runs[run].source)
        self._injector = RunInfo.Runs[run].injector
        self._pmt_df = self._build_pmt_df(tree, wvar, fit)

        if fit and self._diffuser is not 'diffuser':
            self._loc_signal()

        self._plot(fit, rot)

    def _build_pmt_df(self, tree, wvar, fit):
        print '\tBinning...'
        map_df = pd.read_csv('pmt_id_pos_map.csv')

        x_min = 0
        x_max = 11147

        hist = ROOT.TH1F('hist', 'hist', x_max, x_min-0.5, x_max-0.5)

        tree.Draw('cable_vec>>hist', wvar, 'goff')

        if fit:
            hist_barrel = ROOT.TH2F('barrel', 'barrel', sk_constants.WCBarrelNumPMTHorizontal, -sk_constants.WCIDCircumference/2.0, sk_constants.WCIDCircumference/2.0, int(sk_constants.WCBarrelNRings*sk_constants.WCPMTperCellVertical), -sk_constants.WCIDHeight/2.0, sk_constants.WCIDHeight/2.0)

            tree.Draw('(pmtz_vec/100.):((((pmtx_vec/100.)^2 + (pmty_vec/100.)^2  )^0.5)*TMath::ATan2(pmtx_vec/100.,pmty_vec/100.))>>barrel', 'pmtz_vec<1809 && pmtz_vec>-1809'+wvar, 'goff')

            self._barrel_root_hist = hist_barrel

        self._no_events = tree.GetEntries()
        array, edges = root_numpy.hist2array(hist, return_edges=True)

        data = pd.DataFrame({'val': array, 'pmtid': (np.array(edges[0][0:-1])+0.5).astype(int)})

        merged_df = pd.merge(data, map_df, on='pmtid')

        merged_df['pmtx'] = merged_df['pmtx']*cm
        merged_df['pmty'] = merged_df['pmty']*cm
        merged_df['pmtz'] = merged_df['pmtz']*cm

        return merged_df

    def _loc_signal(self):
        print '\tFitting...'

        injector_pos = self._injector.Pos
        target_pos = self._injector.Tar

        beam_l = ((injector_pos.X*cm - target_pos.X*cm)**2 + (injector_pos.Y*cm - target_pos.Y*cm)**2 + (injector_pos.Z*cm - target_pos.Z*cm)**2)**0.5

        barrel_hist = self._barrel_root_hist
        source = self._diffuser
        if source is 'barefibre':

            fwhm_theta_air = 22.80075328628417
            beam_init_r = 0.1*mm

        elif source is 'collimator':

            fwhm_theta_air = 1.80
            beam_init_r = 0.9*mm

        fwhm_r = beam_init_r + beam_l*np.tan(np.radians(fwhm_theta_air*1.0003/1.333))

        target_pos = self._injector.Tar
        target_pos_s = np.arctan2(target_pos.X*cm, target_pos.Y*cm)*((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5



        if source is 'barefibre':

            gaus = ROOT.TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", target_pos_s - fwhm_r, target_pos_s+fwhm_r, target_pos.Z*cm - fwhm_r, target_pos.Z*cm + fwhm_r)
            sigma = 20.0
            gaus.SetParameters(1.0,target_pos_s,sigma,target_pos.Z*cm,sigma)
            gaus.SetParLimits(1, target_pos_s - fwhm_r, target_pos_s+fwhm_r)
            gaus.SetParLimits(2, 0.0, 50.0)
            gaus.SetParLimits(3, target_pos.Z*cm - fwhm_r, target_pos.Z*cm + fwhm_r)
            gaus.SetParLimits(4, 0.0, 50.0)

        elif source is 'collimator':

            gaus = ROOT.TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", target_pos_s - fwhm_r*5.0, target_pos_s+5.0*fwhm_r, target_pos.Z*cm - 5.0*fwhm_r, target_pos.Z*cm + 5.*fwhm_r)
            sigma = 1.25
            gaus.SetParameters(1e4,target_pos_s+2.5,sigma,target_pos.Z*cm -2.0,sigma)
            gaus.SetParLimits(1, target_pos_s - 5.*fwhm_r, target_pos_s+5.*fwhm_r)
            gaus.SetParLimits(2, 0.0, 3.0)
            gaus.SetParLimits(3, target_pos.Z*cm - 5.*fwhm_r, target_pos.Z*cm + 5.*fwhm_r)
            gaus.SetParLimits(4, 0.0, 3.0)

        self._barrel_root_hist.Fit(gaus, 'SURF')
        s, z = gaus.GetParameter(1), gaus.GetParameter(3)
        s_err, z_err = gaus.GetParError(1), gaus.GetParError(3)

        self._signal = (s,z)
        self._signal_err = (s_err,z_err)

        target_pos_r = ((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5

        self._signal_cart = (target_pos_r*np.sin(self._signal[0]/target_pos_r), target_pos_r*np.cos(self._signal[0]/target_pos_r), self._signal[1])

        self._signal_err_cart = (abs(s_err*np.cos(self._signal[0]/target_pos_r)), abs(s_err*np.sin(self._signal[0]/target_pos_r)), z_err)

        return

    def _draw_inj_tar(self, ax, fit):

        source = self._diffuser

        injector_pos = self._injector.Pos
        target_pos = self._injector.Tar

        injector_pos_s = np.arctan2(injector_pos.X*cm, injector_pos.Y*cm)*((injector_pos.X*cm)**2 + (injector_pos.Y*cm)**2)**0.5
        target_pos_s = np.arctan2(target_pos.X*cm, target_pos.Y*cm)*((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5

        beam_l = ((injector_pos.X*cm - target_pos.X*cm)**2 + (injector_pos.Y*cm - target_pos.Y*cm)**2 + (injector_pos.Z*cm - target_pos.Z*cm)**2)**0.5

        ax.plot(injector_pos_s, injector_pos.Z*cm, 'wx', markerfacecolor='blue', alpha=0.35)

        if source == 'barefibre' or source == 'collimator':
            print 'target'
            if source is 'barefibre':

                fwhm_theta_air = 22.80075328628417
                beam_init_r = 0.1*mm

            elif source is 'collimator':

                fwhm_theta_air = 1.80
                beam_init_r = 0.9*mm

            fwhm_r = beam_init_r + beam_l*np.tan(np.radians(fwhm_theta_air*1.0003/1.333))

            fwhm_circle = mpl.patches.Circle((target_pos_s, target_pos.Z*cm), radius=fwhm_r, fill=False, edgecolor='blue', linewidth=0.5, alpha=0.3)
            ax.add_artist(fwhm_circle)

        if fit:

            s_sig, z_sig = self._signal

            ax.plot(s_sig, z_sig, 'rx', markerfacecolor='green', alpha=0.5)

        ax.plot(target_pos_s, target_pos.Z*cm, 'bx', markerfacecolor='red', alpha=0.5)

        return

    def _plot(self, fit=False, rot=0.):

        print 'Plotting...'
        self._setup_pyplot()

        df = self._pmt_df

        if rot is not 0.:
            print 'rot'
            df = self._rotate_detector(df, rot)
        
        df = self._segment_detector(df)

        fig, ax = plt.subplots()

        fig.patch.set_facecolor('xkcd:black')

        det_frame, det_geom = self._draw_detector_frame(ax)
        self._draw_hits(ax, df, det_geom)
        self._draw_inj_tar(ax, fit)
        self._add_text(ax)

        for ext in ['.png', '.pdf']:

            plt.savefig('%s/%s_%s_occ%s' % (self._diffuser, Injector.tostr(self._injector), self._diffuser, ext), dpi=800)

        return

    def _add_text(self, ax):

        ax.text(0.10, 0.98,
        '''
        KOREAN LASER FEB'19

        OCCUPANCIES

        RUN %s
        %s %s
        
        %s EVENTS SCANNED

            ''' % (self._run, Injector.tostr(self._injector), self._diffuser.upper(), str(self._no_events)),
            color='w', fontsize=4, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

        tar_r_vec = (self._injector.Tar.X*cm - self._injector.Pos.X*cm, 
                    self._injector.Tar.Y*cm - self._injector.Pos.Y*cm,
                    self._injector.Tar.Z*cm - self._injector.Pos.Z*cm)

        sig_r_vec = (self._signal_cart[0] - self._injector.Pos.X*cm, 
                    self._signal_cart[1] - self._injector.Pos.Y*cm,
                    self._signal_cart[2] - self._injector.Pos.Z*cm)

        sig_err = self._signal_err_cart

        sig_r = (sig_r_vec[0]**2 + sig_r_vec[1]**2 + sig_r_vec[2]**2)**0.5
        tar_r = (tar_r_vec[0]**2 + tar_r_vec[1]**2 + tar_r_vec[2]**2)**0.5

        sig_r_err = (((sig_r_vec[0]*sig_err[0])**2 + (sig_r_vec[1]*sig_err[1])**2 + (sig_r_vec[2]*sig_err[2])**2)**0.5)/sig_r

        tar_theta = np.arccos(tar_r_vec[2]/tar_r)
        tar_phi = np.arctan(tar_r_vec[1]/tar_r_vec[0])

        sig_theta = np.arccos(sig_r_vec[2]/sig_r)
        sig_phi = np.arctan(sig_r_vec[1]/sig_r_vec[0])

        sig_phi_err = (abs(sig_r_vec[1]/sig_r_vec[0])*((sig_err[1]/sig_r_vec[1])**2 + (sig_err[0]/sig_r_vec[0])**2 )**0.5)/(1. + (sig_r_vec[1]/sig_r_vec[0])**2)

        sig_theta_err = ( abs(sig_r_vec[2]/sig_r)*( (sig_err[2]/sig_r_vec[2])**2 + (sig_r_err/sig_r)**2 )**0.5 )/(  (1. - (sig_r_vec[2]/sig_r)**2)**0.5 )

        off_axis_theta = abs(np.degrees(sig_theta - tar_theta))
        off_axis_phi = abs(np.degrees(sig_phi - tar_phi))

        ax.text(0.85, 0.98, 
            u'''
            INJECTOR COORDS
            [%.2f, %.2f, %.2f] m

            TARGET CENTROID COORDS
            [%.2f, %.2f, %.2f] m

            SIGNAL CENTROID COORDS
            [%.2f, %.2f, %.2f] m 

            SIGNAL OFF-AXIS BY:
            \u03d1 = %.2f \u00b1 %.2f\u00b0
            \u03c6 = %.2f \u00b1 %.2f\u00b0

            ''' % (self._injector.Pos.X*cm, self._injector.Pos.Y*cm, self._injector.Pos.Z*cm, 
                    self._injector.Tar.X*cm, self._injector.Tar.Y*cm, self._injector.Tar.Z*cm,
                    self._signal_cart[0], self._signal_cart[1], self._signal_cart[2],
                    off_axis_theta, abs(np.degrees(sig_theta_err)), off_axis_phi, abs(np.degrees(sig_phi_err)) ),
            color='w', fontsize=4, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

        return

    def _setup_pyplot(self):

        mpl.rcParams['savefig.facecolor'] = 'black'
        mpl.rcParams['savefig.edgecolor'] = 'black'
        mpl.rcParams['font.family'] = 'monospace'
        
        return

    def _rotate_detector(self, data, rot):

        pmtx = data['pmtx'].values
        pmty = data['pmty'].values

        data['pmtx'] = pmtx*np.cos(rot) - pmty*np.sin(rot)
        data['pmty'] = pmtx*np.sin(rot) + pmty*np.cos(rot)

        return data

    def _segment_detector(self, data):

        data['pmtphi'] = np.arctan2(data['pmtx'], data['pmty'])
        data['pmtr'] = np.sqrt(np.power(data['pmtx'], 2.) + np.power(data['pmty'], 2.))
        data['pmts'] = data['pmtr']*data['pmtphi']

        conditions = [ (data['pmtz'] == 18.1*m), 
                        (data['pmtz'] == -18.1*m),
                         ( 18.1*m > data['pmtz']) & (data['pmtz'] > -18.1*m)]

        choices = ['top', 'bottom', 'barrel']
        data['pmt_det_region'] = np.select(conditions, choices, default='none')

        return data

    def _draw_detector_frame(self, ax):

        plt_barrel_width = sk_constants.WCIDCircumference
        plt_barrel_height = sk_constants.WCIDHeight
        plt_id_radius = sk_constants.WCIDDiameter/2.0

        top_c = (0.0, 0.0 + plt_barrel_height)
        bottom_c = (0.0, 0.0 - plt_barrel_height)
        barrel_c = (0.0 - plt_barrel_width/2.0, 0.0 - plt_barrel_height/2.0)

        barrel = mpl.patches.Rectangle(barrel_c, width=plt_barrel_width, height=plt_barrel_height, fill=False, edgecolor='black')
        top_cap = mpl.patches.Circle(top_c, radius=plt_id_radius, fill=False, edgecolor='black')
        bottom_cap = mpl.patches.Circle(bottom_c, radius=plt_id_radius, fill=False, edgecolor='black')

        patches = [barrel, top_cap, bottom_cap]

        collection = mpl.collections.PatchCollection(patches, match_original=True)

        #ax.add_collection(collection)
        plt.axis('equal')
        plt.xlim(-plt_barrel_width/1.8, plt_barrel_width/1.8)
        plt.ylim(-1.1*(2.*plt_id_radius+(plt_barrel_height/2.0)), 1.1*(2.*plt_id_radius+(plt_barrel_height/2.0)))

        plt.axis('off')

        return collection, (barrel_c, top_c, bottom_c)

    def _draw_barrel_hits(self, ax, data, det_geom, cmap):

        data_barrel = data.query('pmt_det_region == \'barrel\'')

        hits = []

        for hit in data_barrel.itertuples(index=True, name='Pandas'):

            s_pos = getattr(hit, 'pmts')
            z_pos = getattr(hit, 'pmtz')
            z_val = getattr(hit, 'val')

            pmt_r = sk_constants.IDPMTRadius
            pmt_patch = mpl.patches.Circle((s_pos, z_pos), pmt_r, fill=True, linewidth=None, facecolor=cmap(z_val))

            hits.append(pmt_patch)

        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)


        ax.add_collection(hit_patches)

        return

    def _draw_top_hits(self, ax, data, det_geom, cmap):

        data_top = data.query('pmt_det_region == \'top\'')

        pmt_r = sk_constants.IDPMTRadius
        hits = []

        for hit in data_top.itertuples(index=True, name='Pandas'):

            x_pos = getattr(hit, 'pmtx')
            y_pos = getattr(hit, 'pmty')
            z_val = getattr(hit, 'val')
            pmt_patch = mpl.patches.Circle((x_pos + det_geom[1][0], y_pos + det_geom[1][1]), pmt_r, fill=True, facecolor=cmap(z_val), linewidth=None)

            hits.append(pmt_patch)


        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)

        ax.add_collection(hit_patches)

        return 

    def _draw_bottom_hits(self, ax, data, det_geom, cmap):

        data_bottom = data.query('pmt_det_region == \'bottom\'')
        pmt_r = sk_constants.IDPMTRadius
        hits = []

        for hit in data_bottom.itertuples(index=True, name='Pandas'):

            x_pos = getattr(hit, 'pmtx')
            y_pos = getattr(hit, 'pmty')
            z_val = getattr(hit, 'val')

            pmt_patch = mpl.patches.Circle((x_pos + det_geom[2][0], y_pos + det_geom[2][1]), pmt_r, fill=True, facecolor=cmap(z_val), linewidth=None)

            hits.append(pmt_patch)

        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)

        ax.add_collection(hit_patches)

        return 

    def _draw_hits(self, ax, data, det_geom):

        zmin = data['val'].values.min()
        zmax = data['val'].values.max()

        cmap = mpl.cm.get_cmap('viridis') 

        data['val'] = data['val']/zmax

        self._draw_barrel_hits(ax, data, det_geom, cmap)
        self._draw_top_hits(ax, data, det_geom, cmap)
        self._draw_bottom_hits(ax, data, det_geom, cmap)

        return

if __name__ == "__main__":
    main()
