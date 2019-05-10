import os
import ROOT
import root_numpy
import numpy as np
import pandas as pd
from glob import glob
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from collections import namedtuple
from argparse import ArgumentParser

from constants import *

class EventDisplay():
    # Recommended cmaps: viridis, inferno, magma, afmhot, kindlmann, kindlmannext
    def __init__(self, tree, run, wvar='occ', fit=True, norm='area', rot=0., cmap='plasma', invert=False):

        self._run = run
        self._diffuser = Source.tostr(RunInfo.Runs[run].source)
        self._injector = RunInfo.Runs[run].injector

        print '\nPreparing run %s (%s %s) event display...' % (self._run, Injector.tostr(self._injector), self._diffuser)

        self._pmt_df = self._build_pmt_df(tree, wvar, fit)

        if fit and self._diffuser is not 'diffuser':
            self._loc_signal()

        self._plot(cmap, fit, rot, invert, wvar)

    def _build_pmt_df(self, tree, wvar, fit):

        print '\tBinning...'

        map_df = pd.read_csv('%s/pmt_id_pos_map.csv' % os.path.dirname(__file__))

        x_min = 0
        x_max = 11147

        hist = ROOT.TH1F('hist', 'hist', x_max, x_min-0.5, x_max-0.5)

        if wvar is 'occ':
            wstr = '1/Entries$'
            self._plot_name = 'OCCUPANCY'
        elif wvar is 'charge':
            wstr = 'charge_vec/Entries$'
            self._plot_name = 'CHARGE DISPLAY'
        tree.Draw('cable_vec>>hist', wstr, 'goff')

        if fit:
            hist_barrel = ROOT.TH2F('barrel', 'barrel', sk_constants.WCBarrelNumPMTHorizontal, -sk_constants.WCIDCircumference/2.0, sk_constants.WCIDCircumference/2.0, int(sk_constants.WCBarrelNRings*sk_constants.WCPMTperCellVertical), -sk_constants.WCIDHeight/2.0, sk_constants.WCIDHeight/2.0)

            tree.Draw('(pmtz_vec/100.):((((pmtx_vec/100.)^2 + (pmty_vec/100.)^2 )^0.5)*TMath::ATan2(pmtx_vec/100.,pmty_vec/100.))>>barrel', '(pmtz_vec<1809 && pmtz_vec>-1809)*(%s)' % wstr, 'goff')

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

        self._barrel_root_hist.Fit(gaus, 'QRSMN')
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

    def _plot(self, cmap, fit=False, rot=0., invert=False, wvar='occ', draw_frame=False):

        print '\tPlotting...'
        self._setup_pyplot(invert)

        df = self._pmt_df

        if int(rot) is not 0:
            df = self._rotate_detector(df, rot)
        
        df = self._segment_detector(df)

        fig, ax = plt.subplots()

        det_frame, det_geom = self._draw_detector_frame(ax, draw_frame)
        self._draw_hits(ax, df, det_geom, cmap, invert)
        self._draw_inj_tar(ax, fit)
        self._add_text(ax)

        print '\tSaving figures...'
        for ext in ['.png', '.pdf']:
            fname = '%s/%s_%s_%s%s' % (self._diffuser, Injector.tostr(self._injector), self._diffuser, wvar, ext)

            plt.savefig(fname, dpi=800)
            print '\t\t%s saved!' % fname
        return

    def _add_text(self, ax):

        ax.text(0.13, 0.98,
        '''
        KOREAN LASER FEB'19

        %s

        RUN %s
        %s %s
        
        %s EVENTS SCANNED

            ''' % (self._plot_name, self._run, Injector.tostr(self._injector), self._diffuser.upper(), str(self._no_events)),
            fontsize=4, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

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
            \u03b8 = %.2f \u00b1 %.2f\u00b0
            \u03c6 = %.2f \u00b1 %.2f\u00b0

            ''' % (self._injector.Pos.X*cm, self._injector.Pos.Y*cm, self._injector.Pos.Z*cm, 
                    self._injector.Tar.X*cm, self._injector.Tar.Y*cm, self._injector.Tar.Z*cm,
                    self._signal_cart[0], self._signal_cart[1], self._signal_cart[2],
                    off_axis_theta, abs(np.degrees(sig_theta_err)), off_axis_phi, abs(np.degrees(sig_phi_err)) ),
            fontsize=4, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

        return

    def _setup_pyplot(self, invert):

        mpl.rcParams['font.monospace'] = "DIN Alternate"
        mpl.rcParams['font.family'] = "monospace"

        if not invert:
            mpl.rcParams['savefig.facecolor'] = 'black'
            mpl.rcParams['savefig.edgecolor'] = 'black'
            mpl.rcParams['text.color'] = 'w'
            mpl.rcParams['patch.edgecolor'] = 'w'
        else:
            mpl.rcParams['savefig.facecolor'] = 'white'
            mpl.rcParams['savefig.edgecolor'] = 'white'
            mpl.rcParams['text.color'] = 'k'
            mpl.rcParams['patch.edgecolor'] = 'k'
        
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

    def _draw_detector_frame(self, ax, draw_frame):

        plt_barrel_width = sk_constants.WCIDCircumference
        plt_barrel_height = sk_constants.WCIDHeight
        plt_id_radius = sk_constants.WCIDDiameter/2.0

        top_c = (0.0, 0.0 + plt_barrel_height)
        bottom_c = (0.0, 0.0 - plt_barrel_height)
        barrel_c = (0.0 - plt_barrel_width/2.0 - 0.5*sk_constants.WCCapPMTSpacing, 0.0 - plt_barrel_height/2.0)

        barrel = mpl.patches.Rectangle(barrel_c, width=plt_barrel_width+sk_constants.WCCapPMTSpacing, height=plt_barrel_height, fill=False, linewidth=0.1, alpha=0.5)
        top_cap = mpl.patches.Circle(top_c, radius=plt_id_radius+0.5*sk_constants.WCCapPMTSpacing, fill=False, linewidth=0.1, alpha=0.5)
        bottom_cap = mpl.patches.Circle(bottom_c, radius=plt_id_radius+0.5*sk_constants.WCCapPMTSpacing, fill=False, linewidth=0.1, alpha=0.5)

        patches = [barrel, top_cap, bottom_cap]

        collection = mpl.collections.PatchCollection(patches, match_original=True)
        if draw_frame:
            ax.add_collection(collection)
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

    def _draw_hits(self, ax, data, det_geom, cmap_name, invert):

        zmin = data['val'].values.min()
        zmax = data['val'].values.max()

        if cmap_name is not 'kindlmann' and cmap_name is not 'kindlmannext':

            cmap = mpl.cm.get_cmap(cmap_name)

        elif cmap_name is 'kindlmann':

            cmap_csv = pd.read_csv('%s/kindlmann-table-float-1024.csv' % os.path.dirname(__file__))
            cmap = mpl.colors.ListedColormap(cmap_csv[['RGB_r', 'RGB_g', 'RGB_b']].values)

        elif cmap_name is 'kindlmannext':

            cmap_csv = pd.read_csv('%s/extended-kindlmann-table-float-1024.csv' % os.path.dirname(__file__))
            cmap = mpl.colors.ListedColormap(cmap_csv[['RGB_r', 'RGB_g', 'RGB_b']].values)

        if invert:
            cmap = cmap.reversed()

        data['val'] = data['val']/zmax

        self._draw_barrel_hits(ax, data, det_geom, cmap)
        self._draw_top_hits(ax, data, det_geom, cmap)
        self._draw_bottom_hits(ax, data, det_geom, cmap)

        return
