import os
import ROOT
import time
import root_numpy
import numpy as np
import pandas as pd
from glob import glob
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from collections import namedtuple
from argparse import ArgumentParser
import datetime
from hepunits.constants import c_light
from constants import *

class EventDisplay():
    
    def __init__(self, tree, run, wvar='occupancy', fit=True, norm='area', rot=0., cmap='plasma', invert=False, draw_frame=False, tof_cut=True, event_id_cut=True, draw_timing=True, correct=False, logz=False, tof_cut_override=None):
        """Generates an Super-Kamiokande event display in Matplotlib from a ROOT tree.

            
        Note:
            We recommend perceptually uniform colourmaps for displaying physics data for clarity, such viridis, inferno, plasma, magma, afmhot, kindlmann and kindlmannext.

        Args:
            tree (ROOT:TTree/ROOT:TChain): The tree or chain to read data from.
            run (int): The run number of the data loaded.
            wvar (str, optional): The plot type to generate, choices are `occ` (occupancy) or `charge` (charge weighted), defaults to `occ`.
            fit (bool, optional): Whether to make a fit to the wall data in order to locate the signal, defaults to True, overriden to False later if the source is a diffuser.
            norm (str, optional): The normalisation mode, not implemented properly, defaults to `area`.
            rot (float, optional): The angle to rotate the detector around the z-axis, doesn't work properly, defaults to 0.0.
            cmap (str, optional): The default colourmap to use, this can be chosen from ones built-in to matplotlib or which have been interfaced from .csv in `_get_cmap`. 
            invert (bool, optional): Inverts the colours of the plot (from dark mode to light mode), defaults to False.
            draw_frame (bool, optional): Draws a detector frame around the detector regions, defaults to False.
            tof_cut (bool, optional): Whether to cut on the time-of-flight corrected timing distribution, defaults to True.
            event_id_cut (bool, optional): Whether to cut on the time-of-flight corrected timing distribution, defaults to True.
            draw_timing (bool, optional): Whether to draw the time-of-flight corrected timing distribution in lower right of the event display, defaults to True.
            correct (bool, optional): Whether to apply corrections to the data, the types of correction to be applied are selected in `_do_corrections`, defaults to False.
            logz (bool, optional): Whether to make the colourbar scale logarithmic, defaults to False.
            tof_cut_override (tuple(float, float), optional): Overrides the default time-of-flight time cut in `_calc_hit_in_time`, defaults to None.

        """

        self._tree = tree 
        self._run = run 

        self._diffuser = Source.tostr(RunInfo.Runs[run].source) #:str: The name of the diffuser type (eg. diffuser, collimator, barefibre).
        self._injector = RunInfo.Runs[run].injector #:namedtuple: The injector position and target position contained within a namedtuple.
        self._monitor = RunInfo.Runs[run].monitor #:monitor: The ID of the monitor available for this run.

        print '\nPreparing run %s (%s %s) %s event display...' % (self._run, Injector.tostr(self._injector), self._diffuser, wvar)

        # A dataframe storing information for every ID-PMT integrated over the run.
        self._pmt_df = self._build_pmt_df(tree=tree, 
                                            wvar=wvar, 
                                            fit=fit, 
                                            tof_cut=tof_cut, 
                                            event_id_cut=event_id_cut, 
                                            draw_timing=draw_timing, 
                                            tof_cut_override=tof_cut_override)

        if fit and self._diffuser is not 'diffuser': 
            self._loc_signal() # Attempt to find the wall signal with a fit.

        if correct: 
            self._do_corrections() # Apply some corrections.

        self._plot(cmap=cmap, 
                    fit=fit, 
                    rot=rot, 
                    invert=invert, 
                    wvar=wvar, 
                    draw_frame=draw_frame, 
                    draw_timing=draw_timing, 
                    correct=correct, 
                    logz=logz) # Draw the plot.

    def _do_corrections(self, gain=False, solid_angle=True, angular=False, attenuation=False):
        """Applies selected corrections to the PMT DataFrame.

        Note:
            The gain correction is not implemented and maybe the angular correction is broken, I don't remember.

        Args:
            gain (bool, optional): Whether to apply the gain correction.
            solid_angle (bool, optional): Whether to apply to solid angle correction.
            angular (bool, optional): Whether to apply the angular correction.
            attenuation (bool, optional): Whether to apply the attenuation correction.

        Returns:
            None

        """

        self._pmt_df['tot_cor'] = 1.0

        self._add_to_df()

        if gain:
            self._do_gain_correction()
        if solid_angle:
            self._do_solid_angle_correction()
        if angular:
            self._do_angular_correction()
        if attenuation:
            self._do_attenuation_correction()

        corrections = {'gain': gain, 'solid angle': solid_angle, 'pmt angular acceptance' : angular, 'water attenuation' : attenuation}

        correction_str = ''

        for correction in corrections: # Build string to indicate corrections on plot.

            if corrections[correction]:

                correction_str += '''
                %s ''' % correction.upper()

        self._corrections = correction_str

    def _do_gain_correction(self):
        """Applies the gain correction to the PMT DataFrame.

        Note:
            Not implemented.

        Returns:
            None

        """
        raise NotImplementedError('Gain correction not implemented.')

    def _add_to_df(self):
        """Adds some useful variables characterising the source directionality to the PMT DataFrame.

        Returns:
            None

        """
        injector = self._injector.Pos
        target = self._injector.Tar

        self._pmt_df['pmt_inj_vec_x'] = self._pmt_df['pmtx'] - injector.X*cm
        self._pmt_df['pmt_inj_vec_y'] = self._pmt_df['pmty'] - injector.Y*cm
        self._pmt_df['pmt_inj_vec_z'] = self._pmt_df['pmtz'] - injector.Z*cm

        self._pmt_df['pmt_inj_vec_mag'] = np.sqrt(np.square(self._pmt_df['pmt_inj_vec_x']) + np.square(self._pmt_df['pmt_inj_vec_y']) + np.square(self._pmt_df['pmt_inj_vec_z'])) # Distance from PMT to injector.

        self._pmt_df['tar_inj_vec_x'] = (target.X - injector.X)*cm
        self._pmt_df['tar_inj_vec_y'] = (target.Y - injector.Y)*cm
        self._pmt_df['tar_inj_vec_z'] = (target.Z - injector.Z)*cm

        self._pmt_df['tar_inj_vec_mag'] = np.sqrt(np.square(self._pmt_df['tar_inj_vec_x']) + np.square(self._pmt_df['tar_inj_vec_y']) + np.square(self._pmt_df['tar_inj_vec_z'])) # Distance from target to injector.

        self._pmt_df['theta_dir'] = np.arccos((self._pmt_df['tar_inj_vec_x']*self._pmt_df['pmt_inj_vec_x'] + self._pmt_df['tar_inj_vec_y']*self._pmt_df['pmt_inj_vec_y'] + self._pmt_df['tar_inj_vec_z']*self._pmt_df['pmt_inj_vec_z'])/(self._pmt_df['pmt_inj_vec_mag']*self._pmt_df['tar_inj_vec_mag'])) # Angle subtending vector from injector position to PMT position and vector from injector position to target position.

    def _do_solid_angle_correction(self):
        """Applies the solid angle correction to the PMT DataFrame.

        Returns:
            None

        """
        self._pmt_df['solid_angle_corr'] = 1./(1.-(np.abs(self._pmt_df['theta_dir'])/np.pi))
        self._pmt_df['tot_cor'] = self._pmt_df['tot_cor']*self._pmt_df['solid_angle_corr']

    def _do_angular_correction(self):
        """Applies the angular correction to the PMT DataFrame.

        Note:
            I don't remember if this works or not.

        Returns:
            None

        """
        conditions = [ self._pmt_df['pmt_det_region'] == 'top', 
                        self._pmt_df['pmt_det_region'] == 'bottom',
                        self._pmt_df['pmt_det_region'] == 'barrel']

        choices = [1./np.square(np.cos(np.pi/2. - self._pmt_df['theta_dir'])),
                    1./np.square(np.cos(np.pi/2. - self._pmt_df['theta_dir'])), 
                    1./np.square(np.cos(self._pmt_df['theta_dir']))]

        self._pmt_df['angular_correction'] = np.select(conditions, choices, default=1.0)

        self._pmt_df['tot_cor'] = self._pmt_df['tot_cor']*self._pmt_df['angular_correction']

    def _do_attenuation_correction(self):
        """Applies an attenuation correction to the PMT DataFrame.

        Returns:
            None

        """

        l = 9800*cm # attenuation length in cm

        self._pmt_df['attenuation_correction'] = np.exp(-(np.abs(self._pmt_df['pmtz'] - self._injector.Pos.Z*cm))/l)/np.exp(-(self._pmt_df['pmt_inj_vec_mag'])/l)
        self._pmt_df['tot_cor'] = self._pmt_df['tot_cor']*self._pmt_df['attenuation_correction']

    def _get_tof_time_exp(self, monitor=True):
        """Gets the time-of-flight corrected detector timing varexp string to use with ROOT:TTree.Draw().

        Args:
            monitor (bool, optional): Selects whether detector timing should be relative to the monitor timing, defaults to True.

        Returns:
            str: The time-of-flight corrected varexp.

        """

        if monitor: # Use monitor information.

            tof_time_exp = '(time_vec - 1.0e9*sqrt(pow(pmtx_vec/100. - %s/100., 2)+pow(pmty_vec/100. - %s/100., 2)+pow(pmtz_vec/100. - %s/100., 2))/%s - (Sum$(mon_time_vec*(mon_cable_vec==%s))/Sum$((mon_cable_vec==%s))))' % (self._injector.Pos.X, self._injector.Pos.Y, self._injector.Pos.Z, c_light*1e6/1.333, self._monitor, self._monitor)
        else: # Exclude monitor information.
            tof_time_exp = '(time_vec - 1.0e9*sqrt(pow(pmtx_vec/100. - %s/100., 2)+pow(pmty_vec/100. - %s/100., 2)+pow(pmtz_vec/100. - %s/100., 2))/%s)' % (self._injector.Pos.X, self._injector.Pos.Y, self._injector.Pos.Z, c_light*1e6/1.333)

        return tof_time_exp

    def _calc_hit_in_time(self, override=None):
        """Gets the time-of-flight corrected timing cut expressions to use as selection in ROOT:TTree:Draw(varexp, selection), selects the signal region.

        Args:
            override (tuple(float, float)): These timing markers will override the default TOF timing markers f not None, defaults to None.

        Returns:
            tuple(str, tuple(float, float)): The selection expression and the TOF timing markers applied.

        """
        injector = self._injector
        source = self._diffuser

        #CORRECTED_TIME = "(time_vec - 1.0e9*sqrt(pow(pmtx_vec/100. - %s/100., 2)+pow(pmty_vec/100. - %s/100., 2)+pow(pmtz_vec/100. - %s/100., 2))/%s)" % (injector.Pos.X, injector.Pos.Y, injector.Pos.Z, c_light*1e6/1.333)
        #CORRECTED_TIME = "(time_vec - 1.0e9*sqrt(pow(pmtx_vec/100. - %s/100., 2)+pow(pmty_vec/100. - %s/100., 2)+pow(pmtz_vec/100. - %s/100., 2))/%s - (Sum$(mon_time_vec*(mon_cable_vec==11256))/Sum$((mon_cable_vec==11256))))" % (injector.Pos.X, injector.Pos.Y, injector.Pos.Z, c_light*1e6/1.333)

        if self._monitor is not None:
            CORRECTED_TIME = self._get_tof_time_exp(self._monitor)
        else:
            CORRECTED_TIME = self._get_tof_time_exp()

        if override is not None:
            print '\t\tTiming selection overriden with markers: %s.' % str(override)
            return ('(%s > %s)&&(%s < %s)' % (CORRECTED_TIME, override[0], CORRECTED_TIME, override[1]), override)
        #diffuser_t = (1060, 1095)
        #collimator_t = (1065, 1090)
        #bare_t = (1060, 1090)

        diffuser_t = (600, 850)
        collimator_t = (600, 1000)
        bare_t = (600, 780)

        if Injector.tostr(injector) == 'OLD_TOP':
          if source == 'barefibre':
            hit_in_time = "(%s > 700.) && (%s < 730.)" % (CORRECTED_TIME, CORRECTED_TIME)
          elif source == 'collimator':
            hit_in_time = "(%s > 700.) && (%s < 730.)" % (CORRECTED_TIME, CORRECTED_TIME)
          elif source == 'diffuser':
            hit_in_time = "(%s > 700.) && (%s < 820.)" % (CORRECTED_TIME, CORRECTED_TIME)
        else:
          if source == 'barefibre':
            hit_in_time = "(%s > %s) && (%s < %s)" % (CORRECTED_TIME, bare_t[0], CORRECTED_TIME, bare_t[1])
          elif source == 'collimator':
            hit_in_time = "(%s > %s) && (%s < %s)" % (CORRECTED_TIME, collimator_t[0], CORRECTED_TIME, collimator_t[1])
          elif source == 'diffuser':
            hit_in_time = "(%s > %s) && (%s < %s)" % (CORRECTED_TIME, diffuser_t[0], CORRECTED_TIME, diffuser_t[1])

        if source == 'barefibre':
            return (hit_in_time, bare_t)
        elif source == 'collimator':
            return (hit_in_time, collimator_t)
        elif source == 'diffuser':
            return (hit_in_time, diffuser_t)

    def _get_event_id_cut_exp(self):
        """Gets a selection expression for cutting on event ID, this is used when sometimes runs don't end properly. Cut tables stored in the cuts/ folder, if theres a match to the current run the appropriate event ID cut is loaded, otherwise no cut is applied.

        Returns:
            tuple(str, tuple(float, float)): The event ID selection expression.

        """

        flist = glob(os.path.join(os.path.dirname(__file__), 'cuts', '*.csv')) # Load all cut tables

        cut_df = pd.concat([pd.read_csv(f) for f in flist])

        matched_df =  cut_df.loc[cut_df['Run No.'] == self._run] # Find matching entry for this run

        if not matched_df.empty:
            event_cut = matched_df['Event No.'].values[0]
            varexp = '(nev < %s)' % event_cut
        else:
            varexp = '(1)'

        return varexp
        
    def _build_pmt_df(self, tree, wvar, fit, tof_cut, event_id_cut, draw_timing, tof_cut_override):
        """Builds a pandas DataFrame of information (position, id, z value...) for each ID-PMT. 

        Args:
            tree (ROOT:TTree/ROOT:TChain): The tree or chain to read data from.
            wvar (str, optional): The plot type to generate, choices are `occupancy` (occupancy) or `charge` (charge weighted), defaults to `occupancy`.
            fit (bool, optional): Whether to make a fit to the wall data in order to locate the signal, defaults to True, overriden to False later if the source is a diffuser.
            tof_cut (bool, optional): Whether to cut on the time-of-flight corrected timing distribution, defaults to True.
            event_id_cut (bool, optional): Whether to cut on the time-of-flight corrected timing distribution, defaults to True.
            draw_timing (bool, optional): Whether to draw the time-of-flight corrected timing distribution in lower right of the event display, defaults to True.
            tof_cut_override (tuple(float, float), optional): Overrides the default time-of-flight time cut in `_calc_hit_in_time`, defaults to None.

        Returns:
            pandas:DataFrame: The PMT information DataFrame.

        """

        print '\tBinning...'

        map_df = pd.read_csv('%s/tables/pmt_id_pos_map.csv' % os.path.dirname(__file__)) # Reads a file which maps PMT ID to position in the tank.

        x_min = 0
        x_max = 11147 

        hist = ROOT.TH1F('hist', 'hist', x_max, x_min-0.5, x_max-0.5) # Init histogram of PMT ID.

        if wvar is 'occupancy':
            wstr = '1/Entries$'
            self._plot_name = 'OCCUPANCY'
        elif wvar is 'charge':
            wstr = 'charge_vec/Entries$'
            self._plot_name = 'CHARGE DISPLAY'


        if not tof_cut and not event_id_cut:
            varexp = '(%s)' % wstr

        elif tof_cut and not event_id_cut:
            time_str, self._time_markers = self._calc_hit_in_time(tof_cut_override)
            varexp = '(%s)*(%s)' % (wstr, time_str)

        elif tof_cut and event_id_cut:
            time_str, self._time_markers = self._calc_hit_in_time(tof_cut_override)
            event_id_sel_exp = self._get_event_id_cut_exp()
            varexp = '(%s)*(%s)*(%s)' % (wstr, time_str, event_id_sel_exp)

        tree.Draw('cable_vec>>hist', varexp, 'goff') # Fill PMT ID histogram.

        if draw_timing: # Generate the TOF corrected timing histogram.

            self._timing_data = root_numpy.tree2array(tree, '(%s)*(%s)' % (self._get_tof_time_exp(self._monitor), self._get_event_id_cut_exp()))

        if fit: # Generate a ROOT histogram of wall hits.

            hist_barrel = ROOT.TH2F('barrel', 'barrel', sk_constants.WCBarrelNumPMTHorizontal, -sk_constants.WCIDCircumference/2.0, sk_constants.WCIDCircumference/2.0, int(sk_constants.WCBarrelNRings*sk_constants.WCPMTperCellVertical), -sk_constants.WCIDHeight/2.0, sk_constants.WCIDHeight/2.0)

            tree.Draw('(pmtz_vec/100.):((((pmtx_vec/100.)^2 + (pmty_vec/100.)^2 )^0.5)*TMath::ATan2(pmtx_vec/100.,pmty_vec/100.))>>barrel', '(pmtz_vec<1809 && pmtz_vec>-1809)*(%s)' % wstr, 'goff')

            self._barrel_root_hist = hist_barrel

        self._no_events = tree.GetEntries()

        for i, event in enumerate(tree): # Get run start and end timestamps.

            if i == 0: # Gets first event
                year_start = event.year
                month_start = event.month
                day_start = event.day
                hour_start = event.hour
                minute_start = event.minute
                second_start = event.second
                millisecond_start = event.millisecond

            if i == tree.GetEntries()-1: # Gets last event
                year_end = event.year
                month_end = event.month
                day_end = event.day
                hour_end = event.hour
                minute_end = event.minute
                second_end = event.second
                millisecond_end = event.millisecond
        
        self._run_start_tree = datetime.datetime(year_start, month_start, day_start, hour_start, minute_start, second_start, millisecond_start)
        self._run_end_tree = datetime.datetime(year_end, month_end, day_end, hour_end, minute_end, second_end, millisecond_end)

        array, edges = root_numpy.hist2array(hist, return_edges=True) # Convert PMT ID histogram to numpy array.

        data = pd.DataFrame({'val': array, 'pmtid': (np.array(edges[0][0:-1])+0.5).astype(int)}) # Fill DataFrame with PMT ID and value of bin. 

        merged_df = pd.merge(data, map_df, on='pmtid') # Match PMT positions to PMT ID.

        merged_df['pmtx'] = merged_df['pmtx']*cm
        merged_df['pmty'] = merged_df['pmty']*cm
        merged_df['pmtz'] = merged_df['pmtz']*cm

        merged_df = self._segment_detector(merged_df) # Add some additional columns to DataFrame.

        return merged_df

    def _loc_signal(self):
        """Locates the signal by attempting a 2D Gaussian to the wall. 

        Returns:
            None

        """
        print '\tFitting...'

        injector_pos = self._injector.Pos
        target_pos = self._injector.Tar

        beam_l = ((injector_pos.X*cm - target_pos.X*cm)**2 + (injector_pos.Y*cm - target_pos.Y*cm)**2 + (injector_pos.Z*cm - target_pos.Z*cm)**2)**0.5

        barrel_hist = self._barrel_root_hist
        source = self._diffuser

        if source is 'barefibre':

            fwhm_theta_air = 22.80075328628417 # The FWHM half-opening angle of the bare fibre as measured in air.
            beam_init_r = 0.1*mm # Assume the beam is initially as large as the core radius of the fibre.

        elif source is 'collimator':

            fwhm_theta_air = 1.80 # The FWHM half-opening angle of the collimator as measured in air.
            beam_init_r = 0.9*mm # Assume the beam is initially as large as the secondary aperture of the collimator.

        fwhm_r = beam_init_r + beam_l*np.tan(np.radians(fwhm_theta_air*1.0003/1.333))

        target_pos = self._injector.Tar
        target_pos_s = np.arctan2(target_pos.X*cm, target_pos.Y*cm)*((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5

        if source is 'barefibre': # Do fit to a bare fibre.

            gaus = ROOT.TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", target_pos_s - fwhm_r, target_pos_s+fwhm_r, target_pos.Z*cm - fwhm_r, target_pos.Z*cm + fwhm_r)
            sigma = 10.0
            gaus.SetParameters(1.0,target_pos_s,sigma,target_pos.Z*cm,sigma)
            gaus.SetParLimits(1, target_pos_s - 5.0*fwhm_r, target_pos_s+5.0*fwhm_r)
            gaus.SetParLimits(2, 0.0, 50.0)
            gaus.SetParLimits(3, target_pos.Z*cm - 5.*fwhm_r, target_pos.Z*cm + 5.*fwhm_r)
            gaus.SetParLimits(4, 0.0, 50.0)

        elif source is 'collimator': # Do fit to a collimator.

            gaus = ROOT.TF2("f2","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", target_pos_s - fwhm_r*7.0, target_pos_s+7.0*fwhm_r, target_pos.Z*cm - 7.0*fwhm_r, target_pos.Z*cm + 7.*fwhm_r)
            sigma = 1.
            gaus.SetParameters(1e4,target_pos_s,sigma,target_pos.Z*cm,sigma)
            gaus.SetParLimits(1, target_pos_s - 4.*fwhm_r, target_pos_s+4.*fwhm_r)
            gaus.SetParLimits(2, 0.0, 3.0)
            gaus.SetParLimits(3, target_pos.Z*cm - 4.*fwhm_r, target_pos.Z*cm + 4.*fwhm_r)
            gaus.SetParLimits(4, 0.0, 3.0)

        self._barrel_root_hist.Fit(gaus, 'QRSMN')
        s, z = gaus.GetParameter(1), gaus.GetParameter(3)
        s_err, z_err = gaus.GetParError(1), gaus.GetParError(3)

        self._signal = (s,z)
        self._signal_err = (s_err,z_err)

        target_pos_r = ((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5

        self._signal_cart = (target_pos_r*np.sin(self._signal[0]/target_pos_r), target_pos_r*np.cos(self._signal[0]/target_pos_r), self._signal[1]) # Signal position in the tank in cartesian coordinates.

        self._signal_err_cart = (abs(s_err*np.cos(self._signal[0]/target_pos_r)), abs(s_err*np.sin(self._signal[0]/target_pos_r)), z_err) # Signal position error in the tank in cartesian coordinates.

    def _draw_inj_tar(self, ax, fit, det_geom):
        """Draws the injector, target, signal positions onto the plot. 

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            fit (bool): Whether there was fit made to the data and whether to draw it.
            det_geom (tuple): Some storage of geometry variables for the detector plot.

        Returns:
            None

        """
        top_c, bottom_c = det_geom[1], det_geom[2]

        source = self._diffuser

        injector_pos = self._injector.Pos
        target_pos = self._injector.Tar

        injector_pos_s = np.arctan2(injector_pos.X*cm, injector_pos.Y*cm)*((injector_pos.X*cm)**2 + (injector_pos.Y*cm)**2)**0.5
        target_pos_s = np.arctan2(target_pos.X*cm, target_pos.Y*cm)*((target_pos.X*cm)**2 + (target_pos.Y*cm)**2)**0.5

        beam_l = ((injector_pos.X*cm - target_pos.X*cm)**2 + (injector_pos.Y*cm - target_pos.Y*cm)**2 + (injector_pos.Z*cm - target_pos.Z*cm)**2)**0.5

        ax.plot(injector_pos_s, injector_pos.Z*cm, 'wx', alpha=0.35, label='INJECTOR')
        ax.plot(injector_pos.X*cm+top_c[0], injector_pos.Y*cm+top_c[1], 'wx', markersize=3, alpha=0.35)
        ax.plot(injector_pos.X*cm+bottom_c[0], injector_pos.Y*cm+bottom_c[1], 'wx', markersize=3, alpha=0.35)

        if source == 'barefibre' or source == 'collimator':  # Draw circle representing the expected FWHM angle of either the bare fibre or collimator.

            if source is 'barefibre':

                fwhm_theta_air = 22.80075328628417
                full_theta_air = np.degrees(0.5)
                beam_init_r = 0.1*mm

            elif source is 'collimator':

                fwhm_theta_air = 1.80
                beam_init_r = 0.9*mm

            fwhm_r = beam_init_r + beam_l*np.tan(np.radians(fwhm_theta_air*1.0003/1.333))
            #full_r = beam_init_r + beam_l*np.tan(np.radians(full_theta_air*1.0003/1.333))

            fwhm_circle = mpl.patches.Circle((target_pos_s, target_pos.Z*cm), radius=fwhm_r, fill=False, edgecolor='cyan', linestyle='--', linewidth=0.5, alpha=0.3)

            #full_circle = mpl.patches.Circle((target_pos_s, target_pos.Z*cm), radius=full_r, fill=False, edgecolor='cyan', linewidth=0.5, alpha=0.3)
            #ax.add_artist(fwhm_circle)
            #ax.add_artist(full_circle)

        if fit and source != 'diffuser': # Draw the fit results onto the plot.

            s_sig, z_sig = self._signal
            x_sig, y_sig, z_sig = self._signal_cart

            ax.plot(s_sig, z_sig, 'rx', alpha=0.5, label='SIGNAL')
            ax.plot(x_sig+top_c[0], y_sig+top_c[1], 'rx', markersize=3, alpha=0.5)
            ax.plot(x_sig+bottom_c[0], y_sig+bottom_c[1], 'rx', alpha=0.5, markersize=3)

        ax.plot(target_pos_s, target_pos.Z*cm, 'bx', alpha=0.5, label='TARGET')
        ax.plot(target_pos.X*cm+top_c[0], target_pos.Y*cm+top_c[1], 'bx', markersize=3, alpha=0.5)
        ax.plot(target_pos.X*cm+bottom_c[0], target_pos.Y*cm+bottom_c[1], 'bx', markersize=3, alpha=0.5)
        ax.legend(loc=(0.795, 0.670), markerscale=0.7)._legend_box.align='right'

    def _plot(self, cmap, fit, rot, invert, wvar, draw_frame, draw_timing, correct, logz, save_dir=os.getcwd()):
        """Draws the event display, with all details. 

        Args:
            cmap (str): The colourmap to use, this can be chosen from ones built-in to matplotlib or which have been interfaced from .csv in `_get_cmap`. 
            fit (bool): Whether to make a fit to the wall data in order to locate the signal.
            rot (float): The angle to rotate the detector around the z-axis, doesn't work properly.
            invert (bool): Inverts the colours of the plot (from dark mode to light mode).
            wvar (str): The plot type to generate, choices are `occupancy` (occupancy) or `charge` (charge weighted).
            draw_frame (bool): Draws a detector frame around the detector regions.
            draw_timing (bool): Whether to draw the time-of-flight corrected timing distribution in lower right of the event display.
            correct (bool): Whether to apply corrections to the data, the types of correction to be applied are selected in `_do_corrections`.
            logz (bool): Whether to make the colourbar scale logarithmic.

        Returns:
            None

        """
        print '\tPlotting...'

        df = self._pmt_df

        self._setup_pyplot(invert, usetex=False) # Setup matplotlib rcParams for event display plot.

        if int(rot) is not 0:
            df = self._rotate_detector(df, rot)

        fig, ax = plt.subplots(figsize=(5.2, 4.8))

        det_frame, det_geom = self._draw_detector_frame(ax, draw_frame) # Draw frames around each detector region.
        self._draw_hits(ax, df, det_geom, cmap, invert, logz) # Draw the hits.
        self._draw_inj_tar(ax, fit, det_geom) # Draw markers for the injector, target and signal positions.
        self._add_text(ax, fit, correct) # Add run information to the top of the plot.
        if draw_timing:
            self._add_timing_plot(fig) # Add the TOF corrected timing distribution to the plot.
        self._add_colourbar(fig, cmap, invert) # Add the colourbar.

        plt_barrel_width = sk_constants.WCIDCircumference
        plt_barrel_height = sk_constants.WCIDHeight
        plt_id_radius = sk_constants.WCIDDiameter/2.0
        
        ax.set_xlim(-plt_barrel_width/1.8, plt_barrel_width/1.8)

        print '\tSaving figures...'
        for ext in ['.png', '.pdf']:
            if logz:
                fname = '%s/%s_%s_%s_log%s' % (self._diffuser, Injector.tostr(self._injector), self._diffuser, wvar, ext)
            else:
                fname = '%s/%s_%s_%s%s' % (self._diffuser, Injector.tostr(self._injector), self._diffuser, wvar, ext)
            
            par_dir = os.path.abspath(os.path.join(save_dir, fname, '..'))

            if not os.path.exists(par_dir):
                os.makedirs(par_dir)

            plt.savefig(fname, dpi=1400, bbox_inches='tight', pad_inches=0)
            print '\t\t%s saved!' % os.path.join(save_dir, fname)

    def _add_colourbar(self, fig, cmap_name, invert):
        """Draws the event display, with all details. 

        Args:
            fig (matplotlib.Figure): The Matplotlib figure being used to draw the event display.
            cmap_name (str): The name of the colourmap to use, this can be chosen from ones built-in to matplotlib or which have been interfaced from .csv in `_get_cmap`. 
            invert (bool): Inverts the colours of the plot (from dark mode to light mode), defaults to False.

        Returns:
            None

        """
        ax = fig.add_axes([0.85, 0.15, 0.005, 0.7], frameon=True, facecolor='w')

        cmap = self._get_cmap(cmap_name, invert)

        norm = self._cmap_norm
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)
        cb1.outline.set_visible(False)

    def _add_timing_plot(self, fig):
        """Draws the time-of-flight corrected timing distribution onto the plot.

        Args:
            fig (matplotlib.Figure): The Matplotlib figure being used to draw the event display.

        Returns:
            None

        """
        timing_data = np.concatenate(self._timing_data).ravel()

        ax = fig.add_axes([0.635, 0.17, 0.195, 0.18], facecolor='k')

        n, bins, patches = ax.hist(timing_data, 5000, edgecolor='w', facecolor='k', linewidth=0.1, histtype='step')
        ax.set_xlabel('TOF CORRECTED TIME W.R.T MONITOR (ns)', fontsize=5)
        
        #ax.set_ylim(0.0, 0.2e6)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        #ax.set_ylim(0.0, 0.1e6)
        #ax.set_xlim(1.025e3, 1.2e3)
        ax.set_xlim(0, 1200)

        t1, t2 = self._time_markers

        ax.axvline(t1, color='w', ls='--', dashes=(5, 10), linewidth=0.1)
        ax.axvline(t2, color='w', ls='--', dashes=(5, 10), linewidth=0.1)

    def _add_text(self, ax, fit, correct):
        """Add various information about the run to the top of the plot. 

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            fit (bool): Whether to make a fit to the wall data in order to locate the signal.
            correct (bool): Whether to apply corrections to the data, the types of correction to be applied are selected in `_do_corrections`.

        Returns:
            None

        """
        ax.text(0.077, 0.97,
        '''
        LIVERPOOL LASER JULY'19

        RUN %s
        %s %s
        %s
        
        RUN START %s 
        RUN END %s
        %s EVENTS IN RUN

        TOF CORRECTED TIME CUT:
        %d - %d ns

            ''' % (self._run, Injector.tostr(self._injector), self._diffuser.upper(), self._plot_name, str(self._run_start_tree), str(self._run_end_tree), str(self._no_events), self._time_markers[0], self._time_markers[1]),
            fontsize=5, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
        if correct:
            ax.text(0.077, 0.76,
            '''
            CORRECTIONS APPLIED: %s
            
            
            ''' % (self._corrections),
                fontsize=4, horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)

        if fit and self._diffuser is not 'diffuser':

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

            ax.text(0.90, 0.97, 
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

        else:

            ax.text(0.90, 0.97, 
                u'''
                INJECTOR COORDS
                [%.2f, %.2f, %.2f] m

                TARGET CENTROID COORDS
                [%.2f, %.2f, %.2f] m

                ''' % (self._injector.Pos.X*cm, self._injector.Pos.Y*cm, self._injector.Pos.Z*cm, 
                        self._injector.Tar.X*cm, self._injector.Tar.Y*cm, self._injector.Tar.Z*cm),
                fontsize=4, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)

    def _setup_pyplot(self, invert, usetex=False):
        """Sets up rcParams for the event display. 

        Note:
            Turning on usetex won't work due to the way the text is currently written out in `_add_text`.

        Args:
            invert (bool): Inverts the colours of the plot (from dark mode to light mode).
            usetex (bool, optional): Uses LaTeX for displaying text, defaults to False. 

        Returns:
            None

        """
        if usetex:
            mpl.rcParams['text.usetex'] = True
            mpl.rcParams['font.family'] = 'serif' 
        else:
            mpl.rcParams['font.monospace'] = "DIN Alternate"
            mpl.rcParams['font.family'] = "monospace"

        mpl.rcParams['font.size'] = 10
        mpl.rcParams['axes.labelsize'] = 5
        mpl.rcParams['axes.linewidth'] = 0.1

        mpl.rcParams['legend.frameon'] = False
        mpl.rcParams['legend.fontsize'] = 4
        mpl.rcParams['xtick.minor.width'] = 0.1
        mpl.rcParams['ytick.minor.width'] = 0.1        
        mpl.rcParams['xtick.major.width'] = 0.1
        mpl.rcParams['ytick.major.width'] = 0.1

        mpl.rcParams['xtick.major.size'] = 1.0
        mpl.rcParams['xtick.minor.size'] = 0.3        
        mpl.rcParams['ytick.major.size'] = 1.0
        mpl.rcParams['ytick.minor.size'] = 0.3             

        mpl.rcParams['xtick.major.pad'] = 1.0
        mpl.rcParams['xtick.minor.pad'] = 1.0       
        mpl.rcParams['ytick.major.pad'] = 1.0
        mpl.rcParams['ytick.minor.pad'] = 1.0    

        mpl.rcParams['xtick.labelsize'] = 5
        mpl.rcParams['ytick.labelsize'] = 5

        if not invert:

            mpl.rcParams['savefig.facecolor'] = 'k'
            mpl.rcParams['savefig.edgecolor'] = 'k'
            mpl.rcParams['text.color'] = 'w'
            mpl.rcParams['patch.edgecolor'] = 'w'
            mpl.rcParams['axes.edgecolor'] = 'w'
            mpl.rcParams['axes.labelcolor'] = 'w'
            mpl.rcParams['xtick.color'] = 'w'
            mpl.rcParams['ytick.color'] = 'w'

        else:
            mpl.rcParams['savefig.facecolor'] = 'w'
            mpl.rcParams['savefig.edgecolor'] = 'w'
            mpl.rcParams['text.color'] = 'k'
            mpl.rcParams['patch.edgecolor'] = 'k'
            mpl.rcParams['axes.edgecolor'] = 'k'
            mpl.rcParams['axes.labelcolor'] = 'k'
            mpl.rcParams['xtick.color'] = 'w'
            mpl.rcParams['ytick.color'] = 'w'

    def _rotate_detector(self, data, rot):
        """Rotates the detector around the z-axis

        Note:
            This doesn't work it just fucks up the plot if rot is not 0.0.

        Args:
            data (pandas.DataFrame): The DataFrame of PMT information.
            rot (float): The rotation angle in degrees. 

        Returns:
            pandas.DataFrame: The new rotated PMT information DataFrame.

        """
        pmtx = data['pmtx'].values
        pmty = data['pmty'].values

        data['pmtx'] = pmtx*np.cos(rot) - pmty*np.sin(rot)
        data['pmty'] = pmtx*np.sin(rot) + pmty*np.cos(rot)

        return data

    def _segment_detector(self, data):
        """Adds a few useful columns to the PMT information DataFrame (cylindrical coordinates, detector regions...).

        Args:
            data (pandas.DataFrame): The DataFrame of PMT information.

        Returns:
            pandas.DataFrame: The new PMT information DataFrame with additional columns.

        """
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
        """Draws a frame around the detector regions.

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            draw_frame (bool): Whether to draw the detector frame.
        Returns:
            pandas.DataFrame: The new PMT information DataFrame with additional columns.

        """

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

        ax.axis('equal')
        ax.axis('off')

        return collection, (barrel_c, top_c, bottom_c)

    def _draw_barrel_hits(self, ax, data, det_geom, cmap):
        """Draws the hits on the barrel of the detector.

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            data (pandas.DataFrame): The DataFrame of PMT information.
            det_geom (tuple): Some storage of geometry variables for the detector plot.
            cmap (matplotlib.cmap): The cmap object to be applied to the value on each PMT.

        Returns:
            None

        """
        data_barrel = data.query('pmt_det_region == \'barrel\'')

        hits = []

        for hit in data_barrel.itertuples(index=True, name='Pandas'):

            s_pos = getattr(hit, 'pmts')
            z_pos = getattr(hit, 'pmtz')
            z_val = getattr(hit, 'val')

            pmt_r = sk_constants.IDPMTRadius
            pmt_patch = mpl.patches.Circle((s_pos, z_pos), pmt_r, fill=True, linewidth=None, facecolor=cmap(self._cmap_norm(z_val)))

            hits.append(pmt_patch)

        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)

        ax.add_collection(hit_patches)

    def _draw_top_hits(self, ax, data, det_geom, cmap):
        """Draws the hits on the top cap of the detector.

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            data (pandas.DataFrame): The DataFrame of PMT information.
            det_geom (tuple): Some storage of geometry variables for the detector plot.
            cmap (matplotlib.cmap): The cmap object to be applied to the value on each PMT.

        Returns:
            None

        """
        data_top = data.query('pmt_det_region == \'top\'')

        pmt_r = sk_constants.IDPMTRadius
        hits = []

        for hit in data_top.itertuples(index=True, name='Pandas'):

            x_pos = getattr(hit, 'pmtx')
            y_pos = getattr(hit, 'pmty')
            z_val = getattr(hit, 'val')
            pmt_patch = mpl.patches.Circle((x_pos + det_geom[1][0], y_pos + det_geom[1][1]), pmt_r, fill=True, facecolor=cmap(self._cmap_norm(z_val)), linewidth=None)

            hits.append(pmt_patch)

        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)

        ax.add_collection(hit_patches)

    def _draw_bottom_hits(self, ax, data, det_geom, cmap):
        """Draws the hits on the bottom cap of the detector.

        Args:
            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            data (pandas.DataFrame): The DataFrame of PMT information.
            det_geom (tuple): Some storage of geometry variables for the detector plot.
            cmap (matplotlib.cmap): The cmap object to be applied to the value on each PMT.

        Returns:
            None

        """
        data_bottom = data.query('pmt_det_region == \'bottom\'')
        pmt_r = sk_constants.IDPMTRadius
        hits = []

        for hit in data_bottom.itertuples(index=True, name='Pandas'):

            x_pos = getattr(hit, 'pmtx')
            y_pos = getattr(hit, 'pmty')
            z_val = getattr(hit, 'val')

            pmt_patch = mpl.patches.Circle((x_pos + det_geom[2][0], y_pos + det_geom[2][1]), pmt_r, fill=True, facecolor=cmap(self._cmap_norm(z_val)), linewidth=None)

            hits.append(pmt_patch)

        hit_patches = mpl.collections.PatchCollection(hits, match_original=True)

        ax.add_collection(hit_patches)

    def _get_cmap(self, cmap_name, invert):
        """Loads the colourmap to be applied to the value at each PMT.

        Args:
            cmap_name (str) The name of the cmap to load.
            invert (bool) Whether to reverse the cmap colours.

        Returns:
            mpl.colors.Colormap: The colourmap to be applied.

        """
        if cmap_name is not 'kindlmann' and cmap_name is not 'kindlmannext':

            cmap = mpl.cm.get_cmap(cmap_name)

        elif cmap_name is 'kindlmann':

            cmap_csv = pd.read_csv('%s/tables/kindlmann-table-float-1024.csv' % os.path.dirname(__file__))
            cmap = mpl.colors.ListedColormap(cmap_csv[['RGB_r', 'RGB_g', 'RGB_b']].values)

        elif cmap_name is 'kindlmannext':

            cmap_csv = pd.read_csv('%s/tables/extended-kindlmann-table-float-1024.csv' % os.path.dirname(__file__))
            cmap = mpl.colors.ListedColormap(cmap_csv[['RGB_r', 'RGB_g', 'RGB_b']].values)

        if invert:
            cmap = cmap.reversed()

        return cmap

    def _draw_hits(self, ax, data, det_geom, cmap_name, invert, logz, norm_det=True, mask_injector=True):
        """Draws the detector hits.

        Args:

            ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axis being used to draw the EventDisplay.
            data (pandas.DataFrame): The DataFrame of PMT information.
            det_geom (tuple): Some storage of geometry variables for the detector plot.
            cmap_name (str) The name of the cmap to load.
            invert (bool) Whether to reverse the cmap colours.
            logz (bool): Whether to make the colourbar scale logarithmic.
            norm_det (bool, optional): If True the detector is normalised as a whole, if False each region is normalised individually, defaults to True.
            mask_injector (bool, optional): If True PMTs immediately around the injector position are removed to avoid saturation, defaults to True
        Returns:
            None

        """
        if mask_injector:

            mask_r = 0.707*3

            injector_pos = self._injector.Pos

            injector_pos_s = np.arctan2(injector_pos.X*cm, injector_pos.Y*cm)*((injector_pos.X*cm)**2 + (injector_pos.Y*cm)**2)**0.5

            barrel_mask = data.query('pmt_det_region == \'barrel\' and (%f < pmts < %f) and (%f < pmtz < %f)' % (injector_pos_s - mask_r, injector_pos_s + mask_r, injector_pos.Z*cm - mask_r, injector_pos.Z*cm + mask_r)).index

            data.iloc[barrel_mask, data.columns.get_loc("val")] = np.nan
            data = data.dropna()

        #data['val'] = np.log10(data['val'] + 1e-6)

        if norm_det:

            zmin = data['val'].values.min()
            zmax = data['val'].values.max()

            self._zlims = zmin, zmax

        else:
            z_top = data.query('pmt_det_region == \'top\'').val.values
            z_bottom = data.query('pmt_det_region == \'bottom\'').val.values
            z_barrel = data.query('pmt_det_region == \'barrel\'').val.values

            data.loc[data.pmt_det_region == 'top', 'val'] = (z_top - z_top.min()) /(z_top.max() - z_top.min())
            data.loc[data.pmt_det_region == 'bottom', 'val'] = (z_bottom - z_bottom.min()) /(z_bottom.max() - z_bottom.min())
            data.loc[data.pmt_det_region == 'barrel', 'val'] = (z_barrel - z_barrel.min()) /(z_barrel.max() - z_barrel.min())

        cmap = self._get_cmap(cmap_name, invert)


        if not logz:
            self._cmap_norm = mpl.colors.Normalize(vmin=self._zlims[0],vmax=self._zlims[1], clip=True)
        else:

            data.iloc[data.query('val < 1e-6').index, data.columns.get_loc("val")] = np.nan
            data = data.dropna()

            zmin = data['val'].values.min()
            zmax = data['val'].values.max()

            self._zlims = zmin, zmax
            self._cmap_norm = mpl.colors.LogNorm(vmin=1e-4,vmax=self._zlims[1], clip=True)

        self._draw_barrel_hits(ax, data, det_geom, cmap)
        self._draw_top_hits(ax, data, det_geom, cmap)
        self._draw_bottom_hits(ax, data, det_geom, cmap)
