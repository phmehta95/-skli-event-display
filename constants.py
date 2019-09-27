import os
import numpy as np
from collections import namedtuple
from argparse import ArgumentParser

m = 1.
cm = 1e-2
mm = 1e-3

Pos = namedtuple("Pos", ["X", "Y", "Z"])
Inj = namedtuple("Inj", ["Pos", "Tar"])

class RunPeriod:

    TOP_DEPLOYMENT = 1
    KOREAN_LASER = 2
    LIVERPOOL_LASER = 3

    _enums = {'top deployment': TOP_DEPLOYMENT,
                'korean laser': KOREAN_LASER,
                'liverpool laser': LIVERPOOL_LASER,
    }    

    _names = {TOP_DEPLOYMENT: 'top deployment',
                KOREAN_LASER: 'korean laser',
                LIVERPOOL_LASER: 'liverpool laser',
    }    

    _dates = {TOP_DEPLOYMENT: 'Jan. \'18',
                KOREAN_LASER: 'Feb. \'19',
                LIVERPOOL_LASER: 'liverpool laser',
    }

    ALL = [TOP_DEPLOYMENT, KOREAN_LASER, LIVERPOOL_LASER]

    @classmethod
    def tostr(cls, code):
        return cls._names[code]

    @classmethod
    def datetostr(cls, code):
        return cls._dates[code]

class Source:

    BARE = 1
    COLLIMATOR = 2
    DIFFUSER = 3
    MONITOR = 4

    _names = {BARE:"barefibre",
              COLLIMATOR:"collimator",
              DIFFUSER:"diffuser",
              MONITOR:"monitor",
    }

    ALL = [BARE, COLLIMATOR, DIFFUSER, MONITOR]

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
    RunInfo = namedtuple("RunInfo", ["runperiod", "runnum", "injector", "source", "intensity", "quality", "monitor", "time_sel"])
    RI = RunInfo
    top_deployment_mon_id = 11161
    liverpool_mon_id = 11256

    Runs = {r.runnum:r for r in [
        # # Tuesday 23rd January 2018
        # RI(RunPeriod.TOP_DEPLOYMENT, 77480, Injector.OLD_TOP, Source.BARE, 5, False, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77481, Injector.OLD_TOP, Source.BARE, 0, False, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77483, Injector.OLD_TOP, Source.BARE, 6, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77484, Injector.OLD_TOP, Source.BARE, 7, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77485, Injector.OLD_TOP, Source.BARE, 5, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77486, Injector.OLD_TOP, Source.BARE, 8, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77488, Injector.OLD_TOP, Source.COLLIMATOR, 5, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77489, Injector.OLD_TOP, Source.COLLIMATOR, 6, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77490, Injector.OLD_TOP, Source.COLLIMATOR, 7, True, top_deployment_mon_id, None),
        # # Wednesday 24th January 2018
        # RI(RunPeriod.TOP_DEPLOYMENT, 77496, Injector.OLD_TOP, Source.COLLIMATOR, 8, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77497, Injector.OLD_TOP, Source.DIFFUSER, 10, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77498, Injector.OLD_TOP, Source.DIFFUSER, 15, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77499, Injector.OLD_TOP, Source.DIFFUSER, 20, True, top_deployment_mon_id, None),
        # RI(RunPeriod.TOP_DEPLOYMENT, 77500, Injector.OLD_TOP, Source.DIFFUSER, 25, True, top_deployment_mon_id, None),
        # # Tuesday 5th February 2019
        # RI(RunPeriod.KOREAN_LASER, 80174, Injector.B1, Source.COLLIMATOR, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80175, Injector.B2, Source.COLLIMATOR, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80176, Injector.B3, Source.COLLIMATOR, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80177, Injector.B4, Source.COLLIMATOR, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80178, Injector.B5, Source.COLLIMATOR, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80180, Injector.B1, Source.DIFFUSER, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80181, Injector.B2, Source.DIFFUSER, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80182, Injector.B3, Source.DIFFUSER, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80183, Injector.B4, Source.DIFFUSER, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80184, Injector.B5, Source.DIFFUSER, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80186, Injector.B1, Source.BARE, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80187, Injector.B2, Source.BARE, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80188, Injector.B3, Source.BARE, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80189, Injector.B4, Source.BARE, None, True, None, None),
        # RI(RunPeriod.KOREAN_LASER, 80190, Injector.B5, Source.BARE, None, True, None, None),
        # # July 2019
        RI(RunPeriod.LIVERPOOL_LASER, 81390, Injector.B1, Source.COLLIMATOR, None, True, liverpool_mon_id, (350, 700)),
        RI(RunPeriod.LIVERPOOL_LASER, 81391, Injector.B2, Source.COLLIMATOR, None, True, liverpool_mon_id, (350, 750)),
        RI(RunPeriod.LIVERPOOL_LASER, 81392, Injector.B3, Source.COLLIMATOR, None, True, liverpool_mon_id, (350, 700)),
        RI(RunPeriod.LIVERPOOL_LASER, 81393, Injector.B4, Source.COLLIMATOR, None, True, liverpool_mon_id, (350, 750)),
        RI(RunPeriod.LIVERPOOL_LASER, 81394, Injector.B5, Source.COLLIMATOR, None, True, liverpool_mon_id, (350, 500)),
        RI(RunPeriod.LIVERPOOL_LASER, 81397, Injector.B1, Source.DIFFUSER, None, True, liverpool_mon_id, (350, 600)),
        RI(RunPeriod.LIVERPOOL_LASER, 81398, Injector.B2, Source.DIFFUSER, None, True, liverpool_mon_id, (350, 575)),
        RI(RunPeriod.LIVERPOOL_LASER, 81399, Injector.B3, Source.DIFFUSER, None, True, liverpool_mon_id, (350, 575)),
        RI(RunPeriod.LIVERPOOL_LASER, 81400, Injector.B4, Source.DIFFUSER, None, True, liverpool_mon_id, (350, 550)),
        RI(RunPeriod.LIVERPOOL_LASER, 81401, Injector.B5, Source.DIFFUSER, None, True, liverpool_mon_id, (350, 600)),
        RI(RunPeriod.LIVERPOOL_LASER, 81395, Injector.B5, Source.BARE, None, True, liverpool_mon_id, (250, 500)),
        RI(RunPeriod.LIVERPOOL_LASER, 81402, Injector.B1, Source.BARE, None, True, liverpool_mon_id, (325, 550)),
        RI(RunPeriod.LIVERPOOL_LASER, 81403, Injector.B2, Source.BARE, None, True, liverpool_mon_id, (350, 550)),
        RI(RunPeriod.LIVERPOOL_LASER, 81404, Injector.B3, Source.BARE, None, True, liverpool_mon_id, (330, 530)),
        RI(RunPeriod.LIVERPOOL_LASER, 81405, Injector.B4, Source.BARE, None, True, liverpool_mon_id, (600, 800)),
        RI(RunPeriod.LIVERPOOL_LASER, 81409, Injector.B1, Source.MONITOR, None, True, liverpool_mon_id, (360, 400)),
        # Varied intensity runs
        #RI(81410, Injector.B1, Source.COLLIMATOR, None, True, True),
        #RI(81411, Injector.B3, Source.DIFFUSER, None, True, True),
    ]}

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