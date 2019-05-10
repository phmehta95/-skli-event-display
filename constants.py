import os
import numpy as np
from argparse import ArgumentParser
from collections import namedtuple

m = 1.
cm = 1e-2
mm = 1e-3

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