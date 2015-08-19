#!/usr/bin/env python
# Original filename: showRerun.py
#
# Author: Steve Bickerton
# Email: 
# Date: Wed 2013-08-21 14:19:46
# 
# Summary: 
# 

import sys
import os
import re
import math
import argparse
import numpy
import pyfits

import hsc.pipe.base.butler          as hscButler

import lsst.afw.coord                as afwCoord
import hsc.tools.bick.utils as util
import hsc.tools.bick.meeus as meeus

class Visit(object):

    def __init__(self, visit, coord='cel'):
        self.visit = visit
        self.coord = coord
        
    def addCalexp(self, calexp):
        pixelScale = 0.168  # arcsec/pix on HSC
        self.exptime = calexp.get("EXPTIME")
        self.filter  = calexp.get("FILTER")
        self.obj     = calexp.get("OBJECT")
        self.dateobs = calexp.get("DATE-OBS")
        self.boresight = afwCoord.Fk5Coord(calexp.get("RA2000"), calexp.get("DEC2000"))
        self.hst     = calexp.get("HST")
        self.seeing  = float(calexp.get("SEEING_MODE"))*pixelScale
        self.rot     = 0.5*(calexp.get("INR-STR") + calexp.get("INR-END"))
        yy, mm, dd = map(int, self.dateobs.split("-"))
        jd = meeus.calendar2JD(yy, mm, dd)
        self.moon_phase = meeus.moonIllumFrac(jd, isJD=True)
        

    @staticmethod
    def head(coord='cel'):
        ra, dec = "RA", "Dec"
        if coord == 'ecl':
            ra, dec = "Alpha", "Beta"
        if coord == 'gal':
            ra,dec = "l", "b"
        heads = ("Visit", "F", "Date", "HST", "Target", ra, dec, "ExpTim", "FWHM", "Moon", "Rot")
        headStr = "%-7s %-2s  %-10s %-8s  %-16s  %-10s %-9s  %-6s  %-4s  %-5s  %-4s" % (heads)
        return headStr
    
    def __str__(self):
        bs = self.boresight.toFk5()
        if self.coord == 'ecl':
            bs = self.boresight.toEcliptic()
        if self.coord == 'gal':
            bs = self.boresight.toGalactic()
        ra, dec = bs.getLongitude(), bs.getLatitude()
        values = (self.visit, self.filter, self.dateobs, self.hst[0:8],
                  self.obj, ra.asDegrees(), dec.asDegrees(),
                  self.exptime, self.seeing, self.moon_phase, self.rot)
        return "%7d %-2s  %10s %8s  %-16s  %10.6f %9.6f  %6s  %4.2f  %5.3f  %4.0f" % values

    
#############################################################
#
# Main body of code
#
#############################################################

def main(rerun, root=None, nshow=0, alt=1, coord='cel'):

    butler = util.getButler(rerun, root=root)
    dataIds = butler.queryMetadata("raw", ["visit"], format=["visit"]) #, dataId={'field': rerun})
    visits = []

    dataIds = dataIds[::alt]

    i_show = 0
    outputType = "None"
    for dataIdTmp in dataIds:

        visit= dataIdTmp
        ccd = 1
        dataId = {'visit':visit, "ccd": ccd}
        try:
            dataRef = hscButler.getDataRef(butler, dataId)
        except Exception, e:
            #print "getDataRef() failed for %d %d" % (visit, ccd), str(e)
            continue

        haveSomething = False
        srcfile = dataRef.get("src_filename", immediate=True)[0]
        postfile = dataRef.get("postISRCCD_filename", immediate=True)[0]

        # hopefully have a calexp
        if os.path.exists(srcfile):
            haveSomething = True
            outputType = "calexp"
            md = dataRef.get("calexp_md", immediate=True)

        # might find a postISRCCD
        elif os.path.exists(postfile):
            haveSomething = True
            outputType = "postISRCCD"
            md = dataRef.get('postISRCCD_md', immediate=True)
            
        else:
            #print "both " + srcfile + " and " + postfile + " missing.  continuing."
            continue
        
        
        v = Visit(visit, coord=coord)
        v.addCalexp(md)
        visits.append(v)

        i_show += 1
        if nshow and i_show >= nshow:
            break
        

    print "\nOutput: ", outputType
    print Visit.head(coord=coord)
    for v in visits:
        print v
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", action='version', version="0.0")
    parser.add_argument("rerun", help="Rerun to check")
    parser.add_argument('-r', '--root',   default=None, help="Root directory (overrides DATA_PATH)")
    parser.add_argument("-a", "--alt", default=1, type=int, help="Alternate ... show every nth frame")
    parser.add_argument("-c", "--coord", default='cel', choices=('cel', 'ecl', 'gal'), help="Coord system")
    parser.add_argument("-n", "--nshow", default=0, type=int, help="Number to show")
    args = parser.parse_args()

    main(args.rerun, root=args.root, nshow=args.nshow, alt=args.alt, coord=args.coord)
