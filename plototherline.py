#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import math
import argparse
import plotline_fctns as plf
from pyproj import Proj
from scipy.interpolate import griddata
from netCDF4 import Dataset

def round_to_nths(num, n):
    return int(num*n)/n

def add_to_legend(legns, legvars, legcolors, legn, legendname, varname, psxyspec):
    legns.append(legn)
    legvars.append(varname)
    legcolors.append(psxyspec)
    
    return legns, legvars, legcolors

def plot_legend_var(legn, legendname, varname, psxyspec):
    legenddat = pd.DataFrame({'x':[1,2], 'y':[1+legn,2+legn]})
    legenddat.to_csv('tempinterp1.txt',header=False,index=False,sep=' ')
    os.system("gmt psxy tempinterp1.txt -R -J %s -O -K >> %s" %(psxyspec, legendname))
    os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f9p,Helvetica,black >> %s'%(5, 1+legn, varname, legendname))

    print ('adding to legend ....')
    print("gmt psxy tempinterp.txt -R -J %s -O -K >> %s" %(psxyspec, legendname))
    print('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f9p,Helvetica,black >> %s'%(5, 1+legn, varname, legendname))

def plot_legend(legns, legvars, legcolors, legendname):
    legmax = np.max(legns) + 5
    os.system("gmt psbasemap -JX10 -R0.1/11.1/0.1/%s -B100 -P -K --MAP_GRID_PEN_PRIMARY=black > %s" % (legmax,legendname))
    for i in range(len(legns)):
        plot_legend_var(legns[i], legendname, legvars[i], legcolors[i])

def main(args):
    
    os.system("gmt set FONT 10p,Helvetica")
    if args.vertexg is not None:
        VE = args.vertexg
    else:
        VE = 2

    MorC = args.MorC

    myProj = Proj("+proj=utm +zone=10T, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    if args.utmzone is not None:
        myProj = Proj("+proj=utm +zone=%sT, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs" % args.utmzone)

    if (args.tapename is not None) and (args.tapedir is not None):
        tapedir = args.tapedir
        tapename = args.tapename
        nameandtail = tapename.split('.')
        linename = nameandtail[0]
        navname = '%s.nav' % linename
        mapname = '%s_map.ps' % linename
        legendname = 'Output/%s_legend.ps' % linename
    else:
        print ('must enter -p flag or -d and -n flags. See help desc')
        print ('exiting ... ')
        exit()

    print ('tapedir , tapename:',tapedir,',',tapename)
    navdir = 'Output/navigationfiles'
    os.system("mkdir %s"%navdir)
    
    if args.usetracr is not None:
        if args.usetracr == 'y':
            usetracr = True
        else:
            usetracr = False
    else:
        usetracr = False

    if args.justcontours is not None:
        if args.justcontours == 'y':
            justcontours = True
        else:
            justcontours = False
    else:
        justcontours = False

    if args.justcurve is not None:
        if args.justcurve == 'y':
            justcurve = True
        else:
            justcurve = False
    else:
        justcurve = False

    if MorC == 'C' or usetracr:
        if args.sx is not None:
            try:
                try:
                    print ('running plf.getNavTracrWBDelrtFromHeaderSX')
                    navline = plf.getNavTracrWBDelrtFromHeaderSX(tapedir, navdir, tapename, navname)
                except:
                    print ('running plf.getNavTracrWBFromHeaderSX')
                    navline = plf.getNavTracrWBFromHeaderSX(tapedir, navdir, tapename, navname)
                print (navline)
            except:
                print ('running plf.getNavTracrFromHeaderSX')
                navline = plf.getNavTracrFromHeaderSX(tapedir, navdir, tapename, navname)
        else:
            try:
                print ('running plf.getNavTracrWBFromHeader')
                navline = plf.getNavTracrWBFromHeader(tapedir, navdir, tapename, navname)
            except:
                print ('running plf.getNavTracrFromHeader')
                navline = plf.getNavTracrFromHeader(tapedir, navdir, tapename, navname)
    else:
        print (tapedir, navdir, tapename, navname)
        if args.sx is not None:
            try:
                navline = plf.getNavCdpWBFromHeaderSX(tapedir, navdir, tapename, navname)
            except:
                navline = plf.getNavCdpFromHeaderSX(tapedir, navdir, tapename, navname)
        else:
            try:
                navline = plf.getNavCdpWBFromHeader(tapedir, navdir, tapename, navname)
            except:
                navline = plf.getNavCdpFromHeader(tapedir, navdir, tapename, navname)


    navline['lon'],navline['lat'] = myProj(navline['gx'].values,
                                        navline['gy'].values, inverse=True)

    if args.arcseconds is not None:
        if args.arcseconds == 'y':
            usearcseconds = True
        else:
            usearcseconds = False
    else:
        usearcseconds = False

    if navline['lon'].max() > 360 or navline['lon'].min() < -360 or usearcseconds:
        print ('arcseconds ')
        navline['lon'] = navline['gx'].values / 60 / 60 / 10
        navline['lat'] = navline['gy'].values / 60 / 60 / 10
        if navline['lon'].max() > 360 or navline['lon'].min() < -360:
            print ('arcseconds again')
            navline['lon'] = navline['lon'].values / 10
            navline['lat'] = navline['lat'].values / 10
        if navline['lon'].max() > 360 or navline['lon'].min() < -360:
            print ('arcseconds again again')
            navline['lon'] = navline['lon'].values / 10
            navline['lat'] = navline['lat'].values / 10

    #if MorC == 'C':
    #    navline['lon'] = navline['gx'].values/60/10000
    #    navline['lat'] = navline['gy'].values/60/10000

    print ('meanx, meany', navline['gx'].mean(), navline['gy'].mean())
    print ('meanlon, meanlat', navline['lon'].mean(), navline['lat'].mean())

    ns, dt = plf.get_ns_dt(tapedir,tapename)
    print ('ns, dt', ns, dt, ns*dt)
    '''
    lon1 = navline['lon'].values[0]
    lat1 = navline['lat'].values[0]
    lon2 = navline['lon'].values[-1]
    lat2 = navline['lat'].values[-1]
    maxoffset, az = plf.cosine(lon1, lat1, lon2, lat2)
    lonslats1 = np.zeros((len(navline)-1, 3))
    lonslats2 = np.zeros((len(navline)-1, 3))
    lonslats1[:,0] = navline['lon'].values[1:]
    lonslats1[:,1] = navline['lat'].values[1:]
    lonslats2[:,0] = navline['lon'].values[:-1]
    lonslats2[:,1] = navline['lat'].values[:-1]
    print ('looping through ',len(navline),'lines ... ')
    for i in range(len(lonslats1)):
        lon1 = lonslats1[i,0]
        lat1 = lonslats1[i,1]
        lon2 = lonslats2[i,0]
        lat2 = lonslats2[i,1]
        lonslats1[i,2], az = plf.cosine(lon1, lat1, lon2, lat2)

    middist = np.mean(lonslats1[:,2])
    maxoffset = middist*float(len(set(list(navline['cdp'].values))))
    maxoffset *= 1000
    print (middist, maxoffset, float(len(set(list(navline['cdp'].values)))))
    '''

    lon1 = navline['lon'].values[0]
    lat1 = navline['lat'].values[0]
    lon2 = navline['lon'].values[-1]
    lat2 = navline['lat'].values[-1]
    print (lon1, lat1, lon2, lat2)
    maxoffset, az = plf.cosine(lon1, lat1, lon2, lat2)
    maxoffset *= 1000
    offsets = np.linspace(0,maxoffset,num=len(navline))
    navline['offset'] = offsets

    wOmin = 0
    wOmax = maxoffset
    xscal=1
    centlon = navline['lon'].values[int(len(navline)/2)]
    centlat = navline['lat'].values[int(len(navline)/2)]
    centcdp = navline['cdp'].values[int(len(navline)/2)]

    print (navline)
    print (centlon, centlat, centcdp)
    if args.window is not None:
        if args.intrange is not None:
            print (' can only use -w OR -i option OR -it')
            print ('exiting ... ')
            exit()
        if args.intrange_t is not None:
            print (' can only use -w OR -i option OR -it')
            print ('exiting ... ')
            exit()
        wCmin = args.window[0]
        wCmax = args.window[1]
        wTmin = args.window[2]
        wTmax = args.window[3]
        windnav = navline[(navline.cdp > wCmin) & (navline.cdp < wCmax)]
        wOmax = windnav['offset'].max()
        wOmin = windnav['offset'].min()
        centcdp = int((wCmax+wCmin)/2)
        nlon = navline[navline.cdp == centcdp]['lon'].values[0]
        nlat = navline[navline.cdp == centcdp]['lat'].values[0]
    elif args.intrange is not None:
        if args.window is not None:
            print (' can only use -w OR -i option OR -it')
            print ('exiting ... ')
            exit()
        if args.intrange_t is not None:
            print (' can only use -w OR -i option OR -it')
            print ('exiting ... ')
            exit()
        nlon = args.intrange[0]
        nlat = args.intrange[1]
        bego = args.intrange[2]
        endo = args.intrange[3]
        wTmin = args.intrange[4]
        wTmax = args.intrange[5]
        if navline['lon'].mean() < 0:
            if nlon > 180:
                nlon -= 360
        if navline['lon'].mean() > 0:
            if nlon < 0:
                nlon += 360
        diflons = np.abs(navline['lon'].values - nlon)
        diflats = np.abs(navline['lat'].values - nlat)
        alldiff = (diflons*diflons + diflats*diflats)
        navline['diffs'] = alldiff
        mindiff = navline['diffs'].min()
        thispt = navline[navline.diffs == mindiff]
        centcdp = thispt['cdp'].values[0]
        centlon = thispt['lon'].values[0]
        centlat = thispt['lat'].values[0]
        centoff = thispt['offset'].values[0]
        thischunk = navline[(navline.offset > centoff+bego) &
                            (navline.offset < centoff+endo)]
        minoffset = thischunk['offset'].min()
        maxoffset = thischunk['offset'].max()
        wOmin = minoffset-centoff
        wOmax = maxoffset-centoff
        wCmin = thischunk['cdp'].min()
        wCmax = thischunk['cdp'].max()
        odiff = endo-bego
        ndiff = wOmax-wOmin
        if ndiff < odiff:
            xscal = ndiff/odiff
    elif args.intrange_t is not None:
        wCmin = np.float64(navline['cdp'].min())
        wCmax = np.float64(navline['cdp'].max())
        wTmin = args.intrange_t[0]
        wTmax = args.intrange_t[1]
        centcdp = int((wCmax+wCmin)/2)
        nlon = navline[navline.cdp == centcdp]['lon'].values[0]
        nlat = navline[navline.cdp == centcdp]['lat'].values[0]
    else:
        wCmin = np.float64(navline['cdp'].min())
        wCmax = np.float64(navline['cdp'].max())
        wTmin = 0
        wTmax = float(ns)*dt
        centcdp = int((wCmax+wCmin)/2)
        nlon = navline[navline.cdp == centcdp]['lon'].values[0]
        nlat = navline[navline.cdp == centcdp]['lat'].values[0]

    if 'delrt' in navline.columns:
        delrtmax = navline['delrt'].max() * 1e-3
        delrtmin = navline['delrt'].min() * 1e-3
        wTmin = wTmin + delrtmin
        wTmax = wTmax + delrtmin

    if args.cv is not None:
        cvelocity = args.cv
    else:
        cvelocity = 1500
    
    wDmin = wTmin * cvelocity / 2
    wDmax = wTmax * cvelocity / 2

    xintL = 100
    xintM = 50
    tintL = 1
    tintM = 0.5
    if wOmax - wOmin > 500:
        xintL = 200
        xintM = 100
    if wOmax - wOmin > 2000:
        xintL = 500
        xintM = 250
    if wOmax - wOmin > 3000:
        xintL = 1000
        xintM = 500
    if wOmax - wOmin >5000:
        xintL = 2000
        xintM = 1000
    if wOmax - wOmin >10000:
        xintL = 4000
        xintM = 2000
    if wOmax - wOmin >20000:
        xintL = 10000
        xintM = 10000
    if wOmax - wOmin >80000:
        xintL = 20000
        xintM = 20000
    if wTmax - wTmin < 2.0:
        tintL = 0.4
        tintM = 0.2
    if wTmax - wTmin < 1.1:
        tintL = 0.2
        tintM = 0.1
    if wTmax - wTmin < 1.1:
        tintL = 0.2
        tintM = 0.1
    if wTmax - wTmin < 0.5:
        tintL = 0.1
        tintM = 0.05
    if wTmax - wTmin < 0.2:
        tintL = 0.04
        tintM = 0.02
    if wTmax - wTmin < 0.1:
        tintL = 0.02
        tintM = 0.01
    
    ointL = xintL
    ointM = xintM
    dintL = 400*tintL
    dintM = 400*tintM

    if VE < 4:
        dintL = 1000*tintL
        dintM = 1000*tintM

    xdist = (wOmax-wOmin)

    ydist = cvelocity*(wTmax-wTmin)/2

    xyratio = xdist/ydist
    xplotmax=15*xscal
    yplotmax=-1*ydist*xplotmax/xdist*VE
    shiftmap = -1*yplotmax+4
    title = '%s'%(tapename)
    Rl = '-R%f/%f/%f/%f' % (wCmin,wCmax,wTmin,wTmax)
    Rl2 = '-R%f/%f/%f/%f' % (wOmin,wOmax,wDmin,wDmax)
    Jl = '-JX%s/%s' % (xplotmax,yplotmax)
    if args.skew is not None:
        Pl = '-py140/60'
    else:
        Pl = ' '

    print ('PLLLLLLLLLLL!!!!!!!!!',Pl)

    cdpint = int((navline['cdp'].max() - navline['cdp'].min()) / 5)
    offint = int((navline['offset'].max() - navline['offset'].min()) / 5)
    if args.bs is not None:
        Bl2 = '-B%fg%f:" ":/%fg%f:"~Depth(m)":%s:." ":' % (offint, offint,dintL,dintM, args.bs[:2])
        Bl = '-B%fg%f:" ":/%fg%f:"Time(s)":%s:." ":' % (cdpint, cdpint, tintL,tintM, args.bs[-2:])
    else:
        Bl2 = '-B%fg%f:" ":/%fg%f:"~Depth(m)":EN:." ":' % (offint, offint,dintL,dintM)
        Bl = '-B%fg%f:" ":/%fg%f:"Time(s)":WS:." ":' % (cdpint, cdpint, tintL,tintM)

    sf = 'Output/%s' % mapname
    if MorC == 'C':
        Sl = '-Sb'
    else:
        Sl = '-Sc'

    navline.to_csv('temp.txt',header=True,index=False)

    mXmin, mXmax = np.float64(navline['lon'].min()), np.float64(navline['lon'].max())
    mYmin, mYmax = np.float64(navline['lat'].min()), np.float64(navline['lat'].max())
    cmin, cmax = np.float64(navline['cdp'].min()), np.float64(navline['cdp'].max())
    if cmin > wCmin:
        wCmin += cmin
        wCmax += cmin

    firstx = navline['lon'].values[0]
    firsty = navline['lat'].values[0]
    lastx = navline['lon'].values[-1]
    lasty = navline['lat'].values[-1]
    if firstx < 0:
        firstx += 360
    if lastx < 0:
        lastx += 360
    londiff = lastx - firstx
    latdiff = lasty - firsty
    revdir = False
    if abs(londiff) > 3*abs(latdiff):
        '''
        pmXmin = navline['lon'].mean() - londiff/2.0 - 0.2
        pmXmax = navline['lon'].mean() + londiff/2.0 + 0.2
        pmYmin = navline['lat'].mean() - londiff/2.0 - 0.1
        pmYmax = navline['lat'].mean() + londiff/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.025
        pmXmax = navline['lon'].max() + 0.025
        pmYmin = navline['lat'].min() - 0.05
        pmYmax = navline['lat'].max() + 0.05
        rlabel = 'E'
        llabel = 'W'
        if londiff < 0:
            revdir = True
    elif abs(latdiff) > 3*abs(londiff):
        '''
        pmXmin = navline['lon'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmXmax = navline['lon'].mean() + max(latdiff, londiff)/2.0 + 0.1
        pmYmin = navline['lat'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmYmax = navline['lat'].mean() + max(latdiff, londiff)/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.075
        pmXmax = navline['lon'].max() + 0.075
        pmYmin = navline['lat'].min() - 0.025
        pmYmax = navline['lat'].max() + 0.025
        
        rlabel = 'N'
        llabel = 'S'
        if latdiff < 0:
            revdir = True
    elif londiff>0 and latdiff>0:
        '''
        pmXmin = navline['lon'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmXmax = navline['lon'].mean() + max(latdiff, londiff)/2.0 + 0.1
        pmYmin = navline['lat'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmYmax = navline['lat'].mean() + max(latdiff, londiff)/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.05
        pmXmax = navline['lon'].max() + 0.05
        pmYmin = navline['lat'].min() - 0.025
        pmYmax = navline['lat'].max() + 0.025
        rlabel = 'NE'
        llabel = 'SW'
    elif londiff<0 and latdiff<0:
        '''
        pmXmin = navline['lon'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmXmax = navline['lon'].mean() + max(latdiff, londiff)/2.0 + 0.1
        pmYmin = navline['lat'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmYmax = navline['lat'].mean() + max(latdiff, londiff)/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.05
        pmXmax = navline['lon'].max() + 0.05
        pmYmin = navline['lat'].min() - 0.025
        pmYmax = navline['lat'].max() + 0.025
        rlabel = 'NE'
        llabel = 'SW'
        revdir = True
    elif londiff<0 and latdiff>0:
        '''
        pmXmin = navline['lon'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmXmax = navline['lon'].mean() + max(latdiff, londiff)/2.0 + 0.1
        pmYmin = navline['lat'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmYmax = navline['lat'].mean() + max(latdiff, londiff)/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.05
        pmXmax = navline['lon'].max() + 0.05
        pmYmin = navline['lat'].min() - 0.025
        pmYmax = navline['lat'].max() + 0.025
        rlabel = 'NW'
        llabel = 'SE'
    elif londiff>0 and latdiff<0:
        '''
        pmXmin = navline['lon'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmXmax = navline['lon'].mean() + max(latdiff, londiff)/2.0 + 0.1
        pmYmin = navline['lat'].mean() - max(latdiff, londiff)/2.0 - 0.1
        pmYmax = navline['lat'].mean() + max(latdiff, londiff)/2.0 + 0.1
        '''
        pmXmin = navline['lon'].min() - 0.05
        pmXmax = navline['lon'].max() + 0.05
        pmYmin = navline['lat'].min() - 0.025
        pmYmax = navline['lat'].max() + 0.025
        rlabel = 'NW'
        llabel = 'SE'
        revdir = True
    else:
        pmXmin = mXmin
        pmXmax = mXmax
        pmYmin = mYmin
        pmYmax = mYmax
        print ('first and last x are nan, check nav in headers')
        print ('firstx: ',firstx)
        print ('lastx: ',lastx)
        print ('firsty: ',firsty)
        print ('lasty: ',lasty)
        rlabel = 'NW'
        llabel = 'SE'

    if revdir == True:
        xplotmax2 = -1*xplotmax
        Jl = '-JX%s/%s' % (xplotmax2,yplotmax)

    ogmXmin, ogmXmax, ogmYmin, ogmYmax = mXmin, mXmax, mYmin, mYmax

    print ('mXmin, mXmax, mYmin, mYmax',mXmin, mXmax, mYmin, mYmax)
    mXmin = round_to_nths(mXmin,100)
    mXmax = round_to_nths(mXmax,100)
    mYmin = round_to_nths(mYmin,100)
    mYmax = round_to_nths(mYmax,100)
    print ('mXmin, mXmax, mYmin, mYmax',mXmin, mXmax, mYmin, mYmax)


    morex = (pmXmax - pmXmin) / 3.0
    clipmXmin = pmXmin - morex
    clipmXmax = pmXmax + morex
    clipmYmin = pmYmin
    clipmYmax = pmYmax

    plotboxX = [clipmXmin, clipmXmax, clipmXmax, clipmXmin, clipmXmin]
    plotboxY = [clipmYmin, clipmYmin, clipmYmax, clipmYmax, clipmYmin]
    plotclip = pd.DataFrame({'x':plotboxX, 'y':plotboxY})
    plotclip.to_csv('plotclip.dat', header=False, index=False, sep=' ')
    mXmin = pmXmin - 3
    mXmax = pmXmax + 3
    mYmin = pmYmin - 1
    mYmax = pmYmax + 1
    print ('mXmin, mXmax, mYmin, mYmax',mXmin, mXmax, mYmin, mYmax)
    
    if args.bounds is not None:
        mXmin = args.bounds[0]
        mXmax = args.bounds[1]
        mYmin = args.bounds[2]
        mYmax = args.bounds[3]
    print ('mXmin, mXmax, mYmin, mYmax',mXmin, mXmax, mYmin, mYmax)

    UTMxmin, UTMymin = myProj(mXmin, mYmin, inverse=False)
    UTMxmax, UTMymax = myProj(mXmax, mYmax, inverse=False)
    Rm = '-R%f/%f/%f/%f' % (UTMxmin,UTMxmax,UTMymin,UTMymax+0.2)
    Rm = '-R%f/%f/%f/%f' % (mXmin,mXmax,mYmin,mYmax)
    Rm2 = '-R%f/%f/%f/%f' % (clipmXmin,clipmXmax,clipmYmin,clipmYmax)
    Jx = '-Jx1:160000'
    JXu = '-Ju10T/1:160000'
    mapwidth = 7
    Jm = '-JM%i' % mapwidth
    if mXmax - mXmin < 0.1 and mYmax - mYmin < 0.1:
        Bm = '-Bx0.04 -By0.02 -BNWse'
    elif mXmax - mXmin < 0.5 and mYmax - mYmin < 0.5:
        Bm = '-Bx0.2 -By0.1 -BNWse'
    elif mXmax - mXmin < 1 and mYmax - mYmin < 1:
        Bm = '-Bx0.4 -By0.2 -BNWse'
    elif mXmax - mXmin < 5 and mYmax - mYmin < 5:
        Bm = '-Bx1 -By0.5 -BNWse'
    elif mXmax - mXmin < 10 and mYmax - mYmin < 10:
        Bm = '-Bx2 -By1 -BNWse'
    elif mXmax - mXmin < 50 and mYmax - mYmin < 50:
        Bm = '-Bx10 -By5 -BNWse'
    else:
        Bm = '-Bx30 -By30 -BNWse'

    if clipmXmax - clipmXmin < 0.1 and clipmYmax - clipmYmin < 0.1:
        Bm2 = '-Bx0.04 -By0.02 -BNWse'
    elif clipmXmax - clipmXmin < 0.5 and clipmYmax - clipmYmin < 0.5:
        Bm2 = '-Bx0.2 -By0.1 -BNWse'
    elif clipmXmax - clipmXmin < 1 and clipmYmax - clipmYmin < 1:
        Bm2 = '-Bx0.4 -By0.2 -BNWse'
    elif clipmXmax - clipmXmin < 5 and clipmYmax - clipmYmin < 5:
        Bm2 = '-Bx1 -By0.5 -BNWse'
    elif clipmXmax - clipmXmin < 10 and clipmYmax - clipmYmin < 10:
        Bm2 = '-Bx1 -By0.5 -BNWse'
    elif clipmXmax - clipmXmin < 50 and clipmYmax - clipmYmin < 50:
        Bm2 = '-Bx10 -By5 -BNWse'
    else:
        Bm2 = '-Bx30 -By30 -BNWse'

    Xm = '-Y%s' % shiftmap
    Xm2 = '-X%s' % int(mapwidth + 2)
    DEM = 'library/puget_sound_13_navd88_2014.nc'
    DEM2 = '/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc'
    if args.dem is not None:
        DEM = args.dem

    gmcpt = 'library/gmtopo.cpt'
    PSSnav = 'PSSnav.txt'
    KNUnav = '/Users/ginevramoore/Documents/research/AACSE/Knudson/alllinenav.dat'
    KNUnavb = '/Users/ginevramoore/Documents/research/AACSE/Knudson/alllinenav_beg.dat'
    KNUnave = '/Users/ginevramoore/Documents/research/AACSE/Knudson/alllinenav_end.dat'
    LWBnav = 'LWBnav.txt'
    sproulnav = 'sproulnav.txt'
    PSSnavtxt = 'PSSnav_txt.txt'
    KNUnavtxt = '/Users/ginevramoore/Documents/research/AACSE/Knudson/alllinenav_lab.dat'
    LWBnavtxt = 'LWBnav_txt.txt'
    sproulnavtxt = 'sproulnav_txt.txt'
    upliftpts = 'library/uplift/tenbrink2006_uplift.txt'
    prattfiles = 'library/lines/allprattlines.txt'
    oligocene = '/Users/ginevramoore/Documents/research/segymanip/library/geology/olig.dat'
    miocene = '/Users/ginevramoore/Documents/research/segymanip/library/geology/miocene.dat'
    
    namelist = tapename.split('.')
    gridname = '%s.grd'%linename

    legns = []
    legvars = []
    legcolors = []
    legn = 0
    legshift = 3
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'blah', '-W0.5,blue')

    if 'delrt' in navline.columns:
        wTmin = wTmin - delrtmin
        wTmax = wTmax - delrtmin

    Rg = '-R%i/%i/%.2f/%.2f'%(cmin,cmax,wTmin,wTmax)
    Ig = '-I1/%f'%dt
    Gg = '-GOutput/%s'%(gridname)
    Tg = '%s/%s'%(tapedir, tapename)
    Mg = '-M%i'%(wCmax)

    if 'delrt' in navline.columns:
        wTmin = wTmin + delrtmin
        wTmax = wTmax + delrtmin

    if args.append is not None:
        linfile = '%s_lin.txt' % args.append[:-4]
        dotfile = '%s_dot.txt' % args.append[:-4]
        segfile = '%s_seg.txt' % args.append[:-4]
        labfile = '%s_lab.txt' % args.append[:-4]

        with open(dotfile, 'a') as file:
            file.write('%s %s %s \n' % (str(centlon), str(centlat), linename))

    if args.nc is not None:
        if args.nc == 'y':
            nocontours = True
        else:
            nocontours = False
    else:
        nocontours = False

    if args.thresh is not None:
        thresh = args.thresh
    else:
        thresh = 0.01
    if 'swdep' in navline.columns:
        navline['swdep'] = navline['swdep'].values * 1e-6
        if navline['lat'].max() > 50:
            if 'delrt' in navline.columns:
                try:
                    navline['swdep'] = navline['swdep'].values * 2 / cvelocity * 1e4
                    writeWB = navline[['cdp', 'swdep']]
                    writeWB.to_csv('wbplottempunfilt.txt', header=False, index=False, sep=' ')
                    iii = 0
                    while iii < 15:
                        iii += 1
                        mvarr = np.zeros((len(navline)-1, 7))
                        mvarr[:, 0] = navline['cdp'].values[1:]
                        mvarr[:, 1] = navline['cdp'].values[:-1]
                        mvarr[:, 2] = navline['swdep'].values[0:-1]
                        mvarr[:, 3] = navline['swdep'].values[1:]
                        mvarr[:, 4] = np.abs(mvarr[:, 2] - mvarr[:, 3])
                        mvarr[:, 5] = np.abs(mvarr[:, 2] - mvarr[:, 3])
                        mvarr[:, 5][mvarr[:, 5] > thresh] = -9999
                        mvarr[:, 5][mvarr[:, 2] <= 0.0001] = -9999
                        mvarr[:, 5][mvarr[:, 3] <= 0.0001] = -9999
                        mvarr[:, 6] = griddata(mvarr[:, 0][mvarr[:, 5] != -9999], mvarr[:, 3][mvarr[:, 5] != -9999], mvarr[:, 0], method='linear')
                        mvarr = mvarr[np.isfinite(mvarr[:, 6])]

                        temp = pd.DataFrame({'index1': mvarr[:, 0], 'index2':mvarr[:, 1],
                                             'swdep1': mvarr[:, 2], 'swdep2':mvarr[:, 3],
                                             'diff': mvarr[:, 4], 'diffmod':mvarr[:, 5],
                                             'swdepfilt': mvarr[:, 6]})

                        navline['swdep'] = griddata(mvarr[:, 0][mvarr[:, 6] != 0], mvarr[:, 6][mvarr[:, 6] != 0], navline['cdp'].values, method='nearest')

                    AKbathy = '/Users/ginevramoore/Documents/research/AACSE/bathymetry_grids/bathygrid_09.11.2019.txt'
                    AKbathyd = pd.read_csv(AKbathy, delim_whitespace=True, names=['x','y','z'])
                    AKbathyd = AKbathyd[(AKbathyd.x >= mXmin) & (AKbathyd.x <= mXmax) & \
                                        (AKbathyd.y >= mYmin) & (AKbathyd.y <= mYmax)]
                          
                    '''
                    multibeam = np.zeros((len(AKbathyd), 3))
                    multibeam[:, 0] = AKbathyd['x'].values
                    multibeam[:, 1] = AKbathyd['y'].values
                    multibeam[:, 2] = AKbathyd['z'].values
                    interpWB = np.zeros((len(navline), 3))
                    interpWB[:, 0] = navline['lon'].values
                    interpWB[:, 1] = navline['lat'].values
                    interpWB[:, 2] = griddata(multibeam[:, 0:2], multibeam[:, 2], interpWB[:, 0:2], method='linear')
                    navline['swdep'] = interpWB[:, 2] * -2.0 * cvelocity
                    '''
                except:
                    print ('delrt in columns but not AK dataset')

        writeWB = navline[['cdp', 'swdep']]
        writeWB.to_csv('wbplottemp.txt', header=False, index=False, sep=' ')
        writeWB.to_csv('Output/%s_WB.txt'%linename, header=False, index=False, sep=' ')
        plotWB = True
    else:
        plotWB = False

    if args.nwb is not None:
        if args.nwb == 'y':
            plotWB = False

    if args.tapedir is not None:
        if 'mute' in args.tapedir:
            plotWB = False

    if 'delrt' in navline.columns:
        navline['delrt'] = navline['delrt'].values * 1e-3
        writeWB = navline[['cdp', 'delrt']]
        writeWB.to_csv('delrtplottemp.txt', header=False, index=False, sep=' ')
        plotdelrt = True
    else:
        plotdelrt = False

    navline.to_csv('Output/%s_nav.txt'%linename,header=True,index=False,na_rep=np.nan)

    os.system("gmt segy2grd %s %s %s %s %s -V" % (Tg, Gg, Ig, Rg, Mg))

    if 'delrt' in navline.columns:
        if navline['lat'].max() > 50:
            if delrtmin > 0:
                nc_fid = Dataset('Output/%s'%(gridname), 'r')
                x = nc_fid.variables['x'][:]
                y = nc_fid.variables['y'][:]
                z = nc_fid.variables['z'][:]
                y += navline['delrt'].min()
                xx, yy = np.meshgrid(x, y)
                semb2data = pd.DataFrame({'x':xx.flatten(),'y':yy.flatten(),'z':z.flatten()})
                semb2data = semb2data[['x','y','z']]
                dx = x[1] - x[0]
                dy = y[1] - y[0]
                xmin = np.min(x)
                xmax = np.max(x)
                ymin = np.min(y)
                ymax = np.max(y)
                semb2GridFile = 'Output/%s'%(gridname)
                semb2TextFile = '%s_delrtgrid.txt'%(gridname[:-4])
                print ('writing file to shift')
                semb2data.to_csv(semb2TextFile, header=False, index=False, sep=' ', na_rep=0.0)
                gflag="-G%s" % (semb2GridFile)
                g2flag="%s" % (semb2GridFile)
                rflag="-R%s/%s/%s/%s" %(np.floor(xmin),np.ceil(xmax),np.floor(ymin),np.ceil(ymax))
                iflag="-I%s/%s" %(dx,dy)
                print ('converting to grid ... ')
                os.system("gmt xyz2grd %s %s %s %s" % (semb2TextFile,rflag,iflag,gflag))
                os.system("rm %s" % semb2TextFile)

    phc = False
    if args.phc is not None:
        if args.phc == 'y':
            phc = True

    os.system("gmt grdgradient %s -A90/60 -Ghillshade.grd"%(DEM))
    os.system("gmt grdgradient %s -A90/60 -Ghillshade2.grd"%(DEM2))
    os.system("gmt makecpt -Ctopo -D -Fr -T-300/300/0.01 -Z > dep.cpt")
    os.system("gmt makecpt -Cocean -D -Fr -T-1/1/0.01 -Z > hill.cpt")
    os.system("gmt makecpt -Cpolar -D -Fr -T-100000/100000/1000 -Z > seisJC.cpt" )
    if args.saturation is not None:
        if not phc:
            os.system("gmt makecpt -Cpolar -D -Fr -T-%f/%f/0.01 -Z > seis.cpt" % (args.saturation, args.saturation))
        else:
            os.system("gmt makecpt -Cgray -D -Fr -T-%f/%f/0.01 -Z > seis.cpt" % (args.saturation, args.saturation))

        os.system("gmt makecpt -Cgray -D -Fr -T0/%s/1 -I -Z > chirp.cpt" % args.saturation)
    else:
        if not phc:
            os.system("gmt makecpt -Cpolar -D -Fr -T-1/1/0.01 -Z > seis.cpt")
        else:
            os.system("gmt makecpt -Cgray -D -Fr -T-1/1/0.01 -Z > seis.cpt")
        os.system("gmt makecpt -Cgray -D -Fr -T0/200/1 -I -Z > chirp.cpt")

    trackline = navline[['lon','lat']]
    trackline.to_csv('track_4.xyg',header=False,index=False,sep=' ')
    os.system("gmt grdtrack track_4.xyg -GupliftE_0.001.grd > trackE_4.xygt")
    os.system("gmt grdtrack track_4.xyg -GupliftS_0.001.grd > trackS_4.xygt")

    trackdat = pd.read_csv('trackE_4.xygt', names=['lon','lat','elev'], delim_whitespace=True)
    trackdatS = pd.read_csv('trackS_4.xygt', names=['lon','lat','str'], delim_whitespace=True)

    try:
        trackdat['cdp'] = navline['cdp'].values
        trackdat['strike'] = trackdatS['str'].values
        deriv = np.gradient(trackdat['elev'].values)
        trackdat['grad'] = deriv
        gtS = 'up'
        n = 0
        sides = []
        for index,row in trackdat.iterrows():
            if row['grad'] > 0:
                gtI = 'up'
            elif row['grad'] < 0:
                gtI = 'down'
            else:
                gtI = gtS
            if gtI != gtS:
                n += 1
                gtS = gtI
            sides.append(n)

        trackdat['side'] = sides
    except:
        print ('line not within uplift contour grid!')
        print ('setting all values to 0.0')
        trackdat = navline[['lon','lat','cdp']]
        trackdat['str'] = 0.0
        trackdat['elev'] = 0.0
        trackdat['grad'] = 0.0
        trackdat['side'] = 1

    contz = np.arange(-1.5, 8.5, 0.5)
    if args.ci is not None:
        contz = np.arange(-1.5, 8.5, args.ci)

    clons, clats, cdists, cconts, ccdps, cgrad = [],[],[],[],[],[]
    
    for side in set(list(trackdat['side'].values)):
        trackdati = trackdat[trackdat.side == side]
        mincontour = trackdati['elev'].min()
        maxcontour = trackdati['elev'].max()
        for c in contz:
            if c > maxcontour or c < mincontour:
                continue
            else:
                dists = np.abs(trackdati['elev'].values - c)
                trackdati['dist'] = dists
                mindist = trackdati['dist'].min()
                mintrack = trackdati[trackdati.dist == mindist]
                clons.append(mintrack['lon'].values[0])
                clats.append(mintrack['lat'].values[0])
                ccdps.append(mintrack['cdp'].values[0])
                cdists.append(mindist)
                cconts.append(c)
                cgrad.append(mintrack['grad'].values[0])

    crosscontours = pd.DataFrame({'lon':clons, 'lat':clats, 'dist':cdists,
                                    'contour':cconts, 'cdp':ccdps, 'grad':cgrad})
    crosscontours = crosscontours[crosscontours.dist < 0.009]
    crosscontours.to_csv('Output/%s_crosscontours.csv' % tapename[:-4], header=True,index=False)
    
    colors = ['lightgoldenrod4','mediumpurple4','purple4','darkorchid4',
            'mediumorchid4','plum4','orchid4','magenta4','violetred4',
            'maroon4','palevioletred4','red4','orangered4','tomato4','coral4',
            'darkorange4','orange4','brown4','firebrick4','tan4','wheat4',
            'burlywood4','sienna4','indianred4','rosybrown4','darkgoldenrod4',
            'goldenrod4','gold4','yellow4','lightyellow4',
            'khaki4','darkolivegreen4']

    if args.horizonfile is not None:
        horizondat = pd.read_csv(args.horizonfile)
        dipmax = horizondat['slope'].max()
        dipmin = horizondat['slope'].min()
        sminmax = max(dipmax, abs(dipmin))
        sminmax = 30
        os.system("gmt makecpt -Cpolar -D -Fr -T-%f/%f/0.01 -Z > seis2.cpt" % (sminmax, sminmax))
        sinterval = int(sminmax / 2)

    if justcurve:
        if args.horizonfile is None:
            print ('if using -ju argument, must include contour file associated with this curve')
            print ('exiting ... ')
            exit()
        wDmin = -15
        wDmax = 15
        wTmin = 0
        wTmax = 100
        Rlc = '-R%f/%f/%f/%f' % (wCmin,wCmax,wTmin,wTmax)
        Rl2c = '-R%f/%f/%f/%f' % (wCmin,wCmax,wDmin,wDmax)
        yplotmax = -1*yplotmax
        Jlc = '-JX%s/%s' % (xplotmax,yplotmax/2)
        if revdir == True:
            xplotmax2 = -1*xplotmax
            Jlc = '-JX%s/%s' % (xplotmax2,yplotmax/2)
        if args.bs is not None:
            Bl2c = '-B%fg%f:" ":/%fg%f:"Dip(deg)":%s:." ":' % (cdpint, cdpint,10,10, args.bs[:2])
            Blc = '-B%fg%f:" ":/%fg%f:"N Horizons":%s:." ":' % (cdpint, cdpint, 20,20, args.bs[-2:])
        else:
            Bl2c = '-B%fg%f:" ":/%fg%f:"Dip(deg)":EN:." ":' % (cdpint, cdpint,10,10)
            Blc = '-B%fg%f:" ":/%fg%f:"N Horizons":WS:." ":' % (cdpint, cdpint, 20,20)

        curvefile = '%s_valdips.csv' % args.horizonfile[:-13]
        print (curvefile)
        curvedat = pd.read_csv(curvefile)
        curvedat2 = curvedat[curvedat.nseg > 10]
        print (curvedat)
        tempinterp = curvedat[['cdp','dip','dip']]
        tempinterp.to_csv('tempinterp.txt',header=True,index=False,sep=' ')

        os.system("gmt psbasemap %s %s %s -P -K --MAP_GRID_PEN_PRIMARY=darkorange %s > %s" % (Jlc, Rlc, Blc, Pl, sf))
        tempinterp = pd.DataFrame({'cdp':[wCmin, wCmax], 'nseg':[10,10]})
        tempinterp.to_csv('tempinterp.txt',header=True,index=False,sep=' ')
        os.system("gmt psxy tempinterp.txt -R -J -W0.3,red,- -O -K %s  >> %s" % (Pl, sf))
        legn += legshift
        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, '10 segment threshold', '-W0.3,red')
        tempinterp = curvedat[['cdp','nseg']]
        tempinterp.to_csv('tempinterp.txt',header=True,index=False,sep=' ')
        os.system("gmt psxy tempinterp.txt -R -J -W0.3,darkorange -O -K %s  >> %s" % (Pl, sf))
        legn += legshift
        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'Number of segments per CDP', '-W0.3,darkorange')
        os.system("gmt psbasemap %s %s %s -P -K -O --MAP_GRID_PEN_PRIMARY=lightgrey %s >> %s" % (Jlc, Rl2c, Bl2c, Pl, sf))
        tempinterp = curvedat2[['cdp','dip','dip']]
        tempinterp.to_csv('tempinterp.txt',header=True,index=False,sep=' ')
        os.system("gmt psxy tempinterp.txt -R -J -Gdarkgray -Sc0.075 -O -K %s  >> %s" % (Pl, sf))
        os.system("gmt psxy tempinterp.txt -R -J -Cseis2.cpt -Sc0.055 -O -K %s  >> %s" % (Pl, sf))
        if args.bs is not None:
            if 'S' in args.bs:
                os.system('printf "CDP" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "CDP" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        if args.bs is not None:
            if 'n' not in args.bs:
                os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        if args.bs is not None:
            if 'n' in args.bs:
                os.system('printf "Characteristic Apparent Dip Along Profile" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
            
        os.system('printf "- = N Segments" | gmt pstext -R -J -O -K -F+cTL+f10p,Helvetica,orange %s >> %s'%(Pl, sf))
        os.system('printf "-- = 10 Segement Threshold" | gmt pstext -R -J -O -K -F+cBR+f10p,Helvetica,red %s >> %s'%(Pl, sf))

    elif MorC != 'C':
        os.system("gmt psbasemap %s %s %s -P -K --MAP_GRID_PEN_PRIMARY=lightgrey %s > %s" % (Jl, Rl2, Bl2, Pl, sf))
        os.system("gmt psbasemap %s %s %s -P -K -O -K --MAP_GRID_PEN_PRIMARY=lightgrey %s >> %s" % (Jl, Rl, Bl, Pl, sf))
        if justcontours:
            os.system("gmt grdimage Output/%s -R -J -O -K -CseisJC.cpt %s >> %s"%(gridname, Pl, sf))
        else:
            os.system("gmt grdimage Output/%s -R -J -O -K -Cseis.cpt %s >> %s"%(gridname, Pl, sf))
        if args.interp is not None:
            if args.interp == 'y':
                ci = 0
                for interp in os.listdir('library/interpretations'):
                    if args.onlyapril is not None:
                        if args.onlyapril == 'y':
                            continue
                    if 'blak' in interp:
                        continue
                    if 'BL' in interp:
                        continue
                    n = -1
                    n_horizons = 0
                    horizons = ['line','cdp','x','y']
                    parfile = 'library/interpretations/%s' % interp
                    for line in open(parfile):
                        n += 1
                        plist = line.split()
                        if n == 0:
                            n_horizons = int(plist[0])
                            continue
                        elif n <= n_horizons and n_horizons > 0:
                            horizons.append(plist[0])
                            continue
                        elif n > n_horizons:
                            break
                    n_skip = n_horizons + 1
                    dat = pd.read_csv(parfile, names=horizons, header=None, skiprows=n_skip, delim_whitespace=True)
                    linename = tapename[:-4]
                    idat1 = dat[dat.line == linename]
                    if len(idat1) > 0:
                        i = -1
                        for k in horizons:
                            i += 1
                            if i < 4:
                                continue
                            idat = idat1[idat1.ix[:,i] >0 ]
                            interpxy = np.zeros((len(idat),3))
                            navxy = np.zeros((len(navline),3))
                            interpxy[:,0] = idat['x'].values
                            interpxy[:,1] = idat['y'].values
                            navxy[:,0] = navline['gx'].values
                            navxy[:,1] = navline['gy'].values
                            navxy[:,2] = navline['cdp'].values
                            interpxy[:, 2] = griddata(navxy[:, 0:2], navxy[:, 2], interpxy[:, 0:2], method='nearest')
                            idat['cdp'] = interpxy[:,2]
                            idat = idat[['cdp',k]]
                            idat[k] = idat[k].values / 1000.0
                            idat.to_csv('tempinterp.txt', header=False, index=False)
                            plotcolor = colors[i]
                            plotcolor = 'cyan'
                            if 'dip' in k:
                                plotcolor = 'yellow'
                            if 'seafloor' in k:
                                plotcolor = 'red'
                            if 'shallow' in k:
                                plotcolor = 'orange'
                            if 'OPFb' in k or 'sproulomps' in k:
                                print ('not plotting OPFb fault ...')
                                continue

                            os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.03 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                            legn += legshift
                            legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, k, '-W0.8,%s'%plotcolor)
                            os.system("rm tempinterp.txt")
                            with open('Output/%s_legend.txt'%linename, 'a') as file:
                                file.write('%s: %s \n' % (plotcolor,k))
                            ci += 1
            
                # going through april 2019 interpretations
                ci = 0
                i = 0
                for interp in os.listdir('library/april19_interpretations'):
                    if not interp.endswith('.dat'):
                        continue
                    if 'fault' in interp:
                        continue
                    n = -1
                    n_horizons = 0
                    horizons = ['line','cdp','gx','gy','time']
                    parfile = 'library/april19_interpretations/%s' % interp
                    dat = pd.read_csv(parfile, names=horizons, header=None, delim_whitespace=True)
                    linename = tapename[:-4]
                    idat1 = dat[dat.line == linename]
                    if len(idat1) > 0:
                        k = interp
                        print ('plotting %s ...' %k)
                        
                        interpxy = np.zeros((len(idat1),3))
                        navxy = np.zeros((len(navline),3))
                        interpxy[:,0] = idat1['gx'].values
                        interpxy[:,1] = idat1['gy'].values
                        navxy[:,0] = navline['gx'].values
                        navxy[:,1] = navline['gy'].values
                        navxy[:,2] = navline['cdp'].values
                        interpxy[:, 2] = griddata(navxy[:, 0:2], navxy[:, 2], interpxy[:, 0:2], method='nearest')
                        idat1['cdp'] = interpxy[:,2]
                        #idat = idat[['cdp',k]]
                        #idat[k] = idat[k].values / 1000.0
                        
                        # to save interpretations
                        if args.keepinterp is not None:
                            if args.keepinterp == 'y':
                                idat1.to_csv('tempinterp_%s.txt' % interp[:-4], header=False, index=False)
                        
                        idat1 = idat1[['cdp','time']]
                        idat1.to_csv('tempinterp.txt', header=False, index=False)

                        plotcolor = colors[i]
                        if 'dip' in k and 'main' in k:
                            plotcolor = 'red'
                        if 'lumpys' in k:
                            plotcolor = 'blue'
                        if 'lumpyn' in k:
                            plotcolor = 'black'
                        if 'seafloor' in k:
                            plotcolor = 'red'
                        if 'shallow' in k:
                            plotcolor = 'orange'
                        if 'blak' in k:
                            plotcolor = 'black'
                        if 'OPFb' in k or 'sproulomps' in k:
                            print ('not plotting OPFb ... ')
                            continue
                            
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.03 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        legn += legshift
                        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, k, '-W0.8,%s'%plotcolor)
                        os.system("rm tempinterp.txt")
                        with open('Output/%s_legend.txt'%linename, 'a') as file:
                            file.write('%s: %s \n' % (colors[ci],k))
                        ci += 1
                        i += 1

        if args.horizonfile is not None:
            horizonfile = args.horizonfile
            horzdata = pd.read_csv(horizonfile)
            horzdata = horzdata[horzdata.g > 0]
            for s in set(list(horzdata['s'].values)):
                idat1 = horzdata[horzdata.s == s]
                if len(idat1) > 0:
                    idat1 = idat1[['x','y','slope']]
                    idat1.to_csv('tempinterp.txt', header=False, index=False)
                    plotcolor = 'black'
                    if not phc:
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.02 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                    else:
                        plotcolor = 'gray'
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.025 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        os.system("gmt psxy tempinterp.txt -R -J -Cseis2.cpt -Sc0.02 -O -K %s  >> %s" % (Pl, sf))
                    os.system("rm tempinterp.txt")
        if args.horizonfile2 is not None:
            horizonfile = args.horizonfile2
            horzdata = pd.read_csv(horizonfile)
            horzdata = horzdata[horzdata.g > 0]
            for s in set(list(horzdata['s'].values)):
                idat1 = horzdata[horzdata.s == s]
                if len(idat1) > 0:
                    idat1 = idat1[['x','y','slope']]
                    idat1.to_csv('tempinterp.txt', header=False, index=False)
                    plotcolor = 'black'
                    if not phc:
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.02 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                    else:
                        plotcolor = 'gray'
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.025 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        os.system("gmt psxy tempinterp.txt -R -J -Cseis2.cpt -Sc0.02 -O -K %s  >> %s" % (Pl, sf))
                    os.system("rm tempinterp.txt")
        if args.horizonfile3 is not None:
            horizonfile = args.horizonfile3
            horzdata = pd.read_csv(horizonfile)
            horzdata = horzdata[horzdata.g > 0]
            for s in set(list(horzdata['s'].values)):
                idat1 = horzdata[horzdata.s == s]
                if len(idat1) > 0:
                    idat1 = idat1[['x','y','slope']]
                    idat1.to_csv('tempinterp.txt', header=False, index=False)
                    plotcolor = 'black'
                    if not phc:
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.02 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                    else:
                        plotcolor = 'gray'
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.025 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        os.system("gmt psxy tempinterp.txt -R -J -Cseis2.cpt -Sc0.02 -O -K %s  >> %s" % (Pl, sf))
                    os.system("rm tempinterp.txt")
        if args.horizonfile4 is not None:
            horizonfile = args.horizonfile4
            horzdata = pd.read_csv(horizonfile)
            horzdata = horzdata[horzdata.g > 0]
            for s in set(list(horzdata['s'].values)):
                idat1 = horzdata[horzdata.s == s]
                if len(idat1) > 0:
                    idat1 = idat1[['x','y','slope']]
                    idat1.to_csv('tempinterp.txt', header=False, index=False)
                    plotcolor = 'black'
                    if not phc:
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.02 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                    else:
                        plotcolor = 'gray'
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.025 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        os.system("gmt psxy tempinterp.txt -R -J -Cseis2.cpt -Sc0.02 -O -K %s  >> %s" % (Pl, sf))
                    os.system("rm tempinterp.txt")

        if plotWB:
            os.system("gmt psxy wbplottemp.txt -R -J -Wblack -O -K %s  >> %s" % (Pl, sf))
            legn += legshift
            legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'water bottom', '-Wblack')

        if args.addbox1 is not None:
            ax1 = args.addbox1[0]
            ax2 = args.addbox1[1]
            ax3 = args.addbox1[1]
            ax4 = args.addbox1[0]
            ax5 = args.addbox1[0]
            ay1 = args.addbox1[2]
            ay2 = args.addbox1[2]
            ay3 = args.addbox1[3]
            ay4 = args.addbox1[3]
            ay5 = args.addbox1[2]
            
            plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
            plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
            os.system("gmt psxy temp.txt -R -J -W0.03,black -O -K %s >> %s" % (Pl, sf))
            os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black %s >> %s'%(ax2, ay3, args.addbox1[4],Pl, sf))

        if args.addbox2 is not None:
            ax1 = args.addbox2[0]
            ax2 = args.addbox2[1]
            ax3 = args.addbox2[1]
            ax4 = args.addbox2[0]
            ax5 = args.addbox2[0]
            ay1 = args.addbox2[2]
            ay2 = args.addbox2[2]
            ay3 = args.addbox2[3]
            ay4 = args.addbox2[3]
            ay5 = args.addbox2[2]
            
            plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
            plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
            os.system("gmt psxy temp.txt -R -J -W0.03,black -O -K %s >> %s" % (Pl, sf))
            os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black %s >> %s'%(ax2, ay3, args.addbox2[4],Pl, sf))

        os.system("echo VE=%sx | gmt pstext -R -J -O -K -F+cTL+f10p,Helvetica,black %s >> %s"%(VE, Pl, sf))
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBL+f10p,Helvetica,black %s >> %s"%(llabel, Pl, sf))
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBR+f10p,Helvetica,black %s >> %s"%(rlabel, Pl, sf))
        #os.system('printf "%s %s\n%s %s" | gmt psxy -R -J -W0.7,blue,- -O -K %s >> %s' % (centcdp, wTmin, centcdp, wTmax, Pl, sf))
        if args.bs is not None:
            if 'n' not in args.bs:
                os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))

        if phc:
            nocontours = True
        else:
            green = 'green4'
        if not nocontours:
            for index,row in crosscontours.iterrows():
                ccdp, cont = row['cdp'], row['contour']
                if cont < 0:
                    continue
                linewidth = 0.5
                cval = 9.0-cont
                color = 'grey%i'%int(cval*10)
                color = '-W%0.2f,%s'%(cont/10*cont/10*cont/10*5,green)
                os.system('printf "%s %s\n%s %s" | gmt psxy -R -J %s -O -K %s >> %s' % (ccdp, wTmin, ccdp, wTmax, color, Pl, sf))
                legn += legshift
                legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, '%s meters of uplift'%cont, color)
                if float(cont) != 6.0 and float(cont) != 4.0 and float(cont) != 2.0 and float(cont) != 0.0:
                    continue
                os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -F+f11p,Helvetica,%s %s >> %s'%(ccdp, (wTmax-wTmin)/1.2, cont,green,Pl, sf))
        if args.bs is not None:
            if 'S' in args.bs:
                os.system('printf "CDP" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "CDP" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        linelist = linename.split('_')
        if len(linelist) > 1:
            plotname = linelist[0]
        else:
            plotname = linename
        os.system('printf "%s" | gmt pstext -R -J -O -K -F+cTR+f10p,Helvetica,blue %s >> %s'%(plotname,Pl, sf))
    else:
        os.system("gmt grdinfo Output/%s > temp.txt" % (gridname))
        zmax = -900000
        for line in open('temp.txt'):
            plist = line.split()
            if len(plist)==7 and 'zmax' in plist:
                zmax = float(plist[4])
        os.system("gmt grdmath Output/%s %s DIV = Output/%s_norm.grd" % (gridname, zmax, gridname[:-4]))
        normgridname = '%s_norm.grd' % gridname[:-4]
        pngridname = '%s_pn.grd' % gridname[:-4]
        os.system("gmt grdgradient Output/%s -A0 -GOutput/%s -V" % (normgridname, pngridname))
        os.system("gmt psbasemap %s %s %s -P -K --MAP_GRID_PEN_PRIMARY=lightgrey  %s > %s" % (Jl, Rl2, Bl2, Pl, sf))
        print ("gmt psbasemap %s %s %s -P -K --MAP_GRID_PEN_PRIMARY=lightgrey  %s > %s" % (Jl, Rl2, Bl2, Pl, sf))
        print (' ^^^ ' )
        os.system("gmt psbasemap %s %s %s -P -K -O -K --MAP_GRID_PEN_PRIMARY=lightgrey %s >> %s" % (Jl, Rl, Bl, Pl, sf))
        print ("gmt psbasemap %s %s %s -P -K -O -K --MAP_GRID_PEN_PRIMARY=lightgrey %s >> %s" % (Jl, Rl, Bl, Pl, sf))
        print (' ^^^ ' )
        if not justcontours:
            os.system("gmt grdimage Output/%s -R -J -O -K -Cchirp.cpt  %s >> %s"%(gridname, Pl, sf))
        #os.system("gmt grdimage Output/%s -R -J -O -K -Cseis.cpt  %s >> %s"%(pngridname, Pl, sf))
        #if plotWB:
            #os.system("gmt psxy wbplottempunfilt.txt -R -J -Wyellow -O -K %s  >> %s" % (Pl, sf))
            #os.system("gmt psxy wbplottemp.txt -R -J -Wred -O -K %s  >> %s" % (Pl, sf))
        #os.system("gmt grdcontour Output/%s -C10 -L9/10 -R -J -O -K -W0.04,blue >> %s"%(gridname,sf))
        #os.system("gmt grdcontour Output/%s -C10 -L19/20 -R -J -O -K -W0.03,blue >> %s"%(gridname,sf))
        #os.system("gmt grdcontour Output/%s -C10 -L49/50 -R -J -O -K -W0.02,yellow >> %s"%(gridname,sf))
        #os.system("gmt grdcontour Output/%s -C10 -L99/100 -R -J -O -K -W0.01,red >> %s"%(gridname,sf))
        #os.system("gmt grdcontour gradient.nc -C10 -L-1/1 -R -J -O -K -W0.04,black >> %s"%(sf))
        #os.system("gmt grdcontour gradient2.nc -C10 -L-1/1 -R -J -O -K -W0.02,blue >> %s"%(sf))
        if args.interp is not None:
            if args.interp == 'y':
                ci = 0
                for interp in os.listdir('library/interpretations'):
                    n = -1
                    n_horizons = 0
                    horizons = ['line','cdp','x','y']
                    parfile = 'library/interpretations/%s' % interp
                    for line in open(parfile):
                        n += 1
                        plist = line.split()
                        if n == 0:
                            n_horizons = int(plist[0])
                            continue
                        elif n <= n_horizons and n_horizons > 0:
                            horizons.append(plist[0])
                            continue
                        elif n > n_horizons:
                            break
                    n_skip = n_horizons + 1
                    dat = pd.read_csv(parfile, names=horizons, header=None, skiprows=n_skip, delim_whitespace=True)
                    linename = tapename[:-4]
                    idat1 = dat[dat.line == linename]
                    if len(idat1) > 0:
                        i = -1
                        for k in horizons:
                            i += 1
                            if i < 4:
                                continue
                            idat = idat1[idat1.ix[:,i] >0 ]
                            interpxy = np.zeros((len(idat),3))
                            navxy = np.zeros((len(navline),3))
                            interpxy[:,0] = idat['x'].values
                            interpxy[:,1] = idat['y'].values
                            navxy[:,0] = navline['gx'].values
                            navxy[:,1] = navline['gy'].values
                            navxy[:,2] = navline['cdp'].values
                            interpxy[:, 2] = griddata(navxy[:, 0:2], navxy[:, 2], interpxy[:, 0:2], method='nearest')
                            idat['cdp'] = interpxy[:,2]
                            idat = idat[['cdp',k]]
                            idat[k] = idat[k].values / 1000.0
                            idat.to_csv('tempinterp.txt', header=False, index=False)
                            plotcolor = colors[i]
                            if 'dip' in k:
                                plotcolor = 'yellow'
                            if 'seafloor' in k:
                                plotcolor = 'red'
                            if 'shallow' in k:
                                plotcolor = 'orange'
                            os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.03 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                            legn += legshift
                            legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, k, '-W0.8,%s'%plotcolor)
                            os.system("rm tempinterp.txt")
                            with open('Output/%s_legend.txt'%linename, 'a') as file:
                                file.write('%s: %s \n' % (colors[ci],k))
                            ci += 1

                # going through april 2019 interpretations
                ci = 0
                i = 0
                for interp in os.listdir('library/april19_interpretations'):
                    if not interp.endswith('.dat'):
                        continue
                    if 'fault' in interp:
                        continue
                    n = -1
                    n_horizons = 0
                    horizons = ['line','cdp','gx','gy','time']
                    parfile = 'library/april19_interpretations/%s' % interp
                    dat = pd.read_csv(parfile, names=horizons, header=None, delim_whitespace=True)
                    linename = tapename[:-4]
                    idat1 = dat[dat.line == linename]
                    if len(idat1) > 0:
                        k = interp
                        print ('plotting %s ...' %k)
                        
                        interpxy = np.zeros((len(idat1),3))
                        navxy = np.zeros((len(navline),3))
                        interpxy[:,0] = idat1['gx'].values
                        interpxy[:,1] = idat1['gy'].values
                        navxy[:,0] = navline['gx'].values
                        navxy[:,1] = navline['gy'].values
                        navxy[:,2] = navline['cdp'].values
                        interpxy[:, 2] = griddata(navxy[:, 0:2], navxy[:, 2], interpxy[:, 0:2], method='nearest')
                        idat1['cdp'] = interpxy[:,2]
                        #idat = idat[['cdp',k]]
                        #idat[k] = idat[k].values / 1000.0
                        
                        # to save interpretations
                        if args.keepinterp is not None:
                            if args.keepinterp == 'y':
                                idat1.to_csv('tempinterp_%s.txt' % interp[:-4], header=False, index=False)
                        
                        idat1 = idat1[['cdp','time']]
                        idat1.to_csv('tempinterp.txt', header=False, index=False)
                        plotcolor = colors[i]
                        if 'dip' in k and 'main' in k:
                            plotcolor = 'red'
                        if 'lumpys' in k:
                            plotcolor = 'blue'
                        if 'lumpyn' in k:
                            plotcolor = 'black'
                        if 'seafloor' in k:
                            plotcolor = 'red'
                        if 'shallow' in k:
                            plotcolor = 'orange'
                        if 'blak' in k:
                            plotcolor = 'black'
                        if 'OPFb' in k or 'sproulomps' in k:
                            print ('not plotting OPFb fault ... ')
                            continue
                        
                        os.system("gmt psxy tempinterp.txt -R -J -G%s -Sc0.03 -O -K %s  >> %s" % (plotcolor,Pl, sf))
                        legn += legshift
                        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, k, '-W0.8,%s'%plotcolor)
                        os.system("rm tempinterp.txt")
                        with open('Output/%s_legend.txt'%linename, 'a') as file:
                            file.write('%s: %s \n' % (colors[ci],k))
                        ci += 1
                        i += 1

        if args.addbox1 is not None:
            ax1 = args.addbox1[0]
            ax2 = args.addbox1[1]
            ax3 = args.addbox1[1]
            ax4 = args.addbox1[0]
            ax5 = args.addbox1[0]
            ay1 = args.addbox1[2]
            ay2 = args.addbox1[2]
            ay3 = args.addbox1[3]
            ay4 = args.addbox1[3]
            ay5 = args.addbox1[2]
            plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
            plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
            os.system("gmt psxy temp.txt -R -J -W0.03,red -O -K %s >> %s" % (Pl, sf))
            os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,red %s >> %s'%(ax2, ay3, args.addbox1[4],Pl, sf))

        if args.addbox2 is not None:
            ax1 = args.addbox2[0]
            ax2 = args.addbox2[1]
            ax3 = args.addbox2[1]
            ax4 = args.addbox2[0]
            ax5 = args.addbox2[0]
            ay1 = args.addbox2[2]
            ay2 = args.addbox2[2]
            ay3 = args.addbox2[3]
            ay4 = args.addbox2[3]
            ay5 = args.addbox2[2]
            plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
            plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
            os.system("gmt psxy temp.txt -R -J -W0.03,red -O -K %s >> %s" % (Pl, sf))
            os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,red %s >> %s'%(ax2, ay3, args.addbox2[4],Pl, sf))
        

        os.system("echo VE=%sx | gmt pstext -R -J -O -K -F+cTL+f10p,Helvetica,red %s >> %s"%(VE, Pl, sf))
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBL+f10p,Helvetica,red %s >> %s"%(llabel, Pl, sf))
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBR+f10p,Helvetica,red %s >> %s"%(rlabel, Pl, sf))
        #os.system('printf "%s %s\n%s %s" | gmt psxy -R -J -W0.7,blue,- -O -K %s >> %s' % (centcdp, wTmin, centcdp, wTmax, Pl, sf))
        if args.bs is not None:
            if 'n' not in args.bs:
                os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        else:
            os.system('printf "Distance(m)" | gmt pstext -R -J -O -K -F+cTC+f10p,Helvetica,black %s >> %s'%(Pl, sf))
        if phc:
            nocontours = True
        else:
            green = 'green4'
        
        if not nocontours:
            for index,row in crosscontours.iterrows():
                ccdp, cont = row['cdp'], row['contour']
                if cont < 0:
                    continue
                linewidth = 0.5
                cval = 9.0-cont
                color = 'grey%i'%int(cval*10)
                color = '-W%0.2f,%s'%(cont/10*cont/10*cont/10*5,green)
                os.system('printf "%s %s\n%s %s" | gmt psxy -R -J %s -O -K %s >> %s' % (ccdp, wTmin, ccdp, wTmax, color, Pl, sf))
                legn += legshift
                legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, '%s meters of uplift'%cont, color)
                if float(cont) != 6.0 and float(cont) != 4.0 and float(cont) != 2.0 and float(cont) != 0.0:
                    continue
                os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -F+f11p,Helvetica,%s %s >> %s'%(ccdp, (wTmax-wTmin)/1.2, cont,green,Pl, sf))
        if args.bs is not None:
            if 'S' in args.bs:
                os.system('printf "RecNum" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black  %s >> %s'%(Pl, sf))
        else:
                os.system('printf "RecNum" | gmt pstext -R -J -O -K -F+cBC+f10p,Helvetica,black  %s >> %s'%(Pl, sf))
        linelist = linename.split('_')
        if len(linelist) > 1:
            plotname = linelist[0]
        else:
            plotname = linename
        os.system('printf "%s" | gmt pstext -R -J -O -K -F+cTR+f10p,Helvetica,blue %s >> %s'%(plotname,Pl, sf))

    if 'delrt' in navline.columns:
        if navline['lat'].max() > 50:
            AKbathy = '/Users/ginevramoore/Documents/research/AACSE/bathymetry_grids/bathygrid_09.11.2019.txt'
            AKbathyd = pd.read_csv(AKbathy, delim_whitespace=True, names=['x','y','z'])
            AKbathyd = AKbathyd[(AKbathyd.x >= mXmin) & (AKbathyd.x <= mXmax) & \
                                (AKbathyd.y >= mYmin) & (AKbathyd.y <= mYmax)]

            print (AKbathyd['z'].min(), AKbathyd['z'].max())
            os.system("gmt makecpt -Crainbow -D -Fr -T%f/%f/1 -Z > multibeam.cpt" %
                                        (AKbathyd['z'].min(), AKbathyd['z'].max()))

    if args.nomap is not None:
        if phc:
            os.system("gmt psscale -D0.4/12/2/0.8 -Cseis2.cpt -B%i -O -K >> %s"%(sinterval, sf))
        if args.nomap == 'y':
            print ('indicated to not plot map using -nm flag')
            print ('plotting legend ... ')
            plot_legend(legns, legvars, legcolors, legendname)
            print ('exiting ... ')
            exit()

    '''                          zoomed out map!                             '''

    DEM2 = '/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc'
    os.system("gmt psbasemap %s %s %s %s -O -K --FORMAT_GEO_MAP=D >> %s" % (Jm, Rm, Bm, Xm, sf))
    print ("gmt psbasemap %s %s %s %s -O -K --FORMAT_GEO_MAP=D >> %s" % (Jm, Rm, Bm, Xm, sf))
    print (' ^^^ ')
    #os.system("gmt psbasemap %s %s %s %s -O -K --FORMAT_GEO_MAP=D >> %s" % (Jx, Rm, Bm, Xm, sf))
    os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade2.grd >> %s"%(DEM2, gmcpt, sf))
    #os.system(" gmt psxy -R -J %s -O -K -Sc0.005 -Cmultibeam.cpt >> %s" % (AKbathy, sf))
    os.system("gmt psxy plotclip.dat -R -J -W1,red -O -K >> %s" % (sf))
    
    os.system("gmt grdcontour %s -C100000 -R -J -O -K -W0.1,black >> %s"%(DEM2,sf))
    '''
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.3,green4 -GN-4 -L-2/1 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.4,green4 -GN-4 -L1/2 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.5,green4 -GN-4 -L2/3 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.6,green4 -GN-4 -L3/4 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.7,green4 -GN-4 -L4/5 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.8,green4 -GN-4 -L5/6 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.9,green4 -GN-4 -L6/7 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.10,green4 -GN-4 -L7/8 >> %s"%(Rm, Jm, sf))
    '''
    linefile = '%s/temp.txt'%navdir
    tempnav = navline[['lon','lat']]
    tempnav.to_csv(linefile,header=False,index=False, sep = ' ')
    #os.system("gmt psxy %s -R -J -W0.1,darkblue -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt psxy %s -R -J -Sc0.05 -Gblack -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt pstext %s -R -J -Wblack -O -K -F+f6p,Helvetica,darkblue >> %s" % (prattfiles, sf))
    
    '''
    os.system("gmt psxy %s -R -J -W0.3,blue -O -K >> %s" % (PSSnav, sf))
    print (KNUnav)
    '''
    os.system("gmt psxy %s -R -J -Sc0.02 -Gdarkblue -O -K >> %s" % (KNUnav, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'profile locations', '-W0.4,darkblue')
    os.system("gmt psxy %s -R -J -Sc0.1 -Gcyan -O -K >> %s" % (KNUnave, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'profile ends', '-Sc0.1 -Gcyan')
    os.system("gmt psxy %s -R -J -Sc0.1 -Gcyan -O -K >> %s" % (KNUnavb, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f5p,Helvetica,cyan >> %s" % (KNUnavtxt, sf))
    '''
    os.system("gmt psxy %s -R -J -W0.3,blue -O -K >> %s" % (LWBnav, sf))
    os.system("gmt psxy %s -R -J -W0.3,cyan -O -K >> %s" % (sproulnav, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,blue >> %s" % (PSSnavtxt, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,blue >> %s" % (LWBnavtxt, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,cyan >> %s" % (sproulnavtxt, sf))
    os.system("gmt psxy %s -R -J -W1,yellow -O -K >> %s" % (linefile, sf))
    '''
    windnav = navline[(navline.cdp > wCmin) & (navline.cdp < wCmax)]
    windnav = windnav[['lon','lat']]
    
    windnav.to_csv(linefile,header=False,index=False, sep=' ')

    firstx = windnav['lon'].values[0]
    firsty = windnav['lat'].values[0]
    lastx = windnav['lon'].values[-1]
    lasty = navline['lat'].values[-1]

    os.system("gmt psxy %s -R -J -W0.4,red -O -K >> %s" % (linefile, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'profile shown', '-W0.4,red')
    #os.system("echo %s %s | gmt psxy -R -J -Sc0.1 -Gblue -O -K >> %s" % (centlon,centlat, sf))
    #os.system("echo B | gmt pstext -R -J -F+f5p,Helvetica,red -O -K >> %s" % (firstx,firsty,linename,sf))
    '''
    for faultname in os.listdir('library/faults/ndipping'):
        print('library/faults/ndipping/%s'%faultname)
        os.system("gmt psxy %s -R -J -W1.5,navy -O -K >> %s" % ('library/faults/ndipping/%s'%faultname, sf))
    for faultname in os.listdir('library/faults/sdipping'):
        os.system("gmt psxy %s -R -J -W1.5,darkred -O -K >> %s" % ('library/faults/sdipping/%s'%faultname, sf))
        print ('library/faults/sdipping/%s'%faultname)
    '''
    '''
    os.system("gmt psxy %s -R -J -W1,red -O -K >> %s" % (oligocene, sf))
    os.system("gmt psxy %s -R -J -W1,darkorange -O -K >> %s" % (miocene, sf))
    os.system("gmt pstext library/labels/sfzmap.txt -R -J -O -K >> %s" % (sf))
    os.system("gmt psxy %s -R -J -Sc0.1 -Ggreen2 -O -K >> %s" % (upliftpts, sf))
    os.system("gmt pstext %s -R -J -W0.5,green2 -X0.2 -O -K >> %s" % (upliftpts, sf))
    '''
    #os.system("gmt psscale -D0.4/2.5/4/0.4 -C%s -B1000 -O -K >> %s"%(gmcpt, sf))
    #os.system("gmt psxy %s -R -J -Sc0.1 -Gwhite -O -K >> %s" % (upliftpts, sf))
    #os.system("gmt pstext %s -R -J -Wwhite -X0.2 -O -K >> %s" % (upliftpts, sf))


    '''                          zoomed in map!                             '''

    if 'delrt' in navline.columns:
        if navline['lat'].max() > 50:
            AKbathy = '/Users/ginevramoore/Documents/research/AACSE/bathymetry_grids/bathygrid_09.11.2019.txt'
            AKbathyd = pd.read_csv(AKbathy, delim_whitespace=True, names=['x','y','z'])
            AKbathyd = AKbathyd[(AKbathyd.x >= clipmXmin) & (AKbathyd.x <= clipmXmax) & \
                                (AKbathyd.y >= clipmYmin) & (AKbathyd.y <= clipmYmax)]

            print (AKbathyd['z'].min(), AKbathyd['z'].max())
            os.system("gmt makecpt -Crainbow -D -Fr -T%f/%f/1 -Z > multibeam.cpt" %
                                        (AKbathyd['z'].min(), AKbathyd['z'].max()))

    os.system("gmt psbasemap %s %s %s %s -O -K --FORMAT_GEO_MAP=D >> %s" % (Jm, Rm2, Bm2, Xm2, sf))
    print ("gmt psbasemap %s %s %s %s -O -K --FORMAT_GEO_MAP=D >> %s" % (Jm, Rm2, Bm2, Xm2, sf))
    print (' ^^^ ')
    os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade.grd >> %s"%(DEM, gmcpt, sf))
    print ("gmt grdimage %s -R -J -O -K -C%s -Ihillshade.grd >> %s"%(DEM, gmcpt, sf))
    print ( ' ^^^ ' )
    if 'delrt' in navline.columns:
        if navline['lat'].max() > 50:
            os.system(" gmt psxy -R -J %s -O -K -Sc0.02 -Cmultibeam.cpt >> %s" % (AKbathy, sf))
    else:
        print ('no AK bathy')

    os.system("gmt grdcontour %s -C1000 -R -J -O -K -W0.1,black >> %s"%(DEM,sf))

    if phc:
        green = 'green4'
    else:
        green = 'green4'
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.3,%s -GN-4 -L-2/1 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.4,%s -GN-4 -L1/2 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.5,%s -GN-4 -L2/3 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.6,%s -GN-4 -L3/4 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.7,%s -GN-4 -L4/5 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.8,%s -GN-4 -L5/6 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.9,%s -GN-4 -L6/7 >> %s"%(Rm2, Jm, green, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.10,%s -GN-4 -L7/8 >> %s"%(Rm2, Jm, green, sf))

    linefile = '%s/temp.txt'%navdir
    tempnav = navline[['lon','lat']]
    tempnav.to_csv(linefile,header=False,index=False, sep = ' ')
    #os.system("gmt psxy %s -R -J -W0.1,darkblue -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt psxy %s -R -J -Sc0.05 -Gblack -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt pstext %s -R -J -Wblack -O -K -F+f6p,Helvetica,darkblue >> %s" % (prattfiles, sf))

    os.system("gmt psxy %s -R -J -W0.3,blue -O -K >> %s" % (PSSnav, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, '2017 survey locations', '-W0.3,blue')
    print (KNUnav)
    os.system("gmt psxy %s -R -J -Sc0.02 -Gdarkblue -O -K >> %s" % (KNUnav, sf))
    os.system("gmt psxy %s -R -J -Sc0.1 -Gcyan -O -K >> %s" % (KNUnave, sf))
    os.system("gmt psxy %s -R -J -Sc0.1 -Gcyan -O -K >> %s" % (KNUnavb, sf))
    os.system("gmt psxy %s -R -J -W0.3,blue -O -K >> %s" % (LWBnav, sf))
    os.system("gmt psxy %s -R -J -W0.3,cyan -O -K >> %s" % (sproulnav, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, '2011 survey locations', '-W0.3,cyan')
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,blue >> %s" % (PSSnavtxt, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f5p,Helvetica,cyan >> %s" % (KNUnavtxt, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,blue >> %s" % (LWBnavtxt, sf))
    os.system("gmt pstext %s -R -J -O -K -F+f2p,Helvetica,cyan >> %s" % (sproulnavtxt, sf))
    os.system("gmt psxy %s -R -J -W1,yellow -O -K >> %s" % (linefile, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'profile shown', '-W1,yellow')
    windnav = navline[(navline.cdp > wCmin) & (navline.cdp < wCmax)]
    windnav = windnav[['lon','lat']]

    windnav.to_csv(linefile,header=False,index=False, sep=' ')

    firstx = windnav['lon'].values[0]
    firsty = windnav['lat'].values[0]
    lastx = windnav['lon'].values[-1]
    lasty = navline['lat'].values[-1]

    os.system("gmt psxy %s -R -J -W0.4,red -O -K >> %s" % (linefile, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'profile segment shown', '-W0.4,red')
    #os.system("echo %s %s | gmt psxy -R -J -Sc0.1 -Gblue -O -K >> %s" % (centlon,centlat, sf))
    os.system("echo %s %s %s | gmt pstext -R -J -F+f5p,Helvetica,red -O -K >> %s" % (firstx,firsty,linename,sf))

    dontplotfaults = False
    if args.nofaults is not None:
        if args.nofaults == 'y':
            dontplotfaults = True
    
    if not dontplotfaults:
        for faultname in os.listdir('library/faults/ndipping'):
            print('library/faults/ndipping/%s'%faultname)
            os.system("gmt psxy %s -R -J -W1.5,navy -O -K >> %s" % ('library/faults/ndipping/%s'%faultname, sf))
        legn += legshift
        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'north dipping faults', '-W1.5,navy')
        for faultname in os.listdir('library/faults/sdipping'):
            os.system("gmt psxy %s -R -J -W1.5,darkred -O -K >> %s" % ('library/faults/sdipping/%s'%faultname, sf))
            print ('library/faults/sdipping/%s'%faultname)
        legn += legshift
        legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'south dipping faults', '-W1.5,darkred')

    if args.ppf is not None:
        nomag = False
        if args.oma is not None:
            if args.oma == 'y':
                nomag = True

        justdf = False
        if args.jdf is not None:
            if args.jdf == 'y':
                justdf = True

        if args.ppf == 'y':
            colors = ['lightgoldenrod4','mediumpurple4','purple4','darkorchid4',
            'mediumorchid4','plum4','orchid4','magenta4','violetred4',
            'maroon4','palevioletred4','red4','orangered4','tomato4','coral4',
            'darkorange4','orange4','brown4','firebrick4','tan4','wheat4',
            'burlywood4','sienna4','indianred4','rosybrown4','darkgoldenrod4',
            'goldenrod4','gold4','yellow4','lightyellow4',
            'khaki4','darkolivegreen4']
            ci = 0
            os.system("rm %s_legend.txt "%sf[:-3])
            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/blakely2002'):
                print('/Users/ginevramoore/Documents/research/segymanip/library/faults/blakely2002/%s'%faultname)
                color = colors[ci]
                if justdf and 'iDEF' not in faultname:
                    continue
                
                if 'GEO1' in faultname:
                    if nomag:
                        continue
                    else:
                        color = 'purple4,-'
                if 'GEO2' in faultname:
                    if nomag:
                        continue
                    else:
                        color = 'darkorange3,-'
                if 'GEO3' in faultname:
                    if nomag:
                        continue
                    else:
                        color = 'darkorange3,-'
                if 'GEO4' in faultname:
                    if nomag:
                        continue
                    else:
                        color = 'darkorange3,-'
                if 'iDEF' in faultname and not justdf:
                    continue # deformation front
                if 'iDEF' in faultname and justdf:
                    color = 'gray30'
                if 'DEF' in faultname and 'iDEF' not in faultname:
                    color = 'red4'
                if 'FLT1' in faultname:
                    color = 'red4'
                if 'FLT2' in faultname:
                    continue # frontal fault
                if 'FLT3' in faultname:
                    continue # frontal fault
                if 'FLT4' in faultname:
                    continue # deformation front
                if 'FLT5' in faultname:
                    continue
                if 'FLT6' in faultname:
                    continue
                if 'FLT7' in faultname:
                    color = 'red4'# west extension of orchard point fault
                if 'FLT8_' in faultname:
                    continue
                    #color = 'red4' # mid (s2?)
                if 'FLT9' in faultname:
                    color = 'red4' # N crescent fm contact west of LW
                if 'FLT10' in faultname:
                    color = 'red4'
                if 'FLT11_' in faultname:
                    color = 'red4' # orchard point fault
                if 'FLT12' in faultname:
                    color = 'red4' # blakely harbor fault
                if 'FLT13_' in faultname:
                    color='red4'
                if color == 'red4' or color == 'gray30':
                    thickness = 1.2
                else:
                    thickness = 0.7
                os.system("gmt psxy %s -R -J -W%s,%s -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/blakely2002/%s'%faultname, thickness, color, sf))
                with open('%s_legend.txt'%sf[:-3], 'a') as file:
                    file.write('%s: %s \n' % (color,faultname))
                ci += 1

            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/ndipping'):
                if 'pratt' not in faultname and 'moore' not in faultname and 'kelsey' not in faultname:
                    continue
                if justdf:
                    continue
                print('/Users/ginevramoore/Documents/research/segymanip/library/faults/ndipping/%s'%faultname)
                if args.pbl is not None:
                    if args.pbl == 'y':
                        pbl = True
                    else:
                        pbl = False
                else:
                    pbl = False

                if pbl:
                    os.system("gmt psxy %s -R -J -W1.0,blue2 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/ndipping/%s'%faultname, sf))
                else:
                    dat = pd.read_csv('/Users/ginevramoore/Documents/research/segymanip/library/faults/ndipping/%s'%faultname, delim_whitespace=True)
                    print(dat)
                    dat['PLOT_ANG'] = 270
                    try:
                        dat = dat[dat.lon != '>']
                    except:
                        print ("'>' not in file")
                    dat['lat'] = dat['lat'].values.astype(float) + 0.001
                    dat['len'] = '0.01i'
                    dat = dat[['lon','lat','PLOT_ANG','len']]
                    dat.to_csv('aztest0.txt', header=False, index=False, sep=' ', na_rep=np.nan)

                    os.system('gmt psxy aztest0.txt -Sv0.05i+e -Gblue -R -J -O -K  >> %s' % (sf))
                    
    os.system("gmt psxy %s -R -J -W1,red -O -K >> %s" % (oligocene, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'Oligocene contacts', '-W1,red')
    os.system("gmt psxy %s -R -J -W1,darkorange -O -K >> %s" % (miocene, sf))
    legn += legshift
    legns, legvars, legcolors  = add_to_legend(legns, legvars, legcolors, legn, legendname, 'Miocene contacts', '-W1,darkorange')
    os.system("gmt pstext library/labels/sfzmap.txt -R -J -O -K >> %s" % (sf))
    #os.system("gmt psxy %s -R -J -Sc0.1 -Ggreen2 -O -K >> %s" % (upliftpts, sf))
    #os.system("gmt pstext %s -R -J -W0.5,green2 -X0.2 -O -K >> %s" % (upliftpts, sf))
    if 'delrt' in navline.columns:
        if navline['lat'].max() > 50:
            interval = int((AKbathyd['z'].max() - AKbathyd['z'].min()) / 5)
            os.system("gmt psscale -D0.4/2.5/4/0.4 -Cmultibeam.cpt -B%i -O -K >> %s"%(interval,sf))
    else:
        print ('no AKbathyd')

    if args.append is not None:
        # write whole line to line file
        with open(linfile, 'a') as file:
            file.write('> %s \n' % linename)
            for i,r in tempnav.iterrows():
                file.write('%s %s \n' % (str(r['lon']), str(r['lat'])))
        
        # write shown segment to seg file
        with open(segfile, 'a') as file:
            file.write('> %s \n' % linename)
            for i,r in windnav.iterrows():
                file.write('%s %s \n' % (str(r['lon']), str(r['lat'])))

        # write line label to file
        with open(labfile, 'a') as file:
            file.write('%s %s %s \n' % (str(firstx), str(firsty), str(linename)))

    #os.system("rm Output/%s"%gridname)
    if phc:
        os.system("gmt psscale -D0.4/10/2/0.8 -Cseis2.cpt -B%i -O -K >> %s"%(sinterval, sf))

    '''
    os.system("gmt psbasemap %s %s %s -X10 -O -K --FORMAT_GEO_MAP=D >> %s" % (Jm, Rm, Bm, sf))
    os.system("gmt grdimage hillshade.grd -R -J -O -K -Chill.cpt >> %s"%(sf))
    os.system("gmt grdcontour %s -C1000 -R -J -O -K -W0.1,black >> %s"%(DEM,sf))
    for faultname in os.listdir('library/faults/ndipping'):
        os.system("gmt psxy %s -R -J -W0.3,white,- -O -K >> %s" % ('library/faults/ndipping/%s'%faultname, sf))
    for faultname in os.listdir('library/faults/sdipping'):
        os.system("gmt psxy %s -R -J -W0.3,white -O -K >> %s" % ('library/faults/sdipping/%s'%faultname, sf))
    #os.system("gmt psxy %s -R -J -Sc0.1 -Gwhite -O -K >> %s" % (upliftpts, sf))
    #os.system("gmt pstext %s -R -J -Wwhite -X0.2 -O -K >> %s" % (upliftpts, sf))
    os.system("gmt psxy %s -R -J -W0.1,black -O -K >> %s" % (PSSnav, sf))
    os.system("gmt psxy %s -R -J -W0.1,black -O -K >> %s" % (LWBnav, sf))
    linefile = '%s/temp.txt'%navdir
    tempnav = navline[['lon','lat']]

    tempnav.to_csv(linefile,header=False,index=False, sep = ' ')
    os.system("gmt psxy %s -R -J -W0.5,yellow -O -K >> %s" % (linefile, sf))
    windnav = navline[(navline.cdp > wCmin) & (navline.cdp < wCmax)]
    windnav = windnav[['lon','lat']]
    windnav.to_csv(linefile,header=False,index=False, sep=' ')
    os.system("gmt psxy %s -R -J -W0.5,red -O -K >> %s" % (linefile, sf))
    os.system("echo %s %s %s | gmt pstext -R -J -F+f7p,Helvetica,red -O -K >> %s" % (firstx,firsty,linename,sf))
    os.system("echo %s %s | gmt psxy -R -J -Sc0.1 -Gblue -O >> %s" % (centlon,centlat, sf))
    print (nlon, nlat)
    '''

    plot_legend(legns, legvars, legcolors, legendname)

# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        Windows and plots a line with specified range and VE
        Makes grid file of the line
        
        Arguments:
            -l line number
            -z PSS for Puget sound Sparker or LWB for lake washington boomer. PSC or LWC for chirp
            -p M for MCS line C for chirp line, required if -d and -t are not provided
            -d directory containing tape, required if -p is not provided
            -t name of tape, required if -p is not provided
            -w cdpmin cdpmax tmin tmax (optional, assumes no windowing)
            -v vertical exaggeration (optional, assumes 10x VE)
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-p', '--MorC', dest='MorC', type=str, required=True,
                        help='flag M or C to plot MCS or chirp profile, required if d and t are null')
    parser.add_argument('-d', '--tapedir', dest='tapedir', type=str,
                        help='directory containing tape, required if MorC flag is null')
    parser.add_argument('-n', '--tapename', dest='tapename', type=str,
                        help='name of tape, required if MorC flag is null')
    parser.add_argument('-w', '--window', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax'),
                        dest='window', type=float, nargs=4,
                        help='window [cdpmin cdpmax tmin tmax]')
    parser.add_argument('-b', '--bounds', metavar=('lonmin', 'lonmax', 'latmin', 'latmax'),
                        dest='bounds', type=float, nargs=4,
                        help='map bounds [lonmin lonmax latmin latmax]')
    parser.add_argument('-a', '--append', dest='append', type=str,
                        help='name of existing map to append to')
    parser.add_argument('-sk', '--skew', dest='skew', type=str,
                        help='skew plot?')
    parser.add_argument('-i', '--intrange', metavar=('nearlon', 'nearlat', 'bego', 'endo','tmin','tmax'),
                        dest='intrange', type=float, nargs=6,
                        help='range [bego,endo] centered on nearest point to (nlon,nlat) [nlon nlat bego endo tmin tmax]')
    parser.add_argument('-it', '--intrange_t', metavar=('tmin','tmax'),
                        dest='intrange_t', type=float, nargs=2,
                        help='min/max times to plot - entire profile [tmin tmax]')
    parser.add_argument('-v', '--vertexg', dest='vertexg', type=int,
                        help='integer value by which to vertically exaggerate section')
    parser.add_argument('-s', '--saturation', dest='saturation', type=float,
                        help='value at which all amplitudes above are saturated')
    parser.add_argument('-ci', '--ci', dest='ci', type=float,
                        help='interval to plot crossing contours')
                        
    parser.add_argument('-thresh', '--thresh', dest='thresh', type=float,
                        help='threshold to clip water bottom times by. Default is 0.01')
    parser.add_argument('-in', '--interp', dest='interp', type=str,
                        help='enter "-in y" to add interpretations to profile')
    parser.add_argument('-ki', '--keepinterp', dest='keepinterp', type=str,
                        help='enter "-ki y" to write line interpretation points to file')

    parser.add_argument('-sx', '--sxorgx', dest='sx', type=str,
                        help='enter "-sx sx" to read sx, sy as navigation input from segy file rather than gx, gy')
    parser.add_argument('-dem', '--dem', dest='dem', type=str,
                        help='enter "-dem otherdem.grd" to plot DEM other than puget sound')
    parser.add_argument('-utmzone', '--utmzone', dest='utmzone', type=str,
                        help='enter "-utmzone 10" for Puget Sound, "-utmzone 11" for Catalina')
    parser.add_argument('-bs', '--bs', dest='bs', type=str,
                        help='enter e.g. "-bs ENws" for labels on top and right and ticks on bottom and right of map. second two correspond to cdp and time, first two correspond to distance and depth')
    parser.add_argument('-nm', '--nomap', dest='nomap', type=str,
                        help='enter "-nm y" to not include map')
                    
    parser.add_argument('-nf', '--nofaults', dest='nofaults', type=str,
                        help='enter "-nf y" to not include map')
                        
    parser.add_argument('-hf', '--horizonfile', dest='horizonfile', type=str,
                        help='name of csv with minimum columns intcdp, twtt, s that were calculated with this grid')
    parser.add_argument('-hf2', '--horizonfile2', dest='horizonfile2', type=str,
                        help='name of csv with minimum columns intcdp, twtt, s that were calculated with this grid')
    parser.add_argument('-hf3', '--horizonfile3', dest='horizonfile3', type=str,
                        help='name of csv with minimum columns intcdp, twtt, s that were calculated with this grid')
    parser.add_argument('-hf4', '--horizonfile4', dest='horizonfile4', type=str,
                        help='name of csv with minimum columns intcdp, twtt, s that were calculated with this grid')
    parser.add_argument('-phc', '--phc', dest='phc', type=str,
                        help='enter "-phc y" to plot horizons in color and segy in gray')
    parser.add_argument('-as', '--arcseconds', dest='arcseconds', type=str,
                        help='enter "-ac y" to force arcseconds instead of UTM conversion')
    parser.add_argument('-jc', '--justcontours', dest='justcontours', type=str,
                        help='enter "-jc y" to omit plotting grid file below contours')
    parser.add_argument('-nwb', '--nwb', dest='nwb', type=str,
                        help='enter "-nwb y" to not plot water bottom')
    parser.add_argument('-ju', '--justcurve', dest='justcurve', type=str,
                        help='enter "-jc y" to just plot the dip curve')
    parser.add_argument('-onlyapril', '--onlyapril', dest='onlyapril', type=str,
                        help='enter "-onlyapril y" to skip interpretations outside of april19interpretations folder')
    parser.add_argument('-cv', '--cv', dest='cv', type=float,
                        help='enter "-cv [m/s]" to scale by constant velocity. default is 1500 m/s')
    parser.add_argument('-ut', '--usetracr', dest='usetracr', type=str,
                        help='enter "-ut y" to force tracr as rec num')
    parser.add_argument('-nc', '--nc', dest='nc', type=str,
                        help='enter "-nc y" to not include cross contours on profile')
    parser.add_argument('-ab1', '--addbox1', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax', 'trabel'),
                        dest='addbox1', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab2', '--addbox2', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox2', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ppf', '--ppf', dest='ppf', type=str,
                        help='enter "-ppf y" to plot pratt faults')
    parser.add_argument('-oma', '--oma', dest='oma', type=str,
                        help='enter "-oma y" to not plot magnetic anomaly concacts with fault zone')
                        
    parser.add_argument('-jdf', '--jdf', dest='jdf', type=str,
                        help='enter "-jdf y" to only plot deformation front. "-ppf y" flag must also be included')
    parser.add_argument('-pbl', '--pbl', dest='pbl', type=str,
                        help='enter "-pbl y" to plot back thrusts as lines instead of arrows')
    pargs = parser.parse_args()
    
    #cProfile.run('main(pargs)')
    main(pargs)
