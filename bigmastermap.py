#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import math
import argparse
import plotline_fctns as plf
from pyproj import Proj

def round_to_nths(num, n):
    return int(num*n)/n

def main(args):

    sf = args.savefile
    mXmin = args.bounds[0]
    mXmax = args.bounds[1]
    mYmin = args.bounds[2]
    mYmax = args.bounds[3]

    linc = 0.2
    jm2val = 13.0 - linc/13.0
    if args.jm is not None:
        Jm = '-JM%f' % args.jm
    else:
        Jm = '-JM10'
    Jm2 = '-JM%0.4f' % jm2val
    Rm = '-R%f/%f/%f/%f' % (mXmin,mXmax,mYmin,mYmax)
    Rm2 = '-R%f/%f/%f/%f' % (mXmin,mXmax-linc,mYmin,mYmax)
    if args.bs is not None:
        if 'n' in args.bs and 'w' in args.bs and 's' in args.bs and 'e' in args.bs:
            Bm = '-Bx30 -By1 -B%s' % args.bs
        else:
            if mXmax - mXmin > 1:
                Bm = '-Bx0.5 -By0.5 -B%s' % args.bs
            else:
                Bm = '-Bx0.2 -By0.1 -B%s' % args.bs
    else:
        if mXmax - mXmin > 1:
            Bm = '-Bx0.5 -By0.5 -BNWse'
        else:
            Bm = '-Bx0.2 -By0.1 -BNWse'
    if args.dem is not None:
        DEM = args.dem
    else:
        if args.res == 'high':
            DEM = '/Users/ginevramoore/Documents/research/segymanip/library/puget_sound_13_navd88_2014.nc'
        elif args.res == 'low':
            DEM = '/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc'
        else:
            print ('resolution must be specified as "high" or "low"')
            print ('exiting ... ')
            exit()

    DEM2 = '/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc'
    GEO = '/Users/ginevramoore/Documents/research/segymanip/library/geology/geo.grd'
    GEOc = '/Users/ginevramoore/Documents/research/segymanip/library/geology/geo.dat'
    GEOt = '/Users/ginevramoore/Documents/research/segymanip/library/geology/geo.txt'
    PSSnav = 'PSScdpnav.txt'
    LWBnav = 'LWBcdpnav.txt'
    upliftpts = '/Users/ginevramoore/Documents/research/segymanip/library/uplift/tenbrink2006_uplift_supp.txt'
    prattfiles = '/Users/ginevramoore/Documents/research/segymanip/library/lines/allprattlines.txt'
    prattlabels = '/Users/ginevramoore/Documents/research/segymanip/library/lines/allprattlabels.txt'
    oligocene = '/Users/ginevramoore/Documents/research/segymanip/library/geology/olig.dat'
    miocene = '/Users/ginevramoore/Documents/research/segymanip/library/geology/miocene.dat'
    karlincdir = '/Users/ginevramoore/Documents/research/segymanip/library/cores'
    karlinsdir = '/Users/ginevramoore/Documents/research/segymanip/library/sed_ages'

    os.system("gmt gmtset MAP_FRAME_TYPE = plain")
    if args.plot3D is not None:
        os.system("gmt grdgradient %s -A60/60 -Ghillshade.grd"%(DEM))
    else:
        os.system("gmt grdgradient %s -A110/20 -Ghillshade.grd"%(DEM))
    if args.res == 'high':
        os.system("gmt grdgradient %s -A110/20 -Ghillshade2.grd"%('/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc'))
    os.system("gmt makecpt -Cdrywet -I -D -Fr -T-800/600/1 -Z > dep.cpt")
    os.system("gmt makecpt -Cgray -D -Fr -T-20000/20000/10000 -Z > dep3.cpt")
    os.system("gmt makecpt -Cinferno -D -Fr -T-300/0/1 -Z > dep2.cpt")
    os.system("gmt makecpt -Cocean -D -Fr -T-1/1/0.01 -Z > hill.cpt")
    os.system("gmt makecpt -Ctopo -D -Fr -T0/1100/1 > geo.cpt")
    ghayescpt = '/Users/ginevramoore/Documents/GitHub/slab2/slab2code/library/forplotting/ghayes2.cpt'
    gmcpt = '/Users/ginevramoore/Documents/research/segymanip/library/gmtopo.cpt'
    #gmcpt = 'dep2.cpt'
    
    os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
    if args.geology is not None:
        os.system("gmt grdimage %s -R -J -O -K -Cdep2.cpt -Ihillshade.grd >> %s"%(DEM, sf))
        os.system("gmt grdimage %s -R -J -O -K -Cgeo.cpt >> %s"%(GEO, sf))
        os.system("gmt psxy %s -R -J -W0.5,cyan -O -K >> %s" % (GEOc, sf))
        os.system("gmt pstext %s -R -J -O -K -F+f4p,Helvetica,black >> %s" % (GEOt, sf))
    else:
        if args.pic is not None:
            if args.pic == 'y':
                plotincolor = True
            else:
                plotincolor = False
        else:
            plotincolor = False
            
        if plotincolor:
            os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade2.grd >> %s"%(DEM2, gmcpt, sf))
            os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade.grd >> %s"%(DEM, gmcpt, sf))
        else:
            os.system("gmt grdimage %s -R -J -O -K -Cdep3.cpt -Ihillshade2.grd >> %s"%(DEM2, sf))
            os.system("gmt grdimage %s -R -J -O -K -Cdep3.cpt -Ihillshade.grd >> %s"%(DEM, sf))

    if args.pls is not None:
        if args.pls == 'y':
            plotseis = True
        else:
            plotseis = False
    else:
        plotseis = False

    if plotseis:
        dat = pd.read_csv('library/PNSN_1970-2020_seis.csv')
        dat = dat[['Lon','Lat','Depth Km','Magnitude','Evid']]
        dat['Magnitude'] = (dat['Magnitude'].values / 8) * (dat['Magnitude'].values / 8)
        maxval = dat['Depth Km'].max()
        minval = dat['Depth Km'].min()
        dat = dat.sort_values(by=['Magnitude'])
        os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/1 -Z > mag.cpt" % ('purple4,yellow1', minval, maxval))
        dat.to_csv('library/sies_toplot.dat', header=True,index=False,sep=' ')
        os.system("gmt psxy library/sies_toplot.dat -R -J -Sc -Cmag.cpt -O -K >> %s" % (sf))
        intval=10
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/5/2/0.4 -Cmag.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/5/2/0.4 -Cmag.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))

    if args.pma is not None:
        MAG = '/Users/ginevramoore/Documents/research/magnetics/aeromag/pugetsound1.tif'
        #os.system("gmt grdimage %s -R -J -O -K -Cdep.cpt >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L-700/-200 -R -J -O -K -W0.3,red >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L-200/0 -R -J -O -K -W0.3,yellow >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L0/200 -R -J -O -K -W0.3,cyan >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L200/700 -R -J -O -K -W0.3,blue >> %s"%(MAG, sf))

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
        os.system("gmt psxy temp.txt -R -J -W2,gray60 -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox1[4], sf))

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
        os.system("gmt psxy temp.txt -R -J -W2,gray60 -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox2[4], sf))

    if args.addbox3 is not None:
        ax1 = args.addbox3[0]
        ax2 = args.addbox3[1]
        ax3 = args.addbox3[1]
        ax4 = args.addbox3[0]
        ax5 = args.addbox3[0]
        ay1 = args.addbox3[2]
        ay2 = args.addbox3[2]
        ay3 = args.addbox3[3]
        ay4 = args.addbox3[3]
        ay5 = args.addbox3[2]
        
        plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
        plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
        os.system("gmt psxy temp.txt -R -J -W2,gray60 -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox3[4], sf))

    if args.addbox4 is not None:
        ax1 = args.addbox4[0]
        ax2 = args.addbox4[1]
        ax3 = args.addbox4[1]
        ax4 = args.addbox4[0]
        ax5 = args.addbox4[0]
        ay1 = args.addbox4[2]
        ay2 = args.addbox4[2]
        ay3 = args.addbox4[3]
        ay4 = args.addbox4[3]
        ay5 = args.addbox4[2]
        
        plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
        plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
        os.system("gmt psxy temp.txt -R -J -W2,gray60 -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox4[4], sf))

    if args.addbox5 is not None:
        ax1 = args.addbox5[0]
        ax2 = args.addbox5[1]
        ax3 = args.addbox5[1]
        ax4 = args.addbox5[0]
        ax5 = args.addbox5[0]
        ay1 = args.addbox5[2]
        ay2 = args.addbox5[2]
        ay3 = args.addbox5[3]
        ay4 = args.addbox5[3]
        ay5 = args.addbox5[2]
        
        plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
        plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
        os.system("gmt psxy temp.txt -R -J -W2,gray60 -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox5[4], sf))

    if args.addbox6 is not None:
        ax1 = args.addbox6[0]
        ax2 = args.addbox6[1]
        ax3 = args.addbox6[1]
        ax4 = args.addbox6[0]
        ax5 = args.addbox6[0]
        ay1 = args.addbox6[2]
        ay2 = args.addbox6[2]
        ay3 = args.addbox6[3]
        ay4 = args.addbox6[3]
        ay5 = args.addbox6[2]
        
        plotbox = pd.DataFrame({'x':[ax1, ax2, ax3, ax4, ax5], 'y':[ay1, ay2, ay3, ay4, ay5]})
        plotbox.to_csv('temp.txt', header=True, index=False, sep=' ')
        os.system("gmt psxy temp.txt -R -J -W2,black -O -K >> %s" % (sf))
        os.system('printf "%s %s %s" | gmt pstext -R -J -O -K -D-D0.15/0.15 -F+f11p,Helvetica,black >> %s'%(ax2, ay3, args.addbox6[4], sf))
    
    if args.p17 is not None:
        if args.p17 == 'y':
            os.system("gmt psxy %s -R -J -W0.9,black -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/WP17nav.dat', sf))
            os.system("gmt psxy %s -R -J -W0.9,black -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/WL17nav.dat', sf))
    if args.p11 is not None:
        if args.p11 == 'y':
            os.system("gmt psxy %s -R -J -W0.5,yellow -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/WSnav.dat', sf))
    if args.p70 is not None:
        if args.p70 == 'y':
            os.system("gmt psxy %s -R -J -W0.5,black -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/W70nav.dat', sf))
    if args.p95 is not None:
        if args.p95 == 'y':
            os.system("gmt psxy %s -R -J -W0.5,blue4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/W95nav.dat', sf))
    if args.p97 is not None:
        if args.p97 == 'y':
            os.system("gmt psxy %s -R -J -W0.5,purple4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/W97nav.dat', sf))

    if args.puc is not None:
        if args.puc == 'y':
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.3,green4 -GN-4 -L-2/1 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.4,green4 -GN-4 -L1/2 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.5,green4 -GN-4 -L2/3 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.6,green4 -GN-4 -L3/4 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.7,green4 -GN-4 -L4/5 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.8,green4 -GN-4 -L5/6 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.9,green4 -GN-4 -L6/7 >> %s"%(Rm, Jm, sf))
            os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.10,green4 -GN-4 -L7/8 >> %s"%(Rm, Jm, sf))

    if args.oqf is not None:
        if args.oqf != 'y':
            plotquaternaryfaults = True
        else:
            plotquaternaryfaults = False
    else:
        plotquaternaryfaults = True

    if args.jm is not None:
        if int(args.jm) > 15:
            faultwidth = '0.04'
        else:
            faultwidth = '0.04'
    else:
        faultwidth = '0.02'

    if plotquaternaryfaults:
        if args.osf is not None:
            if args.osf == 'y':
                dat = pd.read_csv('/Users/ginevramoore/Documents/research/segymanip/library/pugetlowlandfaults_lonlat.dat', delim_whitespace=True)
                dat = dat[(dat.lat > 47.7) | (dat.lat<47.5)]
                dat.to_csv('temp.txt',header=True,index=False,na_rep=np.nan)
                os.system("gmt psxy %s -R -J -Sc%s -Gred4 -O -K >> %s" % ('temp.txt', faultwidth, sf))
        else:
            os.system("gmt psxy %s -R -J -Sc%s -Gred4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/pugetlowlandfaults_lonlat.dat', faultwidth, sf))
        os.system("gmt psxy %s -R -J -Sc%s -Gred4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/catalinafaults.txt', faultwidth, sf))
    
    nodash = False
    if args.odl is not None:
        if args.odl == 'y':
            nodash = True
    if not nodash:
        os.system("gmt grdcontour %s -L-150/0 -C150 %s %s -O -K -W0.34,black,- >> %s"%(DEM,Rm,Jm, sf))
        os.system("gmt grdcontour %s -C1000000 %s %s -O -K -W0.3,black >> %s"%(DEM,Rm,Jm, sf))
    else:
        os.system("gmt grdcontour %s -C1000000 %s %s -O -K -W0.01,black >> %s"%(DEM,Rm,Jm, sf))

    #os.system("gmt psxy %s -R -J -W0.1,blue2 -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt psxy %s -R -J -Sc0.05 -Gblack -O -K >> %s" % (prattfiles, sf))
    #os.system("gmt pstext %s -R -J -O -K -F+f6p,Helvetica,blue2 >> %s" % (prattlabels, sf))
    
    if args.linfile is not None:
        os.system("gmt psxy %s -R -J -W0.9,red -O -K >> %s" % (args.linfile, sf))
    if args.segfile is not None:
        os.system("gmt psxy %s -R -J -W0.9,red -O -K >> %s" % (args.segfile, sf))
    if args.dotfile is not None:
        os.system("gmt psxy %s -R -J -Sc0.1 -Gblue -O -K >> %s" % (args.dotfile, sf))
    if args.labfile is not None:
        os.system("gmt pstext %s -R -J -F+f12p,Helvetica,red -O -K >> %s" % (args.labfile,sf))

    if args.pqf is not None:
        if args.pqf == 'y':
            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/folds'):
                if 'brocher' not in faultname:
                    continue
                print('/Users/ginevramoore/Documents/research/segymanip/library/faults/folds/%s'%faultname)
                os.system("gmt psxy %s -R -J -W0.7,violet -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/folds/%s'%faultname, sf))

            
    if args.valfile is not None:
        tcpt = args.valfile[3]
        minval = args.valfile[1]
        maxval = args.valfile[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        if args.plotvec is not None:
            if args.plotvec == 'y':
                plotvec = True
            else:
                plotvec = False
        else:
            plotvec = False

        if plotvec:
            strdata = pd.read_csv(args.valfile[0], delim_whitespace=True, names=['x','y','val','str'],skiprows=1)

            strdata['len'] = '0.015i'
            cptdat = pd.read_csv('dep.cpt', names=['v1', 'rgb1', 'v2', 'rgb2'], delim_whitespace=True)
            for index, row in cptdat.iterrows():
                v1 = row['v1']
                v2 = row['v2']
                rgb1 = row['rgb1']
                try:
                    if float(v1) >= float(args.valfile[2]) or float(v2) >= float(args.valfile[2]):
                        plotdat = strdata[(strdata.val >= float(v1))]
                        print ('more than',v1, v2, args.valfile[1], args.valfile[2],len(plotdat))
                    elif float(v1) <= float(args.valfile[1]) or float(v2) <= float(args.valfile[1]):
                        plotdat = strdata[(strdata.val < float(v2))]
                        print ('less than',v1, v2, args.valfile[1], args.valfile[2],len(plotdat))
                    else:
                        plotdat = strdata[(strdata.val >= float(v1)) & (strdata.val < float(v2))]
                        print ('between',v1, v2, args.valfile[1], args.valfile[2],len(plotdat))
                    plotdat = plotdat[['x','y','str','len']]
                    plotdat['str'] = plotdat['str'].values-90
                    plotdat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    plotdat['str'] = plotdat['str'].values+180
                    plotdat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    #print (plotdat)
                    os.system("gmt psxy aztest1.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                    os.system("gmt psxy aztest2.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                except:
                    print (row)

            '''
            strdata['str'] = strdata['str'].values-90
            strdata.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            strdata['str'] = strdata['str'].values+180
            strdata.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)

            os.system("gmt psxy aztest1.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            os.system("gmt psxy aztest2.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            '''

        else:
            if mXmax - mXmin > 1:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile[0], sf)) # for sediment accumulation plots
            else:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile[0], sf)) # for sediment accumulation plots

        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile[0], sf))
        intval = (float(maxval)-float(minval))/5.0

        plothorizontal = False
        if args.valfile2 is None and args.valfile3 is None and args.valfile4 is None and args.valfile5 is None:
            plothorizontal = True
        if args.forcehorizontal is not None:
            if args.forcehorizontal == 'y':
                plothorizontal = True

        if args.ns is not None:
            if args.ns != 'y':
                if plothorizontal:
                    os.system('gmt psscale -D0/-0.5+w%i/0.5+h -Cdep.cpt -B%f -O -K >> %s'%(float(Jm[3:]), intval, sf))
                else:
                    os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
        else:
            if plothorizontal:
                os.system('gmt psscale -D0/-0.5+w%i/0.5+h -Cdep.cpt -B%f -O -K >> %s'%(float(Jm[3:]), intval, sf))
            else:
                os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))

    if args.valfile2 is not None:
        tcpt = args.valfile2[3]
        minval = args.valfile2[1]
        maxval = args.valfile2[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        if args.plotvec is not None:
            if args.plotvec == 'y':
                plotvec = True
            else:
                plotvec = False
        else:
            plotvec = False

        if plotvec:
            strdata = pd.read_csv(args.valfile2[0], delim_whitespace=True, names=['x','y','val','str'],skiprows=1)

            strdata['len'] = '0.015i'
            cptdat = pd.read_csv('dep.cpt', names=['v1', 'rgb1', 'v2', 'rgb2'], delim_whitespace=True)
            for index, row in cptdat.iterrows():
                v1 = row['v1']
                v2 = row['v2']
                rgb1 = row['rgb1']
                try:
                    if float(v1) >= float(args.valfile2[2]) or float(v2) >= float(args.valfile2[2]):
                        plotdat = strdata[(strdata.val >= float(v1))]
                        print ('more than',v1, v2, args.valfile2[1], args.valfile2[2],len(plotdat))
                    elif float(v1) <= float(args.valfile2[1]) or float(v2) <= float(args.valfile2[1]):
                        plotdat = strdata[(strdata.val < float(v2))]
                        print ('less than',v1, v2, args.valfile2[1], args.valfile2[2],len(plotdat))
                    else:
                        plotdat = strdata[(strdata.val >= float(v1)) & (strdata.val < float(v2))]
                        print ('between',v1, v2, args.valfile2[1], args.valfile2[2],len(plotdat))
                    plotdat = plotdat[['x','y','str','len']]
                    plotdat['str'] = plotdat['str'].values-90
                    plotdat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    plotdat['str'] = plotdat['str'].values+180
                    plotdat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    #print (plotdat)
                    os.system("gmt psxy aztest1.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                    os.system("gmt psxy aztest2.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                except:
                    print (row)

            '''
            strdata['str'] = strdata['str'].values-90
            strdata.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            strdata['str'] = strdata['str'].values+180
            strdata.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)

            os.system("gmt psxy aztest1.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            os.system("gmt psxy aztest2.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            '''

        else:
            if mXmax - mXmin > 1:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile2[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile2[0], sf)) # for sediment accumulation plots
            else:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile2[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile2[0], sf)) # for sediment accumulation plots

        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile2[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+1.95, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+1.95, intval, sf))
            
    if args.valfile3 is not None:
        tcpt = args.valfile3[3]
        minval = args.valfile3[1]
        maxval = args.valfile3[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        if args.plotvec is not None:
            if args.plotvec == 'y':
                plotvec = True
            else:
                plotvec = False
        else:
            plotvec = False

        if plotvec:
            strdata = pd.read_csv(args.valfile3[0], delim_whitespace=True, names=['x','y','val','str'],skiprows=1)

            strdata['len'] = '0.015i'
            cptdat = pd.read_csv('dep.cpt', names=['v1', 'rgb1', 'v2', 'rgb2'], delim_whitespace=True)
            for index, row in cptdat.iterrows():
                v1 = row['v1']
                v2 = row['v2']
                rgb1 = row['rgb1']
                try:
                    if float(v1) >= float(args.valfile3[2]) or float(v2) >= float(args.valfile3[2]):
                        plotdat = strdata[(strdata.val >= float(v1))]
                        print ('more than',v1, v2, args.valfile3[1], args.valfile3[2],len(plotdat))
                    elif float(v1) <= float(args.valfile3[1]) or float(v2) <= float(args.valfile3[1]):
                        plotdat = strdata[(strdata.val < float(v2))]
                        print ('less than',v1, v2, args.valfile3[1], args.valfile3[2],len(plotdat))
                    else:
                        plotdat = strdata[(strdata.val >= float(v1)) & (strdata.val < float(v2))]
                        print ('between',v1, v2, args.valfile3[1], args.valfile3[2],len(plotdat))
                    plotdat = plotdat[['x','y','str','len']]
                    plotdat['str'] = plotdat['str'].values-90
                    plotdat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    plotdat['str'] = plotdat['str'].values+180
                    plotdat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    #print (plotdat)
                    os.system("gmt psxy aztest1.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                    os.system("gmt psxy aztest2.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                except:
                    print (row)

            '''
            strdata['str'] = strdata['str'].values-90
            strdata.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            strdata['str'] = strdata['str'].values+180
            strdata.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)

            os.system("gmt psxy aztest1.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            os.system("gmt psxy aztest2.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            '''

        else:
            if mXmax - mXmin > 1:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile3[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile3[0], sf)) # for sediment accumulation plots
            else:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile3[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile3[0], sf)) # for sediment accumulation plots

        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile3[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+3.55, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+3.55, intval, sf))
            
    if args.valfile4 is not None:
        tcpt = args.valfile4[3]
        minval = args.valfile4[1]
        maxval = args.valfile4[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        if args.plotvec is not None:
            if args.plotvec == 'y':
                plotvec = True
            else:
                plotvec = False
        else:
            plotvec = False

        if plotvec:
            strdata = pd.read_csv(args.valfile4[0], delim_whitespace=True, names=['x','y','val','str'],skiprows=1)

            strdata['len'] = '0.015i'
            cptdat = pd.read_csv('dep.cpt', names=['v1', 'rgb1', 'v2', 'rgb2'], delim_whitespace=True)
            for index, row in cptdat.iterrows():
                v1 = row['v1']
                v2 = row['v2']
                rgb1 = row['rgb1']
                try:
                    if float(v1) >= float(args.valfile4[2]) or float(v2) >= float(args.valfile4[2]):
                        plotdat = strdata[(strdata.val >= float(v1))]
                        print ('more than',v1, v2, args.valfile4[1], args.valfile4[2],len(plotdat))
                    elif float(v1) <= float(args.valfile4[1]) or float(v2) <= float(args.valfile4[1]):
                        plotdat = strdata[(strdata.val < float(v2))]
                        print ('less than',v1, v2, args.valfile4[1], args.valfile4[2],len(plotdat))
                    else:
                        plotdat = strdata[(strdata.val >= float(v1)) & (strdata.val < float(v2))]
                        print ('between',v1, v2, args.valfile4[1], args.valfile4[2],len(plotdat))
                    plotdat = plotdat[['x','y','str','len']]
                    plotdat['str'] = plotdat['str'].values-90
                    plotdat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    plotdat['str'] = plotdat['str'].values+180
                    plotdat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    #print (plotdat)
                    os.system("gmt psxy aztest1.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                    os.system("gmt psxy aztest2.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                except:
                    print (row)

            '''
            strdata['str'] = strdata['str'].values-90
            strdata.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            strdata['str'] = strdata['str'].values+180
            strdata.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)

            os.system("gmt psxy aztest1.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            os.system("gmt psxy aztest2.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            '''

        else:
            if mXmax - mXmin > 1:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile4[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile4[0], sf)) # for sediment accumulation plots
            else:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile4[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile4[0], sf)) # for sediment accumulation plots

        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile4[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+5.17, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+5.17, intval, sf))
            
    if args.valfile5 is not None:
        tcpt = args.valfile5[3]
        minval = args.valfile5[1]
        maxval = args.valfile5[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        
        
        if args.plotvec is not None:
            if args.plotvec == 'y':
                plotvec = True
            else:
                plotvec = False
        else:
            plotvec = False

        if plotvec:
            strdata = pd.read_csv(args.valfile5[0], delim_whitespace=True, names=['x','y','val','str'],skiprows=1)

            strdata['len'] = '0.015i'
            cptdat = pd.read_csv('dep.cpt', names=['v1', 'rgb1', 'v2', 'rgb2'], delim_whitespace=True)
            for index, row in cptdat.iterrows():
                v1 = row['v1']
                v2 = row['v2']
                rgb1 = row['rgb1']
                try:
                    if float(v1) >= float(args.valfile5[2]) or float(v2) >= float(args.valfile5[2]):
                        plotdat = strdata[(strdata.val >= float(v1))]
                        print ('more than',v1, v2, args.valfile5[1], args.valfile5[2],len(plotdat))
                    elif float(v1) <= float(args.valfile5[1]) or float(v2) <= float(args.valfile5[1]):
                        plotdat = strdata[(strdata.val < float(v2))]
                        print ('less than',v1, v2, args.valfile5[1], args.valfile5[2],len(plotdat))
                    else:
                        plotdat = strdata[(strdata.val >= float(v1)) & (strdata.val < float(v2))]
                        print ('between',v1, v2, args.valfile5[1], args.valfile5[2],len(plotdat))
                    plotdat = plotdat[['x','y','str','len']]
                    plotdat['str'] = plotdat['str'].values-90
                    plotdat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    plotdat['str'] = plotdat['str'].values+180
                    plotdat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
                    #print (plotdat)
                    os.system("gmt psxy aztest1.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                    os.system("gmt psxy aztest2.txt -SV0.02i -W0.03,%s -R -J -O -K >> %s"%(rgb1,sf))
                except:
                    print (row)

            '''
            strdata['str'] = strdata['str'].values-90
            strdata.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            strdata['str'] = strdata['str'].values+180
            strdata.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)

            os.system("gmt psxy aztest1.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            os.system("gmt psxy aztest2.txt -SV0.01i -Cdep.cpt -R -J -O -K >> %s"%(sf))
            '''

        else:
            if mXmax - mXmin > 1:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile5[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile5[0], sf)) # for sediment accumulation plots
            else:
                os.system("gmt psxy %s -R -J -Sc0.09 -Cdep.cpt -O -K >> %s" % (args.valfile5[0], sf)) # for apparent dip plots
                #os.system("gmt psxy %s -R -J -Sc0.2 -Cdep.cpt -O -K >> %s" % (args.valfile5[0], sf)) # for sediment accumulation plots

        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile5[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+6.77, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/2/0.4 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+6.77, intval, sf))
            
    # uncomment to plot segment and line files over the val files
    
    #if args.linfile is not None:
    #    os.system("gmt psxy %s -R -J -W0.9,black -O -K >> %s" % (args.linfile, sf))
    #if args.segfile is not None:
    #    os.system("gmt psxy %s -R -J -W0.9,black -O -K >> %s" % (args.segfile, sf))
    if args.dotfile is not None:
        os.system("gmt psxy %s -R -J -Sc0.1 -Gblue -O -K >> %s" % (args.dotfile, sf))
    if args.labfile is not None:
        os.system("gmt pstext %s -R -J -F+f12p,Helvetica,red -O -K >> %s" % (args.labfile,sf))
           
    if args.jm is not None:
        if int(args.jm) > 15:
            geosize = '0.9'
        else:
            geosize = '0.5'
    else:
        geosize = '0.5'

    if args.omo is not None:
        if args.omo == 'y':
            omitgeology = True
        else:
            omitgeology = False
    else:
        omitgeology = False

    if not omitgeology:
        os.system("gmt psxy %s -R -J -W%s,darkred -O -K >> %s" % (oligocene, geosize, sf))
        os.system("gmt psxy %s -R -J -W%s,darkorange -O -K >> %s" % (miocene, geosize, sf))
        geodir = 'library/geology/contacts'
        for geo in os.listdir(geodir):
            if geo[0] == 'E':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,purple3 -O -K >> %s" % (geodir, geo, geosize, sf))
            if geo[0] == 'M':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,darkorange2 -O -K >> %s" % (geodir, geo, geosize, sf))
            if geo[0] == 'O':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,red -O -K >> %s" % (geodir, geo, geosize, sf))
            if geo[0] == 'T':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,purple -O -K >> %s" % (geodir, geo, geosize, sf))
            if geo[0] == 'K':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,pink -O -K >> %s" % (geodir, geo, geosize, sf))
            if geo[0] == 'J':
                dat = pd.read_csv('%s/%s' % (geodir, geo), names=['lon','lat'])
                os.system("gmt psxy %s/%s -R -J -W%s,green -O -K >> %s" % (geodir, geo, geosize, sf))

    if args.llabel is not None:
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBL+f14p,Helvetica,black >> %s"%(args.llabel, sf))
    if args.lrabel is not None:
        os.system("echo %s | gmt pstext -R -J -O -K -F+cBR+f14p,Helvetica,black >> %s"%(args.lrabel, sf))
    if args.tlabel is not None:
        os.system("echo %s | gmt pstext -R -J -O -K -F+cTL+f14p,Helvetica,black >> %s"%(args.tlabel, sf))
    if args.trabel is not None:
        os.system("echo %s | gmt pstext -R -J -O -K -F+cTR+f14p,Helvetica,black >> %s"%(args.trabel, sf))

    if args.ppn is not None:
        if args.ppn == 'y':
            os.system("gmt pstext library/labels/sfzmap.txt -R -J -O -K >> %s" % (sf))

    if args.pup is not None:
        if args.pup == 'y':
            os.system("gmt psxy %s -R -J -Sc0.1 -Ggreen4 -O -K >> %s" % (upliftpts, sf))
    if args.psd100 is not None:
        if args.psd100 == 'y':
            dat = pd.read_csv('/Users/ginevramoore/Documents/research/segymanip/library/geology/attitude/100k_att')

            dat['len'] = '0.05i'
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest0.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            dat['PLOT_ANG'] = dat['PLOT_ANG'].values-90
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            dat['PLOT_ANG'] = dat['PLOT_ANG'].values-90
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
    
            os.system('gmt psxy aztest0.txt -Sv0.05i -W01p,black -R -J -O -K  >> %s' % (sf))
            os.system('gmt psxy aztest1.txt -Sv0.05i -W01p,black -R -J -O -K  >> %s' % (sf))
            os.system('gmt psxy aztest2.txt -Sv0.05i -W01p,black -R -J -O -K  >> %s' % (sf))

    if args.psd24 is not None:
        if args.psd24 == 'y':
            dat = pd.read_csv('/Users/ginevramoore/Documents/research/segymanip/library/geology/attitude/24k_att_nodesc')

            dat = dat[(dat.ATTUD_CD == 1) | (dat.ATTUD_CD == 2) | (dat.ATTUD_CD == 3) | (dat.ATTUD_CD == 8) | (dat.ATTUD_CD == 13) | (dat.ATTUD_CD == 14) | (dat.ATTUD_CD == 15) | (dat.ATTUD_CD == 71) | (dat.ATTUD_CD == 72) | (dat.ATTUD_CD == 75) | (dat.ATTUD_CD == 97) | (dat.ATTUD_CD == 98)]
            dat['len'] = '0.05i'
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest0.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            dat['PLOT_ANG'] = dat['PLOT_ANG'].values-90
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest1.txt', header=False, index=False, sep=' ', na_rep=np.nan)
            dat['PLOT_ANG'] = dat['PLOT_ANG'].values-90
            dat = dat[['lon','lat','PLOT_ANG','len']]
            dat.to_csv('aztest2.txt', header=False, index=False, sep=' ', na_rep=np.nan)
    
            os.system('gmt psxy aztest0.txt -Sv0.03i -W0.02p -R -J -O -K  >> %s' % (sf))
            os.system('gmt psxy aztest1.txt -Sv0.02i -W0.02p -R -J -O -K  >> %s' % (sf))
            os.system('gmt psxy aztest2.txt -Sv0.03i -W0.02p -R -J -O -K  >> %s' % (sf))

    if args.pkc is not None:
        if args.pkc == 'y':
            for filename in os.listdir(karlincdir):
                dat = pd.read_csv('%s/%s' % (karlincdir, filename))
                dat.to_csv('temp.txt', header=False, index=False, sep=' ')
                os.system("gmt psxy %s -R -J -Sc0.1 -W0.01,red4 -O -K >> %s" % ('temp.txt', sf))
    if args.psf is not None:
        if args.pkc == 'y':
            for filename in os.listdir(karlinsdir):
                dat = pd.read_csv('%s/%s' % (karlinsdir, filename))
                dat.to_csv('temp.txt', header=False, index=False, sep=' ')
                os.system("gmt psxy %s -R -J -W0.03,orange4 -O -K >> %s" % ('temp.txt', sf))

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
            '''
            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/sdipping'):
                if 'nelson' not in faultname:
                    continue
                os.system("gmt psxy %s -R -J -W1.0,red4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/sdipping/%s'%faultname, sf))
                print ('/Users/ginevramoore/Documents/research/segymanip/library/faults/sdipping/%s'%faultname)
            
            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/deffronts'):
                if 'prat' not in faultname:
                    continue
                print('/Users/ginevramoore/Documents/research/segymanip/library/faults/deffronts/%s'%faultname)
                os.system("gmt psxy %s -R -J -W1.0,red4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/deffronts/%s'%faultname, sf))
            
            for faultname in os.listdir('/Users/ginevramoore/Documents/research/segymanip/library/faults/strikeslip'):
                if 'ss1' not in faultname:
                    continue
                print('/Users/ginevramoore/Documents/research/segymanip/library/faults/strikeslip/%s'%faultname)
                os.system("gmt psxy %s -R -J -W0.7,red4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/faults/strikeslip/%s'%faultname, sf))

            '''

    #if args.geology is None:
    #    os.system("gmt psscale -D18/2.5/4/0.4 -C%s -B1000 -O -K >> %s"%(gmcpt, sf))
    plotsx = (mXmax+mXmin) / 2.0
    if (mYmax-mYmin) > 1 and (mYmax-mYmin) < 2:
        plotsy = mYmin - ((mYmax-mYmin) / 15)
    elif (mYmax-mYmin) <= 1:
        plotsy = mYmin - ((mYmax-mYmin) / 10)
    else:
        plotsy = mYmin - ((mYmax-mYmin) / 20)
    if (mXmax - mXmin) < 0.1:
        distance = 1
    elif (mXmax - mXmin) < 0.2:
        distance = 2
    elif (mXmax - mXmin) < 0.4:
        distance = 4
    elif (mXmax - mXmin) < 0.8:
        distance = 8
    elif (mXmax - mXmin) < 1:
        distance = 10
    elif (mXmax - mXmin) < 2:
        distance = 20
    else:
        distance = 50

    print ("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -D1/0.25 >> %s " % (plotsx, plotsy, plotsy, distance, sf))

    if args.riv is not None:
        if args.riv == 'y':
            riv = True
        else:
            riv = False
    else:
        riv = False
        
    if riv:
        os.system("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -Ia/3 >> %s " % (plotsx, plotsy, plotsy, distance, sf))
        riverfile = '/Users/ginevramoore/Documents/classes/OCEAN541/Paper_Presentation/NHD_H_Washington_State_GDB/NHD_H_Washington_State_GDB/NHDFlowline_PS_46006.gmt'
        os.system("gmt psxy %s -R -J -W0.5,blue -O -K >> %s" % (riverfile, sf))
    elif not nodash:
        os.system("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -I1/0.25 >> %s " % (plotsx, plotsy, plotsy, distance, sf))
    else:
        print ('not plotting pscoast')

    if args.put is not None:
        if args.put == 'y':
            os.system("gmt pstext %s -R -J -W0.5,green4 -X0.2 -O -K >> %s" % (upliftpts, sf))

    mXmin1, mYmin1, mXmax1, mYmax1 = mXmin, mYmin, mXmax, mYmax
    Jm = '-JM10'

    if args.dem is not None:
        if args.valfile is not None:
            dataval = pd.read_csv(args.valfile[0], delim_whitespace=True)
            Rm = '-R%f/%f/%f/%f' % (dataval['lon'].min()-3,dataval['lon'].max()+3,dataval['lat'].min()-3,dataval['lat'].max()+3)
    else:
        Rm = '-R%f/%f/%f/%f' % (-123.3,-121.7,47,48.5)
    if args.bs is not None:
        if 'n' in args.bs and 'w' in args.bs and 's' in args.bs and 'e' in args.bs:
            Bm = '-Bx30 -By1 -B%s' % args.bs
        else:
            Bm = '-Bx1 -By1 -B%s' % args.bs
    else:
        Bm = '-Bx1 -By1 -BNWse'
    #DEM2 = '/Users/ginevramoore/Documents/research/segymanip/library/GMRTv3_7_20200421topo.grd'
    os.system("gmt grdgradient %s -A60/60 -Ghillshade2.grd"%(DEM2))

    sf = '%s_inset.ps'%args.savefile[:-3]
    lons = [mXmin, mXmin, mXmax, mXmax, mXmin]
    lats = [mYmin, mYmax, mYmax, mYmin, mYmin]
    dat = pd.DataFrame({'lon':lons, 'lat':lats})
    dat = dat[['lon','lat']]
    dat.to_csv('temp.txt', header=False,index=False,sep=' ')
    os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
    
    '''
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.1,green4 -L-2/1 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.2,green4 -L1/2 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.3,green4 -L2/3 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.4,green4 -L3/4 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.5,green4 -L4/5 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.6,green4 -L5/6 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.7,green4 -L6/7 >> %s"%(Rm, Jm, sf))
    os.system("gmt grdcontour upliftE_0.001.grd -C1 -A2+an %s %s -O -K -W0.8,green4 -L7/8 >> %s"%(Rm, Jm, sf))
    '''
    #os.system("gmt grdcontour %s -C1000 %s %s -O -K -W0.6,black >> %s"%(DEM,Rm,Jm, sf))
    #os.system("gmt psxy library/pnwcoasts.con -R -J -W0.5,black -O -K >> %s" % (sf))
    os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade2.grd >> %s"%(DEM2,gmcpt, sf))
    os.system("gmt grdcontour %s -C10000 %s %s -O -K -W0.1,black >> %s"%('/Users/ginevramoore/Documents/research/segymanip/library/nw_pacific_crm_v1.nc',Rm,Jm, sf))
    if args.valfile is not None:
        tcpt = args.valfile[3]
        minval = args.valfile[1]
        maxval = args.valfile[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        os.system("gmt psxy %s -R -J -Sc0.01 -Cdep.cpt -O -K >> %s" % (args.valfile[0], sf))
        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
                if plotseis:
                    os.system("gmt psscale -D%0.1f/5/3/0.5 -Cmag.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
            if plotseis:
                os.system("gmt psscale -D%0.1f/5/3/0.5 -Cmag.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
                    
    if args.valfile2 is not None:
        tcpt = args.valfile2[3]
        minval = args.valfile2[1]
        maxval = args.valfile2[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        os.system("gmt psxy %s -R -J -Sc0.01 -Cdep.cpt -O -K >> %s" % (args.valfile2[0], sf))
        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile2[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
            
    if args.valfile3 is not None:
        tcpt = args.valfile3[3]
        minval = args.valfile3[1]
        maxval = args.valfile3[2]
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        os.system("gmt psxy %s -R -J -Sc0.01 -Cdep.cpt -O -K >> %s" % (args.valfile3[0], sf))
        #os.system("gmt psxy %s -R -J -Sc0.05 -Gred@10 -O -K >> %s" % (args.valfile3[0], sf))
        intval = (float(maxval)-float(minval))/5.0
        if args.ns is not None:
            if args.ns != 'y':
                os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
        else:
            os.system("gmt psscale -D%0.1f/1.5/3/0.5 -Cdep.cpt -B%f -O -K >> %s"%(float(Jm[3:])+0.15, intval, sf))
            
    if args.pma is not None:
        MAG = '/Users/ginevramoore/Documents/research/magnetics/aeromag/pugetsound1.tif'
        #os.system("gmt grdimage %s -R -J -O -K -Cdep.cpt >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L-700/-200 -R -J -O -K -W0.3,red4 >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L-200/0 -R -J -O -K -W0.3,orange4 >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L0/200 -R -J -O -K -W0.3,yellow4 >> %s"%(MAG, sf))
        os.system("gmt grdcontour %s -C50 -A100 -L200/700 -R -J -O -K -W0.3,green4 >> %s"%(MAG, sf))

    os.system("gmt psxy temp.txt -R -J -W3.5,red -O -K >> %s" % (sf))
    os.system("gmt psxy %s -R -J -Sc0.04 -Gred4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/pugetlowlandfaults_lonlat.dat', sf))
    os.system("gmt psxy %s -R -J -Sc0.04 -Gred4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/catalinafaults.txt', sf))
    os.system("gmt psxy temp.txt -R -J -W4,red -O -K >> %s" % (sf))

    '''
    if args.valfile is not None:
        if '.grd' not in args.valfile:
            dataval = pd.read_csv(args.valfile[0], delim_whitespace=True)
            mXmin, mXmax, mYmin, mYmax = dataval['lon'].min()-3,dataval['lon'].max()+3,dataval['lat'].min()-3,dataval['lat'].max()+3
    else:
        mXmin, mXmax, mYmin, mYmax = -123.5,-121.5,46.8,49
    '''

    mXmin, mXmax, mYmin, mYmax = -123.5,-121.5,46.8,49

    plotsx = (mXmax+mXmin) / 2.0
    if (mYmax-mYmin) > 1 and (mYmax-mYmin) < 2:
        plotsy = mYmin - ((mYmax-mYmin) / 350)
    elif (mYmax-mYmin) <= 1:
        plotsy = mYmin - ((mYmax-mYmin) / 300)
    else:
        plotsy = mYmin - ((mYmax-mYmin) / 400)

    if (mXmax - mXmin) < 0.1:
        distance = 1
    elif (mXmax - mXmin) < 0.2:
        distance = 2
    elif (mXmax - mXmin) < 0.4:
        distance = 4
    elif (mXmax - mXmin) < 0.8:
        distance = 8
    elif (mXmax - mXmin) < 1:
        distance = 10
    elif (mXmax - mXmin) < 2:
        distance = 20
    elif (mXmax - mXmin) < 3:
        distance = 30
    elif (mXmax - mXmin) < 4:
        distance = 40
    elif (mXmax - mXmin) < 5:
        distance = 50
    elif (mXmax - mXmin) < 6:
        distance = 60
    elif (mXmax - mXmin) < 7:
        distance = 70
    elif (mXmax - mXmin) < 8:
        distance = 80
    elif (mXmax - mXmin) < 9:
        distance = 90
    elif (mXmax - mXmin) < 10:
        distance = 100
    else:
        distance = 200

    print ("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -D1/0.25 >> %s " % (plotsx, plotsy, plotsy, distance, sf))
    os.system("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -I1/0.25 >> %s " % (plotsx, plotsy, plotsy, distance, sf))


    # bigger inset
    DEM2 = '/Users/ginevramoore/Documents/research/segymanip/library/GMRTv3_7_20200421topo.grd'
    os.system("gmt grdgradient %s -A60/60 -Ghillshade2.grd"%(DEM2))
    Jm = '-JM5'

    bigin = False
    if args.bigin is not None:
        if args.bigin == 'y':
            bigin = True

    #Rm = '-R%f/%f/%f/%f' % (-155,-100,25,65)
    mXXmin, mXXmax, mYYmin, mYYmax = -131,  -117.2, 42.1, 51.4
    Rm = '-R%f/%f/%f/%f' % (mXXmin, mXXmax, mYYmin, mYYmax)
    if args.bs is not None:
        if 'n' in args.bs and 'w' in args.bs and 's' in args.bs and 'e' in args.bs:
            Bm = '-Bx30 -By1 -B%s' % args.bs
        else:
            Bm = '-Bx1 -By1 -B%s' % args.bs
    else:
        Bm = '-Bx1 -By1 -BNWse'
    Bm = '-Bx200 -By200 -Bnwse'
    sf = '%s_inset_big.ps'%args.savefile[:-3]

    plotsx = (mXXmax+mXXmin) / 2.0
    if (mYYmax-mYYmin) > 1 and (mYYmax-mYYmin) < 2:
        plotsy = mYYmin - ((mYYmax-mYYmin) / 25)
    elif (mYYmax-mYYmin) <= 1:
        plotsy = mYYmin - ((mYYmax-mYYmin) / 20)
    else:
        plotsy = mYYmin - ((mYYmax-mYYmin) / 30)
    distance = 200

    if not bigin:
        lons = [mXmin, mXmin, mXmax, mXmax, mXmin]
        lats = [mYmin, mYmax, mYmax, mYmin, mYmin]
    else:
        lons = [args.bounds[0], args.bounds[0], args.bounds[1], args.bounds[1], args.bounds[0]]
        lats = [args.bounds[2], args.bounds[3], args.bounds[3], args.bounds[2], args.bounds[2]]

    dat = pd.DataFrame({'lon':lons, 'lat':lats})
    dat = dat[['lon','lat']]
    dat.to_csv('temp.txt', header=False,index=False,sep=' ')
    os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
    #os.system("gmt grdimage %s -R -J -O -K -C%s >> %s"%('/Users/ginevramoore/Documents/research/Slab2/slab2code/library/forplotting/world_lr.grd', '/Users/ginevramoore/Documents/research/Slab2/slab2code/library/forplotting/ghayes2.cpt', sf))
    os.system("gmt grdimage %s -R -J -O -K -C%s -Ihillshade2.grd >> %s"%(DEM2,gmcpt, sf))

    os.system("gmt pscoast -Dc -R -J -O -K -ENA+p0.25p,red+r -EUS.%s+p0.05p,navy+r >> %s" % ('OR', sf))
    os.system("gmt pscoast -Dc -R -J -O -K -ENA+p0.25p,red+r -EUS.%s+p0.05p,navy+r >> %s" % ('WA', sf))
    os.system("gmt pscoast -Dc -R -J -O -K -ENA+p0.25p,red+r -ECA.%s+p0.05p,navy+r >> %s" % ('BC', sf))
    #os.system("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -O -K -Df -I3/0.025 >> %s " % (plotsx, plotsy, plotsy, distance, sf))
    gpsdat = pd.read_csv('library/gps_mccaffrey_2013.txt', delim_whitespace=True)
    gpsdat = gpsdat[['Longitude','Latitude','Ve','Vn','Se','Sn','NEcor.']]
    gpsdat = gpsdat[['Longitude','Latitude','Ve','Vn']]
    gpsdat.to_csv('gpstemp.txt',header=False,index=False,sep=' ')
    os.system("gmt psvelo gpstemp.txt -R -J -Gblack -Se0.02/1/11 -V -O -K >> %s" % (sf))
    os.system("gmt psxy temp.txt -R -J -W2,red -O -K >> %s" % (sf))
    '''
    os.system("gmt psxy %s -R -J -Sc0.02 -Gyellow4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/pugetlowlandfaults_lonlat.dat', sf))
    os.system("gmt psxy %s -R -J -Sc0.02 -Gyellow4 -O -K >> %s" % ('/Users/ginevramoore/Documents/research/segymanip/library/catalinafaults.txt', sf))
    '''
    #os.system("gmt pscoast -R -J -Df -O -K -I1/0.25 >> %s " % (sf))

    #os.system("gmt pscoast -R -J -Di -O -K -W0.1,black >> %s " % (sf))
    os.system("gmt grdcontour %s -C1000000 %s %s -O -K -W0.1,black >> %s"%(DEM2,Rm,Jm, sf))

    #os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))

    if args.plot3D is not None:
        pval = args.plot3D
        mXmin = args.bounds[0]
        mXmax = args.bounds[1]
        mYmin = args.bounds[2]
        mYmax = args.bounds[3]
        minval = float(args.valfile[2]) * -750.0
        maxval = float(args.valfile[1]) * -750.0
        if args.jm is not None:
            Jm = '-JM%f' % args.jm
            jm2 = args.jm / 2.0
            Jz = '-JZ%fc' % jm2
        else:
            Jm = '-JM10'
            Jz = 'JZ5c'
        Jm2 = '-JM%0.4f' % jm2val
        Rm = '-R%f/%f/%f/%f/%f/%f' % (mXmin,mXmax,mYmin,mYmax,float(minval),float(maxval))
        Rm2 = '-R%f/%f/%f/%f' % (mXmin,mXmax-linc,mYmin,mYmax)
        '''
        if args.bs is not None:
            if 'n' in args.bs and 'w' in args.bs and 's' in args.bs and 'e' in args.bs:
                Bm = '-Bx30 -By1 -Bz0.1 -B%s' % args.bs
            else:
                if mXmax - mXmin > 1:
                    Bm = '-Bx0.5 -By0.5 -Bz100 -B%sZ' % args.bs
                else:
                    Bm = '-Bx0.2 -By0.1 -Bz100 -B%sZ' % args.bs
        else:
            if mXmax - mXmin > 1:
                Bm = '-Bx0.5 -By0.5 -Bz100 -BNWseZ'
            else:
                Bm = '-Bx0.2 -By0.1 -Bz100 -BNWseZ'
        '''
        if mXmax - mXmin > 1:
            Bm = '-Bx0.5 -By0.5 -Bz100 -BNWSEZ'
        else:
            Bm = '-Bx0.05 -By0.05 -Bz100 -BNWSEZ'


        sf = '%s_3D.ps'%args.savefile[:-3]
        os.system("gmt psbasemap %s %s %s %s -p%s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Jz, Rm, Bm, pval, sf))
        #pscoast -Rd -JX8id/5id -Dc -Gblack -E200/40 -K -U"Example 10 in Cookbook" > $ps
        #os.system("gmt grdview %s -R -J -JZ -O -K -C%s -Ihillshade2.grd -p%s >> %s"%(DEM2, gmcpt, sf))
        if '.grdz' not in args.valfile[0]:
            os.system("gmt grdview %s -R -J -JZ -O -K -C%s -Ihillshade.grd -Qs -p%s>> %s"%(DEM, gmcpt, pval, sf))
        if not omitgeology:
            os.system("gmt grdtrack %s -GupliftE_0.001.grd > trackE_4.xygt" % oligocene)
            os.system("gmt psxyz trackE_4.xygt -R -J -JZ -W%s,red -O -K -p%s >> %s" % (geosize, pval, sf))
            os.system("gmt grdtrack %s -GupliftE_0.001.grd > trackE_4.xygt" % miocene)
            os.system("gmt psxyz trackE_4.xygt -R -J -JZ -W%s,orange -O -K -p%s >> %s" % (geosize, pval, sf))

        tcpt = args.valfile[3]
        minval = float(args.valfile[2]) * -750.0
        maxval = float(args.valfile[1]) * -750.0
        print ('minval, maxval',minval, maxval)
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        if '.grd' not in args.valfile[0]:

            dat = pd.read_csv(args.valfile[0], delim_whitespace=True)
            dat['time2'] = dat['time'].values * -750
            dat['time'] = dat['time'].values * -750
            dat.to_csv('temp.txt', header=True, index=False, sep=' ')
            os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
        else:
            sgrdfile = args.valfile[0]
            os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
            #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
            os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if args.valfile2 is not None:
            tcpt = args.valfile2[3]
            minval = float(args.valfile2[2]) * -750.0
            maxval = float(args.valfile2[1]) * -750.0
            if '*' in tcpt:
                tcpt = tcpt[:-1]
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

            if '.grd' not in args.valfile2[0]:

                dat = pd.read_csv(args.valfile2[0], delim_whitespace=True)
                dat['time2'] = dat['time'].values * -750
                dat['time'] = dat['time'].values * -750
                dat.to_csv('temp.txt', header=True, index=False, sep=' ')
                os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
            else:
                sgrdfile = args.valfile2[0]
                os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
                #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
                os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if args.valfile3 is not None:
            tcpt = args.valfile3[3]
            minval = float(args.valfile3[2]) * -750.0
            maxval = float(args.valfile3[1]) * -750.0
            if '*' in tcpt:
                tcpt = tcpt[:-1]
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            if '.grd' not in args.valfile3[0]:

                dat = pd.read_csv(args.valfile3[0], delim_whitespace=True)
                dat['time2'] = dat['time'].values * -750
                dat['time'] = dat['time'].values * -750
                dat.to_csv('temp.txt', header=True, index=False, sep=' ')
                os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
            else:
                sgrdfile = args.valfile3[0]
                os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
                #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
                os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if args.valfile4 is not None:
            tcpt = args.valfile4[3]
            minval = float(args.valfile4[2]) * -750.0
            maxval = float(args.valfile4[1]) * -750.0
            if '*' in tcpt:
                tcpt = tcpt[:-1]
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            if '.grd' not in args.valfile4[0]:

                dat = pd.read_csv(args.valfile4[0], delim_whitespace=True)
                dat['time2'] = dat['time'].values * -750
                dat['time'] = dat['time'].values * -750
                dat.to_csv('temp.txt', header=True, index=False, sep=' ')
                os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
            else:
                sgrdfile = args.valfile4[0]
                os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
                #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
                os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if args.valfile5 is not None:
            tcpt = args.valfile5[3]
            minval = float(args.valfile5[2]) * -750.0
            maxval = float(args.valfile5[1]) * -750.0
            if '*' in tcpt:
                tcpt = tcpt[:-1]
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            if '.grd' not in args.valfile5[0]:

                dat = pd.read_csv(args.valfile5[0], delim_whitespace=True)
                dat['time2'] = dat['time'].values * -750
                dat['time'] = dat['time'].values * -750
                dat.to_csv('temp.txt', header=True, index=False, sep=' ')
                os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
            else:
                sgrdfile = args.valfile5[0]
                os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
                #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
                os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if args.valfile6 is not None:
            tcpt = args.valfile6[3]
            minval = float(args.valfile6[2]) * -750.0
            maxval = float(args.valfile6[1]) * -750.0
            if '*' in tcpt:
                tcpt = tcpt[:-1]
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                if float(maxval) - float(minval) < 0.1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                elif float(maxval) - float(minval) < 1:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                else:
                    os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            if '.grd' not in args.valfile6[0]:

                dat = pd.read_csv(args.valfile6[0], delim_whitespace=True)
                dat['time2'] = dat['time'].values * -750
                dat['time'] = dat['time'].values * -750
                dat.to_csv('temp.txt', header=True, index=False, sep=' ')
                os.system("gmt psxyz temp.txt -R -J -JZ -Sc0.1 -Cdep.cpt -O -K -p%s >> %s" % (pval,sf))
            else:
                sgrdfile = args.valfile6[0]
                os.system("gmt grdgradient %s -A60/60 -GhillshadeG.grd"%(sgrdfile))
                #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -IhillshadeG.grd -Qs -p%s>> %s"%(sgrdfile, pval, sf))
                os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs -p%s>> %s"%(sgrdfile, pval, sf))

        if '.grdz' in args.valfile[0]:
            os.system("gmt grdview %s -R -J -JZ -O -K -C%s -Ihillshade.grd -Qs -p%s>> %s"%(DEM, gmcpt, pval, sf))

# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''
        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-d', '--dotfile', dest='dotfile', type=str,
                        help='list of dots')
    parser.add_argument('-l', '--linfile', dest='linfile', type=str,
                        help='list of lines')
    parser.add_argument('-a', '--labfile', dest='labfile', type=str,
                        help='list of labels')
    parser.add_argument('-s', '--segfile', dest='segfile', type=str,
                        help='list of segments')
    parser.add_argument('-b', '--bounds', metavar=('lonmin', 'lonmax', 'latmin', 'latmax'),
                        dest='bounds', type=float, nargs=4,required=True,
                        help='map bounds [lonmin lonmax latmin latmax]')
    parser.add_argument('-f', '--savefile', dest='savefile', type=str,
                        required=True, help='file to save to')
    parser.add_argument('-g', '--geology', dest='geology', type=str,
                        help='enter -g y to plot geology instead of elevation')
    parser.add_argument('-r', '--res', dest='res', type=str,
                        help='enter -r high or -r low ')
    parser.add_argument('-v', '--valfile', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-v2', '--valfile2', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile2', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-v3', '--valfile3', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile3', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-v4', '--valfile4', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile4', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-v5', '--valfile5', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile5', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-v6', '--valfile6', metavar=('file', 'min', 'max', 'cpt'),
                        dest='valfile6', type=str, nargs=4,
                        help='values [file minval maxval cpt]')
    parser.add_argument('-dem', '--dem', dest='dem', type=str,
                        help='enter "-dem otherdem.grd" to plot DEM other than puget sound')
    parser.add_argument('-ll', '--llabel', dest='llabel', type=str,
                        help='string to print in lower-left corner of map')
    parser.add_argument('-lr', '--lrabel', dest='lrabel', type=str,
                        help='string to print in lower-right corner of map')
    parser.add_argument('-tl', '--tlabel', dest='tlabel', type=str,
                        help='string to print in top-left corner of map')
    parser.add_argument('-tr', '--trabel', dest='trabel', type=str,
                        help='string to print in top-right corner of map')
                        
    parser.add_argument('-jm', '--jm', dest='jm', type=float,
                        help='enter e.g. "-jm 10" for 10 inch wide figure')
                        
    parser.add_argument('-bs', '--bs', dest='bs', type=str,
                        help='enter e.g. "-bs NWse" for labels on top and left and ticks on bottom and right of map')
                        
    parser.add_argument('-pv', '--plotvec', dest='plotvec', type=str,
                        help='enter "-pv y" to plot as vectors')
                        
    parser.add_argument('-ns', '--ns', dest='ns', type=str,
                        help='enter "-ns y" to not include scale bar')
                        
    parser.add_argument('-fh', '--fh', dest='forcehorizontal', type=str,
                        help='enter "-fh y" to force horizontal scale bar for the first value file')
                        
    parser.add_argument('-puc', '--puc', dest='puc', type=str,
                        help='enter "-puc y" to plot uplift contours')
                        
    parser.add_argument('-pup', '--pup', dest='pup', type=str,
                        help='enter "-pup y" to plot uplift points')
                        
    parser.add_argument('-put', '--put', dest='put', type=str,
                        help='enter "-pup y" to plot uplift text')
                        
    parser.add_argument('-psd100', '--psd100', dest='psd100', type=str,
                        help='enter "-psd100 y" to plot strike and dip values for 100k')
                        
    parser.add_argument('-psd24', '--psd24', dest='psd24', type=str,
                        help='enter "-psd24 y" to plot strike and dip values for 24k')
                        
    parser.add_argument('-ppf', '--ppf', dest='ppf', type=str,
                        help='enter "-ppf y" to plot pratt faults')
                        
    parser.add_argument('-pqf', '--pqf', dest='pqf', type=str,
                        help='enter "-pqf y" to plot quaternary folds')
                        
    parser.add_argument('-oqf', '--oqf', dest='oqf', type=str,
                        help='enter "-oqf y" to omit quaternary fault database')
                        
    parser.add_argument('-odl', '--odl', dest='odl', type=str,
                        help='enter "-odl y" to omit dashed line marking basin edge')
                        
    parser.add_argument('-omo', '--omo', dest='omo', type=str,
                        help='enter "-omo y" to omit miocene and oligocene contacts')
                        
    parser.add_argument('-ppn', '--ppn', dest='ppn', type=str,
                        help='enter "-ppn y" to plot place names')
                        
    parser.add_argument('-pic', '--pic', dest='pic', type=str,
                        help='enter "-pic y" to plot topography in color')
                        
    parser.add_argument('-p17', '--p17', dest='p17', type=str,
                        help='enter "-p17 y" to plot pss17 and lwb17 nav')

    parser.add_argument('-p70', '--p70', dest='p70', type=str,
                        help='enter "-p70 y" to plot 1970 profiles')
                        
    parser.add_argument('-p95', '--p95', dest='p95', type=str,
                        help='enter "-p95 y" to plot 1995 profiles')
                        
    parser.add_argument('-p97', '--p97', dest='p97', type=str,
                        help='enter "-p97 y" to plot 1997 profiles')
                        
    parser.add_argument('-p11', '--p11', dest='p11', type=str,
                        help='enter "-p11 y" to plot 1911 profiles')
                        
    parser.add_argument('-pma', '--pma', dest='pma', type=str,
                        help='enter "-pma y" to plot magnetic anomaly contours')
                        
    parser.add_argument('-pkc', '--pkc', dest='pkc', type=str,
                        help='enter "-pkc y" to plot the locations of the cores in Karlin et al 2004')
                        
    parser.add_argument('-psf', '--psf', dest='psf', type=str,
                        help='enter "-psf y" to plot the locations of slope failure in Karlin et al 2004')
                        
    parser.add_argument('-p3D', '--plot3D', dest='plot3D', type=str,
                        help='enter "130/20" to view from an azimuth of 130 degrees and an elevation of 20 meters')
                        
    parser.add_argument('-pbl', '--pbl', dest='pbl', type=str,
                        help='enter "-pbl y" to plot back thrusts as lines instead of arrows')
                        
    parser.add_argument('-pls', '--pls', dest='pls', type=str,
                        help='enter "-pls y" to plot seismicity')
                        
    parser.add_argument('-osf', '--osf', dest='osf', type=str,
                        help='enter "-osf y" to exclude the USGS Quaternary fault and fold database from the map')
                        
    parser.add_argument('-oma', '--oma', dest='oma', type=str,
                        help='enter "-oma y" to not plot magnetic anomaly concacts with fault zone')
                        
    parser.add_argument('-jdf', '--jdf', dest='jdf', type=str,
                        help='enter "-jdf y" to only plot deformation front. "-ppf y" flag must also be included')
                        
    parser.add_argument('-bigin', '--bigin', dest='bigin', type=str,
                        help='enter "-bigin y" to plot inset of original map on big outset map')
                        
    parser.add_argument('-riv', '--riv', dest='riv', type=str,
                        help='enter "-riv y" to plot major rivers')
                        
    parser.add_argument('-ab1', '--addbox1', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax', 'trabel'),
                        dest='addbox1', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab2', '--addbox2', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox2', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab3', '--addbox3', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox3', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab4', '--addbox4', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox4', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab5', '--addbox5', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox5', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
    parser.add_argument('-ab6', '--addbox6', metavar=('cdpmin', 'cdpmax', 'tmin', 'tmax','trlabel'),
                        dest='addbox6', type=str, nargs=5,
                        help='bounds of another profile inset [cdpmin cdpmax tmin tmax trlabel]')
                        
    pargs = parser.parse_args()
    
    #cProfile.run('main(pargs)')
    main(pargs)


