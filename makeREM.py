import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from netCDF4 import Dataset
from pyproj import Proj

def main(args):
    
    logfile = '%s_log.txt' % args.og[:-4]
    os.system("rm %s" % logfile)
    outF = open(logfile, 'w')
    outF.write('Input arguments to makeREM.py: ')
    outF.write('\n')
    outF.write('\n')
    print (args, file=outF)
    outF.close()
    
    homedir = '/Users/ginevramoore/Documents/rivers/NHD'
    smoothgrid = '%s_relevs.grd' % args.og[:-4]
    smoothgrid2 = '%s_relevs2.grd' % args.og[:-4]
    remgrid = '%s_rem.grd' % args.og[:-4]
    remgrid2 = '%s_rem2.grd' % args.og[:-4]
    cutgrid = '%s_cut.grd' % args.og[:-4]
    cutgrid2 = '%s_cut2.grd' % args.og[:-4]
    maskgrid = '%s_mask.grd' % args.og[:-4]
    hillshade = '%s_HS.grd' % args.og[:-4]
    hillshade2 = '%s_HS2.grd' % args.og[:-4]
    
    os.system("rm %s" % smoothgrid)
    os.system("rm %s" % remgrid)
    os.system("rm %s" % remgrid2)
    os.system("rm %s" % cutgrid)
    os.system("ogr2ogr -f CSV %s.csv %s -lco GEOMETRY=AS_WKT" % (args.nhdfile[:-4], args.nhdfile))
    dat = pd.read_csv('%s.csv' % args.nhdfile[:-4])
    
    if args.nhdfile2 is not None:
        os.system("ogr2ogr -f CSV %s.csv %s -lco GEOMETRY=AS_WKT" % (args.nhdfile2[:-4], args.nhdfile2))
        dat2 = pd.read_csv('%s.csv' % args.nhdfile2[:-4])
        dat = pd.concat([dat, dat2])
    
    print (len(dat))

    nofcodecol = False
    if 'IX' in dat.columns:
        dat['GNIS_Name'] = dat['IX'].values.astype('str')
        nofcodecol = True
    
    if args.rn4 is not None:
        dat3 = dat[dat.GNIS_Name == args.rn4]
        if len(dat3) == 0:
            print ('not a valid river name for this NHD file')
            print ('exiting ...')
            exit()

    if args.rn3 is not None:
        dat2 = dat[dat.GNIS_Name == args.rn3]
        if len(dat2) == 0:
            print ('not a valid river name for this NHD file')
            print ('exiting ...')
            exit()

    if args.rn2 is not None:
        dat1 = dat[dat.GNIS_Name == args.rn2]
        if len(dat1) == 0:
            print ('not a valid river name for this NHD file')
            print ('exiting ...')
            exit()

    if args.rn is not None:
        dat = dat[dat.GNIS_Name == args.rn]
        if len(dat) == 0:
            print ('not a valid river name for this NHD file')
            print ('exiting ...')
            exit()
        if args.rn2 is not None:
            dat = pd.concat([dat, dat1])
        if args.rn3 is not None:
            dat = pd.concat([dat, dat2])
        if args.rn4 is not None:
            dat = pd.concat([dat, dat3])


    if nofcodecol and args.rn is None and args.rn2 is None:
        print ('importing shapefile from topotoolbox')
        print ('use info tool in QGIS to find "IX" number for desired river')
        print ('use ^^^ that number as "-rn" flag')
        print ('exiting ...')
        exit()

    if args.rn is None and args.rn2 is None:
        dat = dat[dat.FCode == 55800]

    print (len(dat))
    #exit()
    i = 0
    lons = []
    lats = []
    if nofcodecol:
        lonsp = latsp = 1
    else:
        lonsp = latsp = 0.0001
    for index,row in dat.iterrows():
        i += 1
        line = row['WKT']
        [pre, linelist] = line.split('(')
        linelist = linelist[:-1]
        pts = linelist.split(',')
        if i%5000 == 0:
            print (i, len(dat), len(pts))
        lons0 = []
        lats0 = []
        for k in range(len(pts)):
            coords = pts[k].split(' ')
            lons0.append(float(coords[0]))
            lats0.append(float(coords[1]))

        if not nofcodecol:
            xdiff = np.max(lons0) - np.min(lons0)
            ydiff = np.max(lats0) - np.min(lats0)
            lons0 = np.array(lons0)
            lats0 = np.array(lats0)
            if xdiff < ydiff:
                latx = np.arange(np.min(lats0), np.max(lats0), latsp)
                lons0a = griddata(lats0, lons0, latx, method='linear')
                lats0a = griddata(lons0, lats0, lons0a, method='linear')
            else:
                lonx = np.arange(np.min(lons0), np.max(lons0), lonsp)
                lats0a = griddata(lons0, lats0, lonx, method='linear')
                lons0a = griddata(lats0, lons0, lats0a, method='linear')
            lons.extend(list(lons0a))
            lats.extend(list(lats0a))

        else:
            lons.extend(list(lons0))
            lats.extend(list(lats0))
            
    os.system("gmt grdinfo %s > temp1.txt -C" % args.dem)
    file1 = open('temp1.txt', 'r')
    lines = file1.readlines()
    line = lines[0]
    (name, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nrows, ncols) = line.split()
    demProj = Proj("+proj=utm +zone=%s, +north +ellps=GRS80 +datum=WGS84 +units=m +no_defs" % (args.utm))
    os.system("gmt grdinfo %s > temp1.txt" % args.dem)
    file1 = open('temp1.txt', 'r')
    lines = file1.readlines()
    for iline in lines:
        if '+proj' in iline:
            print (iline)
            demProj = Proj("%s" % iline)
            break

    deminutm = False
    if float(xmin) > 1000 or float(ymin) > 1000 or float(xmax) < -1000 or float(ymax) < -1000:
        print ('DEM in UTM')
        print ('run grdinfo on DEM to check utmzone and rerun if you have not already')
        print ('continuing ... ')
        deminutm = True

    nhdinutm = False
    if min(lons) > 1000 or min(lats) > 1000 or max(lons) < -1000 or max(lats) < -1000:
        print ('NHD in UTM')
        print ('run grdinfo on DEM to check utmzone and rerun if you have not already')
        print ('continuing ... ')
        nhdinutm = True

    myProj = Proj("+proj=utm +zone=%s, +north +ellps=GRS80 +datum=WGS84 +units=m +no_defs" % (args.utm))
    if deminutm and not nhdinutm:
        print ('converting nhd to utm .... ')
        lons, lats = myProj(lons, lats, inverse=False)

    xmin1 = min(lons)
    xmax1 = max(lons)
    ymin1 = min(lats)
    ymax1 = max(lats)

    print (name, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nrows, ncols)

    if not args.bounds:
        if xmin1 > float(xmin):
            xmin = xmin1
        if xmax1 < float(xmax):
            xmax = xmax1
        if ymin1 > float(ymin):
            ymin = ymin1
        if ymax1 < float(ymax):
            ymax = ymax1

    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)

    forcebounds = False
    if args.fb:
        if args.fb == 'y':
            forcebounds = True

    if args.bounds and forcebounds:
        mxmin = args.bounds[0]
        mxmax = args.bounds[1]
        mymin = args.bounds[2]
        mymax = args.bounds[3]
        xmin = args.bounds[0]
        xmax = args.bounds[1]
        ymin = args.bounds[2]
        ymax = args.bounds[3]
    elif args.bounds:
        mxmin = args.bounds[0]
        mxmax = args.bounds[1]
        mymin = args.bounds[2]
        mymax = args.bounds[3]
        if deminutm and not (mxmin > 1000 or mxmax > 1000 or mymin > 1000 or mymax > 1000):
            [mxmin, mxmax], [mymin,mymax] = myProj([mxmin, mxmax], [mymin,mymax], inverse=False)
        if xmin < mxmin:
            xmin = mxmin
        if xmax > mxmax:
            xmax = mxmax
        if ymin < mymin:
            ymin = mymin
        if ymax > mymax:
            ymax = mymax
    else:
        print ('no bounds entered ...')
        mxmin = min(lons)
        mxmax = max(lons)
        mymin = min(lats)
        mymax = max(lats)

    print ("gmt grdcut -G%s -R%s/%s/%s/%s %s" %(cutgrid, xmin, xmax, ymin, ymax, args.dem))
    os.system("gmt grdcut -G%s -R%s/%s/%s/%s %s" %(cutgrid, xmin, xmax, ymin, ymax, args.dem))

    os.system("gmt grdgradient %s -A315/45 -G%s"%(cutgrid, hillshade))

    os.system("gmt grdinfo %s > temp1.txt -C" % cutgrid)
    file1 = open('temp1.txt', 'r')
    lines = file1.readlines()
    line = lines[0]
    (name, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nrows, ncols) = line.split()

    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)

    print (xmin, xmax, ymin, ymax)
    
    rflag="-R%s/%s/%s/%s" %(xmin,xmax,ymin,ymax)
    iflag="-I%s/%s" %(dx, dy)
    iflagb="-I%s/%s" %(float(dx)*args.dxm, float(dy)*args.dxm)

    newdat = pd.DataFrame({'lon':lons, 'lat':lats})
    newdat = newdat[['lon','lat']]
    newdat.to_csv('tempriv.txt', header=False, index=False, sep=' ')
    newdat = newdat[newdat.lon >= float(xmin)]
    newdat = newdat[newdat.lon <= float(xmax)]
    newdat = newdat[newdat.lat >= float(ymin)]
    newdat = newdat[newdat.lat <= float(ymax)]
    newdat.to_csv('temp.txt', header=False, index=False, sep=' ')

    if len(newdat) < 1:
        print ('            no river locations within bounds of DEM')
        print ('            check DEM and NHD inputs and bounds and re-run')
        print ('            exiting ... ')
        exit()

    print ("gmt grdtrack temp.txt -G%s > temp2.txt" % (args.dem))
    os.system("gmt grdtrack temp.txt -G%s > temp2.txt" % (args.dem))
    
    
    print ("gmt blockmean %s %s %s | gmt surface -G%s %s %s -T0.0i -T0.0b -V -r" % ('temp2.txt', rflag, iflagb, smoothgrid, iflag, rflag))
    os.system("gmt blockmean %s %s %s | gmt surface -G%s %s %s -T0.0i -T0.0b -V -r" % ('temp2.txt', rflag, iflagb, smoothgrid, iflag, rflag))

    print ("gmt grdmath %s %s SUB = %s" % (cutgrid, smoothgrid, remgrid))
    os.system("gmt grdmath %s %s SUB = %s" % (cutgrid, smoothgrid, remgrid))
    os.system("gmt grdmath %s 10 LT 0 NAN = %s" % (remgrid, maskgrid))
    os.system("gmt grdmath %s %s MUL = %s" % (remgrid, maskgrid, remgrid2))


    if args.plot is not None:
        tcpt = args.plot[2]
        minval = float(args.plot[0])
        maxval = float(args.plot[1])
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.0001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.0001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))

        os.system("gmt gmtset MAP_FRAME_TYPE = plain")
        Rm = rflag
        Bm = '-Bx11000002 -By11000002 -BNWse'
        sf = '%s_map.ps' % args.og[:-4]
        if (float(xmax) - float(xmin)) > (float(ymax) - float(ymin)):
            if deminutm:
                Jm = '-JX20/%f' % ((ymax-ymin) * 20 / (xmax-xmin))
                mapwidth = 20
                Jz = '-JZ20c'
            else:
                Jm = '-JM20'
                mapwidth = 20
                Jz = '-JZ20c'
            os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
        else:
            if deminutm:
                Jm = '-JX10/%f' % ((ymax-ymin) * 10 / (xmax-xmin))
                mapwidth = 10
                Jz = '-JZ10c'
            else:
                mapwidth = 10
                Jm = '-JM10'
                Jz = '-JZ10c'
            os.system("gmt psbasemap %s %s %s %s -K -P --FORMAT_GEO_MAP=D > %s" % (Jm, Jz, Rm, Bm, sf))

        plotsx = (xmax+xmin) / 2
        plotsy = (ymax+ymin) / 2

        if (xmax - xmin) < 0.1:
            distance = 1
        elif (xmax - xmin) < 0.2:
            distance = 2
        elif (xmax - xmin) < 0.4:
            distance = 4
        elif (xmax - xmin) < 0.8:
            distance = 8
        elif (xmax - xmin) < 1:
            distance = 10
        elif (xmax - xmin) < 2:
            distance = 20
        else:
            distance = 50

        #os.system("gmt grdview %s -R -J -JZ -O -K -Cdep.cpt -Qs >> %s "%(remgrid2, sf))
        os.system("gmt grdimage %s -R -Cdep.cpt -J -O -K -Q >> %s"%(remgrid, sf))
        os.system("gmt grdimage %s -R -Cgray -J -O -K >> %s"%(hillshade, sf))
        os.system("gmt grdimage %s -R -Cdep.cpt -J -O -K -Q >> %s"%(remgrid, sf))
        if not deminutm:
            os.system("gmt pscoast -R -J -L%s/%s+c%s+w%sk+u -Df -O -K -I1/0.25 >> %s " % (plotsx, plotsy, plotsy, distance, sf))
        #os.system("gmt psscale -D10/5/3/0.5 -Cdep.cpt -B5 -O -K >> %s"%(sf))
        os.system('gmt psscale -D0/-0.9+w%i/0.5+h -Cdep.cpt -B%f+l"meters above river surface" -O -K >> %s'%(mapwidth, int(maxval-minval) / 5, sf))


    plotinsets = False
    if args.pi is not None:
        if args.pi == 'y':
            plotinsets = True

    plotincon = False
    if args.pic is not None:
        plotincon = True

    if plotinsets:
        if deminutm:
            [xmin, xmax], [ymin, ymax] = myProj([xmin, xmax], [ymin, ymax], inverse=True)

        Rm = '-R%f/%f/%f/%f' % (xmin-10,xmax+10,ymin-10,ymax+10)
        Jm = '-JM4'
        mapwidth = 4
        Bm = '-Bx1 -By1 -BNWse'

        DEM2 = '/Users/ginevramoore/Documents/rivers/DEMs/world_lr.grd'
        gmcpt = '/Users/ginevramoore/Documents/research/segymanip/library/gmtopo.cpt'
        #os.system("gmt grdgradient %s -A60/60 -Ghillshade2.grd"%(DEM2))

        sf = '%s_inset_big.ps'%args.og[:-4]
        lons = [xmin, xmin, xmax, xmax, xmin]
        lats = [ymin, ymax, ymax, ymin, ymin]
        dat = pd.DataFrame({'lon':lons, 'lat':lats})
        dat = dat[['lon','lat']]
        dat.to_csv('temp.txt', header=False,index=False,sep=' ')
        os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
        #os.system("gmt grdimage %s -R -C%s -J -O -K >> %s" % (DEM2, gmcpt, sf))
        os.system("gmt pscoast -Dc -R -J -O -K -ENA+p0.25p,red+r -EUS.%s+gwhite+p0.25p,black+r >> %s" % (args.state, sf))
        os.system("gmt psxy temp.txt -R -J -W2,red -O -K >> %s" % (sf))

        os.system("gmt grdinfo %s > temp1.txt -C" % smoothgrid)
        file1 = open('temp1.txt', 'r')
        lines = file1.readlines()
        line = lines[0]
        (name, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nrows, ncols) = line.split()
        zmin = float(zmin)
        zmax = float(zmax)
        os.system("gmt grdtrack tempriv.txt -G%s > tempriv2.txt" % (args.dem))
        dat = pd.read_csv('tempriv2.txt', names=['x','y','z'], delim_whitespace=True)
        zmin = dat['z'].min()
        zmax = dat['z'].max()
        os.system("gmt makecpt -C%s -D -Fr -T%f/%f/10 -Z > riv.cpt" % (tcpt, zmin, zmax))

        os.system("gmt grdinfo %s > temp1.txt -C" % args.dem)
        file1 = open('temp1.txt', 'r')
        lines = file1.readlines()
        line = lines[0]
        (name, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, nrows, ncols) = line.split()
        xmin = float(xmin)
        xmax = float(xmax)
        ymin = float(ymin)
        ymax = float(ymax)
        zmin = float(zmin)
        zmax = float(zmax)
        os.system("gmt makecpt -C%s -D -Fr -T%f/%f/10 -Z > elev.cpt" % ('gray', zmin, zmax))

        if deminutm:
            dat = pd.DataFrame({'lon':lons, 'lat':lats})
            dat['lon'], dat['lat'] = myProj(dat['lon'].values, dat['lat'].values, inverse=False)
            dat = dat[['lon','lat']]
            dat.to_csv('temp.txt', header=False,index=False,sep=' ')

        Rm = '-R%f/%f/%f/%f' % (xmin,xmax,ymin,ymax)
        if deminutm:
            Jm = '-JX4'
            mapwidth = 4
        else:
            Jm = '-JM4'
            mapwidth = 4
        Bm = '-Bx1000000 -By10000000 -BNWse'

        sf = '%s_inset.ps'%args.og[:-4]
        print ("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
        os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
        os.system("gmt grdimage %s -R -Celev.cpt -J -O -K >> %s" % (args.dem, sf))
        os.system("gmt psxy temp.txt -R -J -W2,red -O -K >> %s" % (sf))
        if plotincon:
            #os.system("gmt grdcontour %s -L%f/%f -C%f %s %s -O -K -W0.34,white,- >> %s"%(args.dem, args.pic, args.pic*2, args.pic, Rm,Jm, sf))
            os.system("gmt grdcontour %s -C1 %s %s -O -K -W0.1,white -L%s/%s >> %s"%(args.dem, Rm, Jm, args.pic, args.pic+1, sf))
        os.system("gmt psxy tempriv2.txt -R -J -Sc0.04 -Criv.cpt -O -K >> %s" % (sf))

        os.system("gmt psscale -D10/5/3/0.5 -Cdep.cpt -B5 -O -K >> %s"%(sf))
        os.system("gmt psscale -D12/5/3/0.5 -Celev.cpt -B500 -O -K >> %s"%(sf))
        os.system("gmt psscale -D14/5/3/0.5 -Criv.cpt -B100 -O -K >> %s"%(sf))


    if args.plot2 is not None:
        tcpt = args.plot2[2]
        minval = float(args.plot2[0])
        maxval = float(args.plot2[1])
        if '*' in tcpt:
            tcpt = tcpt[:-1]
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.0001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -I -D -Fr -T%f/%f/1 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
        else:
            if float(maxval) - float(minval) < 0.1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.0001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            elif float(maxval) - float(minval) < 1:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.001 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
            else:
                os.system("gmt makecpt -C%s -D -Fr -T%f/%f/0.01 -Z > dep.cpt" % (tcpt, float(minval), float(maxval)))
                
        os.system("gmt grdimage %s -R -Cdep.cpt -J -O -K >> %s"%(smoothgrid, sf))
        sf = '%s_map_smooth.ps' % args.og[:-4]
        if (float(xmax) - float(xmin)) > (float(ymax) - float(ymin)):
            Jm = '-JM20'
            mapwidth = 20
            os.system("gmt psbasemap %s %s %s -K --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))
        else:
            Jm = '-JM10'
            mapwdith = 10
            os.system("gmt psbasemap %s %s %s -K -P --FORMAT_GEO_MAP=D > %s" % (Jm, Rm, Bm, sf))


    keeprem = False
    if args.keeprem is not None:
        if args.keeprem == 'y':
            keeprem = True

    keepaf = False
    if args.keepaf is not None:
        if args.keepaf == 'y':
            keepaf = True

    if not keeprem and not keepaf:
        os.system("rm %s" % remgrid)
    if not keepaf:
        os.system("rm %s" % cutgrid)
        os.system("rm %s" % cutgrid2)
        os.system("rm %s" % remgrid2)
        os.system("rm %s" % smoothgrid)
        os.system("rm %s" % smoothgrid2)
        os.system("rm %s" % hillshade)
        os.system("rm %s" % hillshade2)
        os.system("rm %s" % maskgrid)
        os.system("rm tempriv.txt")
        os.system("rm tempriv2.txt")
        os.system("rm temp.txt")
        os.system("rm temp1.txt")
        os.system("rm temp2.txt")
        os.system("rm *.cpt")
        os.system("rm gmt.history")
        os.system("rm gmt.conf")


# Help/description and command line argument parser
if __name__=='__main__':
    desc = '''

        '''
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-b', '--bounds', metavar=('lonmin', 'lonmax', 'latmin', 'latmax'),
                        dest='bounds', type=float, nargs=4,required=False,
                        help='map bounds [lonmin lonmax latmin latmax]')
                        
    parser.add_argument('-nf', '--nhdfile', dest='nhdfile', type=str,
                        required=True, help='nhdfile')
    parser.add_argument('-nf2', '--nhdfile2', dest='nhdfile2', type=str,
                        required=False, help='nhdfile2')
    parser.add_argument('-keeprem', '--keeprem', dest='keeprem', type=str,
                        required=False, help='enter "-keeprem y" to not delete the REM grid')
    parser.add_argument('-keepaf', '--keepaf', dest='keepaf', type=str,
                        required=False, help='enter "-keepaf y" to not delete any grids or text files')
    parser.add_argument('-pi', '--pi', dest='pi', type=str,
                        required=False, help='enter "-pi y" to plot inset maps')
    parser.add_argument('-pic', '--pic', dest='pic', type=float,
                        required=False, help='enter optional float to plot contour at on inset map')
    parser.add_argument('-dem', '--dem', dest='dem', type=str,
                        required=True, help='DEM')
    parser.add_argument('-fb', '--fb', dest='fb', type=str,
                        required=False, help='enter "-fb y" to force use of entered bounds')
    parser.add_argument('-dxm', '--dxm', dest='dxm', type=float,
                        required=True, help='multiplier of dy/dx of DEM for blockmean')
    parser.add_argument('-state', '--state', dest='state', type=str,
                        required=True, help='state - e.g. "-state WA"')
    parser.add_argument('-utm', '--utm', dest='utm', type=str,
                        required=True, help='UTM zone (e.g. "-utm 10T")')
    parser.add_argument('-og', '--og', dest='og', type=str,
                        required=True, help='Out grid')
    parser.add_argument('-p', '--plot', metavar=('min', 'max', 'cpt'),
                        dest='plot', type=str, nargs=3,
                        help='values [minval maxval cpt] for REM grid')
    parser.add_argument('-p2', '--plot2', metavar=('min', 'max', 'cpt'),
                        dest='plot2', type=str, nargs=3,
                        help='values [minval maxval cpt] for smoothed grid')
    parser.add_argument('-rn', '--rn', dest='rn', type=str,
                        required=False, help='river name e.g. "-rn "Roaring Fork River""')
    parser.add_argument('-rn2', '--rn2', dest='rn2', type=str,
                        required=False, help='river name 2 e.g. "-rn2 "Crystal River""')
    parser.add_argument('-rn3', '--rn3', dest='rn3', type=str,
                        required=False, help='river name 3 e.g. "-rn3 "Cattle Creek""')
    parser.add_argument('-rn4', '--rn4', dest='rn4', type=str,
                        required=False, help='river name 4 e.g. "-rn4 "Frying Pan River""')
                        
    pargs = parser.parse_args()
    
    #cProfile.run('main(pargs)')
    main(pargs)
