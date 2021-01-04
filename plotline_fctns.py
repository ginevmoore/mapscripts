import pandas as pd
import numpy as np
import os
import math

def getNavCdpFromHeader(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,cdp > %s/%s"%(tapedir,tapename,navdir,navname))
    print ('%s/%s'%(navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','cdp'])
    
    print (df)
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp']]

    return gxcol

def getNavCdpWBFromHeader(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,cdp,swdep > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','cdp','swdep'])
    
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str
    swdepcol = df[['swdep']]
    swdepcol['string'], swdepcol['swdep'] = swdepcol['swdep'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol['swdep'] = swdepcol['swdep'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp','swdep']]

    return gxcol

def getNavCdpFromHeaderSX(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=sx,sy,cdp > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','cdp'])
    
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp']]

    return gxcol

def getNavCdpWBFromHeaderSX(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=sx,sy,cdp,swdep > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','cdp','swdep'])
    
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str
    swdepcol = df[['swdep']]
    swdepcol['string'], swdepcol['swdep'] = swdepcol['swdep'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol['swdep'] = swdepcol['swdep'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp','swdep']]

    return gxcol

def npcosine(lon1, lat1, lon2, lat2):

    # Written GLM 4.25.17

    ''' Arguments:  lon1 - longitude point that the angle is referenced from
                            (float)[deg]
                    lat1 - latitude point that the angle is referenced from
                            (float)[deg]
                    lon2 - array of longitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                    lat2 - array of latitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                        
        Returns:    dist - array of great circle distance between the two
                            lat/lon points (arr of floats)[km]
                    ang - array of angles between the two points (clockwise from
                            0 degrees from lat1/lon1 point)
                            (arr of floats)[deg] '''

    # Creating degrees/radians conversion constants
    d2r = (math.pi/180.0)
    r2d = (180.0/math.pi)
    ddlon = lon1-lon2
    
    # Ensuring that azimuths are between 0 and 360
    if lon1 < 0.0:
        lon1 += 360.0
    lon2[lon2<0.0] += 360.0

    # Getting distance and angle between the two points (in degrees)
    dist, ang = npcosrule(d2r, lon1, lat1, lon2, lat2)
    ang[(lon1>lon2)&(ddlon<180.0)] = 2*math.pi-ang[(lon1>lon2)&(ddlon<180.0)] # causes runtime
    dist = np.abs(dist * r2d)

    dist[dist > 180.0] = 360-dist[dist > 180] # causes runtime
    ang[dist > 180.0] += math.pi # causes runtime
    ang[ang > 2.0*math.pi] = 2.0*math.pi - ang[ang > 2.0*math.pi] # causes runtime
    dist *= 111.19
    ang *= r2d
    
    lon2[lon2<0]+=360

    return dist, ang


def npcosrule(d2r, lon1, lat1, lon2, lat2):

    # Written GLM 4.25.17

    ''' Arguments:  d2r - degree to radians conversion constant (float)
                    lon1 - longitude point that the angle is referenced from
                            (float)[deg]
                    lat1 - latitude point that the angle is referenced from
                            (float)[deg]
                    lon2 - array of longitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
                    lat2 - array of latitudes that the angle is going to
                            (clockwise from 0 degrees) (arr of floats)[deg]
        
        Returns:    dist2 - array of great circle distance between the two
                            lat/lon points (arr of floats)[km]
                    ang - array of angles between the two points (clockwise from
                            0 degrees from lat1/lon1 point)
                            (arr of floats)[deg] '''

    # breaks when lat1==lat2 or lon1==lon2. Add small difference where needed
    londiff = np.abs(lon2-lon1)
    latdiff = np.abs(lat2-lat1)
    lon2[londiff<0.0001] += 0.0001 # causes runtime
    lat2[latdiff<0.0001] += 0.0001 # causes runtime

    cl1 = (90.0-lat1)*d2r
    cl2 = (90.0-lat2)*d2r
    dlon = (lon2-lon1)*d2r

    coscl2 = np.cos(cl2)
    sincl2 = np.sin(cl2)
    cosdlon = np.cos(dlon)
    
    coscl1 = math.cos(cl1)
    sincl1 = math.sin(cl1)

    dist = (coscl1 * coscl2) + (sincl1 * sincl2 * cosdlon)

    dist[dist < -1] = -1.0 # causes runtime
    dist[dist > 1] = 1.0 # causes runtime

    dist2 = np.arccos(dist)
    dist2[dlon > math.pi] = 2*math.pi - dist2[dlon > math.pi] # causes runtime
    
    ang = np.zeros(len(dist))
    num = np.zeros(len(dist))
    den = np.zeros(len(dist))

    num[dist != 0] = (coscl2[dist != 0] - (dist[dist != 0] * coscl1))
    den[dist != 0] = (np.sin(dist2[dist != 0]) * sincl1)

    ang[dist != 0] = num[dist != 0] / den[dist != 0]
    ang[dist == 0] = 1.0
    ang[ang < -1] = -1.0 # causes runtime
    ang[ang > 1] = 1.0 # causes runtime
    ang2 = np.arccos(ang)
    
    return dist2, ang2


def getMinMaxVal(tapedir,tapename,val):

    navname = 'temp.nav'
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,%s > %s"%(tapedir,tapename,val,navname))
    df = pd.read_csv('%s'%(navname), delim_whitespace=True, names=['gx','gy','cdp'])
    
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp']]
    
    minval = gxcol['cdp'].min()
    maxval = gxcol['cdp'].max()
    #os.system("rm temp.nav")
    
    return minval, maxval

def getHeaderDf(tapename, navname):
    
    print ("segyread tape=%s verbose=0 | segyclean | sugethw key=sx,sy,gx,gy,cdp,tracr,fldr > %s"%(tapename,navname))
    os.system("segyread tape=%s verbose=0 | segyclean | sugethw key=sx,sy,gx,gy,cdp,tracr,fldr > %s"%(tapename,navname))
    df = pd.read_csv('%s'%(navname), delim_whitespace=True, names=['sx','sy','gx','gy','cdp','tracr','fldr'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    sxcol = df[['sx']]
    sxcol['string'], sxcol['sx'] = sxcol['sx'].str.split('=',1).str
    sycol = df[['sy']]
    sycol['string'], sycol['sy'] = sycol['sy'].str.split('=',1).str
    cdpcol = df[['cdp']]
    cdpcol['string'], cdpcol['cdp'] = cdpcol['cdp'].str.split('=',1).str
    fldrcol = df[['fldr']]
    fldrcol['string'], fldrcol['fldr'] = fldrcol['fldr'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['sx'] = sxcol['sx'].values.astype(np.float64)
    gxcol['sy'] = sycol['sy'].values.astype(np.float64)
    gxcol['cdp'] = cdpcol['cdp'].values.astype(int)
    gxcol['fldr'] = fldrcol['fldr'].values.astype(int)
    gxcol['tracr'] = tracrcol['tracr'].values.astype(int)
    gxcol = gxcol[['sx','sy','gx','gy','cdp','tracr','fldr']]

    return gxcol

def getNavTracrFromHeader(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,tracr > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','tracr'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['tracr'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp']]

    return gxcol

def getNavTracrFromHeaderSX(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=sx,sy,tracr > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','tracr'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['tracr'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp']]

    return gxcol

def getNavTracrWBFromHeader(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,tracr,swdep,sdepth > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','tracr','swdep','sdepth'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str
    swdepcol = df[['swdep']]
    swdepcol['string'], swdepcol['swdep'] = swdepcol['swdep'].str.split('=',1).str
    sdepcol = df[['sdepth']]
    sdepcol['string'], sdepcol['sdepth'] = sdepcol['sdepth'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['tracr'].values.astype(int)
    gxcol['swdep'] = swdepcol['swdep'].values.astype(int) + sdepcol['sdepth'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp','swdep']]
    return gxcol

def getNavTracrWBFromHeaderSX(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=sx,sy,tracr,swdep,sdepth > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','tracr','swdep','sdepth'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str
    swdepcol = df[['swdep']]
    swdepcol['string'], swdepcol['swdep'] = swdepcol['swdep'].str.split('=',1).str
    sdepcol = df[['sdepth']]
    sdepcol['string'], sdepcol['sdepth'] = sdepcol['sdepth'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['tracr'].values.astype(int)
    gxcol['swdep'] = swdepcol['swdep'].values.astype(int) + sdepcol['sdepth'].values.astype(int)
    gxcol = gxcol[['gx','gy','cdp','swdep']]
    return gxcol

def getNavTracrWBDelrtFromHeaderSX(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=sx,sy,tracr,swdep,delrt > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','tracr','swdep','delrt'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['tracr']]
    tracrcol['string'], tracrcol['tracr'] = tracrcol['tracr'].str.split('=',1).str
    swdepcol = df[['swdep']]
    swdepcol['string'], swdepcol['swdep'] = swdepcol['swdep'].str.split('=',1).str
    delrtcol = df[['delrt']]
    delrtcol['string'], delrtcol['delrt'] = delrtcol['delrt'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['tracr'].values.astype(int)
    gxcol['swdep'] = swdepcol['swdep'].values.astype(int)# + delrtcol['delrt'].values.astype(int) * 1e3
    gxcol['delrt'] = delrtcol['delrt'].values.astype(int)
    
    gxcol = gxcol[['gx','gy','cdp','swdep','delrt']]
    return gxcol

def getNavFLDRFromHeader(tapedir,navdir,tapename, navname):
    
    os.system("segyread tape=%s/%s verbose=0 | segyclean | sugethw key=gx,gy,fldr > %s/%s"%(tapedir,tapename,navdir,navname))
    df = pd.read_csv('%s/%s'%(navdir,navname), delim_whitespace=True, names=['gx','gy','fldr'])
    gxcol = df[['gx']]
    gxcol['string'], gxcol['gx'] = gxcol['gx'].str.split('=',1).str
    gycol = df[['gy']]
    gycol['string'], gycol['gy'] = gycol['gy'].str.split('=',1).str
    tracrcol = df[['fldr']]
    tracrcol['string'], tracrcol['fldr'] = tracrcol['fldr'].str.split('=',1).str

    gxcol['gx'] = gxcol['gx'].values.astype(np.float64)
    gxcol['gy'] = gycol['gy'].values.astype(np.float64)
    gxcol['cdp'] = tracrcol['fldr'].values.astype(int)
    gxcol = gxcol[['gx','gy','fldr']]

    return gxcol

def cosrule(d2r, lat1, lon1, lat2, lon2):
    
    ''' Arguments:  d2r - degree to radians conversion constant
                    lat1 - latitude point that the angle is referenced from
                    lon1 - longitude point that the angle is referenced from
                    lat2 - latitude point that the angle is going to (clockwise from 0 degrees)
                    lon2 - longitude point that the angle is going to (clockwise from 0 degrees)
        
        Returns:    dist2 - great circle distance between the two lat/lon points
                    ang - the angle between the two points (clockwise from 0 degrees from lat1/lon1 point) '''
    
    if abs(lon1-lon2) < 0.00001 or abs(lat1-lat2) < 0.00001:
        lat2 = lat2+0.0001
        lon2 = lon2+0.0001
    
    cl1 = (90-lat1) * d2r
    cl2 = (90-lat2) * d2r
    dlon = (lon2-lon1) * d2r
    dist = math.cos(cl1) * math.cos(cl2) + math.sin(cl1) * math.sin(cl2) * math.cos(dlon)
    if dist < -1:
        dist = -1.0
    if dist > 1:
        dist = 1.0
    dist2 = math.acos(dist)
    if dlon > math.pi:
        dist2 = 2 * math.pi-dist2
    if dist != 0:
        ang = (math.cos(cl2) - (dist * math.cos(cl1))) / (math.sin(dist2) * math.sin(cl1))
    else:
        ang = 1.0
    if ang < -1:
        ang = -1.0
    if ang > 1:
        ang = 1.0
    ang = math.acos(ang)
    return dist2, ang

def get_ns_dt(tapedir,tapename):
    os.system("segyread tape=%s/%s verbose=0 | segyclean | surange > temp.txt"%(tapedir,tapename))
    for line in open('temp.txt'):
        plist = line.split()
        if len(plist)>1:
            if plist[0] == 'ns':
                ns = plist[1]
            if plist[0] == 'dt':
                dt = plist[1]
    os.system("rm temp.txt")
    ns = int(ns)
    dt = float(dt)/10e5
    return ns, dt


###############################################

### 10 ###

###############################################

## Written GM

def cosine(lon1, lat1, lon2, lat2):
    
    ''' Arguments:  lon1 - latitude point that the angle is referenced from
                    lat1 - longitude point that the angle is referenced from
                    lon2 - latitude point that the angle is going to (clockwise from 0 degrees)
                    lat2 - longitude point that the angle is going to (clockwise from 0 degrees)
        
        Returns:    dist - great circle distance between the two lat/lon points
                    ang - the angle between the two points (clockwise from 0 degrees from lat1/lon1 point)
                    lon1 - same longitude as input argument, used in some applications
                    lat1 - same latitude as input argument, used in some applications  '''
    
    # Ensuring that azimuths are between 0 and 360
    if lon1 < 0:
        lon1 = lon1 + 360
    if lon2 < 0:
        lon2 = lon2 + 360

    # Creating degrees/radians conversion constants
    d2r = (math.pi/180)
    r2d = (180/math.pi)
    ddlon = lon1 - lon2

    # Getting distance and angle between the two points (in degrees)
    dist, ang = cosrule(d2r, lat1, lon1, lat2, lon2)
    if lon1 > lon2 and ddlon < 180:
        ang = 2*math.pi - ang
    dist = abs(dist*r2d)
    if dist > 180:
        dist = 360 - dist
        ang = ang + math.pi
    if ang > 2*math.pi:
        ang = 2*math.pi - ang
    dist = dist * 111.19
    ang = ang * r2d
    
    return dist, ang

def mkSDgrddata(xi, zi, flipornot):
    
    # get dx, dy, and list of lats from zi coordinates (listed in xi)
    xpts, ypts = xi[:, 0], xi[:, 1]
    xpts.shape = zi.shape
    ypts.shape = zi.shape
    dlats = ypts[:, 0]
    dlons = xpts[0, :]
    ny = len(dlats)
    nx = len(dlons)
    dy = abs(dlats[1] - dlats[0])
    dx = abs(dlons[1] - dlons[0])
    
    # flip array over if needed
    if flipornot == 'flip':
        depthgrid = np.flipud(zi)
    else:
        depthgrid = np.copy(zi)

    # initialize grid spacing in km
    alldy = dy * 111.19
    alldx = dx * 111.19
    Xgrad, Ygrad = [],[]

    # loop through lats and get gradient, use different lon spacing for each lat
    for i in range(1, ny - 1):
        thisgrid = depthgrid[i - 1:i + 2,:]
        thisy = math.radians(abs(dlats[i]))
        thisdx = alldx * math.cos(thisy)
        gradyi, gradxi = np.gradient(thisgrid, alldy, thisdx)
        
        # add first two lines to gradient if first loop
        if len(Xgrad) < 1:
            Xgrad = gradxi[0:2, :]
            Ygrad = gradyi[0:2, :]
        
        # otherwise, add just this row to the gradient array
        else:
            Xgrad = np.vstack((Xgrad, gradxi[1, :]))
            Ygrad = np.vstack((Ygrad, gradyi[1, :]))

    # add the last row to the gradient array
    Xgrad = np.vstack((Xgrad, gradxi[2, :]))
    Ygrad = np.vstack((Ygrad, gradyi[2, :]))

    # Get gradient magnitude
    Maggrid = np.sqrt((Ygrad**2)+(Xgrad**2))

    # Define a grid file that is the direction perpendicular to the max gradient
    Strikegrid = np.degrees(np.arctan2(Ygrad, Xgrad))
    Strikegrid = np.where(Strikegrid < 0, Strikegrid+360, Strikegrid)

    # Assign strike and dip arrays to grids with same dimensions as depth grid
    Dipgrid = np.degrees(np.arctan2(Maggrid, 1))

    # flip grids upside down if needed
    if flipornot == 'flip':
        Strikegrid = np.flipud(Strikegrid)
        Dipgrid = np.flipud(Dipgrid)

    return Strikegrid, Dipgrid
