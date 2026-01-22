import numpy as np
import pylab as p
import pandas as pd
import datetime as dt
import matplotlib.dates as mdates
import xarray as xr
import netCDF4 as nc
import time
import dask
import os
import cartopy.crs as ccrs
import cartopy
from glob import glob
import parcels
import sys
from tqdm import tqdm
from parcels import AdvectionRK4_3D
from parcels import FieldSet
from parcels import JITParticle, Variable
from parcels import ParticleFile
from parcels import ParticleSet
from parcels.tools import StatusCode

def getDataNames(firstReleaseDate,lastReleaseDate,trackDuration,dataDir,pad = 1.0):
    #make a constant padRuns of about 1 day to avoid issues with the last time
    padRuns = dt.timedelta(days=pad)
    duration = dt.timedelta(days=trackDuration)
    testPad = dt.timedelta(days=1)
    
    print('making list of filenames for fieldset')
    numDays=((lastReleaseDate-firstReleaseDate)+duration+padRuns).days
    fnameU=[]
    fnameV=[]
    fnameW=[]

    for date in pd.date_range(firstReleaseDate,lastReleaseDate+duration+testPad): #yay pandas date_range!
        thisDay=mdates.date2num(date)- mdates.date2num(np.datetime64('0000-12-31')) #sigh, make into old matlab date format
        #print('add %d to days to look at'%(thisDay,))
        fnameU.append(dataDir+'umerc_phy_%0.0d.nc'%(thisDay,))
        fnameV.append(dataDir+'vmerc_phy_%0.0d.nc'%(thisDay,))
        fnameW.append(dataDir+'wmerc_phy_%0.0d.nc'%(thisDay,))
        assert os.path.exists(fnameU[-1]),fnameU[-1]+' does not exist?'
        assert os.path.exists(fnameV[-1]),fnameV[-1]+' does not exist?'
        assert os.path.exists(fnameW[-1]),fnameW[-1]+' does not exist?'
        
    return(fnameU,fnameV,fnameW,padRuns,duration)


    
def makeFieldset(fnameU,fnameV,fnameW):
    #mask grid location
    maskFile='/data/guppy2/willlush/Mercator/cGrid/MeshFiles/ext-PSY4V3R1_mesh_hgr.nc'
    vmFile = '/data/rip/willlush/parcelsAddFiles/bottomDepth.nc'

    filenames={'U':{'lon':maskFile,'lat':maskFile,'depth':fnameW[0],'data':fnameU},
               'V':{'lon':maskFile,'lat':maskFile,'depth':fnameW[0],'data':fnameV},
               'W':{'lon':maskFile,'lat':maskFile,'depth':fnameW[0],'data':fnameW},}
               #'maxDepth':{'lon':maskFile,'lat':maskFile,'data':vmFile}}
    variables={'U':'vozocrtx','V':'vomecrty','W':'vovecrtz',}#'maxDepth':'bdepth'}
    dimensions = {'lon':'glamf','lat':'gphif','depth':'depthw','time':'time_counter'}
                  
    indices={'depth':list(range(0,35))}

    #cs={'time':('time_counter',1),'depth':('depthw',1),'lat':('gphif',512),'lon':('glamf',512),} #default for now
    #cs={'time':('time_counter',1),'depth':('depthw',1),'lat':('gphif',1024),'lon':('glamf',512),}
    cs={'time':('time_counter',1),'depth':('depthw',8),'lat':('gphif',1024),'lon':('glamf',1024),}
    fset=FieldSet.from_nemo(filenames,variables,dimensions,indices=indices,chunksize=cs)
    print('made fset with',fnameU[0])

    dVars = {'maxDepth':'bdepth'}
    dDims = {'lon':'lon','lat':'lat','depth':'depth'}
    fset.add_field(
        parcels.field.Field.from_netcdf(vmFile, variable=dVars, dimensions=dDims, deferred_load=False))
    lmLoad = np.load('/data/rip/willlush/parcelsOutput/landmask.npz')
    lMask = lmLoad['landmask']#makeLandMask(fnameU[0],fnameV[0])
    coords = xr.open_dataset(maskFile, decode_times=False)
    fset.add_field(
        parcels.Field("landmask",
                      data=lMask,
                      lon=coords["glamf"],
                      lat=coords["gphif"],
                      depth=coords["nav_lev"][0:36],
                      mesh="spherical",))    
    return(fset)

def runFull(lonStarts,latStarts,dateStarts,releaseDepths,
            outFileName,outputFrequencyHours,trackDuration,
            keepDepth,dataDir,is3D = True):

    firstReleaseDate=np.amin(dateStarts)
    lastReleaseDate=np.amax(dateStarts)

    #make lists of filenames to be used:
    nameU,nameV,nameW,padRuns,duration = getDataNames(firstReleaseDate,lastReleaseDate,trackDuration,dataDir=dataDir)
    
    #make fieldset:
    fset = makeFieldset(nameU,nameV,nameW)
    nParts=len(dateStarts)
    timeReleaseVec=[]
    #make lists of release locations and times, all the same length
    #(hence all the stuff with arange(). Make timebase the same as in
    #the model files. NOT SURE IT RELEASES AT RIGHT TIME IN DAY
    firstTime=fset.U.grid.time[0]
    timeReleaseVec=[firstTime+(releaseTime-firstReleaseDate).total_seconds() for releaseTime in dateStarts]
    print('done making timeReleaseVec')
    
    #add a kernel that kills the particle after durationDays
    fset.add_constant('maxage',duration.total_seconds())
    def SampleAge(particle,fieldset,time):
        particle.age=particle.age+math.fabs(particle.dt)
        if particle.age>fieldset.maxage:
            particle.delete()

    #make a new kind of particle that has an age
    class AgeParticle(JITParticle):
        age=Variable('age',initial=0.,dtype=np.float32)
        bd=Variable('botDepth')
        on_land=Variable('on_land')

    #release depth based on depth-keeping rules
    #shortcut, assuming release points are centered in grid cell (should be the case...)
    rdLon = fset.maxDepth.lon.ravel()
    rdLat = fset.maxDepth.lat.ravel()
    rdList = list(zip(rdLon,rdLat))
    rdDict = dict(zip(rdList,fset.maxDepth.data[0].ravel()))

    rDepths = np.array(releaseDepths)
    rLon = np.array(lonStarts)
    rLat = np.array(latStarts)
    rll = list(zip(rLon,rLat))
    bottom = np.array([rdDict[x] for x in rll])
    bMsk = rDepths>=bottom #where particles started on/below bottom
    rDepths[bMsk] = bottom[bMsk]-2.5

    pset=ParticleSet.from_list(fieldset=fset,
                               pclass=AgeParticle,
                               lon=np.array(lonStarts),
                               lat=np.array(latStarts),
                               depth=rDepths, #np.array(releaseDepths),
                               time=np.array(timeReleaseVec),)

    print('made pset with ',len(releaseDepths),'particles')

    #now integrate for some days, save output in outFile
    #function for particle to maintain depth
    fset.add_constant('keepDepth',keepDepth) #add constant for depth keeping
    fset.add_constant('offBottom',2.5)#2.5m from bottom
    fset.add_constant('landVal',1.)
    fset.add_constant('backUnder',.5)
    
    #if keepdepth is too deep, place particle 2.5m above bottom
    #else automatically set depth to 'keepDepth' meters
    def DepthKeeping(particle, fieldset, time):
        if keepDepth>=fieldset.maxDepth[time,
                                        particle.depth,
                                        particle.lat,
                                        particle.lon]:
            particle.depth=particle.depth=fieldset.maxDepth[time,
                                                            particle.depth,
                                                            particle.lat,
                                                            particle.lon] - fieldset.offBottom
                
        else:
            particle.depth = fieldset.keepDepth

    # def Sample_land(particle, fieldset, time):
    #     particle.on_land = fieldset.landmask[time,
    #                                          particle.depth,
    #                                          particle.lat,
    #                                          particle.lon]

    def beached(particle, fieldset, time):
        if particle.on_land == fieldset.landVal:
            particle.delete()
    
    def SampleDepth(particle, fieldset, time):
        particle.botDepth = fieldset.maxDepth[time,particle.depth,particle.lat,particle.lon]

    def DeleteParticle(particle, fieldset, time):
        if (particle.state == StatusCode.ErrorOutOfBounds) or (particle.state == StatusCode.ErrorThroughSurface):
            particle.delete()

    outFile=pset.ParticleFile(name=outFileName,
                              outputdt=dt.timedelta(hours=outputFrequencyHours),chunks=(60000,10)) #chunksize

    print(lastReleaseDate-firstReleaseDate+duration+padRuns)
    #[AdvectionRK4_3D,SampleAge,DepthKeeping,DeleteParticle,SampleDepth,Sample_land,beached]
    if is3D:
        pset.execute([AdvectionRK4_3D,SampleAge,DepthKeeping,beached,DeleteParticle],
                     runtime=lastReleaseDate-firstReleaseDate+duration+padRuns,
                     dt=dt.timedelta(minutes=30.0),output_file=outFile,
                     verbose_progress=False,)
                     #recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    else:
        assert False
    return()

if __name__ == '__main__':
    #generate releases
    assert len(sys.argv)==4, 'Need 3 arguments: year, depth, and machine'
    year = int(sys.argv[1])
    machine = sys.argv[2]
    dpth = int(sys.argv[3])

    #year=2007;machine='rip';dpth=50

    if machine=='rip':
        print('running on rip')
        msDir = '/data/ripBig/pringle/mercator/' #machine-specific directories
        saveLoc = '/data/rip/willlush/parcelsOutput/tempFiles/'
    elif machine=='surge':
        print('running on surge')
        msDir = '/data/surge/pringle/mercator/'
        saveLoc = '/data/surge/willlush/parcels/tempFiles/'
    else:
        assert False,'improper/no machine specified'
        
    print(year,dpth,machine)

    
    rpName = '/data/breakhome/willlush/workfiles/behavior_upwelling/particle_tracking/testing_v1/2coastalReleasePts.npz'
    rLocLoad = np.load(rpName)
    rLon = rLocLoad['lons']
    rLat = rLocLoad['lats']
    
    #lmax = np.logical_and(rLat>65.,rLat<68.25)
    lmax = rLat<68.5
    #lmax = rLat>65.
    rLon = rLon[lmax]
    rLat = rLat[lmax]

    #dates over which to release particles
    timeStart = dt.datetime(year,1,1,0,0)
    timeEnd = dt.datetime(year,12,31,0,0)
    timeArr = pd.date_range(timeStart,timeEnd)
    #timeArr = timeArr[0::3]

    #print(timeStart)
    lonRel = []
    latRel = []
    timeRel = []
    depthRel = []
    for ix in np.arange(len(timeArr)):
        lonRel.extend(rLon)
        latRel.extend(rLat)
        timeRel.extend(np.full(len(rLon),timeArr[ix]))
        depthRel.extend(np.full(len(rLon),dpth))

    lonRel = np.array(lonRel)
    latRel = np.array(latRel)
    timeRel = np.array(timeRel)
    depthRel = np.array(depthRel)

    #saveName = saveLoc+'TEST_year_keepDepth_tempProc_%s_%02d'%(year,dpth)
    saveName = saveLoc+'temp_depthKeep_%s_%02d'%(year,dpth)
    runFull(lonRel,latRel,timeRel,depthRel,saveName,12.0,60,keepDepth=dpth,dataDir=msDir)
    print('done with particle tracking, combine files if parallel')
