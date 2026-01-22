import numpy as np
import pylab as p
import netCDF4 as nc
import pandas as pd
import datetime as dt
import time
import pickle
import zarr
import os
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import glob
import gc
from sklearn.neighbors import NearestNeighbors
from tqdm import tqdm

print('test...')
print('setting up kdtree to find nearest grid cell')
tick = time.time()
gridFile = '/data/guppy2/willlush/Mercator/cGrid/MeshFiles/ext-PSY4V3R1_mesh_hgr.nc'

gDat = nc.Dataset(gridFile,'r') #load data
lonGrid = gDat['glamt'][0,:,:].data #t-point lons (grid cell center)
latGrid = gDat['gphit'][0,:,:].data #t-point lats (grid cell center)
gDat.close()
lonArr = np.ravel(lonGrid)
latArr = np.ravel(latGrid)
assert False
gx = np.array(list(zip(lonArr,latArr)))
nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(gx)
rtime = time.time()-tick
print('Done in %s s'%(rtime))
gx = None

print('Loading trajectory data')
tick = time.time()
#where is data coming from
dataDir = '/data/rip/willlush/parcelsOutput/combined/'

#where to save data
saveDir = '/data/rip/willlush/parcelsOutput/startEnd/'

#parameters (yearlist, depths)
yearList = np.arange(2007,2024)

alreadyRun = [(2007,1),
              (2007,10),
              (2007,20),
              (2007,30),
              (2008,1),
              (2008,10),
              (2008,20),]

#year = 2017
depthList = [1,10,20,30]
for year in yearList:
    for depth in depthList:
        if (year,depth) in alreadyRun:
            print('%s %sm already run'%(year,depth))
            continue
        dataName = dataDir+'traj_depthKeep_%s_%02dm.nc'%(year,depth)
        print('%s, %sm'%(year,depth))

        #load data
        dat = nc.Dataset(dataName,'r')

        age = dat['age'][:].data
        #trim nans at end:
        nMsk = ~np.all(np.isnan(age),axis=0) #where entire column is zero

        tArr = dat['time'][:].data
        lat = dat['lat'][:].data
        lon = dat['lon'][:].data
        dat.close()

        #apply nan mask
        age = age[:,nMsk]
        tArr = tArr[:,nMsk]
        lat = lat[:,nMsk]
        lon = lon[:,nMsk]

        #replace nans with invalid values:
        age[np.isnan(age)] = 1000*24*3600 #total seconds in 1000 days - longer than any runs
        tArr[np.isnan(tArr)] = 1000*24*3600
        lat[np.isnan(lat)] = 0.0
        lon[np.isnan(lon)] = 0.0

        rtime = time.time()-tick
        print('Done in %s s'%(rtime))

        print('Querying kdtree to find grid cells')
        tick = time.time()
        X = np.array(list(zip(np.ravel(lon),np.ravel(lat))))
        distances, indices = nbrs.kneighbors(X)
        rtime = time.time()-tick
        print('Done in %s s (%s min)'%(rtime,round(rtime/60.,2)))
        distances=None
        X = None
        
        indices = np.ravel(indices)
        ind = indices.reshape(np.shape(lat))
        stInd = ind[:,0]
        retInd = np.isin(indices,stInd)
        retInd = retInd.reshape(np.shape(lat))
        rAge = (tArr.T-tArr[:,0]).T

        tick=time.time()

        #saveDir = '/data/break/willlush/behavior/cmTest/'
        saveName = 'startEnd_%s_%02dm.nc'%(year,depth)
        print('creating netcdf...')
        ncout = nc.Dataset(saveDir+saveName,'w')
        ncout.createDimension('trajectory',ind.shape[0])
        ncout.createDimension('obs',ind.shape[1])
        ncout.createDimension('tPoints',len(latArr))
        ncout.createDimension('startYear',1)

        #metadata (trying to follow CF conventions)
        ncout.title = 'Trajectory/return data for %s'%(year)
        ncout.institution = 'University of New Hampshire (UNH)'
        ncout.author = 'William Lush (UNH)'
        ncout.source = 'created using Parcels (https://oceanparcels.org/) within 1/12 deg. velocity fields obtained from the Mercator Ocean global physical model (PSY4V3R1), obtained from the Copernicus Marine Environmental Monitoring Service (https://marine.copernicus.eu/)'
        ncout.history = 'created: %s'%(dt.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

        #grid cell center locations
        gLat = ncout.createVariable('tLat',
                                    'f8',
                                    'tPoints',
                                    compression='zlib',
                                    complevel=8)
        gLat.long_name = 't-point latitude'
        gLat.standard_name = 'latitude_of_grid_cell_center'
        gLat.units = 'degrees_north'

        gLon = ncout.createVariable('tLon',
                                    'f8',
                                    'tPoints',
                                    compression='zlib',
                                    complevel=8)
        gLon.long_name = 't-point longitude'
        gLon.standard_name = 'longitude_of_grid_cell_center'
        gLon.units = 'degrees_east'

        #particle/trajectory data
        pLat = ncout.createVariable('lat',
                                    'f8',
                                    ('trajectory','obs'),
                                    compression='zlib',
                                    complevel=8)
        pLat.long_name = 'particle latitude'
        pLat.standard_name = 'latitude_of_simulated_particle'
        pLat.units = 'degrees_north'
        pLon = ncout.createVariable('lon',
                                    'f8',
                                    ('trajectory','obs'),
                                    compression='zlib',
                                    complevel=8)
        pLon.long_name = 'particle longitude'
        pLon.standard_name = 'longitude_of_simulated_particle'
        pLon.units = 'degrees_east'
        pAge = ncout.createVariable('age','f8',
                                    ('trajectory','obs'),
                                    compression='zlib',
                                    complevel=8)
        pAge.long_name = 'particle age'
        pAge.standard_name = 'seconds_since_particle_release'
        pAge.units = 's'
        pTime = ncout.createVariable('time',
                                     'f8',
                                     ('trajectory','obs'),
                                     compression='zlib',
                                     complevel=8)
        pTime.long_name = 'time at observation'
        pTime.standard_name = 'seconds_since_year_start'
        pTime.units = 'seconds since yyyy-1-1 0:0:0'
        pIndex = ncout.createVariable('index',
                                      'u4',
                                      ('trajectory','obs'),
                                      compression='zlib',
                                      complevel=8)
        pIndex.long_name = 'particle index'
        pIndex.standard_name = 'index_in_flattened_grid_position_arrays'
        pMask = ncout.createVariable('retMask',
                                     'u1',
                                      ('trajectory','obs'),
                                      compression='zlib',
                                      complevel=8)
        pMask.long_name = 'in-domain mask'
        pMask.standard_name = 'mask_to_indicate_if_in_near_shore_grid_cells'

        pSt = ncout.createVariable('startInd',
                                   'u4',
                                   'trajectory',
                                   compression='zlib',
                                   complevel=8)
        pMask.long_name = 'starting indices'
        pMask.standard_name = 'array_of_start_indices'

        gLon[:] = lonArr
        gLat[:] = latArr

        pLon[:] = lon
        pLat[:] = lat
        pAge[:] = rAge
        pTime[:] = tArr
        pIndex[:] = ind
        pMask[:] = retInd.astype(int)
        pSt[:] = stInd

        indices = None
        ind = None
        retInd = None
        rAge = None
        tArr = None
        lat = None
        lon = None
        nMsk = None

        ncout.close()
        rtime = time.time()-tick
        print('Done in %s s (%s min)'%(rtime,round(rtime/60.,2)))
        gc.collect()

# assert False
# #testing...
# pld = 10 #pld in days
# plds = pld*24*3600
# aMsk = rAge==plds
# endInd = np.ravel(ind[aMsk])
# stTime = tArr[:,0]
# endTime = tArr[aMsk]

# retMsk = np.isin(endInd,stInd)
# retSt = stInd[retMsk]
# retStTime = stTime[retMsk]
# tRet = endTime[retMsk]
# retEnd = endInd[retMsk]
# noRet = stInd[~retMsk]

# if False:
#     p.figure(1,figsize=(9.0,9.0))
#     p.clf()
#     ax = p.axes(projection=ccrs.PlateCarree())
#     ax.set_global()
#     #ax.add_feature(cartopy.feature.LAND,zorder=50)
#     p.plot(lonArr[stInd],latArr[stInd],'.')
#     p.plot(lonArr[endInd[retMsk]],latArr[endInd[retMsk]],'.')
#     ax.add_feature(cartopy.feature.COASTLINE)
#     p.show()



