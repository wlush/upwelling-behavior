import numpy as np
import pylab as p
import netCDF4 as nc
import pandas as pd
import datetime as dt
import pickle
import zarr
import os
import xarray as xr
import cartopy.crs as ccrs
import cartopy
import glob
from tqdm import tqdm

#where connectivity matrices are located
cmLoc = '/data/rip/willlush/parcelsOutput/connectivity/'
saveLoc = cmLoc+'monthlyConnectivity/'

#params
depth = 30
pld = 30
#month = 2

moList = np.arange(12)+1

for month in moList:
    #name pattern
    cmList = []
    cmPattern = 'cm_%02d_*_%02dm_%sd.nc'%(month, depth, pld) #month, year, depth, pld
    cmListSub = glob.glob(cmLoc+cmPattern)    
    cmList.extend(cmListSub)


    for ix in tqdm(np.arange(len(cmList))):
        cmName = cmList[ix]
        cmLoad = nc.Dataset(cmName,'r')
        cm = cmLoad['connectivity'][:].data
        nr = cmLoad['numRel'][:].data

        if ix==0:
            cmAvg=cm
            nRel=nr
            cmLon = cmLoad['longitude'][:].data
            cmLat = cmLoad['latitude'][:].data
            cmMonth = month
            cmDepth = depth
            cmPld = pld
        else:
            cmAvg+=cm
            cmAvg=cmAvg/2.0
            nRel+=nr

    saveName = 'cmMonthly_%02d_%02dm_%02dd.nc'%(month, depth, pld)
    print('saving as netcdf...')
    ncout = nc.Dataset(saveLoc+saveName,'w')
    ncout.createDimension('nx',len(cmAvg[0,0,0,:,:]))
    ncout.createDimension('ny',len(cmAvg[0,0,0,:,:]))
    ncout.createDimension('month',1)
    ncout.createDimension('depth',1)
    ncout.createDimension('pld',1)

    moNames = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',
               6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',
               11:'Nov',12:'Dec'}

    #metadata (trying to follow CF conventions)
    ncout.title = '%s-day, %sm mean connectivity matrix for %s, averaged over 2007-2023'%(pld, depth, moNames[month])
    ncout.institution = 'University of New Hampshire'
    ncout.source = 'created using Parcels (https://oceanparcels.org/) within 1/12 deg. velocity fields obtained from the Mercator Ocean global physical model (PSY4V3R1), obtained from the Copernicus Marine Environmental Monitoring Service (https://marine.copernicus.eu/)'
    ncout.history = 'created: %s'%(dt.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))

    cmOut = ncout.createVariable('connectivity',
                                 'f8',
                                 ('month','depth','nx','ny'),
                                 compression='zlib',
                                 complevel=8)
    cmOut.long_name = 'connectivity matrix'

    lonOut = ncout.createVariable('longitude',
                                  'f8',
                                  'nx',compression='zlib',
                                  complevel=8)
    latOut = ncout.createVariable('latitude',
                                  'f8',
                                  'nx',
                                  compression='zlib',
                                  complevel=8)
    numRel = ncout.createVariable('numRel',
                                  'u8',
                                  ('month','depth','nx'),
                                  compression='zlib',
                                  complevel=8)
    numRel.long_name='number of particles released'
    mmonth = ncout.createVariable('mmonth',
                                  'u8',
                                  'month',
                                  compression='zlib',
                                  complevel=8)
    mdepth = ncout.createVariable('mdepth',
                                  'u8',
                                  'depth',
                                  compression='zlib',
                                  complevel=8)
    mdepth.units='m'
    mpld = ncout.createVariable('mpld',
                                'u8',
                                'pld',
                                compression='zlib',
                                complevel=8)

    cmOut[0,0,:,:] = cmAvg[0,0,0,:,:]
    lonOut[:] = cmLon
    latOut[:] = cmLat
    numRel[0,0,:] = nRel[0,0,0,:]
    mmonth[:] = cmMonth
    mdepth[:] = cmDepth
    mpld[:] = cmPld
    ncout.close()
    print('done with %s'%(moNames[month]))

print('done, time to run for next pld/depth (delete this statement if looping...')

