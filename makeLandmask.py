#function to make landmask to load into particle tracking, instead of making it on the fly
import numpy as np
import pylab as p
import netCDF4 as nc
import xarray as xr
import zarr as zr
import matplotlib.dates as mdates
import datetime as dt

#velocity data...
dataLoc = '/data/ripBig/pringle/mercator/'
testdate = dt.datetime(2012,3,15)
testDay = mdates.date2num(testdate)- mdates.date2num(np.datetime64('0000-12-31'))
uDatName = dataLoc+'umerc_phy_%0.0d.nc'%(testDay)
vDatName = dataLoc+'vmerc_phy_%0.0d.nc'%(testDay)

#load data
uDat = nc.Dataset(uDatName,'r')
vDat = nc.Dataset(vDatName,'r')

#extract masks
uMsk = uDat['vozocrtx'][0,0:36,:,:].mask
vMsk = vDat['vomecrty'][0,0:36,:,:].mask

bMsk = np.logical_and(uMsk,vMsk).astype(int)

svMsk = np.savez('landmask.npz',landmask = bMsk)
