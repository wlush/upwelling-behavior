import numpy as np
import pylab as p
import netCDF4 as nc
import xarray as xr
import dask
import zarr
import cartopy.crs as ccrs
import cartopy
from tqdm import tqdm

#load lat/lon for model domain:
ltp = np.load('../tPoints.npz')
pLon = ltp['lon']
pLat = ltp['lat']

fRetDict = {}
nSetDict = {}
monthList = np.arange(1,13)
for pld in tqdm(np.arange(5,65,5)):
    for mo in monthList:
        for depth in [1,10,20,30]:
            dataLoc = '/data/break/willlush/behavior/connectivityMatrices/'
            dn = 'connectivity_pld%02d_mo%02d_%02dm.nc'%(pld,mo,depth)
            dat = xr.open_dataset(dataLoc+dn)
            retNum = dat['cm_raw'].sum(axis=1).values
            retFrac = retNum/dat['cm_counts'].values[1:-1]
            nSet = dat['cm_counts'].sum(axis=0).values
            fRetDict[(mo,pld,depth)]=retFrac
            nSetDict[(mo,pld,depth)]=nSet

np.savez('fRet_and_nSet.npz',nSet=nSetDict,fRet = fRetDict)
