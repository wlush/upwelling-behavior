import numpy as np
import pylab as p
import netCDF4 as nc
import datetime as dt
import xarray as xr
from tqdm import tqdm
assert False
#where are start/end files?
seDir = '/data/rip/willlush/parcelsOutput/startEnd/'

#parameters (PLD, depth)
depth = 30
print('depth is %sm'%(depth))
pldList = np.arange(1,61)
#pldList = [5,10,15,20,25,30,35,40,45,50,55,60]

yearList = np.arange(2007,2024) #current progress...
#yearList = [2007,2008]
#moRange = [1,2,3]
moRange = np.arange(1,13)

#load preliminary data from first file...
tdn = 'startEnd_%s_%02dm.nc'%(2007,depth)
td = nc.Dataset(seDir+tdn,'r')

print('loading t-point locations, building empty cm')
lonArr = td['tLon'][:].data
latArr = td['tLat'][:].data
sInd = td['startInd'][:].data

stLen = len(np.unique(td['startInd'][:].data))
td.close()
bCm = np.zeros((stLen-2,stLen)) #cm size, minus 1st and last
rCt = np.zeros(stLen)
print('done')

for mo in moRange:
    for pld in pldList:
        print('month is %s, depth is %sm, PLD is %s days'%(mo,depth,pld))
        starts = False
        for year in yearList:
            #load data:
            dName = 'startEnd_%s_%02dm.nc'%(year,depth)
            dat = nc.Dataset(seDir+dName,'r')

            #calculate time in seconds for connectivity starts we're exploring:
            stTime = dt.datetime(year,1,1,0,0,0)
            endTime = dt.datetime(year,12,31,0,0,0)+dt.timedelta(hours=12) #time elapsed
            t = np.arange(stTime, endTime, dt.timedelta(hours=12)).astype(dt.datetime)
            months = [mo,]
            tEl_mod = np.array([x.month in months for x in t])
            tEl = t-stTime
            whatTimeStart = np.array([x.total_seconds() for x in tEl[tEl_mod]]) #in seconds

            getTime = dat['time'][:].data
            startMask = np.isin(getTime[:,0],whatTimeStart)

            getTime = getTime[startMask]
            getAge = dat['age'][:].data[startMask]
            getInd = dat['index'][:].data[startMask]
            getMsk = dat['retMask'][:].data[startMask]
            dat.close()
            maxTime = np.amax(getAge,axis=1) #max time (not on 12h when killed)

            tEnd = pld*3600.*24.
            eMsk = getAge==tEnd
            killMask = ~np.all(~eMsk,axis=1)
            endInd = getInd[eMsk]

            #indIt = list(set(getInd))
            stIx,stCts = np.unique(getInd[:,0],return_counts=True)
            st2arr = dict(zip(stIx,np.arange(len(stIx))))
            arr2st = dict(zip(np.arange(len(stIx)),stIx))
            if starts==False:
                indOrd = stIx.copy()
                starts=True

            tvStack = []
            for ix in tqdm(np.arange(1,len(stIx)-1)): #remove 1st and last point for sanity
                tVec = np.zeros(len(stIx)) #width of full domain
                startIx = arr2st[ix]
                stMsk = getInd[killMask,0]==startIx
                endLocs = endInd[stMsk]
                for lc in endLocs:
                    if lc in stIx:
                        tVec[st2arr[lc]]+=1.
                tvStack.append(tVec)
            tvStack = np.array(tvStack)
            bCm += tvStack
            rCt += stCts
            print('done with %s'%(year))

        cmCts = bCm.copy()
        rCts = rCt.copy()
        cms = (bCm/rCt).copy()

        print('done, save to netcdf')
        print(cmCts.shape)
        #assert False
  
        #netcdf:
        #conn. mat., raw conn. mat., counts, lonArr, latArr
        #dims: cm = (pld,from,to)

        #cm axis 0 is from, axis 1 is to
        sDir = '/data/break/willlush/behavior/connectivityMatrices/'
        sName = sDir + 'bigConnectivity.nc'
        
        data_vars=dict(
            cm=([],


        
        ncout = nc.Dataset(sName,'w')
        ncout.createDimension('tPoints',len(lonArr))
        ncout.createDimension('whereFrom',cmCts.shape[1])
        ncout.createDimension('whereTo',cmCts.shape[0])
        ncout.createDimension('indArr',indOrd.shape[0])

        #metadata (trying to follow CF conventions)
        ncout.title = '%s-day trajectory/return data for %s, month %s, depth %s'%(pld, year, mo, depth)
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

        #indices...
        gInd = ncout.createVariable('indices',
                                    'f8',
                                    'indArr',
                                    compression='zlib',
                                    complevel=8)
        gLat.long_name = 'indice for lat-lon arrays in cm order'

        cm_raw = ncout.createVariable('cm_raw',
                                      'f8',
                                      ('whereTo','whereFrom'),
                                      compression='zlib',
                                      complevel=8)

        cm = ncout.createVariable('cm',
                                  'f8',
                                  ('whereTo','whereFrom'),
                                  compression='zlib',
                                  complevel=8)

        cm_counts = ncout.createVariable('cm_counts',
                                         'f8',
                                         'whereFrom',
                                         compression='zlib',
                                         complevel=8)

        
        cm_raw[:] = cmCts
        cm[:] = cms
        cm_counts[:] = rCts
        gInd[:] = indOrd
        gLon[:] = lonArr
        gLat[:] = latArr
        
        ncout.close()
