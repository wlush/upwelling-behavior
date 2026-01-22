import numpy as np
import pylab as p
import cartopy.crs as ccrs
import cartopy
import netCDF4 as nc
from tqdm import tqdm
from scipy.signal import savgol_filter
from shapely.geometry import LineString,MultiPoint, Point

#load lat/lon for model domain:
ltp = np.load('../tPoints.npz')
pLon = ltp['lon']
pLat = ltp['lat']

frd = np.load('fRet_and_nSet.npz',allow_pickle=True)['fRet'].item()

dl = np.load('maxRv.npz',allow_pickle=True)
maxDict = dl['maxDict'].item()

pld=30
#pld = 60
depth = 1

plotTuples = []


from scipy.stats import linregress


avgMos1 = [4,5,6,7,8]
avgMos2 = [10,11,12,1,2]
moArr = [avgMos1, avgMos2]

rcmb = []
fdb = []
for ix in [0,1]:
    avgMos=moArr[ix]
    fd = np.zeros(len(pLon))
    rcMin = np.zeros(len(pLon))
    for mo in tqdm(avgMos):
        fDiff = np.array(frd[(mo,pld,20)])-np.array(frd[(mo,pld,1)])
        modRc = np.array(maxDict[(pld,mo,100,20)])
        modRc[modRc<0]=np.nan

        fd+=fDiff
        rcMin+=modRc

    fd = np.array(fd)/len(avgMos)
    rcMin = np.array(rcMin)/len(avgMos)

    rmNan3 = np.logical_and(~np.isnan(rcMin),~np.isnan(fd))
    #rcf_trend = trd[0]+trd[1]*np.log(np.sort(rcMin[rmNan3]))
    rcmb.extend(rcMin[rmNan3])
    fdb.extend(fd[rmNan3])

lnf3 = np.polyfit(np.log(rcmb),fdb,1.)
slope, intercept, r_value, p_value, std_err = linregress(rcmb, fdb)
rcf_trend = np.log(np.sort(rcmb))*lnf3[0]+lnf3[1]
    
p.figure(figsize=(8.0,6.0))
p.clf()
p.title("$R'_c$ vs. $\Delta F_{ret}$",fontsize='x-large')    
p.plot(rcmb,fdb,'.',markersize=2)
p.plot(np.sort(rcmb),rcf_trend,c='k',alpha=.75,linestyle='dashed')
p.ylabel('$\Delta F_{ret}$',fontsize='large')
p.xlabel("$R'_c$",fontsize='large')
p.tight_layout()
p.show()
p.savefig('rcVsFDiff_defense.png')
