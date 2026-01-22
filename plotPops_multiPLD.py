import numpy as np
import pylab as p
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cft
import netCDF4 as nc
from shapely.geometry import LineString, Point
from tqdm import tqdm

#load lat/lon for model domain:
ltp = np.load('../tPoints.npz')
pLon = ltp['lon']
pLat = ltp['lat']

#square 1
sq1LonMask = np.logical_and(pLon>-169.268,pLon<-156.660)
sq1LatMask = np.logical_and(pLat>57.952,pLat<68.983)
sq1Mask = np.logical_and(sq1LonMask,sq1LatMask)

#square 2
sq2LonMask = np.logical_and(pLon>-80.739,pLon<-76.826)
sq2LatMask = np.logical_and(pLat>3.027,pLat<9.639)
sq2Mask = np.logical_and(sq2LonMask,sq2LatMask)

excMask = np.logical_or(sq1Mask,sq2Mask)

rotated_crs = ccrs.RotatedPole(pole_longitude=120.0, pole_latitude=45.0)
rotated = rotated_crs.transform_points(src_crs=ccrs.PlateCarree(),x=pLon,y=pLat)
rtX,rtY = rotated[:,0],rotated[:,1]

minLat = np.amin(pLat)-2.5
maxLat = np.amax(pLat)
minLon = -179.
maxLon = np.amax(pLon)+2.5

#Blanchette edge dict (upperLat, upperLon, lowerLat, lowerLon)
bbd = {'Aleutian':(58.271,-136.731,57.154,-135.755),
       'Columbian':(54.078,-131.795,48.853,-125.183),
       'Mendocinian':(48.171,-124.736,36.948,-122.065),
       'Montereyan':(36.621,-121.907,34.442,-120.456),
       'Southern Californian':(34.467,-120.278,34.012,-118.793),
       'Ensenadian':(33.715,-118.236,29.949,-115.812),
       'Magdalenian':(27.726,-114.995,26.705,-113.576)}

def blanchette(lmName,lmTuple,axs):
    upperLat = lmTuple[0]
    upperLon = lmTuple[1]
    lowerLat = lmTuple[2]
    lowerLon = lmTuple[3]
    trans = rotated_crs.transform_points(src_crs=ccrs.PlateCarree(),x=np.array([upperLon,lowerLon]),y=np.array([upperLat,lowerLat]))
    t2 = rotated_crs.transform_points(src_crs=ccrs.PlateCarree(),x=np.array([-170.,-170.]),y=np.array([upperLat,lowerLat]))

    upperLon_r = trans[0,0]
    lowerLon_r = trans[1,0]
    upperLat_r = trans[0,1]
    lowerLat_r = trans[1,1]
    blRectLon = [upperLon_r+2.5,upperLon_r-2.5,lowerLon_r-2.5,lowerLon_r+2.5,upperLon_r+2.5]
    blRectLat = [upperLat_r,upperLat_r,lowerLat_r,lowerLat_r,upperLat_r]
    axs.plot(blRectLon,blRectLat,color='k',linestyle='dashed',zorder=100,linewidth=1)
    return()

#load population dictionary:
pdl = np.load('../rValTest.npz',allow_pickle=True)
popDict = pdl['popDict'].item()

pld=30
depth = 1

avgMos1 = [4,5,6,7,8]
am1 = 'upwelling'
avgMos2 = [10,11,12,1,2]
am2 = 'non-upwelling'

moArr = [avgMos1,avgMos2]

from matplotlib.gridspec import GridSpec
fig = p.figure(figsize=(5.5,8),constrained_layout=True)
gs = GridSpec(2, 3, figure=fig)

ax0 = fig.add_subplot(gs[0, 0], projection=rotated_crs)
ax1 = fig.add_subplot(gs[0, 1], projection=rotated_crs)
ax2 = fig.add_subplot(gs[0, 2], projection=rotated_crs)
ax3 = fig.add_subplot(gs[1, 0], projection=rotated_crs)
ax4 = fig.add_subplot(gs[1, 1], projection=rotated_crs)
ax5 = fig.add_subplot(gs[1, 2], projection=rotated_crs)

axList = [ax0,ax1,ax2,ax3,ax4,ax5]
for ix in np.arange(6):
    ax = axList[ix]
    ax.set_extent([minLon+117.5,maxLon+55,minLat-40,maxLat-17.5])
    ax.add_feature(cft.COASTLINE,linewidth=.5,alpha=.5)
    ax.add_feature(cft.BORDERS, linestyle='--',linewidth=.5)
    ax.add_feature(cft.LAND,alpha=.5)
    ax.add_feature(cft.OCEAN,alpha=.5)
    ax.gridlines(draw_labels=False, dms=True,linewidth=.5,x_inline=False)
    ax.plot([-45.13,-43.1],[-8.3,-12.4],'r--')
    ax.plot([-37.5,-41.],[-20.6,-26],'r--')
    ax.plot([-29.7,-34.5],[-30.2,-32.7],'r--')
    ax.plot([-48.8,-51.4],[45.3,40.2],'r--')
    for bb in bbd.keys():
        bbv = bbd[bb]
        blanchette(bb,bbv,ax)

axPairs = [(0,3),(1,4),(2,5)]
pldList = [5,30,60]
for axP in axPairs:
    ax_0 = axList[axP[0]]
    pld = pldList[axP[0]]    
    for ix in [0]:
        avgMos=moArr[ix]
        popAvg = np.zeros(len(pLon))
        for mo in tqdm(avgMos):
            pop = np.array(popDict[(pld, mo, 100, 1, 20)])
            popAvg+=pop

    popAvg = np.array(popAvg)/len(avgMos)
    po = np.argsort(popAvg[~excMask])

    ax_0.set_title('Upwelling (PLD=%s)'%(pld),fontsize='small')
    ax_0.plot(rtX[excMask],rtY[excMask],'.',markersize='1.5',color='grey')
    sct = ax_0.scatter(rtX[~excMask][po],rtY[~excMask][po],c=np.log10(popAvg[~excMask][po]),
                      s=5,alpha=.5,vmin=-4,vmax=-1,cmap='gnuplot2')

    ax_1 = axList[axP[1]]
    for ix in [1]:
        avgMos=moArr[ix]
        popAvg = np.zeros(len(pLon))
        for mo in tqdm(avgMos):
            pop = np.array(popDict[(pld, mo, 100, 1, 20)])
            popAvg+=pop

        popAvg = np.array(popAvg)/len(avgMos)

    po = np.argsort(popAvg[~excMask])

    ax_1.set_title('Non-upwelling (PLD=%s)'%(pld),fontsize='small')
    ax_1.plot(rtX[excMask],rtY[excMask],'.',markersize='1.5',color='grey')
    sct = ax_1.scatter(rtX[~excMask][po],rtY[~excMask][po],c=np.log10(popAvg[~excMask][po]),
                      s=5,alpha=.5,vmin=-4,vmax=-1,cmap='gnuplot2')

p.show()    
#p.savefig('multiFig_popsPLD_noCBar.png')
