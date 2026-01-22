import numpy as np
import pylab as p
import pandas as pd
import cartopy.crs as ccrs
import cartopy
import cartopy.feature as cfeature
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

def blanchette(lmName,lmTuple,axs,lab=False):
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
    if lab==True:
        anName = lmName
        anLon = (lowerLon_r+upperLon_r)/2.
        anLat = (lowerLat_r+upperLat_r)/2.
        axs.annotate(anName,xy=(anLon,anLat), xytext=(anLon-5.,anLat),horizontalalignment='right',va='center',arrowprops=dict(arrowstyle='-',facecolor='black'))
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


fig = p.figure(figsize=(8,7))
p.clf()

for ix in [0]:
    avgMos=moArr[ix]
    popAvg = np.zeros(len(pLon))
    for mo in tqdm(avgMos):
        pop = np.array(popDict[(pld, mo, 100, 1, 20)])
        popAvg+=pop

    popAvg = np.array(popAvg)/len(avgMos)

po = np.argsort(popAvg[~excMask])
    
ax1 = fig.add_subplot(1, 2, 1, projection=rotated_crs)
ax1.set_extent([minLon+117.5,maxLon+55,minLat-40,maxLat-17.5]) #manually figured out borders...
#ax1.set_global()
ax1.add_feature(cfeature.COASTLINE,linewidth=.5,alpha=.5)
ax1.add_feature(cfeature.BORDERS, linestyle='--',linewidth=.5)
ax1.add_feature(cfeature.LAND,alpha=.5)
ax1.add_feature(cfeature.OCEAN,alpha=.5)
ax1.set_title('Upwelling\n(Apr.-Aug.)')
#gl = ax1.gridlines(draw_labels={'right':'y'}, dms=True,linewidth=.5,x_inline=False)
gl = ax1.gridlines(draw_labels=False, dms=True,linewidth=.5,x_inline=False)
ax1.plot(rtX[excMask],rtY[excMask],'.',markersize='1.5',color='grey')
sct = ax1.scatter(rtX[~excMask][po],rtY[~excMask][po],c=np.log10(popAvg[~excMask][po]),
                  s=5,alpha=.5,vmin=-4,vmax=-1,cmap='gnuplot2')
for bb in bbd.keys():
    bbv = bbd[bb]
    blanchette(bb,bbv,ax1)#,lab=True)
    
ax1.plot([-45.13,-43.1],[-8.3,-12.4],'r--')
ax1.plot([-37.5,-41.],[-20.6,-26],'r--')
ax1.plot([-29.7,-34.5],[-30.2,-32.7],'r--')
ax1.plot([-48.8,-51.4],[45.3,40.2],'r--')
# ax1.text(x=-56.2,y=42.5,s='Aleutian Islands\n(Spalding et al.)',ha='right',va='center',c='r',fontsize='small')
# ax1.text(x=-43.7,y=36.2,s='Gulf of Alaska\n(Spalding et al.)',ha='right',va='center',c='r',fontsize='small')
# ax1.text(x=-41.4,y=-7.7,s='Cortezian\n(Spalding et al.)',ha='left',va='center',c='r',fontsize='small')
# ax1.text(x=-43,y=-16.4,s='Mexican Tropical Pacific\n(Spalding et al.)',ha='right',va='center',c='r',fontsize='small')
# ax1.text(x=-34.6,y=-28.5,s='Chiapas-Nicaragua\n(Spalding et al.)',ha='right',va='center',c='r',fontsize='small')


for ix in [1]:
    avgMos=moArr[ix]
    popAvg = np.zeros(len(pLon))
    for mo in tqdm(avgMos):
        pop = np.array(popDict[(pld, mo, 100, 1, 20)])
        popAvg+=pop

    popAvg = np.array(popAvg)/len(avgMos)

po = np.argsort(popAvg[~excMask])

ax2 = fig.add_subplot(1, 2, 2, projection=rotated_crs,)
ax2.set_extent([minLon+117.5,maxLon+55,minLat-40,maxLat-17.5]) #manually figured out borders...
ax2.add_feature(cfeature.COASTLINE,linewidth=.5,alpha=.5)
ax2.add_feature(cfeature.BORDERS, linestyle='--',linewidth=.5)
ax2.add_feature(cfeature.LAND,alpha=.5)
ax2.add_feature(cfeature.OCEAN,alpha=.5)
ax2.set_title('Non-upwelling\n(Oct.-Feb.)')

gl = ax2.gridlines(draw_labels=False, dms=True,linewidth=.5,x_inline=False)
ax2.plot(rtX[excMask],rtY[excMask],'.',markersize='1.5',color='grey')
sct = ax2.scatter(rtX[~excMask][po],rtY[~excMask][po],c=np.log10(popAvg[~excMask][po]),
                  s=5,alpha=.5,vmin=-4,vmax=-1,cmap='gnuplot2')
for bb in bbd.keys():
    bbv = bbd[bb]
    blanchette(bb,bbv,ax2)
ax2.plot([-45.13,-43.1],[-8.3,-12.4],'r--')
ax2.plot([-37.5,-41.],[-20.6,-26],'r--')
ax2.plot([-29.7,-34.5],[-30.2,-32.7],'r--')
ax2.plot([-48.8,-51.4],[45.3,40.2],'r--')

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
divider = make_axes_locatable(ax2)

# Append a new axes for the colorbar, positioned to the right
cax = divider.append_axes("right",size='10%',pad=0.25,axes_class=p.Axes)

# Add the colorbar to the new axes
cbar = p.colorbar(sct,cax=cax,label='Log proportion of 1m species').solids.set(alpha=1)
#fig.subplots_adjust(hspace=0)
p.tight_layout()
p.show()
#p.savefig('blanchette_Pops.png')
