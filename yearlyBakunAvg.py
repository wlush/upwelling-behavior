import numpy as np
import pylab as p
import pandas as pd
import xarray as xr

dataDir = '/home/willlush/workfiles/behavior_upwelling/population_model/'

dat = pd.read_csv(dataDir+'bakunMonthly.txt',header=[1],delim_whitespace=True)
dat = dat[dat.YEAR>2006]
dat = dat[dat.YEAR<2024]

moAvg = dat.groupby(['LAT','LON']).mean()
moAvg.drop(columns='YEAR',inplace=True)

moAvg.to_csv(dataDir+'bakunAvg_2007_2023.csv')
