'''
this skript compares different ways to calculate attenuation, once with areas with large variance masked, once without
Author: Leonie von Terzi
'''

from glob import glob
import numpy as np
import pandas as pd
from sys import argv
import atmFunc
import xarray as xr
#import attenuationLib as attLib
import matplotlib.pyplot as plt
from matplotlib import cm
from mcradar import getHydroAtmAtt
from cmcrameri import cm as colmap
from attenuationFunctions import * 
from matplotlib.colors import ListedColormap
import matplotlib.colors as colors
date2start = '20220206'
date = pd.to_datetime(date2start)
data = xr.open_dataset('../output/{year}/{month}/{day}/{date}_tripex_pol_scan_3fr_lv2_NOTmaskedVar.nc'.format(year=date.strftime('%Y'),
																							month = date.strftime('%m'),
																							day = date.strftime('%d'),
																							date = date.strftime('%Y%m%d')))
																							
dataMasked = xr.open_dataset('../output/{year}/{month}/{day}/{date}_tripex_pol_scan_3fr_lv2_maskedVar.nc'.format(year=date.strftime('%Y'),
																							month = date.strftime('%m'),
																							day = date.strftime('%d'),
																							date = date.strftime('%Y%m%d')))
																							
fig,ax = plt.subplots(nrows=5,figsize=(15,15),constrained_layout=True,sharex=True)
p1 = ax[0].pcolormesh(data.time,data.range,data.Ze.sel(freq=9.6e9).T,cmap = colmap.batlow,vmin=-30,vmax=20)
#p1 = ax[0].pcolormesh(data.time,data.range,data.Ze.sel(freq=9.6e9).where(~data.varFlag.sel(freq=9.6e9)).T,cmap = colmap.batlow,vmin=-30,vmax=20)
cb = plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
cb.set_label('Ze X-Band [dBz]',fontsize=20)
cb.ax.tick_params(labelsize=18)
ax[0].set_title('possible flags to apply, 30 min mean',fontsize=24)
ax[0].set_ylabel('height [m]',fontsize=20)

#ax[1].plot(data.time,data.offsetVar.sel(freq=94e9),c='C0',lw=2,label='W-Band')
#ax[1].plot(data.time,data.offsetVar.sel(freq=9.6e9),c='C1',lw=2,label='X-Band')
ax[1].plot(data.time,data.offsetMean.sel(freq=94e9),c='C0',lw=2,ls=':',label='W-Band')
ax[1].plot(data.time,data.offsetMean.sel(freq=9.6e9),c='C1',lw=2,ls=':',label='X-Band')
ax[1].plot(dataMasked.time,dataMasked.offsetMean.sel(freq=94e9),c='C0',lw=2,label='W-Band var Masked')
ax[1].plot(dataMasked.time,dataMasked.offsetMean.sel(freq=9.6e9),c='C1',lw=2,label='X-Band var Masked')
#ax[1].plot(data.time,data.offsetMean.sel(freq=94e9).where(~data.varFlag.sel(freq=94e9)),c='C0',lw=2,label='W-Band')
#ax[1].plot(data.time,data.offsetMean.sel(freq=9.6e9).where(~data.varFlag.sel(freq=9.6e9)),c='C1',lw=2,label='X-Band')
#ax[1].axhline(y=2,c='r',ls=':',lw=2)
ax[1].set_ylabel('Calculated offset',fontsize=20)

ax[2].plot(data.time,data.offsetVar.sel(freq=94e9),c='C0',lw=2,ls=':',label='W-Band')
ax[2].plot(data.time,data.offsetVar.sel(freq=9.6e9),c='C1',lw=2,ls=':',label='X-Band')
ax[2].plot(dataMasked.time,dataMasked.offsetVar.sel(freq=94e9),c='C0',lw=2,label='W-Band var Masked')
ax[2].plot(dataMasked.time,dataMasked.offsetVar.sel(freq=9.6e9),c='C1',lw=2,label='X-Band var Masked')
ax[2].axhline(y=2,c='r',ls=':',lw=2)
ax[2].set_ylabel('variance of offset',fontsize=20)

ax[3].plot(data.time,data.correlation.sel(freq=94e9),lw=2,c='C0',ls=':',label='W-Band')
ax[3].plot(data.time,data.correlation.sel(freq=9.6e9),lw=2,c='C1',ls=':',label='X-Band')
ax[3].plot(dataMasked.time,dataMasked.correlation.sel(freq=94e9),lw=2,c='C0',label='W-Band var Masked')
ax[3].plot(dataMasked.time,dataMasked.correlation.sel(freq=9.6e9),lw=2,c='C1',label='X-Band var Masked')
ax[3].axhline(y=0.7,c='r',ls=':',lw=2)
ax[3].set_ylabel('correlation \n between Zes',fontsize=20)

ax[4].plot(data.time,data.validPoints.sel(freq=94e9),lw=2,ls=':',c='C0',label='W-Band')
ax[4].plot(data.time,data.validPoints.sel(freq=9.6e9),lw=2,ls=':',c='C1',label='X-Band')
ax[4].plot(dataMasked.time,dataMasked.validPoints.sel(freq=94e9),lw=2,c='C0',label='W-Band var Masked')
ax[4].plot(dataMasked.time,dataMasked.validPoints.sel(freq=9.6e9),lw=2,c='C1',label='X-Band var Masked')
ax[4].axhline(y=300,c='r',ls=':',lw=2)
ax[4].set_ylabel('Number of \n valid points',fontsize=20)
for i,a in enumerate(ax):
	if i > 0:
		a.legend(fontsize=16)
	a.tick_params(labelsize=18)
plt.setp(ax[3].xaxis.get_majorticklabels(), rotation=0,horizontalalignment='left')
#plt.colorbar(p3,ax=ax[0],label='Ka -30 to -10 dB for W-Band')
#ax[1].set_title('Ka-Band between -30 and -10dB for W-Band calibration')
#plt.tight_layout()
plt.savefig('plots/{date}_offset_30minMean_compareMaskvsnoMask.png'.format(date=date.strftime('%Y%m%d')))
plt.close()
quit()
