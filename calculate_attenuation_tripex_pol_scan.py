'''
this skript calculates the attenuation due to atmospheric gases for selected cases of the tripex-pol-scan campaign using the functions implemented in McRadar.
It also matches the reflectivities at the tip of the cloud to calculate attenuation due to liquid and snow using the method described in von Terzi et al 2022
Author: Leonie von Terzi
'''

from glob import glob
import numpy as np
import pandas as pd
from sys import argv
import xarray as xr
import matplotlib.pyplot as plt

from cmcrameri import cm as colmap
from attenuationFunctions import * 
from mcradar import getHydroAtmAtt
import atmFunc

# define all necessary paths:
pathData = '/archive/meteo/external-obs/juelich/' #base path of data
pathCN = pathData+'cloudnet/jue/' # CN path
pathRes = pathData + 'campaign_aux_data/tripex-pol-scan/wband_scan/resampled/' # this contains resampled moments of w-band

date2start = '20220206'
date2end = '20220206'
dateStart = pd.to_datetime(date2start); dateEnd = pd.to_datetime(date2end)
dateList = pd.date_range(dateStart, dateEnd,freq='d')
freq=[9.6e9,35.5e9,94.0e9]


for date in dateList:
	'''
	# calculate attenuation from gases using mcradar and cloudnet product
	CN = xr.open_dataset('{path}{date}_juelich_categorize.nc'.format(path=pathCN,date=date.strftime('%Y%m%d')))
	#results = getAtmAttenuation
	CN['relHum'] = atmFunc.speHumi2RelHum(CN.q, CN.temperature, CN.pressure)
	CN['model_height'] = CN['model_height'] -111
	CN['height'] = CN['height'] -111
	CN1 = CN[['temperature','pressure','relHum']].rename({'model_time':'time','model_height':'range'})
	CN1 = CN1.where(CN1.range < 12000,drop=True)
	CN2 = CN[['lwp','category_bits']].rename({'height':'range'})
	heightRes = np.diff(CN1.range)
	#quit()
	CN1['att_atmo'] = CN1.relHum.copy().expand_dims({'freq':freq})*np.nan
	#quit()
	for t in CN1.time:
		print(t.values)
		for i,h in enumerate(CN1.range[0:-1]):
			#print(h)
			for f in freq:
				rt_kextatmo = getHydroAtmAtt(CN1.temperature.sel(range=h,time=t),CN1.relHum.sel(range=h,time=t),CN1.pressure.sel(range=h,time=t),f*1e-9)
				CN1['att_atmo'].loc[f,t,h] = 10*np.log10(np.exp(rt_kextatmo.values * heightRes[i]))
		
	CN1['atm_att_2way'] = 2*CN1['att_atmo'].cumsum(dim='range')
	CN1.atm_att_2way.attrs['long_name'] = 'path integrated attenuation due to atmospheric gases, 2 way'
	CN1.atm_att_2way.attrs['units'] = 'dBz'
	CN1.to_netcdf('test_files/{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	CN2.to_netcdf('test_files/{date}_cloudnet_attenuation2.nc'.format(date=date.strftime('%Y%m%d')))		
	'''
	CN1 = xr.open_dataset('test_files/{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	CN2 = xr.open_dataset('test_files/{date}_cloudnet_attenuation2.nc'.format(date=date.strftime('%Y%m%d')))		

	#- read in radar Data and reindex to common grid, also, for easier handling: rename all variables to the same name and have freq as coordinate
	ZENW = getWband(pathRes,date,offset=0) 
	ZENKa = getJoyrad(pathData,date,ZENW.range,ZENW.time,'Ka',offset=0) 
	ZENX = getJoyrad(pathData,date,ZENW.range,ZENW.time,'X',offset=0.5) 
	
	#- merge triple-frequency dataset and calculated attenuation as well as other CN variables needed for quality flag generation:
	#CN1 = CN[['atm_att_2way','temperature']].rename({'model_time':'time','model_height':'range'})
	CN2 = CN2.interp({'time':ZENKa.time,'range':ZENKa.range})
	CN1 = CN1.interp({'time':ZENKa.time,'range':ZENKa.range})
	print(CN1,CN2)
	data = xr.merge([ZENKa,ZENX,ZENW,CN1,CN2])
	print(data)
	#quit()
	#data['calibration_offset'] = data.calibration_offset.fillna(0)
	#- add gas attenuation to Ze
	data['Ze'] = data.Ze + data.atm_att_2way + data.calibration_offset
	data.Ze.attrs['units'] = 'dBz'
	print('finished creating combined dataset')
	#- I can potentially masked areas where variance is larger than 2, however that does not make a difference.
	#data = getVarianceFlag(data) # variance flag
	#dataMaskedX = data.Ze.sel(freq=9.6e9).where(data.varRolling.sel(freq=9.6e9) < 2)
	#dataMaskedW = data.Ze.sel(freq=94e9).where(data.varRolling.sel(freq=94e9) < 2)
	#dataMasked = xr.merge([dataMaskedX.expand_dims({'freq':np.array([9.6e9])}),dataMaskedW.expand_dims({'freq':np.array([94e9])}),data.Ze.sel(freq=35.5e9).expand_dims({'freq':np.array([35.5e9])})])
	#print(dataMasked)
	
	#- select only cases where 1km above ML, and within predefined interval 30 to -10dB for Ka-W (Ka-Band reflectivity), -15 to 0 for X-Ka (Ka-Band reflectivity)
	dataXW = getInterval(data,interval=[-35,-10])
	dataXKa = getInterval(data,interval=[-35,-5])
	#quit()
	
	#- calculate offset
	dataOffset = getOffset(dataXKa,dataXW,date)
	#- now reindex offset to dataset 
	dataOffset = dataOffset.reindex({'time':data.time},method='nearest',tolerance='5T')
	data = xr.merge([data,dataOffset])
	data['Ze'] = data.Ze + data.offsetMean
	data.Ze.attrs['comments'] = 'corrected for absolute calibration error, gas attenuation, and attenuation due to hydrometeors. To remove gas and hydrometeor attenuation and calibration, simply subtract atm_att_2way, calibration_offset and offsetMean from Ze'
	print('calculated and applied offset')
	print(data)
	#- now make quality flag to be the same as Joses
	data['quality_flag_offset'] = getFlags(data)
	data = data[['Ze','MDV','WIDTH','SK','atm_att_2way','calibration_offset','temperature','offsetMean','quality_flag_offset','offsetStd','refRadarStd','refRadarVar','validPoints','totalPoints','correlation','offsetVar','lwp']]
	print('generated quality flags')
	data = globalAttr(data)
	encDic = {x: {"zlib": True} for x in data}
	data.to_netcdf('../output/{year}/{month}/{day}/{date}_tripex_pol_scan_3fr_lv2.nc'.format(year=date.strftime('%Y'),
																							month = date.strftime('%m'),
																							day = date.strftime('%d'),
																							date = date.strftime('%Y%m%d')))
	#quit()
	#- plot flags:
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
	ax[1].plot(data.time,data.offsetMean.sel(freq=94e9),c='C0',lw=2,label='W-Band')
	ax[1].plot(data.time,data.offsetMean.sel(freq=35.5e9),c='C1',lw=2,label='Ka-Band')
	#ax[1].plot(data.time,data.offsetMean.sel(freq=94e9).where(~data.varFlag.sel(freq=94e9)),c='C0',lw=2,label='W-Band')
	#ax[1].plot(data.time,data.offsetMean.sel(freq=9.6e9).where(~data.varFlag.sel(freq=9.6e9)),c='C1',lw=2,label='X-Band')
	#ax[1].axhline(y=2,c='r',ls=':',lw=2)
	ax[1].set_ylabel('Calculated offset',fontsize=20)
	
	ax[2].plot(data.time,data.offsetVar.sel(freq=94e9),c='C0',lw=2,label='W-Band')
	ax[2].plot(data.time,data.offsetVar.sel(freq=35.5e9),c='C1',lw=2,label='Ka-Band')
	ax[2].axhline(y=2,c='r',ls=':',lw=2)
	ax[2].set_ylabel('variance of offset',fontsize=20)
	
	ax[3].plot(data.time,data.correlation.sel(freq=94e9),lw=2,label='W-Band')
	ax[3].plot(data.time,data.correlation.sel(freq=35.5e9),lw=2,label='Ka-Band')
	ax[3].axhline(y=0.7,c='r',ls=':',lw=2)
	ax[3].set_ylabel('correlation \n between Zes',fontsize=20)
	
	ax[4].plot(data.time,data.validPoints.sel(freq=94e9),lw=2,label='W-Band')
	ax[4].plot(data.time,data.validPoints.sel(freq=35.5e9),lw=2,label='Ka-Band')
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
	plt.savefig('plots/{date}_Flags_offset_30minMean_maskedVar.png'.format(date=date.strftime('%Y%m%d')))
	plt.close()
	quit()
	
	
	
	
#DWRKaW30 = dataKaW.Ze.sel(freq=35.5e9) - dataKaW.Ze.sel(freq=94e9)+dataOffset30.offsetMean.sel(freq=94e9)
#DWRKaW30 = DWRKaW30.dropna(dim='time',how='all')
#DWRXKa30 = dataXKa.Ze.sel(freq=9.6e9)-dataOffset30.offsetMean.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
#DWRXKa30 = DWRXKa30.dropna(dim='time',how='all')

#DWRKaW10 = dataKaW.Ze.sel(freq=35.5e9) - dataKaW.Ze.sel(freq=94e9)+dataOffset10.offsetMean.sel(freq=94e9)
#DWRKaW10 = DWRKaW10.dropna(dim='time',how='all')
#DWRXKa10 = dataXKa.Ze.sel(freq=9.6e9)-dataOffset10.offsetMean.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
#DWRXKa10 = DWRXKa10.dropna(dim='time',how='all')

#DWRKaW45 = dataKaW.Ze.sel(freq=35.5e9) - dataKaW.Ze.sel(freq=94e9)+dataOffset45.offsetMean.sel(freq=94e9)
#DWRKaW45 = DWRKaW45.dropna(dim='time',how='all')
#DWRXKa45 = dataXKa.Ze.sel(freq=9.6e9)-dataOffset45.offsetMean.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
#DWRXKa45 = DWRXKa45.dropna(dim='time',how='all')

#DWRKaW90 = dataKaW.Ze.sel(freq=35.5e9) - dataKaW.Ze.sel(freq=94e9)+dataOffset90.offsetMean.sel(freq=94e9)
#DWRKaW90 = DWRKaW90.dropna(dim='time',how='all')
#DWRXKa90 = dataXKa.Ze.sel(freq=9.6e9)-dataOffset90.offsetMean.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
#DWRXKa90 = DWRXKa90.dropna(dim='time',how='all')


#bins=np.arange(-5,5.5,0.5)+0.25
#print(bins)
#plt.hist([DWRKaW10.values.flatten(),DWRKaW.values.flatten(),DWRKaW30.values.flatten(),DWRKaW45.values.flatten(),DWRKaW90.values.flatten()],bins=bins,label=['10min','15min','30min','45min','90min'])
#plt.hist(DWRKaW.values.flatten(),bins=bins,alpha=0.5,label='KaW 15min')
#plt.hist(DWRXKa.values.flatten(),bins=bins,alpha=0.5,label='XKa 15min')
#plt.hist(DWRKaW30.values.flatten(),bins=bins,alpha=0.5,label='KaW 30min')
#plt.hist(DWRXKa30.values.flatten(),bins=bins,alpha=0.5,label='XKa 30min')
#plt.legend()
#plt.xlabel('DWR KaW [dB]')
#plt.title('30 min time window')
#plt.savefig('test_attenuation_DWRKaW_hist.png')
#plt.show()
#quit()
#plt.pcolormesh(DWRKaW.time,DWRKaW.range,DWRKaW.T,cmap='jet',vmin=-5,vmax=5)
#plt.colorbar()
#plt.show()
#plt.scatter(dataKaW.Ze.sel(freq=35.5e9),dataKaW.Ze.sel(freq=(94e9))+dataKaW.offsetMean.sel(freq=94e9))
#plt.xlim(-40,0)
#plt.ylim(-40,0)
#plt.show()
#print(dataKaW)
#quit()

#print(dataOffset.correlation.sel(freq=94e9).max())
#print(dataOffset.correlation.sel(freq=94e9).min())
#quit()
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9),c='C0',alpha=0.5)
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9)+dataOffset.offsetStd.sel(freq=94e9),c='C0',ls='--',alpha=0.5)
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9)-dataOffset.offsetStd.sel(freq=94e9),c='C0',ls='--',alpha=0.5)
#dataOffset = dataOffset.where(dataOffset.correlation > 0.7)
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9),c='C0')
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9)+dataOffset.offsetStd.sel(freq=94e9),c='C0',ls='--')
#plt.plot(dataOffset.time,dataOffset.offsetMean.sel(freq=94e9)-dataOffset.offsetStd.sel(freq=94e9),c='C0',ls='--')
#plt.grid()
#plt.ylabel('Offset W-Band relative to Ka-Band [dB]')
#plt.title('Moving window width 15min')
#plt.tight_layout()
#plt.savefig('20220206_test_offset_mowin_15min.png')
#plt.show()
#quit()
#- TODO: vary the width of moving window

	
	
	
