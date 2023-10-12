'''
this skript calculates the attenuation due to atmospheric gases for selected cases of the tripex-pol-scan campaign. Using the functions implemented in McRadar
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

def getclosestTemp(data,T):
	'''
	calculate range closest to a certain Temperature
	input:
	data: xr.dataset with temperature as variable and range as coordinate
	T: temperature to be selected in Kelvin
	returns: 
	closest range to the input temperature
	'''
	closestT = np.abs(T - data.temperature)
	#allna = closestT.isnull().all(dim='Vel')
	return data.range[closestT.argmin(dim='range')]

def maskData(variable,flag):
# apply masks from LV2 dataset:
    maskP = int('1000000000000000',2) #n points
    maskC = int('0100000000000000',2) # correl
    maskV = int('0010000000000000',2) # Variance
    print(maskP)
    print((flag.values & maskV)!= maskV)
    print(flag.values.dtype)
    quit()
    #plt.plot(((flag.values & maskP)!= maskP)*1)
    #plt.show()
    #plt.plot((flag.values >> 15) & 1)
    #plt.show()
    print(~((flag.values >> 13) & 1).astype(bool))
    #quit()
    #plt.close()
    #plt.pcolormesh(variable.time,variable.range,variable.T)
    #plt.show()
    variableMasked = variable.where((flag.values & maskV) != maskV)
    #plt.pcolormesh(variableMasked.time,variableMasked.range,variableMasked.T)
    #plt.show()
    #plt.close()
    #fig,ax=plt.subplots(ncols=2)
    ax[0].pcolormesh(variableMasked.time,variableMasked.range,variableMasked.T)
    variableMasked = variable.where(~((flag.values >> 13) & 1).astype(bool))
    ax[1].pcolormesh(variableMasked.time,variableMasked.range,variableMasked.T)
    
    plt.show()
    variable = variable.where((flag.values & maskC) != maskC)
    variable =variable.where((flag.values & maskV) != maskV)
    
    return variable
pathCN = ''#'/archive/meteo/external-obs/juelich/cloudnet/jue/'
pathPro = '/archive/meteo/external-obs/juelich/tripex-pol/processed/'
pathData = pathPro+'tripex_pol_level_1/' #'/scratch/l/L.Terzi/campaign_aux_data/'
#pathPro = pathData+'campaign_aux_data/tripex-pol-scan/wband_scan/processed/' # this contains the regridded and dealized spectra
#pathRes = pathData + 'campaign_aux_data/tripex-pol-scan/wband_scan/resampled/' # this contains resampled moments of w-band
LV1name = '_tripex_pol_3fr_L1_mom.nc'
date2start = '20181124'
date2end = '20181124'
dateStart = pd.to_datetime(date2start); dateEnd = pd.to_datetime(date2end)
dateList = pd.date_range(dateStart, dateEnd,freq='d')
freq=[9.6e9,35.5e9,94.0e9]
for date in dateList:
	timeWindow = pd.to_timedelta(60, unit='m')
	start = date; end=date+pd.offsets.Day(1)-pd.offsets.Second(1)
	beginTime = pd.date_range(start-timeWindow, end-timeWindow, freq='1min')
	midTime = pd.date_range(start, end, freq='1min')
	endTime = pd.date_range(start+timeWindow, end+timeWindow, freq='1min')
	#print(timesBegin,timesEnd,timesMid)
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(30)-pd.offsets.Second(4),freq='60S')
	#midTime = pd.date_range(date + pd.offsets.Minute(15),date+pd.offsets.Day(1)-pd.offsets.Minute(15),freq='60S')
	#endTime = pd.date_range(date + pd.offsets.Minute(30),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='60S')
	#print(beginTime,endTime)
	#quit()
	dataLV2 = xr.open_dataset('{path}/tripex_pol_level_2/{date}_tripex_pol_3fr_L2_mom.nc'.format(path=pathPro,date=date.strftime('%Y%m%d')))
	dbz_w = dataLV2.W_DBZ_H
	print(dataLV2.quality_flag_offset_w)
	# calculate attenuation from gases using mcradar
	'''
	timeRef = pd.date_range(date,date+pd.offsets.Day(1),freq='30T')
	print(timeRef)
	CN = xr.open_dataset('{path}{date}_juelich_categorize.nc'.format(path=pathCN,date=date.strftime('%Y%m%d')))
	#results = getAtmAttenuation
	CN['relHum'] = atmFunc.speHumi2RelHum(CN.specific_humidity, CN.temperature, CN.pressure)
	CN['model_height'] = CN['model_height'] -105
	CN = CN.where(CN.model_height < 12000,drop=True)
	heightRes = np.diff(CN.model_height)
	print(CN.time.values)
	
	CN = CN.reindex({'time':timeRef},method='nearest',tolerance='100S')
	CN['att_atmo'] = CN.relHum.copy().expand_dims({'freq':freq})*np.nan
	print(CN)
	#quit()
	for t in CN.time:
		print(t.values)
		for i,h in enumerate(CN.model_height[0:-1]):
			#print(h)
			for f in freq:
				rt_kextatmo = getHydroAtmAtt(CN.temperature.sel(model_height=h,time=t),CN.relHum.sel(model_height=h,time=t),CN.pressure.sel(model_height=h,time=t),f*1e-9)
				CN['att_atmo'].loc[f,t,h] = 10*np.log10(np.exp(rt_kextatmo.values * heightRes[i]))
		
		#plt.plot(CN.att_atmo.sel(model_time=t,freq=94).values,CN.model_height)
		#plt.plot(CN.att_atmo.sel(model_time=t,freq=35.5).values,CN.model_height)
		#plt.plot(CN.att_atmo.sel(model_time=t,freq=9.6).values,CN.model_height)
		
		#plt.show()
	CN['atm_att_2way'] = 2*CN['att_atmo'].cumsum(dim='model_height')
	print(CN.atm_att_2way)
	CN.to_netcdf('{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	#quit()
	'''
	'''
	CN = xr.open_dataset('test_files/{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	#- merge triple-frequency dataset and calculated attenuation:
	dataZEN = xr.open_dataset('{path}/{date}{ID}'.format(path=pathData,date=date.strftime('%Y%m%d'),ID = LV1name))
	dataZEN['Ka_DBZ_H'] = dataZEN.Ka_DBZ_H + dataLV2.rain_offset_Ka
	dataZEN['W_DBZ_H'] = dataZEN.W_DBZ_H + dataLV2.rain_offset_W
	dataZEN['X_DBZ_H'] = dataZEN.X_DBZ_H + dataLV2.rain_offset_X
	Ka = dataZEN[['Ka_DBZ_H','Ka_VEL_H']]
	Ka = Ka.expand_dims({'freq':np.array([35.5e9])}).rename({'Ka_DBZ_H':'Ze','Ka_VEL_H':'MDV'})
	X = dataZEN[['X_DBZ_H','X_VEL_H']]
	X = X.expand_dims({'freq':np.array([9.6e9])}).rename({'X_DBZ_H':'Ze','X_VEL_H':'MDV'})
	W = dataZEN[['W_DBZ_H','W_VEL_H']]
	W = W.expand_dims({'freq':np.array([94e9])}).rename({'W_DBZ_H':'Ze','W_VEL_H':'MDV'})
	#dataZEN = xr.merge([X,Ka,W])
	#print(dataZEN)
	
	#quit()
	CN = CN[['atm_att_2way','temperature']].rename({'model_height':'range'})
	#CN['freq'] = CN.freq*1e9
	print(CN)
	CN = CN.interp({'time':dataZEN.time,'range':dataZEN.range})
	data = xr.merge([X,Ka,W,CN])
	print(data)
	
	data['Ze'] = data.Ze + data.atm_att_2way
	data.Ze.attrs['units'] = 'dBz'
	data.Ze.attrs['comments'] = 'corrected for gas attenuation, to remove gas attenuation, simply subtract atm_att_2way from Ze'
	
	data.atm_att_2way.attrs['long_name'] = 'path integrated attenuation due to atmospheric gases, 2 way'
	data.atm_att_2way.attrs['units'] = 'dBz'
	data.to_netcdf('{date}_test_att_correction.nc'.format(date=date.strftime('%Y%m%d')))
	
	data = xr.open_dataset('{date}_test_att_correction.nc'.format(date=date.strftime('%Y%m%d')))
	data = data.where(np.isfinite(data))
	#- select only cases where T<-5°C (one km above ML)
	#- get range of 1km above ML:
	#plt.pcolormesh(data.time,data.range,data.Ze.sel(freq=94e9).T,cmap='Greens')
	#plt.colorbar()
	#plt.show()
	#quit()
	MLrange = getclosestTemp(data,273.15)
	dataCold = data.where(data.range > MLrange+1000)
	#dataCold = data.where(data.range > 5000)
	#- select only cases where Ze within predifend interval: -30 to -10dB for Ka-W (Ka-Band reflectivity), -15 to 0 for X-Ka (Ka-Band reflectivity)
	dataXKa = dataCold.where(dataCold.Ze.sel(freq=35.5e9)>-15)
	dataXKa = dataXKa.where(dataXKa.Ze.sel(freq=35.5e9)<-0)
	
	dataKaW = dataCold.where(dataCold.Ze.sel(freq=35.5e9)>-30)
	dataKaW = dataKaW.where(dataKaW.Ze.sel(freq=35.5e9)<-10)
	print(dataKaW)
	#plt.pcolormesh(dataKaW.time,dataKaW.range,dataKaW.Ze.sel(freq=94e9).T,cmap='Greys')
	#plt.colorbar()
	#offsetW = dataKaW.Ze.sel(freq=94e9) - dataKaW.Ze.sel(freq=35.5e9)
	#offsetX = dataXKa.Ze.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
	
	#- calculate mean of offset (and other statistics) in moving window of 15min width 
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(15)-pd.offsets.Second(4),freq='60S')
	#midTime = pd.date_range(date + pd.offsets.Minute(7)+pd.offsets.Second(30),date+pd.offsets.Day(1)-pd.offsets.Minute(7)-pd.offsets.Second(34),freq='60S')
	#endTime = pd.date_range(date + pd.offsets.Minute(15),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='60S')
	
	#beginTime = pd.date_range(date,date+pd.offsets.Day(1)-pd.offsets.Minute(30)-pd.offsets.Second(4),freq='60S')
	#midTime = pd.date_range(date + pd.offsets.Minute(15),date+pd.offsets.Day(1)-pd.offsets.Minute(15),freq='60S')
	#endTime = pd.date_range(date + pd.offsets.Minute(30),date+pd.offsets.Day(1)-pd.offsets.Second(4),freq='60S')
	
	print(beginTime,midTime,endTime)
	offsetVar = xr.DataArray(dims=('time','freq'),coords={'time': midTime,'freq':[9.6e9,94e9]})
	dataOffset = xr.Dataset({'offsetMean':offsetVar,
							'offsetStd':offsetVar.copy(),
							'offsetVar':offsetVar.copy(),
							'refRadarStd':offsetVar.copy(),
							'refRadarVar':offsetVar.copy(),
							'validPoints':offsetVar.copy(),
							'totalPoints':offsetVar.copy(),
							'correlation':offsetVar.copy(),
							})
	
	dataOffset['offsetMean'].attrs = {'description':'mean of offset calculated between X and Ka-Band (select freq=9.6e9) and W and Ka-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['offsetStd'].attrs = {'description':'standard deviation of offset calculated between X and Ka-Band (select freq=9.6e9) and W and Ka-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['offsetVar'].attrs = {'description':'variance of offset calculated between X and Ka-Band (select freq=9.6e9) and W and Ka-Band (select freq=94e9) in 30 min window','units':'dB'}
	dataOffset['refRadarStd'].attrs = {'description':'standard deviation of Ka-Band Ze (reference radar) in 15 min window','units':'dB'}
	dataOffset['refRadarVar'].attrs = {'description':'variance of Ka-Band Ze (reference radar) in 15 min window','units':'dB'}
	dataOffset['validPoints'].attrs = {'description':'number of valid points in 15 min window','units':''}
	dataOffset['totalPoints'].attrs = {'description':'number of total points in 15 min window','units':''}
	dataOffset['correlation'].attrs = {'description':'correlation between X and Ka-Band Ze (select freq=9.6e9) and between W and Ka-band (select freq=94e9) 30 min window','units':''}
	
	
	for i,tstart,tmid,tend in zip(range(len(beginTime)),beginTime,midTime,endTime):
		#tmid = tstart + pd.offsets.Minute(7)+ pd.offsets.Second(30)
		#tend = tstart + pd.offsets.Minute(15)
		#offsetWsel = offsetW.sel(time=slice(tstart,tend))#.flatten()
		ZeKaWsel = dataKaW.Ze.sel(time=slice(tstart,tend))
		ZeXKasel = dataXKa.Ze.sel(time=slice(tstart,tend))
		offsetW = ZeKaWsel.sel(freq=94e9) - ZeKaWsel.sel(freq=35.5e9)
		#print(offsetW.mean())
		#plt.pcolormesh(offsetW.time,offsetW.range,offsetW.T,cmap='Greens')
		#plt.colorbar()
		#plt.show()
		#quit()
		offsetX = ZeXKasel.sel(freq=9.6e9) - ZeXKasel.sel(freq=35.5e9)
		#offsetXsel = offsetX.sel(time=slice(tstart,tend)).flatten()
		#- calculate statistics (mean, std, count of valid points and total count of points in volume)
		#if not np.isnan(offsetW).all():
		#	plt.hist(offsetW.values.flatten(),bins=50)
		#	plt.axvline(x=offsetW.mean(),c='C1')
		#	plt.axvline(x=offsetW.mean()+offsetW.std(),c='C1',ls='--')
		#	plt.axvline(x=offsetW.mean()-offsetW.std(),ls='--',c='C1')
		#	plt.axvline(x=dataOffset['offsetMean'].loc[tmid,94e9],c='C2')
		#	print(dataOffset['offsetMean'].loc[tmid,94e9])
		#	print(offsetW.mean())
			#plt.axvline(x=offsetW.mean()+offsetW.std(),c='C1',ls='--')
			
			#plt.axvline(x=offsetW.mean()-offsetW.std(),ls='--',c='C1')
			
		#	plt.show()
		dataOffset['offsetMean'].loc[tmid,94e9] = offsetW.mean()
		dataOffset['offsetStd'].loc[tmid,94e9] = offsetW.std()
		dataOffset['offsetVar'].loc[tmid,94e9] = offsetW.var()
		dataOffset['validPoints'].loc[tmid,94e9] = offsetW.count()
		dataOffset['totalPoints'].loc[tmid,94e9] = len(offsetW.values.flatten())
		dataOffset['correlation'].loc[tmid,94e9] = xr.corr(dataKaW.Ze.sel(freq=94e9),dataKaW.Ze.sel(freq=35.5e9))
		dataOffset['refRadarStd'].loc[tmid,94e9] = ZeKaWsel.sel(freq=35.5e9).std()
		dataOffset['offsetMean'].loc[tmid,9.6e9] = offsetX.mean()
		dataOffset['offsetStd'].loc[tmid,9.6e9] = offsetX.std()
		dataOffset['offsetVar'].loc[tmid,9.6e9] = offsetX.var()
		dataOffset['validPoints'].loc[tmid,9.6e9] = offsetX.count()
		dataOffset['totalPoints'].loc[tmid,9.6e9] = len(offsetX.values.flatten())
		dataOffset['correlation'].loc[tmid,9.6e9] = xr.corr(dataXKa.Ze.sel(freq=9.6e9),dataXKa.Ze.sel(freq=35.5e9))
		dataOffset['refRadarStd'].loc[tmid,9.6e9] = ZeXKasel.sel(freq=35.5e9).std()
		print(tend)
		
	print(dataOffset)
	
	dataOffset.to_netcdf('{date}_test_offset_WindowJose.nc'.format(date=date.strftime('%Y%m%d')))
	'''
	data = xr.open_dataset('test_files/{date}_test_att_correction.nc'.format(date=date.strftime('%Y%m%d')))
	#data = xr.open_dataset('{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	#plt.pcolormesh(data.time,data.range,data.Ze.sel(freq=9.6e9).T)
	#plt.show()
	#quit()
	#print(data)
	#quit()
	dataOffset = xr.open_dataset('test_files/{date}_test_offset_WindowJose.nc'.format(date=date.strftime('%Y%m%d')))
	#dataOffset30 = xr.open_dataset('20220206_test_offset_30minWindow.nc')
	
	dataOffset = dataOffset.reindex({'time':data.time},method='nearest',tolerance='60S')
	#dataOffset30 = dataOffset30.reindex({'time':data.time},method='nearest',tolerance='60S')
	data = xr.merge([data,dataOffset])
	
	
	MLrange = getclosestTemp(data,273.15)
	dataCold = data.where(data.range > MLrange+1000)
	
	#- select only cases where Ze within predifend interval: -30 to -10dB for Ka-W (Ka-Band reflectivity), -15 to 0 for X-Ka (Ka-Band reflectivity)
	dataXKa = dataCold.where(dataCold.Ze.sel(freq=35.5e9)>-15)
	dataXKa = dataXKa.where(dataXKa.Ze.sel(freq=35.5e9)<-0)
	
	dataKaW = dataCold.where(dataCold.Ze.sel(freq=35.5e9)>-30)
	dataKaW = dataKaW.where(dataKaW.Ze.sel(freq=35.5e9)<-10)
	DWRKaW = dataKaW.Ze.sel(freq=35.5e9) - dataKaW.Ze.sel(freq=94e9)#+dataKaW.offsetMean.sel(freq=94e9)
	
	
	DWRKaW = DWRKaW.dropna(dim='time',how='all')
	DWRXKa = dataXKa.Ze.sel(freq=9.6e9)-dataXKa.offsetMean.sel(freq=9.6e9) - dataXKa.Ze.sel(freq=35.5e9)
	DWRXKa = DWRXKa.dropna(dim='time',how='all')

	bins=np.arange(-5,5.5,0.5)+0.25
	print(bins)
	#plt.hist([DWRKaW10.values.flatten(),DWRKaW.values.flatten(),DWRKaW30.values.flatten(),DWRKaW45.values.flatten(),DWRKaW90.values.flatten()],bins=bins,label=['10min','15min','30min','45min','90min'])
	#plt.hist(DWRKaW.values.flatten(),bins=bins,alpha=0.5,label='KaW 2h')
	#plt.hist(DWRXKa.values.flatten(),bins=bins,alpha=0.5,label='XKa 2h')
	#plt.hist(DWRKaW30.values.flatten(),bins=bins,alpha=0.5,label='KaW 30min')
	#plt.hist(DWRXKa30.values.flatten(),bins=bins,alpha=0.5,label='XKa 30min')
	#plt.legend()
	#plt.xlabel('DWR KaW [dB]')
	#plt.title('30 min time window')
	#plt.savefig('test_attenuation_tripex-pol.png')
	#plt.show()
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


	fig,ax = plt.subplots(nrows=4,figsize=(15,15),constrained_layout=True,sharex=True)
	#p1 = ax[0].pcolormesh(data.time,data.range,data.Ze.sel(freq=9.6e9).T,cmap = colmap.batlow,vmin=-30,vmax=20)
	#cb = plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
	#cb.set_label('Ze X-Band [dBz]',fontsize=20)
	#p1 = ax[0].pcolormesh(dataLV2.time,dataLV2.range,dbz_w.T,cmap = 'Greys',vmin=-30,vmax=20)
	qFlagW = dataLV2.quality_flag_offset_w
	dbz_w = maskData(dbz_w, qFlagW)
	quit()
	p1 = ax[0].pcolormesh(dataLV2.time,dataLV2.range,dbz_w.T,cmap = colmap.batlow,vmin=-30,vmax=20)
	cb = plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
	cb.set_label('Ze W-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[0].set_title('possible flags to apply, 30 min mean',fontsize=24)
	ax[0].set_ylabel('height [m]',fontsize=20)
	
	ax[1].plot(data.time,data.offsetVar.sel(freq=94e9),c='C0',lw=2,label='W-Band')
	ax[1].plot(data.time,data.offsetVar.sel(freq=9.6e9),c='C1',lw=2,label='X-Band')
	ax[1].axhline(y=2,c='r',ls=':',lw=2)
	ax[1].set_ylabel('Variance of offset',fontsize=20)
	
	ax[2].plot(data.time,data.correlation.sel(freq=94e9),lw=2,label='W-Band')
	ax[2].plot(data.time,data.correlation.sel(freq=9.6e9),lw=2,label='X-Band')
	ax[2].axhline(y=0.7,c='r',ls=':',lw=2)
	ax[2].set_ylabel('correlation between Zes',fontsize=20)
	
	ax[3].plot(data.time,data.validPoints.sel(freq=94e9),lw=2,label='W-Band')
	ax[3].plot(data.time,data.validPoints.sel(freq=9.6e9),lw=2,label='X-Band')
	ax[3].axhline(y=300,c='r',ls=':',lw=2)
	ax[3].set_ylabel('Number of valid points',fontsize=20)
	for i,a in enumerate(ax):
		if i > 0:
			a.legend(fontsize=16)
		a.tick_params(labelsize=18)
	plt.setp(ax[3].xaxis.get_majorticklabels(), rotation=0,horizontalalignment='left')
	#plt.colorbar(p3,ax=ax[0],label='Ka -30 to -10 dB for W-Band')
	#ax[1].set_title('Ka-Band between -30 and -10dB for W-Band calibration')
	#plt.tight_layout()
	plt.savefig('plots/{date}_Flags_2hmean.png'.format(date=date.strftime('%Y%m%d')))
	plt.close()
	quit()

	fig,ax = plt.subplots(nrows=2,figsize=(15,10),constrained_layout=True,sharex=True)
	cmap = cm.get_cmap('coolwarm', 4)
	#p1 = ax[0].pcolormesh(data.time,data.range,(data.Ze.sel(freq=35.5e9)/data.Ze.sel(freq=35.5e9)).T,cmap = cmap,vmin=0.5,vmax=4.5)
	#plt.colorbar(p1,ax=ax[0])
	#p1 = ax[0].pcolormesh(dataCold.time,dataCold.range,(dataCold.Ze.sel(freq=35.5e9)/dataCold.Ze.sel(freq=35.5e9)*2).T,cmap = cmap,vmin=0.5,vmax=4.5)
	#plt.colorbar(p1,ax=ax[0])
	#ax.plot(MLrange.time,MLrange,label='ML height')
	#ax.legend()
	#ax[0].set_title('all values at 1km above ML')
	#p2 = ax[0].pcolormesh(dataXKa.time,dataXKa.range,(dataXKa.Ze.sel(freq=35.5e9)/dataXKa.Ze.sel(freq=35.5e9)*3).T,cmap = cmap,vmin=0.5,vmax=4.5)
	#plt.colorbar(p2,ax=ax[0],label='Ka -15 to 0 dB for X-Band')
	#ax[1].set_title('Ka-Band between -15 and 0dB for X-Band calibration')
	
	#p3 = ax[0].pcolormesh(dataKaW.time,dataKaW.range,(dataKaW.Ze.sel(freq=35.5e9)/dataKaW.Ze.sel(freq=35.5e9)*4).T,cmap = cmap,vmin=0.5,vmax=4.5)
	p1 = ax[0].pcolormesh(data.time,data.range,data.Ze.sel(freq=94e9).T,cmap = colmap.batlow,vmin=-30,vmax=20)
	cb = plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
	cb.set_label('Ze X-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[0].set_ylabel('height [m]',fontsize=20)
	ax[0].tick_params(labelsize=18)
	ax[1].plot(dataOffset.time,-1*dataOffset.offsetMean.sel(freq=94e9),c='C0',label='W-Band, new')
	ax[1].plot(dataOffset.time,-1*dataOffset.offsetMean.sel(freq=9.6e9),c='C1',label='X-Band, new')
	ax[1].plot(dataLV2.time,dataLV2.offset_w.sel(range=1000,method='nearest'),c='C2',label='W-Band, Jose')
	ax[1].plot(dataLV2.time,dataLV2.offset_x.sel(range=1000,method='nearest'),c='C3',label='X-Band, Jose')
	#ax[1].plot(dataOffset30.time,dataOffset30.offsetMean.sel(freq=94e9),c='C2',label='W-Band, 30 min mean')
	#ax[1].plot(dataOffset30.time,dataOffset30.offsetMean.sel(freq=9.6e9),c='C3',label='X-Band,30 min mean')
	ax[1].legend(fontsize=16)
	ax[1].set_ylabel('offset [dB]',fontsize=20)
	ax[1].tick_params(labelsize=18)
	#ax[1].set_xticklabels(rotation=0, ha='left')
	plt.setp(ax[1].xaxis.get_majorticklabels(), rotation=0,horizontalalignment='left')
	#plt.colorbar(p3,ax=ax[0],label='Ka -30 to -10 dB for W-Band')
	#ax[1].set_title('Ka-Band between -30 and -10dB for W-Band calibration')
	#plt.tight_layout()
	plt.savefig('{date}_test_corrrection_tripex_pol_joseHeight.png'.format(date=date.strftime('%y%m%d')))
	plt.show()
	quit()
	
	
	
	
	
	
	
