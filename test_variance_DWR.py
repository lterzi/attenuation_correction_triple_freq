from glob import glob
import numpy as np
import pandas as pd
from sys import argv
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import cm
from cmcrameri import cm as colmap
from attenuationFunctions import * 

pathCN = '/archive/meteo/external-obs/juelich/cloudnet/jue/'
pathData = '/archive/meteo/external-obs/juelich/' #'/scratch/l/L.Terzi/campaign_aux_data/'
pathPro = pathData+'campaign_aux_data/tripex-pol-scan/wband_scan/processed/' # this contains the regridded and dealized spectra
pathRes = pathData + 'campaign_aux_data/tripex-pol-scan/wband_scan/resampled/' # this contains resampled moments of w-band
date2start = '20220206'
date2end = '20220206'
dateStart = pd.to_datetime(date2start); dateEnd = pd.to_datetime(date2end)
dateList = pd.date_range(dateStart, dateEnd,freq='d')
freq=[9.6e9,35.5e9,94.0e9]
for date in dateList:
	#filesKa = glob('../output/{year}/{month}/{day}/*_ZEN_regridded_LV0.nc'.format(year=date.strftime('%Y'),month=date.strftime('%m'),day=date.strftime('%d')))
	#dataZEN = xr.Dataset()
	#for f in filesKa:
	#	d = xr.open_dataset(f)
	#	dataZEN = xr.merge([dataZEN,d[['ZeH_X','ZeH_Ka','ZeH_W','NoisePowH_X','NoisePowH_Ka','NoisePowH_W']]])
	#dataZEN.to_netcdf('test_files/{date}_dataZEN_mom.nc'.format(date=date.strftime('%Y%m%d')))		
	#dataZEN = xr.open_dataset('test_files/{date}_dataZEN_mom.nc'.format(date=date.strftime('%Y%m%d')))		
	dataZEN = xr.open_mfdataset('../output/2022/02/06/20220206_*_ZEN_regridded_LV1.nc')
	#print(dataZEN)
	#quit()
	CN = xr.open_dataset('test_files/{date}_cloudnet_attenuation.nc'.format(date=date.strftime('%Y%m%d')))		
	#- read in radar Data and reindex to common grid, also, for easier handling: rename all variables to the same name and have freq as coordinate
	ZENW = xr.open_dataset('{path}{date}_ZEN_moments_wband_scan.nc'.format(path=pathRes,date=date.strftime('%Y%m%d')))
	ZENW = ZENW.expand_dims({'freq':np.array([94e9])})
	ZENKa = getJoyrad(pathData,date,ZENW.range,ZENW.time,'Ka',offset=3.0) 
	ZENX = getJoyrad(pathData,date,ZENW.range,ZENW.time,'X',offset=1.6) 
	#- merge triple-frequency dataset and calculated attenuation:
	CN = CN[['atm_att_2way','temperature']].rename({'model_time':'time','model_height':'range'})
	CN = CN.interp({'time':ZENW.time,'range':ZENW.range})
	data = xr.merge([ZENKa,ZENW,ZENX,CN])
	data['calibration_offset'] = data.calibration_offset.fillna(0)
	#- add gas attenuation to Ze
	data['Ze'] = data.Ze + data.atm_att_2way + data.calibration_offset
	#dataZEN['ZeH'] = dataZEN.ZeH + data.atm_att_2way + data.calibration_offset
	data['DWRKaW'] = data.Ze.sel(freq=35.5e9) - data.Ze.sel(freq=94e9)
	data['DWRXKa'] = data.Ze.sel(freq=9.6e9) - data.Ze.sel(freq=35.5e9)
	dataKaW = getInterval(data,interval=[-30,-10])
	dataXKa = getInterval(data,interval=[-15,0])
	
	dataAtt = xr.open_dataset('../output/{year}/{month}/{day}/{date}_tripex_pol_scan_3fr_lv2_30min.nc'.format(year=date.strftime('%Y'),
																							month = date.strftime('%m'),
																							day = date.strftime('%d'),
																							date = date.strftime('%Y%m%d')))
	T0 = getclosestTemp(data,273.15)
	SNRKa = 10*np.log10(dataZEN.SNRH.sel(freq=35.5e9))
	#plt.pcolormesh(dataZEN.time,dataZEN.range,10*np.log10(dataZEN.NoisePowH_W).T)
	#plt.pcolormesh(SNRKa.time,SNRKa.range,SNRKa.T)
	#plt.show()
	for t in data.time:
		datasel = data.sel(time=t)
		dataselKaW = dataKaW.sel(time=t)
		if not np.isnan(dataselKaW.DWRKaW).all():
			print(t.values)
			fig,ax=plt.subplots(ncols=2,figsize=(10,7),sharey=True)
			ax[0].plot(datasel.Ze.sel(freq=35.5e9),datasel.range,lw=2,alpha=0.5,c='C0')
			ax[0].plot(datasel.Ze.sel(freq=94e9),datasel.range,lw=2,alpha=0.5,c='C1')
			ax[0].plot(dataselKaW.Ze.sel(freq=35.5e9),dataselKaW.range,lw=2,c='C0',label='Ka-Band')
			ax[0].plot(dataselKaW.Ze.sel(freq=94e9),dataselKaW.range,lw=2,c='C1',label='W-Band')
			ax[0].plot(SNRKa.sel(time=t,method='nearest'),SNRKa.range,c='C2',lw=2,label='SNR Ka')
			
			#ax[0].plot(10*np.log10(dataZEN.ZeH.sel(time=t,freq=35.5e9,method='nearest')),dataZEN.range,c='C3',label='Ka LV1')
			#ax[0].plot(10*np.log10(dataZEN.ZeH.sel(time=t,freq=94e9,method='nearest')),dataZEN.range,c='C4',label='W LV1')
			ax[0].axhline(y=T0.sel(time=t).values,c='r',ls='--',lw=2,label='0°C')
			ax[0].set_xlabel('Ze [dBz]',fontsize=20)
			ax[0].set_ylabel('range [m]',fontsize=20)
			#ax[1].set_xlabel('Ze W-Band [dBz]',fontsize=20)
			ax[1].plot(datasel.DWRKaW,datasel.range,lw=2,alpha=0.5,c='C0')
			ax[1].plot(dataselKaW.DWRKaW,dataselKaW.range,lw=2,label='Variance of DWR:{:.2f}'.format(dataAtt.offsetVar.sel(time=t,freq=94e9,method='nearest')),c='C0')
			#ax[1].plot(10*np.log10(dataZEN.ZeH.sel(time=t,freq=35.5e9,method='nearest')/dataZEN.ZeH.sel(time=t,freq=94e9,method='nearest')),dataZEN.range,c='C1',label='DWR LV1')
			ax[1].axhline(y=T0.sel(time=t).values,c='r',ls='--',lw=2)
			ax[1].set_xlabel('DWR KaW [dB]',fontsize=20)
			
			for a in ax:
				a.legend(fontsize=16)
				a.grid()
				a.tick_params(labelsize=18)
				
			plt.suptitle('{}'.format(pd.to_datetime(t.values).strftime('%Y%m%d %H:%M:%S')),fontsize=24)
			plt.tight_layout()
			plt.savefig('plots/test_variance/KaW/{}_test_variance_KaW.png'.format(pd.to_datetime(t.values).strftime('%Y%m%d_%H%M%S')))
			plt.close()
	'''
	for t in data.time:
		datasel = data.sel(time=t)
		dataselXKa = dataXKa.sel(time=t)
		if not np.isnan(dataselXKa.DWRXKa).all():
			print(t.values)
			fig,ax=plt.subplots(ncols=2,figsize=(10,7),sharey=True)
			ax[0].plot(datasel.Ze.sel(freq=35.5e9),datasel.range,lw=2,alpha=0.5,c='C0')
			ax[0].plot(datasel.Ze.sel(freq=9.6e9),datasel.range,lw=2,alpha=0.5,c='C1')
			ax[0].plot(dataselXKa.Ze.sel(freq=35.5e9),dataselXKa.range,lw=2,c='C0',label='Ka-Band')
			ax[0].plot(dataselXKa.Ze.sel(freq=9.6e9),dataselXKa.range,lw=2,c='C1',label='X-Band')
			ax[0].axhline(y=T0.sel(time=t).values,c='r',ls='--',lw=2,label='0°C')
			ax[0].set_xlabel('Ze [dBz]',fontsize=20)
			ax[0].set_ylabel('range [m]',fontsize=20)
			#ax[1].set_xlabel('Ze W-Band [dBz]',fontsize=20)
			ax[1].plot(datasel.DWRXKa,datasel.range,lw=2,alpha=0.5,c='C0')
			ax[1].plot(dataselXKa.DWRXKa,dataselXKa.range,lw=2,label='Variance of DWR:{:.2f}'.format(dataAtt.offsetVar.sel(time=t,freq=9.6e9,method='nearest')),c='C0')
			ax[1].axhline(y=T0.sel(time=t).values,c='r',ls='--',lw=2,label='0°C')
			ax[1].set_xlabel('DWR XKa [dB]',fontsize=20)
			
			for a in ax:
				a.legend(fontsize=16)
				a.grid()
				a.tick_params(labelsize=18)
				
			plt.suptitle('{}'.format(pd.to_datetime(t.values).strftime('%Y%m%d %H:%M:%S')),fontsize=24)
			plt.tight_layout()
			plt.savefig('plots/test_variance/XKa/{}_test_variance_XKa.png'.format(pd.to_datetime(t.values).strftime('%Y%m%d_%H%M%S')))
			plt.close()
	'''
	quit()		
	fig,ax = plt.subplots(nrows=2,figsize=(15,10))
	p1=ax[0].pcolormesh(data.time,data.range,data.DWRKaW.T,vmin=-5,vmax=5,cmap='jet')
	cb = plt.colorbar(p1,ax=ax[0],pad=0.01,aspect=30)
	cb.set_label('DWR KaW [dB]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[0].tick_params(labelsize=18)
	ax[0].set_ylabel('range [m]',fontsize=20)
	p1=ax[1].pcolormesh(data.time,data.range,data.Ze.sel(freq=94e9).T,vmin=-40,vmax=10,cmap='jet')
	cb = plt.colorbar(p1,ax=ax[1],pad=0.01,aspect=30)
	cb.set_label('Ze W-Band [dBz]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	ax[1].tick_params(labelsize=18)
	ax[1].set_ylabel('range [m]',fontsize=20)
	plt.show()
	print(data)
	quit()
