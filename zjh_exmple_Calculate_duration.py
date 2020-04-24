
import numpy as np
import matplotlib.pyplot as plt
import os
from .zzh_py3_file import readcol,findfile
from .Calculate_duration import Plot,get_txx,save_result
from astropy.io import fits


infor_link = './list.txt'
savedir = './result/'
data_top = '/your_path/trigdata/'


sample_list,dete_list,year_list,el,eu,tl,tu = readcol(infor_link)
start = -20
stop = 10


def ch_to_energy(time,ch,ch_n,e1,e2):
	'''
	
	:param time:
	:param ch:
	:param ch_n:
	:param e1:
	:param e2:
	:return:
	'''
	new_t = np.array([])
	new_energy = np.array([])
	for index,channel in enumerate(ch_n):
		ch_t_index = np.where(ch == channel)
		ch_t = time[ch_t_index]
		energy_array = get_energy_of_ch(ch_t,e1[index],e2[index])
		new_t = np.concatenate((new_t,ch_t))
		new_energy = np.concatenate((new_energy,energy_array))
	index_all = np.argsort(new_t)
	new_t = new_t[index_all]
	new_energy = new_energy[index_all]
	return new_t,new_energy




for i in range(len(sample_list)):
	print(sample_list[i])
	datalink = data_top+str(year_list[i]) + '/'+sample_list[i] + '/'
	filename = findfile(datalink,'glg_tte_'+dete_list[i]+'_'+sample_list[i]+'_v*')[0]
	
	hl = fits.open(datalink + filename)
	trigtime = hl[0].header['TRIGTIME']
	time = hl[2].data.field(0)
	ch = hl[2].data.field(1)
	ch_n = hl[1].data.field(0)
	e1 = hl[1].data.field(1)
	e2 = hl[1].data.field(2)
	t = time - trigtime
	t, e = ch_to_energy(t, ch, ch_n, e1, e2)
	e_index = np.where((e >= el[i]) & (e <= eu[i]))[0]
	t = t[e_index]
	if os.path.exists(savedir) == False:
		os.makedirs(savedir)
	
	result = get_txx(t,binsize = 0.01,time_edges=[start,stop],background_degree=7,sigma = 5,txx = 0.9,it = 300,prior = 5,plot_check=savedir + 'Z_'+sample_list[i]+'_check.png',hardnss=100)
	#----------------------------------------------------------------------------
	
	save_result(result,savedir + 'C_'+sample_list[i]+'_T90.csv',float_format='%.3f')
	#----------------------------------------------------------------------------
	
	
	myplt = Plot(result)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_light_curve(sigma=4.5)
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'A_'+sample_list[i]+'_lightcurve.png')
	plt.close()
	if result['good']:
		for ij in range(len(result['txx'])):
			plt.title(sample_list[i])
			myplt.plot_distribution('90',num = ij)
			plt.savefig(savedir + 'D_'+sample_list[i]+'_distribution_'+str(ij)+'.png')
			plt.close()
	#--------------------------------------------------------------------------
	#
	plt.figure(figsize = (10,10))
	plt.subplot(2,1,1)
	plt.title('GRB'+sample_list[i][2:])
	myplt.plot_Txx1('90')#
	plt.xlim(tl[i],tu[i])
	plt.subplot(2,1,2)
	myplt.plot_Txx2('90')
	plt.xlim(tl[i],tu[i])
	plt.savefig(savedir + 'B_'+sample_list[i]+'_txx.png')
	plt.close()
