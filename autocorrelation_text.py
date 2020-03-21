'''
该工具通过自相关检测来放大脉冲

'''
import numpy as np
from .background_kernel import *
from scipy.interpolate import interp1d


def autocorrelation_text(t,v,step_size = 1,block_n = 50,block_time = None,para = False,
			 time_unified = True):
	'''

	:param t:
	:param v:
	:param step_size:
	:param block_n:
	:param block_time:
	:param para:
	:param time_unified:
	:return:
	'''
	t = np.array(t)
	v = np.array(v)
	index_sort = np.argsort(t)
	t = t[index_sort]
	v = v[index_sort]


	dt = t[1]-t[0]
	if block_time is not None:
		if(block_time<=40*dt):
			print('block_time should bigger than ',40*dt)
		else:
			block_n = int(block_time/dt)
			print('block_n',block_n)
	index_block = np.arange(0,int(block_n),1)

	block_para = []
	block_index = []

	for i in range(int(t.size/step_size)):


		if(index_block[-1]+step_size*i<t.size):

			one_block_index = index_block + step_size*i
		else:

			one_block_index = -1-index_block
		block_index.append(one_block_index)
		para_t = np.mean(t[one_block_index])
		para_n = np.mean(v[one_block_index])
		para_n_std = np.var(v[one_block_index])  #方差
		block_para.append([para_t,para_n,para_n_std])
		
			

	block_para = np.array(block_para).T

	ACC =  block_para[1]*block_para[2]/np.sqrt((block_para[1]**2).sum()*(block_para[2]**2).sum())

	#ACC的背景处理。
	ACC_bf = Baseline_in_time(block_para[0],ACC,case = 'FD',time_unified=time_unified)
	#ACC_bf = AirPLS(ACC)
	ACC_cs = ACC - ACC_bf.bs
	#ACC_cs = ACC_cs - AirPLS(ACC_cs).bottom_airPLS()
	#ACC_cs = ACC_cs - Baseline_in_time(block_para[0],ACC_cs,fitness = 'bottom').bs
	#ACC_cs 均大于等于0
	binsize = np.percentile(ACC_cs,30)

	ACC_edges = np.arange(0,np.max(ACC_cs)+binsize,binsize)

	pn,pe = np.histogram(ACC_cs,bins = ACC_edges)
	pe_c = (pe[1:] + pe[:-1]) * 0.5
	pe_c = np.concatenate((pe[:1], pe_c))
	pn = np.concatenate((pn[:1],pn))

	if(len(pe[pn == 0])>0):
		es = pe[pn ==0][0]
	else:
		es = pe_c[-1]

	new_pe = np.linspace(0, es, 1000)
	new_pn = interp1d(pe_c, pn, kind='cubic')(new_pe)

	threshold = new_pe[new_pn<np.max(new_pn)*0.25][0] #得到阈值

	#index_back_of_block = np.where(ACC_cs<=threshold)[0]
	#index_signl_of_block = np.where(ACC_cs>threshold)[0]

	normallization = np.zeros(t.size)
	for index,value in enumerate(ACC_cs):

		if(value > threshold):

			normallization[block_index[index]] = normallization[block_index[index]] + 1
	#print('normal:\n',list(normallization))
	background_index = np.where(normallization < block_n)[0]

	if para :

		return background_index,normallization/block_n,block_index,ACC_cs-threshold,block_para[0]
	else:
		return background_index






























