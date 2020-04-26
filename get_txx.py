
from astropy.stats import bayesian_blocks
import numpy as np
from .Bayesian_duration import *
from .Baseline import TD_baseline,WhittakerSmooth

def get_txx(t,binsize = 0.064,background_degree = 7,sigma = 5,time_edges = None,txx = 0.9,it = 300,p0 = 0.05,plot_check = None,hardnss=100.):
	'''
	
	:param t: an 1D-array of event times.
	:param binsize: binsize of light curve.
	:param background_degree: the degree to which background floats. lt is used for background assessment.
	:param sigma: The threshold value of bayesian block signal strength used to judge the background area.
	:param time_edges: time_edges = [time_start,time_stop] ,time_start>=t.min time_stop <=t.max
	:param txx: if you want to estimate T90 ,txx = 0.9.
	:param it: Number of mcmc samples.
	:param p0: bayesian_blocks parameter
	:param plot_check:
	:param lamd_: Background line hardness
	:return:
	'''
	t = np.array(t)
	if time_edges is None:
		t_start = t.min()
		t_stop = t.max()
		
		bins = np.arange(t_start,t_stop,binsize)
	else:
		t_start = time_edges[0]
		t_stop = time_edges[1]
		t_index = np.where((t >= t_start)&(t<= t_stop+binsize))[0]
		t = t[t_index]
		bins = np.arange(t_start,t_stop,binsize)
	bin_n,bin_edges = np.histogram(t,bins = bins)
	rate = bin_n/binsize
	t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
	t_c,cs_rate,bs_rate = TD_baseline(t_c,rate)
	rate_sm = cs_rate+bs_rate.mean()
	#bin_n_sm = np.round(rate_sm*binsize)
	if binsize<0.064:
		if len(t)>40000:
			print('the size of t is ',len(t))
			print('binsize < 0.064, so will execute the following command:')
			print("bayesian_blocks(t,fitness='events',p0 = p0)")
			print('but the size of t is over 40000, longer computation times may be required.')
		edges = bayesian_blocks(t,fitness='events',p0 = p0)
		sizes = edges[1:]-edges[:-1]
		gg = np.around(edges/binsize)
		gg = np.unique(gg)
		gg = np.sort(gg)
		edges = gg*binsize
		print('min size :',sizes.min())
	else:
		edges = bayesian_blocks(t_c,bin_n,fitness='events',p0=p0)
	result = background_correction(t_c,rate_sm,edges,degree = background_degree,plot_save=plot_check)
	startedges,stopedges = get_bayesian_duration(result,sigma = sigma)
	w = np.ones(len(t_c))

	for ij in range(len(startedges)):
		index_w = np.where((t_c>=startedges[ij])&(t_c<=stopedges[ij]))[0]
		w[index_w] = 0

	c_rate = result['lc'][1]
	sigma = result['bkg'][2]
	re_rate = result['re_hist'][0]

	result1 = accumulate_counts(t_c,c_rate*binsize,np.sqrt(np.abs(c_rate*binsize)),w,startedges,stopedges,txx = txx,it = it,lamd = hardnss)
	result1['time_edges'] = [startedges,stopedges]
	result1['t_c'] = t_c
	result1['rate'] = c_rate
	result1['bs'] = WhittakerSmooth(c_rate,w,lambda_=hardnss/binsize**1.5)
	result1['good'] = True
	result1['xx'] = str(int(100*txx))
	result1['sigma'] = sigma
	result1['bayesian_edges'] = [edges]
	result1['bayesian_rate'] = [np.concatenate((re_rate[:1], re_rate))]
	return result1
	
def accumulate_counts(t,n,n_err,w,t_start,t_stop,txx = 0.9,it = 1000,lamd = 100.):
	'''

	:param t: lightcurve time
	:param n: lightcurve counts
	:param n_err: lightcurve conts err
	:param sigma: Background region sigma
	:param w: The weight
	:param t_start:
	:param t_stop:
	:param txx: if you want to estimate T90 ,txx = 0.9.
	:param it:
	:return:
	'''
	txx = 1. - txx
	t = np.array(t)
	n = np.array(n)
	dt = t[1]-t[0]
	dif = np.array(t_stop) - np.array(t_start)
	dif[dif > 20] = 20
	dif[dif <5] = 5
	lamd = lamd/dt**2
	tmin_ = t_start[0]-dif[0]
	if tmin_<t[0]:
		tmin_ = t[0]
	tmax_ = t_stop[-1]+dif[-1]
	if tmax_>t[-1]:
		tmax_ = t[-1]

	data_index = np.where((t>=tmin_)&(t<=tmax_))[0]
	
	part_n = len(t_start)

	sigma = np.std(n[w > 0.5])
	bs1 = WhittakerSmooth(n, w, lamd)
	cs1 = n - bs1
	cs_f = np.cumsum(cs1)
	w1 = np.ones(len(cs_f))
	cs_fit = WhittakerSmooth(cs_f, w1, 1)
	t_e = t[data_index]
	cs_fe = cs_f[data_index]
	#durti = t_stop-t_start
	ns = 3*sigma#*durti
	duration = t_stop - t_start
	index_sort = np.argsort(duration)[0]
	

	if len(np.where((t>=t_start[index_sort]) & (t<= t_stop[index_sort]))[0])<100:
		d_t = (t_stop[index_sort]-t_start[index_sort])/100
		print('dt for interp:',d_t)
		t_l = np.arange(tmin_, tmax_+d_t, d_t)
		cs_ff = np.interp(t_l, t, cs_fit)
	else:
		t_l = t[data_index]
		cs_ff = cs_fit[data_index]
	
	csf_fit_list = [np.mean(cs_fe[np.where(t_e < t_start[0])[0]])]
	for i in range(part_n):
		if(i < part_n-1):
			cs1_fit_max = np.mean(cs_fe[np.where((t_e>t_stop[i])&(t_e<t_start[i+1]))[0]])
			csf_fit_list.append(cs1_fit_max)
		else:
			cs1_fit_max = np.mean(cs_fe[np.where(t_e > t_stop[i])[0]])
			csf_fit_list.append(cs1_fit_max)
	csf_fit_list = np.array(csf_fit_list)
	dcsf_fit_list = csf_fit_list[1:]-csf_fit_list[:-1]
	t90 = []
	index_i = []
	t1 = []
	t2 = []
	l1 = []
	l2 = []

	for index,dcsf in enumerate(dcsf_fit_list):
		if dcsf > ns:
			dd = txx * dcsf
			l1i = dd + csf_fit_list[index]
			l2i = csf_fit_list[index+1] - dd
			#------------------------------------------------
			if 3*duration[index]>10:
				bb = 3*duration[index]
			else:
				bb = 10
			t1_range1 = t_start[index]-bb
			t2_range2 = t_stop[index]+bb
			if index <len(t_start)-1:
				
				if t2_range2>t_start[index+1]:
					t2_range2 = t_start[index+1]
			if(index>0):
				if t1_range1<t_stop[index-1]:
					t1_range1 = t_stop[index-1]
				
			#------------------------------------------------
			t90i, t1i, t2i = found_txx(t_l, cs_ff, l1i, l2i,csf_fit_list[index+1],t1_range1,t2_range2)
			t90.append(t90i)
			t1.append(t1i)
			t2.append(t2i)
			l1.append(l1i)
			l2.append(l2i)
			index_i.append(index)
	fit_max = csf_fit_list[index_i+[index_i[-1]+1]]
	if len(t90)<1:
		return {'good':False}

	t90_list = []
	t1_list = []
	t2_list = []
	bs_list = []
	index_list = []
	nnn = 0
	nnnn = 0
	while nnn < it:
		nnnn = nnnn+1
		try:
			
			bin_ratexx = n + n_err * np.random.randn(len(n_err))
			bs11 = WhittakerSmooth(bin_ratexx, w, lamd)

			cs11 = bin_ratexx - bs11
			cs11_f = np.cumsum(cs11)
			cs11_fit = WhittakerSmooth(cs11_f, w1, 1)
			cs11_fe = cs11_f[data_index]

			cs_ff1 = np.interp(t_l, t, cs11_fit)

			csf_fit_list1 = [np.mean(cs11_fe[np.where(t_e < t_start[0])[0]])]

			for i in range(part_n):
				if (i < part_n - 1):
					cs1_fit_max = np.mean(cs11_fe[np.where((t_e > t_stop[i]) & (t_e < t_start[i + 1]))[0]])
					csf_fit_list1.append(cs1_fit_max)
				else:
					cs1_fit_max = np.mean(cs11_fe[np.where(t_e > t_stop[i])[0]])
					csf_fit_list1.append(cs1_fit_max)
			csf_fit_list1 = np.array(csf_fit_list1)

			dcsf_fit_list1 = csf_fit_list1[1:] - csf_fit_list1[:-1]
			pp = 0
			for index, dcsf in enumerate(dcsf_fit_list1):
				if index in index_i:
					
					if dcsf > dcsf_fit_list[index]*0.1:

						#print(dcsf,ns)
						dd = txx * dcsf
						l11 = dd + csf_fit_list1[index]
						l21 = csf_fit_list1[index + 1] - dd
						
						if 3*t90[index]>10:
							bb = 3*t90[index]
						else:
							bb = 10
						t1_range1 = t_start[index]-bb
						t1_range2 = t_start[index]+bb
						t2_range1 = t_stop[index]-bb
						t2_range2 = t_stop[index]+bb
						
						if index <len(t_start)-1:
							if t1_range2>t_start[index+1]:
								t1_range2 = t_start[index+1]
							if t2_range2>t_start[index+1]:
								t2_range2 = t_start[index+1]
						if(index>0):
							if t1_range1<t_stop[index-1]:
								t1_range1 = t_stop[index-1]
							if t2_range1<t_stop[index-1]:
								t2_range1 = t_stop[index-1]
						#print(t1_range1,t1_range2,t2_range1,t2_range2)
						t90i, t1i, t2i = found_txx(t_l, cs_ff1, l11, l21,
						                           csf_fit_list1[index + 1],t1_range1,t2_range2)
						if ((t1i < t1_range2 and t1i > t1_range1) and (t2i > t2_range1 and t2i<t2_range2)):

							t90_list.append(t90i)
							t1_list.append(t1i)
							t2_list.append(t2i)
							index_list.append(index)
							pp = pp + 1

			if pp > 0:
				print(nnn,end = '\r')
				bs_list.append(bs11)
				nnn = nnn + 1

		except:
			continue
		if nnnn > it * 100:
			return {'good': False}
	t90_list = np.array(t90_list)
	t1_list = np.array(t1_list)
	t2_list = np.array(t2_list)
	index_list = np.array(index_list)

	new_t90_list = []
	new_t1_list = []
	new_t2_list = []
	t90_err1 = []
	t90_err2 = []
	t1_err1 = []
	t1_err2 = []
	t2_err1 = []
	t2_err2 = []
	for index,i in enumerate(index_i):
		i_index = np.where(index_list == i)
		t90_mean = np.mean(t90_list[i_index])
		t90_err = np.std(t90_list[i_index],ddof=0)
		t90_err1i = t90[index] - t90_mean + t90_err
		t90_err2i = t90_mean + t90_err - t90[index]
		if t90_err1i<0:
			t90_err1i = 0
		if t90_err2i<0:
			t90_err2i = 0
		t90_err1.append(t90_err1i)
		t90_err2.append(t90_err2i)
		t1_mean = np.mean(t1_list[i_index])
		t2_mean = np.mean(t2_list[i_index])

		t1_err = np.std(t1_list[i_index],ddof=0)
		t2_err = np.std(t2_list[i_index],ddof=0)

		t1_err1i = t1[index] - t1_mean + t1_err
		t1_err2i = t1_mean + t1_err - t1[index]
		if t1_err1i<0:
			t1_err1i = 0
		if t1_err2i<0:
			t1_err2i = 0
		t1_err1.append(t1_err1i)
		t1_err2.append(t1_err2i)
		t2_err1i = t2[index] - t2_mean + t2_err
		t2_err2i = t2_mean + t2_err - t2[index]
		if t2_err1i<0:
			t2_err1i = 0
		if t2_err2i<0:
			t2_err2i = 0
		t2_err1.append(t2_err1i)
		t2_err2.append(t2_err2i)

		new_t90_list.append(t90_list[i_index])
		new_t1_list.append(t1_list[i_index])
		new_t2_list.append(t2_list[i_index])

	result = {'good':True,
		  'txx':t90,'txx_err':[t90_err1,t90_err2],
		  't1':t1,'t2':t2,'t1_err':[t1_err1,t1_err2],'t2_err':[t2_err1,t2_err2],
		  'txx_list':new_t90_list,'t1_list':new_t1_list,'t2_list':new_t2_list,
		  'cs_f_max':fit_max,'cs_f':cs_f,
		  't':t,'n':n,'bs_list':bs_list,'bs1':bs1,
		  'l':[l1,l2]}
	return result



def found_txx(t,v,st1,st2,base2,tmin,tmax):
	index_ = np.where((t>=tmin)&(t<=tmax))
	t = t[index_]
	v = v[index_]
	t1 = []
	for i in range(len(t)):
		if v[i] >= st1:
			if v[i] >= st2:
				if v[i]>base2:
					break
			else:
				t1.append(t[i])

		else:
			t1 = []
	t90 = t1[-1]-t1[0]
	return t90,t1[0],t1[-1]




























