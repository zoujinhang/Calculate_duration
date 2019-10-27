
from .SNR_text import *
from astropy.stats import bayesian_blocks


def get_txx(t,binsize = 0.5,criterion = 1,step_size = 1,block_n = 50,block_time = None,txx = 0.05,it = 1000,
	    bayesian = True,SNR =True):
	'''

	:param t:
	:param binsize:
	:param criterion:
	:param step_size:
	:param block_n:
	:param block_time:
	:param txx:
	:param it:
	:param bayesian:
	:param SNR:
	:return:
	'''
	t = np.array(t)

	edges_bin = np.arange(t[0],t[-1]+binsize,binsize)

	bin_n,bin_edges = np.histogram(t,bins = edges_bin)

	t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
	rate = bin_n/binsize

	SNR_result = SNR_text(t_c,rate,criterion=criterion,step_size = step_size,
			      block_n = block_n,block_time = block_time,
			      SNR = SNR)
	bs = SNR_result['bs']
	good_index = SNR_result['good_index']

	#贝叶斯
	time_edges = []
	max_SNR_list = []
	by_edges_list = []
	by_rate_list = []
	w = np.ones(len(rate))
	for one_index in good_index:

		t_c_in_one = t_c[one_index]
		t_in_one_start = t_c_in_one[0]-5
		t_in_one_stop = t_c_in_one[-1]+5
		if bayesian:
			t_in_one = t[np.where((t>=t_in_one_start)&(t<=t_in_one_stop))[0]]

			if len(t_in_one) <= 10000:

				edges = bayesian_blocks(t_in_one,fitness = 'events')
			else:
				t_b_index = np.where((t_c>=t_in_one_start)&(t_c<=t_in_one_stop))[0]
				t_b = t_c[t_b_index]
				bin_n_b = bin_n[t_b_index]
				edges = bayesian_blocks(t_b,bin_n_b,fitness = 'events')

			if len(edges > 3):#大于3才是有东西

				bs_in_one = bs[np.where((t_c>=t_in_one_start+5)&(t_c<=t_in_one_stop-5))[0]]
				w[np.where((t_c>=t_in_one_start+5)&(t_c<=t_in_one_stop-5))[0]] = 0
				bs_mean = np.mean(bs_in_one)

				bin_n_by ,bin_b_edges = np.histogram(t_in_one,bins = edges)

				bin_b_size = bin_b_edges[1:] - bin_b_edges[:-1]

				bin_b_rate = bin_n_by/bin_b_size
				bin_b_rate = np.concatenate((bin_b_rate[:1],bin_b_rate))

				max_SNR = (np.max(bin_b_rate)-bs_mean)/SNR_result['sigma']
				if SNR :

					if max_SNR >1:#说明有信噪比够好
						time_edges.append([edges[1], edges[-2]])
						max_SNR_list.append(max_SNR)
						by_edges_list.append(bin_b_edges)
						by_rate_list.append(bin_b_rate)
				else:
					time_edges.append([edges[1], edges[-2]])
					max_SNR_list.append(max_SNR)
					by_edges_list.append(bin_b_edges)
					by_rate_list.append(bin_b_rate)
		else:
			time_edges.append([t_c_in_one[0],t_c_in_one[-1]])
	time_start = time_edges[0][0]
	time_stop = time_edges[-1][-1]
	result = accumulate_counts(t_c,bin_n,np.sqrt(bin_n),SNR_result['sigma'],w,time_start,time_stop,txx =txx,it = it)
	result['time_edges'] = time_edges
	result['t_c'] = t_c
	result['rate'] = rate
	if bayesian:
		result['bayesian_edges'] = by_edges_list
		result['bayesian_rate'] = by_rate_list

	return result


def accumulate_counts(t,n,n_err,sigma,w,t_start,t_stop,txx = 0.05,it = 1000):
	'''

	:param t: 时间
	:param n: 计数
	:param n_err: 计数误差
	:param sigma: 背景区平均sigma
	:param w: 权重
	:param t_start: 开始时间
	:param t_stop: 结束时间
	:param txx:
	:param it: 迭代次数
	:return:
	'''
	t = np.array(t)
	n = np.array(n)
	bs1 = WhittakerSmooth(n, w, 100)
	cs1 = n - bs1
	cs_f = np.cumsum(cs1)
	w1 = np.ones(len(cs_f))
	cs_fit = WhittakerSmooth(cs_f, w1, 1)
	#durti = t_stop-t_start
	ns = 0.5*sigma#*durti
	if len(t)<1000:
		t_l = np.linspace(t[0], t[-1], 1000)
		cs_ff = np.interp(t_l, t, cs_fit)
	else:
		t_l = t
		cs_ff = cs_fit

	cs1_fit_max = np.mean(cs_f[np.where(t>t_stop)[0]])
	if cs1_fit_max < ns:
		return {'good':False}
	dd = txx*cs_ff
	l1 = dd
	l2 = cs1_fit_max - dd
	t90, t1, t2 = found_txx(t_l, cs_ff, l1, l2)

	t90_list = []
	t1_list = []
	t2_list = []
	bs_list = []

	i = 0
	while i < it:
		try:
			bin_ratexx = n + n_err * np.random.randn(len(n_err))
			bs11 = WhittakerSmooth(bin_ratexx, w, 100)

			cs11 = bin_ratexx - bs11
			cs11_f = np.cumsum(cs11)
			# cs11_fit = cs11_f
			cs11_fit = WhittakerSmooth(cs11_f, w1, 1)
			# cs11_fit = cs11_f+np.sqrt(cs11_f)*np.random.randn(len(cs11_f))
			# t_l = np.linspace(bin_t[0],bin_t[-1],1000)
			if len(t) < 1000:
				t_l = np.linspace(t[0], t[-1], 1000)
				cs_ff1 = np.interp(t_l, t, cs11_fit)
			else:
				t_l = t
				cs_ff1 = cs_fit
			#cs_ff1 = np.interp(t_l,t, cs11_fit)
			cs11_fit_max = np.mean(cs11_f[np.where(t>t_stop)[0]])
			if cs11_fit_max >= ns:
				dd = 0.05 * cs11_fit_max
				# plt.plot(bin_t, cs11_fit, color='r',alpha=0.1)
				l11 = dd
				l21 = cs1_fit_max - dd
				t901, t11, t21 = found_txx(t_l, cs_ff1, l11, l21)
				if ((t11 <= t_stop )&(t21 >= t_start)):
					bs_list.append(bs11)
					t90_list.append(t901)
					t1_list.append(t11)
					t2_list.append(t21)
					i = i + 1
		except:
			continue
	t90_list = np.array(t90_list)
	t90_err = np.std(t90_list)
	t1_list = np.array(t1_list)
	t2_list = np.array(t2_list)

	t90_mean = t90_list.mean()
	t90_err1 = t90 - t90_mean + t90_err
	t90_err2 = t90_mean + t90_err - t90

	t1_mean = t1_list.mean()
	t2_mean = t2_list.mean()

	t1_err = np.std(t1_list)
	t2_err = np.std(t2_list)

	t1_err1 = t1 - t1_mean + t1_err
	t1_err2 = t1_mean + t1_err - t1

	t2_err1 = t2 - t2_mean + t2_err
	t2_err2 = t2_mean + t2_err - t2

	result = {'good':True,
		  'txx':t90,'txx_err':[t90_err1,t90_err2],
		  't1':t1,'t2':t2,'t1_err':[t1_err1,t1_err2],'t2_err':[t2_err1,t2_err2],
		  'txx_list':t90_list,'t1_list':t1_list,'t2_list':t2_list,
		  'cs_f_max':cs1_fit_max,'cs_f':cs_f,
		  't':t,'n':n,'bs_list':bs_list,'bs':bs1,
		  'l':[l1,l2]}
	return result



def found_txx(t,v,st1,st2):
	t1 = []
	for i in range(len(t)):
		if v[i] >= st1:
			if v[i] >= st2:
				break
			t1.append(t[i])

		else:
			t1 = []
	t90 = t1[-1]-t1[0]
	return t90,t1[0],t1[-1]




























