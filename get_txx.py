
from .SNR_text import *
from astropy.stats import bayesian_blocks


def get_txx(t,binsize = 0.5,sigma = 1,step_size = 1,block_n = 50,block_time = None,txx = 0.9,it = 1000,
	    bayesian = True,SNR =True,
	    time_unified = True,
	    hardness = 100

	    ):
	'''

	:param t: 光子的时间
	:param binsize: 切片的大小
	:param sigma: 信噪比判定标准
	:param step_size: 浮动块移动步长
	:param block_n: 浮动块大小
	:param block_time: or 浮动块时长
	:param txx:
	:param it: MCMC 迭代次数
	:param bayesian: bayesian 判定
	:param SNR: 信噪比判定
	:param time_unified: 背景时间标准化
	:param hardness: 背景硬度
	:return: 字典
	'''
	t = np.array(t)

	edges_bin = np.arange(t[0],t[-1]+binsize,binsize)

	bin_n,bin_edges = np.histogram(t,bins = edges_bin)

	t_c = (bin_edges[1:]+bin_edges[:-1])*0.5
	rate = bin_n/binsize

	SNR_result = SNR_text(t_c,rate,step_size = step_size,
			      block_n = block_n,block_time = block_time,
			      time_unified=time_unified,lambda_= hardness)
	bs = SNR_result['bs']
	good_index = SNR_result['good_index']
	pp = SNR_result['normallization']
	#贝叶斯
	time_edges = []
	max_SNR_list = []
	by_edges_list = []
	by_rate_list = []
	w = np.ones(len(rate))
	print(len(good_index))
	for one_index in good_index:
		#print('one_index:',one_index)
		t_c_in_one = t_c[one_index]
		t_in_one_start = t_c_in_one[0]-15#*binsize
		t_in_one_stop = t_c_in_one[-1]+15#*binsize
		print('range:',t_in_one_stop-t_in_one_start)
		if bayesian:
			t_in_one = t[np.where((t>=t_in_one_start)&(t<=t_in_one_stop))[0]]

			if len(t_in_one) <= 10000:

				edges = bayesian_blocks(t_in_one,fitness = 'events',gamma = np.exp(-12))
			else:
				t_b_index = np.where((t_c>=t_in_one_start)&(t_c<=t_in_one_stop))[0]
				t_b = t_c[t_b_index]
				bin_n_b = bin_n[t_b_index]
				#print(t_b)
				edges = bayesian_blocks(t_b,bin_n_b,fitness = 'events',gamma = np.exp(-12))
			print('edges: ',edges)
			if len(edges) > 3:#大于3才是有东西
				edges_c = (edges[1:] + edges[:-1]) * 0.5


				bin_n_by, bin_b_edges = np.histogram(t_in_one, bins=edges)  # 统计
				bin_b_size = bin_b_edges[1:] - bin_b_edges[:-1]
				bin_b_rate = bin_n_by / bin_b_size


				t_start, t_stop = found_edges(edges, bin_b_rate)
				w[np.where((t_c >= t_start) & (t_c <= t_stop))[0]] = 0
				#w[np.where((t_c >= t_start - 5 * binsize) & (t_c <= t_stop + 5 * binsize))[0]] = 0
				bs = WhittakerSmooth(rate,w,lambda_= hardness) #背景修正

				t_c_in_one1 = t_c[np.where((t_c >= edges[0]) & (t_c <= edges[-1]))[0]]
				bs_in_one = bs[np.where((t_c >= edges[0]) & (t_c <= edges[-1]))[0]]
				bs_mean = np.interp(edges_c, t_c_in_one1, bs_in_one)
				# bs_mean = bined_hist(t_c_in_one1,bs_in_one,bins = edges)


				SNR_block = ((bin_b_rate-bs_mean)/SNR_result['sigma'])
				bin_b_rate = np.concatenate((bin_b_rate[:1], bin_b_rate))
				max_SNR = np.max(SNR_block[1:-1])
				print('max_SNR',max_SNR)
				if SNR :

					if max_SNR > sigma:#说明有信噪比够好
						print('max_good')
						time_edges.append([t_start, t_stop])
						max_SNR_list.append(max_SNR)
						by_edges_list.append(bin_b_edges)
						by_rate_list.append(bin_b_rate)
					else:#信号太弱，忽略脉冲
						print('max_not_good')
						w[np.where((t_c >= t_start) & (t_c <= t_stop))[0]] = 1

						bs = WhittakerSmooth(rate, w, lambda_=hardness)

				else:
					time_edges.append([t_start, t_stop])
					max_SNR_list.append(max_SNR)
					by_edges_list.append(bin_b_edges)
					by_rate_list.append(bin_b_rate)
		else:
			w[np.where((t_c >= t_in_one_start ) & (t_c <= t_in_one_stop ))[
				0]] = 0
			time_edges.append([t_c_in_one[0],t_c_in_one[-1]])
	if(len(time_edges) == 0):
		print('we do not find a pulse !')
		result = {'good':False,
			't_c':t_c,
			'rate':rate,
			'bs':bs,
			'normallization':pp,
			'sigma':SNR_result['sigma']}
		if bayesian:
			print('have bayesian blocks.')
			result['bayesian_edges'] = by_edges_list
			result['bayesian_rate'] = by_rate_list
		return result
	time_edges = np.array(time_edges).T
	
	time_start = time_edges[0]
	time_stop = time_edges[1]
	print('accumulate')
	result = accumulate_counts(t_c,bin_n,np.sqrt(bin_n),w,time_start,time_stop,txx =txx,it = it,lamd=hardness)
	print('accumulate over')
	result['time_edges'] = time_edges
	result['t_c'] = t_c
	result['rate'] = rate
	result['normallization'] = pp
	result['sigma'] = SNR_result['sigma']
	result['bs'] = bs
	if bayesian:
		print('have bayesian blocks.')
		result['bayesian_edges'] = by_edges_list
		result['bayesian_rate'] = by_rate_list

	return result




def bined_hist(t,v,bins):
	t = np.array(t)
	v = np.array(v)
	re = []
	bin_start = bins[:-1]
	bin_stop = bins[1:]
	n = len(bin_start)
	for i in range(n):

		vin = v[np.where((t>=bin_start[i])&(t<=bin_stop[i]))[0]]
		aa = np.mean(vin)
		re.append(aa)
	return np.array(re)

def found_edges(edges,v):
	if len(edges) == 4:
		return 0.5*(edges[0]+edges[1]),0.5*(edges[2]+edges[3])
	edges  =  np.array(edges)
	v = np.array(v)
	edges_start = edges[:-1]
	edges_stop = edges[1:]
	edges_size = edges_stop - edges_start
	size_sort = np.sort(edges_size)[-2:]  #尺寸最大的块
	start_time = edges_start[0]
	stop_time = edges_stop[-1]
	trait = []
	fringe = 1
	cafe = 2
	pulse = 3
	for i in range(len(v)):
		if (i == 0):

			if(v[i] < v[i+1])and(edges_size[i] in size_sort):
				trait.append(cafe)
			else:
				trait.append(fringe)
				v[i] = v[i+1]+1
			if(v[len(v)-1] <= v[len(v)-2])and(edges_size[len(v)-1] not in size_sort):
				v[len(v)-1] = v[len(v)-2]+1
				
		elif(i == len(v)-1):

			if(v[i] < v[i-1])and(edges_size[i] in size_sort):
				trait.append(cafe)
			else:
				trait.append(fringe)

		else:
			if ((v[i] > v[i-1])and(v[i] > v[i+1])):
				trait.append(pulse)
			elif((v[i] < v[i-1])and(v[i]<v[i+1])):
				trait.append(cafe)
			else:
				trait.append(fringe)
	for i in range(len(v)):
		if(trait[i] == cafe):
			start_time = (edges_start[i]+edges_stop[i])*0.5
			break
	for i in range(len(v)):
		if(trait[-1-i] == cafe):
			stop_time = (edges_start[-1-i]+edges_stop[-1-i])*0.5
			break
	return start_time,stop_time

def accumulate_counts(t,n,n_err,w,t_start,t_stop,txx = 0.9,it = 1000,lamd = 100):
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
	txx = 1.-txx
	t = np.array(t)
	n = np.array(n)
	part_n = len(t_start)

	sigma = np.std(n[w > 0.5])
	bs1 = WhittakerSmooth(n, w, lamd)
	cs1 = n - bs1
	cs_f = np.cumsum(cs1)
	w1 = np.ones(len(cs_f))
	cs_fit = WhittakerSmooth(cs_f, w1, 1)
	#durti = t_stop-t_start
	ns = 3*sigma#*durti
	duration = t_stop - t_start
	index_sort = np.argsort(duration)[0]


	if len(np.where((t>=t_start[index_sort]) & (t<= t_stop[index_sort]))[0])<1000:#这里是为了提高精度
		d_t = (t_stop[index_sort]-t_start[index_sort])/1000
		print('dt for interp:',d_t)
		t_l = np.arange(t[0], t[-1]+d_t, d_t)
		cs_ff = np.interp(t_l, t, cs_fit)
	else:
		t_l = t
		cs_ff = cs_fit
	csf_fit_list = [0]
	for i in range(part_n):
		if(i < part_n-1):
			cs1_fit_max = np.mean(cs_f[np.where((t>t_stop[i])&(t<t_start[i+1]))[0]])
			csf_fit_list.append(cs1_fit_max)
		else:
			cs1_fit_max = np.mean(cs_f[np.where(t > t_stop[i])[0]])
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
			t90i, t1i, t2i = found_txx(t_l, cs_ff, l1i, l2i,csf_fit_list[index+1])
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

	while nnn < it:
		try:
			
			bin_ratexx = n + n_err * np.random.randn(len(n_err))
			bs11 = WhittakerSmooth(bin_ratexx, w, lamd)

			cs11 = bin_ratexx - bs11
			cs11_f = np.cumsum(cs11)
			cs11_fit = WhittakerSmooth(cs11_f, w1, 1)


			cs_ff1 = np.interp(t_l, t, cs11_fit)

			csf_fit_list1 = [0]

			for i in range(part_n):
				if (i < part_n - 1):
					cs1_fit_max = np.mean(cs11_f[np.where((t > t_stop[i]) & (t < t_start[i + 1]))[0]])
					csf_fit_list1.append(cs1_fit_max)
				else:
					cs1_fit_max = np.mean(cs11_f[np.where(t > t_stop[i])[0]])
					csf_fit_list1.append(cs1_fit_max)
			csf_fit_list1 = np.array(csf_fit_list1)

			dcsf_fit_list1 = csf_fit_list1[1:] - csf_fit_list1[:-1]
			pp = 0
			for index, dcsf in enumerate(dcsf_fit_list1):
				if index in index_i:
					
					if dcsf > ns:

						#print(dcsf,ns)
						dd = txx * dcsf
						l11 = dd + csf_fit_list1[index]
						l21 = csf_fit_list1[index + 1] - dd
						t90i, t1i, t2i = found_txx(t_l, cs_ff1, l11, l21,csf_fit_list1[index + 1])
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
		t90_err = np.std(t90_list[i_index])
		t90_err1i = t90[index] - t90_mean + t90_err
		t90_err2i = t90_mean + t90_err - t90[index]
		t90_err1.append(t90_err1i)
		t90_err2.append(t90_err2i)
		t1_mean = np.mean(t1_list[i_index])
		t2_mean = np.mean(t2_list[i_index])

		t1_err = np.std(t1_list[i_index])
		t2_err = np.std(t2_list[i_index])

		t1_err1i = t1[index] - t1_mean + t1_err
		t1_err2i = t1_mean + t1_err - t1[index]
		t1_err1.append(t1_err1i)
		t1_err2.append(t1_err2i)
		t2_err1i = t2[index] - t2_mean + t2_err
		t2_err2i = t2_mean + t2_err - t2[index]
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



def found_txx(t,v,st1,st2,base2):
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




























