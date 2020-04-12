import numpy as np
import matplotlib.pyplot as plt


def background_correction(t,rate,edges,backgroundsize = 10,degree = 5,plot_save = None):
	'''
	
	:param t:
	:param rate:
	:param edges:
	:param degree:
	:return:
	'''
	re_rate,re_sigma,index_list = re_histogram(t,rate,edges)
	dt = t[1]-t[0]
	binsize = edges[1:]-edges[:-1]
	sort_index = np.argsort(-binsize)
	sort_binsize = binsize[sort_index]
	sort_sigma = re_sigma[sort_index]
	sort_re_rate = re_rate[sort_index]
	sort_index_list = np.array(index_list)[sort_index]
	background_pool = rate[sort_index_list[0]]
	correction_t = [edges[0],edges[-1]+dt,t[sort_index_list[0]][0],t[sort_index_list[0]][-1]]
	mean1 = sort_re_rate[0]
	correction_rate = [rate[0],rate[-1],mean1,mean1]
	sigma1 = sort_sigma[0]
	binsize1 = sort_binsize[0]
	n = 1
	for i in range(1,len(sort_binsize)):
		if binsize1>backgroundsize:
			binsize1_ = backgroundsize
		else:
			binsize1_ = binsize1
		if sort_binsize[i] > backgroundsize:
			sort_binsizei_ = backgroundsize
		else:
			sort_binsizei_ = sort_binsize[i]
		if confidence_analysis(mean1,sigma1,binsize1_,sort_re_rate[i],sort_sigma[i],sort_binsizei_,dt,degree) and len(t[sort_index_list[i]])>0:
			n = n+1
			correction_t.append(t[sort_index_list[i]][0])
			correction_t.append(t[sort_index_list[i]][-1])
			correction_rate.append(sort_re_rate[i])
			correction_rate.append(sort_re_rate[i])
			background_pool = np.concatenate((background_pool,rate[sort_index_list[i]]))
			mean1 = background_pool.mean()
			sigma1 = background_pool.std(ddof = 1)
			binsize1 = binsize1 + sort_binsize[i]
	print('The number of background blocks is %d.'%n)
	correction_t = np.array(correction_t)
	correction_rate = np.array(correction_rate)
	sort_t_index = np.argsort(correction_t)
	correction_t = correction_t[sort_t_index]
	correction_rate = correction_rate[sort_t_index]
	if plot_save is not None:
	
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.step(edges,np.concatenate((re_rate[:1],re_rate)))
		ax.plot(correction_t,correction_rate)
		fig.savefig(plot_save)
		plt.close(fig)
	
	correction_b = np.interp(t,correction_t,correction_rate)
	new_rate = rate - correction_b + correction_b.mean()
	re_rate,re_sigma,index_list = re_histogram(t,new_rate,edges)
	result = {'lc':(t,new_rate),
	          're_hist':(re_rate,re_sigma,index_list),
	          'bkg':(mean1,sigma1,binsize1),
	          'edges':edges,
	          'dt':dt
	          }
	return result

def get_SNR(edges,re_rate,background_mean,background_sigma,dt,non_negative = True):
	'''
	
	:param edges:
	:param re_rate:
	:param background_mean:
	:param background_sigma:
	:param dt:
	:param non_negative:
	:return:
	'''
	binsize = edges[1:] - edges[:-1]
	sigma = background_sigma * np.sqrt(dt/binsize)
	SNR = (re_rate - background_mean)/sigma
	if non_negative:
		if len(SNR[SNR<0])>0:
			SNR[SNR<0] = 0 #we only care about the part where the SNR is greater than 0.
	return SNR

def re_histogram(t,rate,edges):
	'''
	This is a stable function, which can re-bin the data 't' and 'rate' with the new bin edges 'edges'.
	
	:param t: The time of light curve.
	:param rate: The counts rate of light curve.
	:param edges: the new bin-edges, whose interval can different from each other.
	:return: three array ,The rebined rate and sigma of each bins.
	'''


	re_rate = np.zeros(len(edges)-1)
	re_sigma = np.zeros(len(edges)-1)
	
	index_list = []
	
	strat = edges[:-1]
	stop = edges[1:]
	for i in range(len(strat)):
		indexs = np.where((t>=strat[i])&(t<stop[i]))[0]
		index_list.append(indexs)
	'''
	index = np.where(t >= edges[0])[0]
	t = t[index]
	rate = rate[index]
	edges_num = 1
	ones = []
	for index1,value in enumerate(t):
		
		
		if value > edges[edges_num] or value == t[-1]:
			ones.append(index1)
			index_list.append(ones)
			ones = []
			edges_num = edges_num+1
			if edges_num == len(edges):
				break
		else:
			ones.append(index1)
	'''
	for index1,value in enumerate(index_list):
		one_rate = rate[value]
		index_list[index1] = value#+index[0]
		if len(one_rate) > 0:
			re_rate[index1] = one_rate.mean()
			if len(one_rate) == 1:
				re_sigma[index1] = one_rate.std(ddof=0)
			else:
				re_sigma[index1] = one_rate.std(ddof=1)
		else:
			if index1 == len(index_list)-1:
				re_rate[index1] = re_rate[index1-1]
				re_sigma[index1] = re_sigma[index1-1]
			else:
				re_rate[index1] = 0
				re_sigma[index1] = 1
				
	return re_rate,re_sigma,np.array(index_list)

def confidence_analysis(mean1,sigma1,T1,mean2,sigma2,T2,dt,degree = 3):
	
	'''
	conditional dependence of confidence
	:param mean1:
	:param sigma1:
	:param T1:
	:param mean2:
	:param sigma2:
	:param T2:
	:param dt:
	:param degree:
	:return: True or False
	'''
	
	sigma1_2 = sigma1*np.sqrt(dt/T2)
	sigma2_1 = sigma2*np.sqrt(dt/T1)
	return (mean1>=mean2-degree*sigma2_1)&(mean1<=mean2+degree*sigma2_1)&(mean2>=mean1-degree*sigma1_2)&(mean2<=mean1+degree*sigma1_2)


def get_bayesian_duration(data,sigma = 5):
	'''
	
	:param data: the result from background_correction()
	:param sigma:
	:return:
	'''
	start = 0
	stop = 0
	pulse = 1
	cafe = 2
	fringe_up = 3
	fringe_down = 4
	#t,rate = data['lc']
	edges = data['edges']
	re_rate,re_sigma,index_list = data['re_hist']
	dt = data['dt']
	bkg,bkgsigma,bkgsize = data['bkg']
	SNR = get_SNR(edges,re_rate,bkg,bkgsigma,dt,non_negative=True)
	
	binstart = edges[:-1]
	#binstop = edges[1:]
	trait = []
	for index1,hight in enumerate(re_rate):
		if (index1 == 0):
			trait.append(start)
		elif (index1 == len(re_rate)-1):
			trait.append(stop)
		else:
			if (hight> re_rate[index1-1]) and (hight > re_rate[index1+1]):
				trait.append(pulse)
			elif (hight <= re_rate[index1-1]) and (hight <= re_rate[index1+1]):
				trait.append(cafe)
			elif (hight < re_rate[index1-1]) and (hight >= re_rate[index1+1]):
				trait.append(fringe_down)
			elif (hight >= re_rate[index1-1]) and (hight < re_rate[index1+1]):
				trait.append(fringe_up)
			else:
				trait.append(cafe)
	start_tag = False
	start_edges = []
	stop_edges = []
	
	for index,value in enumerate(trait):

		if (SNR[index]>sigma) and (start_tag == False):
			
			if value != fringe_down:
				start_edges.append(binstart[index])
				start_tag = True
			else:
				stop_edges.pop()
				start_tag = True
				
			
		elif (SNR[index] <= sigma) and start_tag:
			
			if value != fringe_up:
				stop_edges.append(binstart[index])
				start_tag = False
	if start_tag:
		start_edges.pop()
	print(start_edges)
	if start_edges[0] == binstart[0]:
		start_edges = start_edges[1:]
		stop_edges = stop_edges[1:]
	return np.array(start_edges),np.array(stop_edges)

def get_bayesian_flash(data,start_edges,stop_edges):
	'''
	
	:param data: the result from background_correction()
	:param start_edgs: the result from get_bayesian_duration()
	:param stop_edges:
	:param sigma:
	:return:
	'''
	start = 0
	stop = 0
	pulse = 1
	cafe = 2
	fringe = 3
	edges = data['edges']
	edges_c = (edges[1:] + edges[:-1])*0.5
	re_rate,re_sigma,index_list = data['re_hist']
	flash_start = []
	flash_stop = []
	for i in range(len(start_edges)):
		#select blocks and edges between the start edge and stop edges of a burst.
		indexi = np.where((edges_c>=start_edges[i])&(edges_c<=stop_edges[i]))[0]
		
		flash_start.append(start_edges[i])
		if len(indexi)>=3:
			#only the size of block larger than 3 can the burst has the bayesian flash.
			re_ratei = re_rate[indexi]
			edgesi_index = np.where((edges >= start_edges[i]) & (edges <= stop_edges[i]))[0]
			edgesi = edges[edgesi_index]
			edg_star = edgesi[:-1]
			edg_stop = edgesi[1:]
			trait = []
			for index1,hight in enumerate(re_ratei):
				if (index1 == 0):
					trait.append(start)
				elif (index1 == len(re_ratei)-1):
					trait.append(stop)
				else:
					if (hight> re_ratei[index1-1]) and (hight > re_ratei[index1+1]):
						trait.append(pulse)
					elif (hight < re_ratei[index1-1]) and (hight < re_ratei[index1+1]):
						trait.append(cafe)
					else:
						trait.append(fringe)
			trait = np.array(trait)
			indexj = np.where(trait == cafe)[0]
			
			if len(indexj) >= 3:
				for ij in indexj:
					flash_start.append(edg_stop[ij])  #the cafe stop edge is the start edge of a flash.
					flash_stop.append(edg_star[ij])   #the cafe start edge is the stop edge of a flash.
		flash_stop.append(stop_edges[i])
	return np.array(flash_start),np.array(flash_stop)
	
	
	
		

