

from .autocorrelation_text import *
from .background_kernel import WhittakerSmooth

def SNR_text(t,v,criterion = 1,step_size = 1,block_n = 50,block_time = None):

	t = np.array(t)
	v = np.array(v)
	#dt = t[1]-t[0]
	background_index,nornallization,block_index,ACC = autocorrelation_text(t,v,step_size=step_size,block_n = block_n,block_time = block_time)

	w = np.zeros(v.size)
	w[background_index] = 1
	bs = WhittakerSmooth(v,w,lambda_= 200)
	cs = v - bs
	sigma = cs[background_index].std()
	index = np.where((cs > criterion*sigma)&(nornallization >=  1))[0]
	good_background_index = np.where((cs <= criterion*sigma) | (nornallization < 1))[0]
	nsi = (cs/sigma)
	result = {'nsi':nsi,'sigma':sigma,'background_index':background_index,
		  'good_index':piecemeal(index),'bs':bs,'cs':cs,'good_background_index':good_background_index}

	return result


def piecemeal(inputindex):

	blocki = []
	result = []
	N = len(inputindex)
	value_old = 0

	for index,value in enumerate(inputindex):

		if index == 0:
			value_old = value
			blocki.append(value)
		if(index == N-1):
			result.append(blocki)

		elif(value-value_old == 1):
			blocki.append(value)
		elif(value-value_old >1):
			result.append(blocki)
			blocki = [value]

	return np.array(result)














