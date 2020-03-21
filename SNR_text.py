

from .autocorrelation_text import *
from .background_kernel import WhittakerSmooth

def SNR_text(t,v,step_size = 1,block_n = 50,
	     block_time = None,time_unified = True,lambda_ = 200):

	t = np.array(t)
	v = np.array(v)
	#dt = t[1]-t[0]
	background_index,nornallization,block_index,ACC,ACCT = autocorrelation_text(t,v,step_size=step_size,
									       block_n = block_n,
									       block_time = block_time,para = True,
									       time_unified=time_unified)
	w = np.zeros(v.size)
	w[background_index] = 1
	bs = WhittakerSmooth(v,w,lambda_= lambda_)
	cs = v - bs
	sigma = cs[background_index].std()

	index = np.where(nornallization >= 0.95)[0]
	good_background_index = np.where(nornallization < 0.95)[0]
	#print('good_index:',index)
	nsi = (cs/sigma)
	result = {'nsi':nsi,'sigma':sigma,'background_index':background_index,
		  'good_index':piecemeal(index),'bs':bs,'cs':cs,'good_background_index':good_background_index,
		  'normallization':nornallization,
		  'ACCT':ACCT,
		  'ACC':ACC
		  }

	return result


def piecemeal(inputindex):

	blocki = []
	result = []
	N = len(inputindex)
	value_old = 0

	for index,value in enumerate(inputindex):
		#print(blocki)
		if index == 0:

			value_old = value
			blocki.append(value)

		if(index == N-1):

			result.append(blocki)

		elif(value-value_old == 1):

			value_old = value
			blocki.append(value)
		elif(value-value_old >1):

			value_old = value
			result.append(blocki)
			blocki = [value]

	return np.array(result)














