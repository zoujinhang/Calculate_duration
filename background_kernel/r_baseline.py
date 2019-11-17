import rpy2.robjects as robjects
from rpy2.robjects import r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import numpy as np

robjects.r("library(baseline)")
robjects.numpy2ri.activate()

def r_baseline(rate,dt,lam = None,hwi = None,it = None,inti = None):
	r.assign('rrate',rate)
	r('y = matrix(rrate,nrow = 1)')

	#fillpeak_int = str(int(time[-1]-time[0])/dt/10)
	if(lam is None):
		lam = int(0.708711*(len(rate))**0.28228+0.27114)
	if(hwi is None):
		hwi = int(40/dt)
	if(it is None):
		it = 10

	if(inti is None):
		fillpeak_int = str(int(len(rate)/10))
	else:
		fillpeak_int = str(inti)

	if(lam < 1):
		lam = 1
	r("rbase = baseline(y,lam ="+str(lam)+",hwi = "+str(hwi)+",it = "+str(it)+",int = "+fillpeak_int+",method = 'fillPeaks')")
	r("bs = getBaseline(rbase)")
	#r("cs = getCorrected(rbase)")
	bs = np.float32(r('bs'))[0]
	#cs = np.float32(r('cs'))[0]
	return bs
