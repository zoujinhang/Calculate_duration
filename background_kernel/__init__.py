'''
背景相关工具


'''

from .WhittakerSmooth import WhittakerSmooth
import numpy as np
from ..time_unified import Time_transform
from .r_baseline import r_baseline


class Baseline_in_time(object):

	def __init__(self,time,value,case='FD',hardness = False,time_unified = True):
		self.time = time
		self.value = value


		self.t_transform = Time_transform(time,value)
		self.unified_time,self.unified_value = self.t_transform.to_unified_time()
		if time_unified:
			self.AirPLS = AirPLS(self.unified_value,hardness = hardness)
			dt = self.unified_time[1] - self.unified_time[0]
			print('dt',dt)
		else:
			self.AirPLS = AirPLS(self.value, hardness=hardness)
			dt = self.time[1]-self.time[0]
			print('dt',dt)

		self.unified_bs = self.AirPLS.bottom_r(dt = dt,case = case)
		if time_unified:
			self.bs = self.t_transform.to_actual_time(self.unified_time,self.unified_bs)[1]
		else:
			self.bs = self.unified_bs
		self.cs = self.value-self.bs

	def get_value(self):
		return self.time,self.cs,self.bs


	def get_bs(self):
		return self.bs

	def get_cs(self):
		return self.cs

	def get_airPLS(self):
		return self.AirPLS

class AirPLS(object):

	def __init__(self,x,hardness = False):
		self.b_index = np.isnan(x)
		self.x = np.array(x)

		self.m = self.x.shape[0]

		self.w = np.ones(self.m)
		if (True in self.b_index):
			print('数据输入存在存在无效值')
			self.x[self.b_index] = 0
			self.w[self.b_index] = 0

		if hardness:
			self.x = WhittakerSmooth(self.x,self.w,4)


	def bottom_r(self,dt,case = 'FD'):

		return r_baseline(self.x,dt,case = case)

























