

import matplotlib.pyplot as plt
import numpy as np


class Plot(plt):

	def __init__(self,result):
		self.result = result
		super(Plot, self).__init__()

	def plot_light_curve(self,**k):
		t = self.result['t_c']
		rate = self.result['rate']
		bs = self.result['bs']
		self.plot(t,rate,color = 'b',label = 'light curve',**k)
		self.plot(t,bs,color = 'r',label = 'background',**k)

		try:
			by_edges_list = self.result['bayesian_edges']
			by_rate_list = self.result['bayesian_rate']
			for index in range(len(by_edges_list)):

				self.plot(by_edges_list[index],by_rate_list[index],linestyle = 'steps',
					  color = 'k',label = 'bayesian block',**k)
		except:
			pass
	def plot_Txx1(self,txx,**k):

		if self.result['good']:

			t = self.result['t']
			dt = t[1] - t[0]
			n = self.result['n']/dt

			self.plot(t,n,color = 'b',label = 'light curve',**k)
			self.plot(t,self.result['bs']/dt,color = 'r',label = 'background',**k)

			t1_label = r'${T_{'+txx+',1} = ' + str(np.round(self.result['t1'] * 1000) / 1000) + '^{+' + str(
				np.round(self.result['t1_err'][0] * 1000) / 1000) + '}_{-' + str(
				np.round(self.result['t1_err'][1] * 1000) / 1000) + '}s}$'
			t2_label = r'${T_{'+txx+',2} = ' + str(np.round(self.result['t2']* 1000) / 1000) + '^{+' + str(
				np.round(self.result['t2_err'][0] * 1000) / 1000) + '}_{-' + str(
				np.round(self.result['t2_err'][1] * 1000) / 1000) + '}s}$'

			for i in self.result['bs_list']:
				self.plot(t,i,color = 'r',alpha = 0.01)
			self.axvline(x = self.result['t1'],color = 'g',linestyle = '--')
			self.axvline(x = self.result['t2'],color = 'g',linestyle = '--')
			self.plot(0,0,',',color = 'g',label = t1_label)
			self.plot(0,0,',',color = 'g',label = t2_label)
			self.xlim([t[0],t[-1]])
			self.xlabel('time (s)',**k)
			self.ylabel('rate',**k)
			self.legend()
		else:
			print('T'+txx+' is not good!')

	def plot_Txx2(self,txx,**k):

		if self.result['good']:

			t = self.result['t']
			cs_f = self.result['cs_f']
			cs_f_max = self.result['cs_f_max']
			l = self.result['l']
			self.axhline(y = 0,color = 'b')
			self.axhline(y = cs_f_max,color = 'b')
			self.axhline(y = l[0],color = 'g',linestyle = '--')
			self.axhline(y = l[1],color = 'g',linestyle = '--')
			self.plot(t,cs_f,color = 'k')
			label1 = r'${T_{'+txx+'} = ' + str(np.round(self.result['txx'] * 1000) / 1000) + '^{+' + str(
				np.round(self.result['txx_err'][0] * 1000) / 1000) + '}_{-' + str(
				np.round(self.result['txx_err'][1] * 1000) / 1000) + '}s}$'
			self.plot(0,0,',',label = label1)
			self.xlabel('time (s)',**k)
			self.ylabel('Accumulated counts',**k)
			self.xlim([t[0],t[-1]])
			self.legend()
		else:
			print('T'+txx+' is not good!')

	def plot_distribution(self,txx,**k):

		if self.result['good']:
			txx_list = self.result['txx_list']
			t1_list = self.result['t1_list']
			t2_list = self.result['t2_list']
			t90_list_sort = np.sort(txx_list)
			t90_err = np.std(txx_list)

			t90_bin = np.linspace(t90_list_sort[0], t90_list_sort[-1], 100)
			t90_n, t90_edges = np.histogram(txx_list, bins=t90_bin)
			t90_n = np.concatenate((t90_n[:1], t90_n))
			self.subplot(1,3,1)
			self.title('T'+txx+' distribution',**k)
			self.plot(t90_edges, t90_n, linestyle='steps', color='k')
			self.axvline(x=txx_list.mean(), color='r')
			self.axvline(x=txx_list.mean() - t90_err, color='r')
			self.axvline(x=txx_list.mean() + t90_err, color='r')
			self.axvline(x=self.result['txx'], color='g')
			self.xlabel('T'+txx+' (s)',**k)

			t1_list_sort = np.sort(t1_list)
			t90_bin = np.linspace(t1_list_sort[0], t1_list_sort[-1], 100)
			t90_n, t90_edges = np.histogram(t1_list, bins=t90_bin)
			t90_n = np.concatenate((t90_n[:1], t90_n))
			t1_err = np.std(t1_list)

			self.subplot(1,3,2)
			self.title('T'+txx+'1 distribution',**k)
			self.plot(t90_edges, t90_n, linestyle='steps', color='k')
			self.axvline(x=t1_list.mean(), color='r')
			self.axvline(x=t1_list.mean() - t1_err, color='r')
			self.axvline(x=t1_list.mean() + t1_err, color='r')
			self.axvline(x=self.result['t1'], color='g')
			self.xlabel('T'+txx+'1 (s)')

			t90_list_sort = np.sort(t2_list)
			t90_bin = np.linspace(t90_list_sort[0], t90_list_sort[-1], 100)
			t90_n, t90_edges = np.histogram(t2_list, bins=t90_bin)
			t90_n = np.concatenate((t90_n[:1], t90_n))
			t2_err = np.std(t2_list)

			self.subplot(1, 3, 3)
			self.title('T'+txx+'2 distribution')
			self.plot(t90_edges, t90_n, linestyle='steps', color='k')
			self.axvline(x=t2_list.mean(), color='r')
			self.axvline(x=t2_list.mean() - t2_err, color='r')
			self.axvline(x=t2_list.mean() + t2_err, color='r')
			self.axvline(x=self.result['t2'], color='g')
			self.xlabel('T'+txx+'2 (s)')

		else:
			print('T' + txx + ' is not good!')






























