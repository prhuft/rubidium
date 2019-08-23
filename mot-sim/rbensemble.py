"""
RbEnsemble
Preston Huft, Summer 2019

A class for generating Rb atom distributions to be passed to various simulations.
"""

## libraries
from matplotlib import pyplot as plt
import numpy as np
from numpy import linspace,sin,cos,log,exp
import math as m
from math import sqrt,pi,e,exp
from random import random as rand

## local files
from physconsts import *
from rbconsts import *

## classes 

class RbEnsemble:
	""" Ensemble of Rb atoms with a momentum distribution specified
		by temperature. Optional spatial distribution.
		
		'T'= temperature
		'size'= number of atoms
		'xdist': optional parameter specifying the initial
			position distribution
	"""
	global mRb
	
	
	def __init__(self,T,size=None,xdist=None,statedist=None):
		
		# For efficiency, pre-generate a specified number of atoms
		if size is not None:
			self.size = size
			self.temp = T
			self.v = self.sampling_maxboltzv(self.size,[0,1],self.temp) # rms
			self.p = mRb*self.v # rms
			self.x = np.empty(self.size)
			if xdist is None:
				self.x = np.zeros(self.size)
			elif xdist is 'normal':
				self.x = np.random.normal(0,size=self.size)
			if statedist is not None:
				self.amplitudes = self.psi_coeffs(self.size,statedist)
		else:
			self.size = 0
			self.temp = T
			self.v = np.array([]) # rms
			self.p = np.array([]) # rms
			self.x = np.array([])
			self.amplitudes = np.array([])
			
	def phasespace(self):
		""" Plots the ensemble in phase space. 1D x and p only for 
			now.
		"""
		xmax = max(self.x) # like xmas but better
		xmin = min(self.x) # because i said so
		dx = xmax-xmin
		
		pmax = max(self.p)/mRb
		pmin = min(self.p)/mRb
		dp = pmax-pmin
		
		fig, ax = plt.subplots()
		ax.scatter(self.p/mRb,self.x)#,linestyle=None)
		ax.set(xlabel='p [m/(s mRb)]', ylabel='r [arb]',
			   xlim=(pmin-.1*dp,pmax+.1*dp),
			   ylim=(xmin-.1*dx,xmax*+.1*dx))
		plt.show()
		
	def vpt(self):
		""" Return a speed from Maxwell-Boltzmann dist. """
		return sampling_maxboltzv(1,[0,1],self.temp)

	def xpt(self,domain):
		""" Return a position from a flat dist by default. """
		
		x1,x2 = domain
		x = rand()*(x2-x1) # only works for x1,x2 > 0
		return x

	def maxboltzv(self,T,v,normalization=False):
		""" Maxwell-Boltzmann distribution of speeds for 3-dimensional
			gas. Returns f(v) for T. """
		global kB,mRb
		m = mRb

		A = 4*pi*(m/(2*pi*kB*T))**(3/2) # normalization consts
		meanv = sqrt(2*kB*T/m) # the maximum occurs at the mean

		if normalization is True:
			return A
		else:
			return A*v**2*exp(-m*v**2/(2*kB*T))

	def sampling_maxboltzv(self,size,domain,T,vectorial=False,showplot=False):
		""" Sample random speeds with a Maxwell-Boltzmann dist. 
			'size': sample size
			'domain': [v1,v2] the restricted domain of the pdf; e.g.
				a Guassian goes to zero well before inf so we could
				let the domain be a finite region
			'T': temperature
			'vectorial': 
				If False, only return a scalar. 
				Set to True to return velocity vectors with a 
				direction from a flat distribution. 
		"""
		global kB,mRb
		m = mRb

		n = size 
		v1,v2 = domain

		mean = sqrt(2*kB*T/m)
		fmax = self.maxboltzv(T,mean) # the maximum
		y_dist = np.empty(n) 
		f_dist = np.empty(n) 
		v_dist = np.empty(n) # this is the distribution we want
		j = 0 # dist index
		while j < n:
			v = (v2-v1)*rand() # rand val on domain of f(x)
			f = self.maxboltzv(T,v)
			y = rand()*fmax # rand val on range of f(x)
			if y <= f:
				y_dist[j]=y
				f_dist[j]=f
				v_dist[j]=v # x vals with approximate gaussian pdf
				j+=1

		# plot distribution as a check:
		if showplot is not False:
			plt.scatter(v_dist,y_dist,c='red',s=10)
			plt.scatter(v_dist,f_dist,c='blue',s=10)
			plt.show()

		return v_dist
		
	def psi_coeffs(self,size,statedist):
		""" For a state dist {'c_e': w_e,'c_g': w_g } where w_i is the weight in (0,1) 
			for the ith state coeff, return c_e,c_g.
			
			Really, this should be giving rho elements I think. 
		"""
		s = size
		coeffs = np.empty(s,object)
		# c_list = np.zeros(len(statedist))
		# for i,c in enumerate(statedist)[:-1]:
			# x = rand()
			# if x < p:
				# c_list[i] = 1
		
		# probably a better way to initialize the states... maybe like randomly
		# filling a density matrix? It's a littel weird to fill rho for one atom.
		
		# for now, just spit out coeffs for a two-level atom:
		for i in range(s):
			ce = cg = 0
			if rand()>.5:
				ce = 1
			else:
				cg = 1
				
			coeffs[i] = [cg,ce]
	
		return coeffs
		
	# def state_plot(state):
		# TODO:
		# set aspect ratio = 1 for the phase space plot,
		# draw a gaussian waist fit on the the p.s. plot
		# add this state plot function to the rb ensemble class
		# """plot the ensemble in phase space and a level histogram"""
		# x,p,e = state
		# levels = [0,1] # the energy levels. not numpy because i'm lazy
		
		# phase space
		# plt.subplot(221)
		# plt.scatter(x,p,s=1)

		# level population
		# plt.subplot(222)
		# plt.hist(e,bins=2,histtype='step')
		# plt.xticks(levels)

		# plt.show()
			
	
	
	