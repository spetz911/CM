#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : parabolic.py
### Author  : Oleg Baskakov
### Description : parabolic pde solver
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *

class Parabolic_PDE(PDE):
	approximate_init = '1lvl'
	approximate_boundary = '1lvl'
	
	u_x = 0.0
	u = 0.0
	coef_t = [+1, -1, +0]
	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	
	def __init__(self, pde = None):
		super(Parabolic_PDE, self).__init__(pde)
		MetaClass.print(self)

	
	def first_eq_2lvl(self, t):
		"""Find coefficients of first equation"""
		alpha = self.left[0]
		beta  = self.left[1]
		phi0  = self.left[2]
		U = self.grid[-1]
		
		a = self.u_xx
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau
		
		a0 = 0
		b0 = alpha * (2*a*a/h + h/tau - c*h) - beta * (2*a*a - b*h)
		c0 = alpha * (-2*a*a/h)
		d0 = alpha * (U[0] * h/tau) - phi0(t) * (2*a*a - b*h)
		return (a0, b0, c0, d0)
	
	def last_eq_2lvl(self, t):
		"""Find coefficients of last equation"""
		alpha = self.right[0]
		beta  = self.right[1]
		phi1  = self.right[2]
		U = self.grid[-1]
		N = len(U)

		a = self.u_xx
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau

		an = alpha * (-2*a*a/h)
		bn = alpha * (2*a*a/h + h/tau - c*h) + beta * (2*a*a + b*h)
		cn = 0
		dn = alpha * (U[-1] * h/tau) + phi1(t) * (2*a*a + b*h)
		return (an, bn, cn, dn)

	def check_scheme(self):
		N = self.N
		h = self.h
		tau = self.tau
		u = self.orig
		
		for j in range(0, T):
			fu = [u(i*h, tau*j) for i in range(1, N-1)]
			
		
		max(   for i in  )
		
		



#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
