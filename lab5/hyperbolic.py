#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : hyperbolic.py
### Author  : Oleg Baskakov
### Description : hyperbolic pde solver
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *


class Hyperbolic_PDE(PDE):
	approximate_init = '1lvl'
	approximate_boundary = '1lvl'
	method = 'implicit'
	
	u_x = 0.0
	u = 0.0

	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	coef_t = [+1, -2, +1]
	
	def __init__(self, pde = None):
		super(Hyperbolic_PDE, self).__init__(pde)
		MetaClass.print(self)


	def first_eq_2lvl(self, t):
		"""Find coefficients of first equation"""
		alpha = self.left[0]
		beta  = self.left[1]
		phi0  = self.left[2]
		U0 = self.grid[-1]
		U1 = self.grid[-2]
		
		a = self.u_xx #XXX sqrt here
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau
		
		a0 = 0.0
		b0 = alpha*(-2*a / (h*(2*a-b*h)) - h/(2*a-b*h) + c*h / (2*a-b*h)) + beta;
		c0 = alpha*(2*a) / (h*(2*a-b*h))
		d0 = phi0(t) + alpha * h * (U0[0]-2*U1[0]) / (2*a-b*h) / tau**2;

		return (a0, b0, c0, d0)
	
	def last_eq_2lvl(self, t):
		"""Find coefficients of last equation"""
		alpha = self.right[0]
		beta  = self.right[1]
		phi1  = self.right[2]
		U0 = self.grid[-1]
		U1 = self.grid[-2]

		a = self.u_xx #XXX sqrt here
		b = self.u_x
		c = self.u
		h = self.h
		tau = self.tau

		an = -2*alpha*a / (h*(2*a+b*h))
		bn = alpha*(2*a+h/tau**2-c*h) / (2*a+b*h) + beta
		cn = 0.0
		dn = phi1(t) - alpha*h*(U1[-1] - 2*U0[-1]) / (2*a+b*h) / tau**2
		
		return (an, bn, cn, dn)

#=====================================================================

def main():

	f = open("input")
	
	pde = parse_file(f)

#	print(pde.left[2](0))
#	print(pde.right[2](0))
	
	res = pde.solve()
	
	for k in [0,1,2,3]:
		print("iter ", k)
		print_vec(res[k])
		print_vec([pde.fun(x, k*pde.tau) for x in frange(0, pde.l, pde.h)])
	
#	print(pde.u_x)


#=====================================================================

if __name__ == "__main__":
	main()
