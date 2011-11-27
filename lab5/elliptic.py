#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : elliptic.py
### Author  : Oleg Baskakov
### Description : elliptic pde solver
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *

class Elliptic_PDE(PDE):
	approximate_boundary = '1lvl'
	
	u_x = 0.0
	u = 0.0
	coef_x = [+1, -2, +1]
	coef_y = [+1, -2, +1]
	
	def __init__(self, pde = None):
		super(Elliptic_PDE, self).__init__(pde)
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


def Liebmann_method(self):
		"""Just solve equation"""
		U = self.grid[-1]
		N = len(U)
		M = len(U[0])
		h = self.h
		
		U = [[1/4 * (U[i][j-1] + U[i-1][j] + U[i+1][j] + U[i][j+1] - h*h*f(i*h, j*h))
			  for j in range(1, M-1)]
			  for i in range(1, N-1)]
		
		#!solve boundary conditions here! 
		
		
		
		if self.coef_t[2] != 0:  # means that time has 2lvl
			U1 = self.grid[-2]
		else:
			U1 = [0] * N
		
		#TODO f(x,t) != 0

#		print(list(zip(* [self.coef_a, self.coef_b, self.coef_c])))		
#		print_vec(coefficients)
#		print([self.sigma, self.omega, self.eta])
		
		res = [0]*N
		res[1:-1] = self.threads_pool.map(PDE.explicit_fun,
		                        [(coeff,coef_t,U0,U1,i) for i in range(1, N-1)])



		(a0, b0, c0, d0) = self.first_eq(self.tau*N)
		(an, bn, cn, dn) = self.last_eq(self.tau*N)
		
		res[0]  = (d0 - c0*res[1])  / b0
		res[-1] = (dn - an*res[-2]) / bn
	
		return res
	
	
	
	



#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
