#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : adi.py
### Author  : Oleg Baskakov
### Description : Alternating direction implicit method
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *

class ADI_solver(PDE):
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


	def solve_dim():
		"""Solver for one dimension"""
		res_x = []
		
		N = len(field)
		M = len(field[0])
		
		for i in range(1, N-1):
			Eq = []
			Eq.append(self.first_eq(t))
			
			for j in range(1, M-1):
				
				# eta = c*u + f(x,t)
				eta = scalar(self.coef_a, field[i-1:i+2][j]) + self.tau * self.c # TODO add f(x,t)
				coeff = PDE.vec_mat([self.sigma, self.omega, self.eta],
                                     [self.coef_a, self.coef_b, self.coef_c])
				Eq.append((teta * self.coeff[0],
			        teta * self.coeff[1] - coef_t[0],
			        teta * self.coeff[2],
			        coef_t[1] * U[i] + coef_t[2] * U1[i]))
			
			Eq.append(self.last_eq(t))
		
			[M.a,M.b,M.c,M.d] = list(zip(* Eq))
			M.n = N
	
			#	print_mat(Eq)
			x = M.solve()
			res_x.append(x)
			
			res_x.append(implicit_method())
		

##====================================================================================
	def solve(self, method = 'ADI'):
		"""Alternating direction implicit method"""
		self.tau /= 2
		#need division here?
		for t in frange(0, self.t, self.tau):

			field = self.field[-1]

			new_field_x = solve_dim(field)
		
			field = [list(x) for x in zip(*field)]
		
			new_field_y = solve_dim(field)
		
		self.tau *= 2
		return new_field_y


	def solve(self, method = 'FS'):
		"""fractional step method"""
		self.tau /= 2
		#need division here?
		
		for t in frange(0, self.t, self.tau):
			field_x = self.field[-1]
			new_field = []
			for self.grid in field:
				new_field.append(self.implicit_method())
			field_y = zip(*new_field)
			new_field = []
			for self.grid in field:
				new_field.append(self.implicit_method())
			self.field.append(new_field)
		
		self.tau *= 2
		
		return self.field



#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
