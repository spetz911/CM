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

class Parabolic_2d_PDE(PDE):
	"""ADI_solver"""
	approximate_init = '1lvl'
	approximate_boundary = '1lvl'
	method = 'adi'
	
	u_x = 0.0
	u = 0.0
	coef_t = [+1, -1, +0]
	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	
	def __init__(self, pde = None):
		super(Parabolic_2d_PDE, self).__init__(pde)
		tmp = self.grid
		self.grid = None
		MetaClass.print(self)
		self.grid = tmp



	def initial_cond_lvl0(self):
		"""qqq"""
		psi0 = self.initial0
		self.grid = []
		self.grid.append([[psi0(x,y) for x in frange(0, self.lx, self.h)]
		                             for y in frange(0, self.ly, self.h)])

	
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


	def solve_dim(self, field, t):
		"""Solver for one dimension"""
		res_x = []
		
		N = len(field)
		M = len(field[0])
		h = self.h
		
		teta = 1.0
		coef_t = self.coef_t
		
		print("TEST:")
		print_mat(field)
		
		res_x.append(field[0])
		for i in range(1, N-1):
			Eq = []
			Eq.append(self.__first_eq(i*h, t))
			
			for j in range(1, M-1):
				U = field[j]
				# eta = c*u + f(x,t)
				tmp = (field[i-1][j],field[i][j],field[i+1][j])
				eta = PDE.scalar(self.coef_a, tmp) + self.tau * self.c
				# TODO add f(x,t)
				coeff = PDE.vec_mat([self.sigma, self.omega, self.eta],
                                     [self.coef_a, self.coef_b, self.coef_c])
				Eq.append((teta * self.coeff[0],
			        teta * self.coeff[1] - coef_t[0],
			        teta * self.coeff[2],
			        coef_t[1] * U[j] ))
			
			Eq.append(self.__last_eq(i*h, t))
		
			Matrix = Tridiagonal_Matrix()
			[Matrix.a, Matrix.b, Matrix.c, Matrix.d] = list(zip(* Eq))
			Matrix.n = M
	
			#	print_mat(Eq)
			x = Matrix.solve()
			res_x.append(x)

		res_x.append(field[-1])
			
		return res_x
#			res_x.append(implicit_method())
		

##====================================================================================
	def solve_adi(self, t):
		"""Alternating direction implicit method"""

		print("true fuction!")
		
		#need division here?
		
		

		field = self.grid[-1]

		self.left = self.west
		self.right = self.east
		new_field_x = self.solve_dim(field,t+self.tau)

		print("after")
		print_mat(new_field_x)


		#transpose
		field = [list(x) for x in zip(*new_field_x)]

		self.left = self.south
		self.right = self.north
		new_field_y = self.solve_dim(field, t + 2*self.tau)
		
		return [list(x) for x in zip(*new_field_y)]


	def __first_eq(self, z, t):
		"""Find coefficients of first equation"""
		alpha = self.left[0]/self.tau
		beta  = self.left[1]
		phi0  = self.left[2]
		h = self.h
		
		a0 = 0
		b0 = beta - alpha / h
		c0 = alpha / h
		d0 = phi0(z,t)
		print("calltrace")
		print((a0, b0, c0, d0))
		
		return (a0, b0, c0, d0)

	def __last_eq(self, z, t):
		"""Find coefficients of first equation"""
		alpha = self.right[0]/self.tau
		beta  = self.right[1]
		phi1  = self.right[2]
		h = self.h
		
		a0 = alpha / h
		b0 = beta - alpha / h
		c0 = 0
		d0 = phi1(z,t)
		return (a0, b0, c0, d0)


	def solve_fs(self, t):
		"""Fractional step method"""
		#need division here?
		h = self.h
		
		field_x = self.grid[-1]
		w = self.west
		e = self.east
		n = self.north
		s = self.south
		tau = self.tau
		h = self.h
		
		# fun(x,y,t)???
#		print("before")
#		print_mat(field_x)
		
		new_field = []
		for i, zz in enumerate(field_x):
			self.left  = (w[0]/tau, w[1], lambda t: w[2](i*h, t))
			self.right = (e[0], e[1], lambda t: e[2](i*h, t))
			if (i!=0) and (i!=len(field_x)-1):
				new_field.append(self.implicit_method(U = zz, t = t))
			else:
				new_field.append(zz)

#		print("after")
#		print_mat(new_field)
			
		field_y = [list(x) for x in zip(*new_field)]


		new_field = []
		for i, zz in enumerate(field_y):
			self.left  = (s[0]/tau, s[1], lambda t: s[2](i*h, t+tau))
			self.right = (n[0], n[1], lambda t: n[2](i*h, t+tau))
			
			if (i!=0) and (i!=len(field_y)-1):
				new_field.append(self.implicit_method(U = zz, t = t+tau))
			else:
				new_field.append(zz)
		
		res = [list(x) for x in zip(*new_field)]
#		print("after-after")
#		print_mat(res)
		
		return res

	def solve(self, method = 'adi'):
		"""Description"""
		self.M = int(self.lx / self.h)
		self.N = int(self.ly / self.h)
		self.max = 0
		tau = self.tau
		
	#	self.grid = [[0.0]*self.M for i in range(self.N)]
		Us = self.grid

		if method == 'adi':
			one_iteration = self.solve_adi
		elif method == 'fs':
			one_iteration = self.solve_fs
		else:
			print("method %s not supported" % method)
			return
		
		Us.append(one_iteration(0))
		count = 2
		while count * tau > self.t and count < self.max:
			Us.append(one_iteration((count-1) * tau))
			print_mat(Us[-1])
			count += 2
		self.count = count
		print("complete")
		return Us




#=====================================================================

def main():
	pass

if __name__ == "__main__":
	main()
