#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce

from pprint import pprint

from pde import *


class Hyperbolic_PDE(PDE):
	def __init__():
		pde.__init__()
		self.sigma = tau**2 * a**2 / h**2
		self.omega = tau**2 * b / (2*h)
		self.eta = c * tau**2 + 2 # TODO add f(x,t)

		U = []
		psi0 = self.initial0
		psi1 = self.initial1
		
		U.append([ psi0(x) for x in frange(0, T, t)])
		U.append([ psi0(x) + psi1(x)*tau for x in frange(0, T, t)])
		for k in range(1, N-1):
			U[-1][k] += (U[0][k-1] - 2*U[0][k] + U[0][k+1]) * alpha**2 * tau**2 / (2 * h**2)
 		
		self.grid = U


	
	## Description: find coefficients of last equation
	def middle_eq(u):
		"""Find coefficients of middle equation"""
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
	

		ai = sum(koef_a)
		bi = sum(koef_b)
		ci = sum(koef_c)
		di = -u[N][k-1]

		return (ai, bi, ci, di)
		
	def last_eq(u):
		"""Find coefficients of last equation"""
		pass

		
	def solve(self):
		U = self.grid
		print("ololo")

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
