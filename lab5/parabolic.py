#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce

from pprint import pprint

from pde import *

class Parabolic_PDE(PDE):
	u_x = 0.0
	u = 0.0
	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	
	def __init__(self, pde = None):
		super(Parabolic_PDE, self).__init__(pde)
		
		MetaClass.print(self)

		a = sqrt(self.u_xx)
		b = self.u_x
		c = self.u
		h = self.h
		l = self.l

		tau = self.tau
		
		self.sigma = tau * a**2 / h**2
		self.omega = tau * b / (2*h)
		self.eta = tau * c # TODO add f(x,t)
		
		psi0 = self.initial0
		U = []
		U.append([ psi0(x) for x in frange(0, l, h)])
		self.grid = U
		
	
	
	def first_eq(self, t):
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
	
	def last_eq(self, t):
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

	def middle_eq(self, i, teta = 1):
		"""Find coefficients of middle equation"""
		U = self.grid[-1]
		N = len(U)
		
		coeff = PDE.vec_mat([self.sigma, self.omega, self.eta],
		                    [self.coef_a, self.coef_b, self.coef_c])
	
		
		ai = teta * coeff[0]
		bi = teta * coeff[1] - 1
		ci = teta * coeff[2]
		di = -U[i]
		
		

		return (ai, bi, ci, di)
	
	def solve(self):
		Us = self.grid
		
#		for t in frange(0, self.t, self.tau):
#			U.append(self.explicit_method())
	
		for t in frange(0, self.t, self.tau):
			Us.append(self.implicit_method())
		
		
		print("ololo")
		return Us


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
