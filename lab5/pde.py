#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : pde.py
### Author  : Oleg Baskakov
### Description : Partial Differencial Equation
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------


##TODO Successive Over Relaxation!!!

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *
# from parabolic import *
# from hyperbolic import *

from multiprocessing import Pool

from pprint import pprint



eps = 0.0001

def print_mat(mat):
	print("matrix:")
	for row in mat:
		for elem in row:
			print(round(elem, 3), end = "\t")
		print()
	print("--------")

def print_vec(row, k = 3):
#	print("vector:")
	for elem in row:
		print(round(elem, k), end = "\t")
	print()
#	print("--------")




class MetaClass:
	def copy(self, src):
		[self.__setattr__(k,v) for k,v in src.__dict__.items()]
	
	def print(self):
		from pprint import pprint
		pprint(dict(self.__dict__))

eps = 0.00001
def frange(x0, x1, d):
	while (x0 < x1 + eps):
		yield x0
		x0 += d



class PDE:
	def __init__(self, pde = None):
		print("PDE constructor", pde)
		if pde:
			MetaClass.copy(self, pde)
		else:
			pass
		self.threads_pool = Pool(processes=4)
		self.set_equation_params()
		self.initial_cond_lvl0()

		# means that time has tau**2
		if self.coef_t[2] != 0:
			self.initial_cond_lvl1()
		
		if self.approximate_init == '2lvl' and self.coef_t[2] != 0:
			self.initial_cond_lvl2()
		
			
		self.coeff = PDE.vec_mat([self.sigma, self.omega, self.eta],
                                 [self.coef_a, self.coef_b, self.coef_c])

	def set_equation_params(self):
		self.a = sqrt(self.u_xx)
		self.b = self.u_x
		self.c = self.u
		
		self.sigma = self.tau * self.a**2 / self.h**2
		self.omega = self.tau * self.b / (2*self.h)
		self.eta = self.tau * self.c # TODO add f(x,t)

	def initial_cond_lvl0(self):
		psi0 = self.initial0
		self.grid = []
		self.grid.append([ psi0(x) for x in frange(0, self.l, self.h)])

	def initial_cond_lvl1(self):
		self.sigma *= self.tau
		self.omega *= self.tau
		self.eta   *= self.tau

		psi0 = self.initial0
		psi1 = self.initial1
		U.append([ psi0(x) + psi1(x)*tau for x in frange(0, self.l, self.h)])
		

	
	def initial_cond_lvl2(self):
		U0 = self.grid[0]
		U1 = self.grid[1]
		N = len(U0)

		U1[1:N-1] = [U[k] + PDE.scalar(coef_t, U0[k-1 : k+2]) * a**2 * tau**2 / (2*h**2)
		             for k in range(1, N-1)]
	

	@staticmethod
	def scalar(v1, v2):
		return sum(x1*x2 for (x1,x2) in zip(v1,v2))

	@staticmethod
	def mat_vec(m, v):
		return [PDE.scalar(vx, v) for vx in m]

	@staticmethod
	def vec_mat(v, m):
		return PDE.mat_vec(zip(*m), v)

	def first_eq(self, t):
		if self.approximate_boundary == '1lvl':
			return self.first_eq_1lvl(t)
		if self.approximate_boundary == '2lvl':
			return self.first_eq_2lvl(t)
		
	def last_eq(self, t):
		if self.approximate_boundary == '1lvl':
			return self.last_eq_1lvl(t)
		if self.approximate_boundary == '2lvl':
			return self.last_eq_2lvl(t)

	def first_eq_1lvl(self, t):
		"""Find coefficients of first equation"""
		alpha = self.left[0]
		beta  = self.left[1]
		phi0  = self.left[2]
		h = self.h
		
		a0 = 0
		b0 = beta - alpha / h
		c0 = alpha / h
		d0 = phi0(t)
		return (a0, b0, c0, d0)

	def last_eq_1lvl(self, t):
		"""Find coefficients of first equation"""
		alpha = self.right[0]
		beta  = self.right[1]
		phi1  = self.right[2]
		h = self.h
		
		a0 = -alpha / h
		b0 = beta - alpha / h
		c0 = 0
		d0 = phi1(t)
		return (a0, b0, c0, d0)



# if __name__ == '__main__':
#    pool = Pool(processes=4)              # start 4 worker processes
#    result = pool.apply_async(f, [10])     # evaluate "f(10)" asynchronously
#    print(result.get(timeout=1))           # prints "100" unless your computer is very slow
#    print())          # prints "[0, 1, 4,..., 81]"



	@staticmethod
	def explicit_fun(coeff, coef_t, U0, U1 ,i): #TODO add threads
		res  = PDE.scalar(coeff, U0[i-1:i+2])
		res -= coef_t[1] * U0[i]
		res -= coef_t[2] * U1[i]
		return res
	
	def explicit_method(self): #TODO add threads
		"""Just solve equation"""
		U0 = self.grid[-1]
		N = len(U0)
		if self.coef_t[2] != 0:  # means that time has 2lvl
			U1 = self.grid[-2]
		else:
			U1 = [0] * N
		
		# TODO f(x,t) != 0

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

	def implicit_method(self, teta = 1):
		"""Solve with method Progonki"""
		U = self.grid[-1]
		N = len(U)
		t = self.tau*N
		coef_t = self.coef_t
		if coef_t[2] != 0:  # means that time has 2lvl
			U1 = self.grid[-2]
		else:
			U1 = [0] * N
		
		
		M = Tridiagonal_Matrix()
		
		Eq = []
		Eq.append(self.first_eq(t))
		Eq.extend([(teta * self.coeff[0],
			        teta * self.coeff[1] - coef_t[0],
			        teta * self.coeff[2],
			        coef_t[1] * U[i] + coef_t[2] * U1[i])
			                       for i in range(1, N-1)])
		Eq.append(self.last_eq(t))
		
		[M.a,M.b,M.c,M.d] = list(zip(* Eq))
		M.n = N
	
	#	print_mat(Eq)

		x = M.solve()
		return x

	def Crank_Nicolson_method(self, teta = 0.5):
		""" Uk+1 - Uk = teta * implicit + (1 - teta) * explicit """
		U = self.grid[-1]
		N = len(self.grid[-1])
		t = self.tau*N
		Eq = []
		Eq.append(self.first_eq(t))
		for i in range(1, N-1):
			Expl = (1 - teta) * PDE.scalar(self.coefficients, U[i-1:i+2])
			Impl = self.middle_eq(i,teta)
			eq = Impl[0], Impl[1], Impl[2], Impl[3] - Expl
			Eq.append(eq)
		Eq.append(self.last_eq(t))

		M = Tridiagonal_Matrix()
		
		[M.a,M.b,M.c,M.d] = list(zip(* Eq))
		M.n = N
		
#		print_mat(Eq)

		x = M.solve()
		return x

	def solve(self, method = 'Crank_Nicolson'):
		"""Description"""
		Us = self.grid
		if method == 'explicit':
			for t in frange(0, self.t, self.tau):
				Us.append(self.explicit_method())
		elif method == 'implicit':
			for t in frange(0, self.t, self.tau):
				Us.append(self.implicit_method())
		elif method == 'Crank_Nicolson':
			for t in frange(0, self.t, self.tau):
				Us.append(self.Crank_Nicolson_method(0.5))
		
		print("complete")
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

