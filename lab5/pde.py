#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *
# from parabolic import *
# from hyperbolic import *

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

class PDE_parser():
	def __init__(self):
		pass
	
	## Description: for parsing one element of summ
	@staticmethod
	def get_token(x):
		tmp = x.split('*')
		if tmp[-1] in ['u', 'u_xx', 'u_x', 'u_yy', 'u_y', 'u_tt', 'u_t']:
			try:
				a = float(tmp[0])
			except:
				a = 1.0
			return (tmp[-1], a)
		else:
			return None

	def parse_pde(self, s):
		"""For parsing partial differential equation"""
		[left, right] = s.lower().rstrip().replace(' ','').split('=')
		self.type = None
		if ('u_xx' in right) and ('u_t' in left):
			self.type = 'parabolic'
		elif ('u_xx' in right) and ('u_tt' in left):
			self.type = 'hyperbolic'
		elif ('u_xx' in right) and ('u_yy' in right):
			self.type = 'elliptic'

		tokens = [ x for x in left.split('+') + right.split('+')]
		rest = []
		for x in tokens:
			t = PDE_parser.get_token(x)
			if t:
				self.__setattr__(t[0], t[1])
			else:
				rest.append(x)
		try: self.fun = eval("lambda x,t: " + "+".join(rest))
		except: self.fun = lambda x,t: 0.0

	def parse_stuff(self, line):
		"""For parsing other coefficients"""
		[left, right] = line.lower().replace(' ','').split('=')
		if left in ['l', 'lx', 'ly', 't', 'w', 'eps', 'tau', 'h']:
			self.__setattr__(left, float(right))
		elif left == 'result':
			try: self.__setattr__('fun', eval("lambda x,t: " + right))
			except: self.fun = lambda x,t: 0.0

	## Description: for solve init eq
	def parse_initial_condition(self, line):
		[left, right] = line.lower().replace(' ','').split('=')
		if left in ['u(x,0)', 'u(x,y,0)', 'u(x,y)']:
			try: self.initial0 = eval("lambda x,y=0: " + right)
			except: self.initial0 = lambda x,y=0: 0.0
		if left in ['u_t(x,0)', 'u_t(x,y,0)', 'u_t(x,y)']:
			try: self.initial1 = eval("lambda x,y=0: " + right)
			except: self.initial1 = lambda x,y=0: 0.0

	## Description: for solve limit eq
	def parse_boundary_condition(self, line):
		[left, right] = line.lower().replace(' ','').split('=')
		#! alpha*u_x + beta*u = f(t)
		alpha = 0.0
		beta = 0.0
		
		for elem in left.split('+'):
			tmp = elem.split('*')
			if ('u_x' in tmp[-1] or 'u_y' in tmp[-1]):
				try: alpha = float(tmp[0])
				except: alpha = 1.0
			elif ('u' in tmp[-1]):
				try: beta = float(tmp[0])
				except: beta = 1.0
		try: result = (alpha, beta, eval("lambda t: " + right))
		except: result = (alpha, beta, lambda t: 0.0)

		if   ('(x,ly' in left):
			self.north = result
		elif ('(0,x' in left):
			self.south = result
		elif ('(l' in left): # (lx,y,0)
			self.east = result
		elif ('(0' in left): # (0,y,0)
			self.west = result
		
		try: self.left = self.west
		except: pass
		try: self.right = self.east
		except: pass



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


	@staticmethod
	def scalar(v1, v2):
		return sum(x1*x2 for (x1,x2) in zip(v1,v2))

	@staticmethod
	def mat_vec(m, v):
		return [PDE.scalar(vx, v) for vx in m]

	@staticmethod
	def vec_mat(v, m):
		return PDE.mat_vec(zip(*m), v)

	
	def explicit_method(self): #TODO add threads
		"""Just solve equation"""
		U = self.grid[-1]
		N = len(U)
		

		# TODO f(x,t) != 0
#		print(list(zip(* [self.coef_a, self.coef_b, self.coef_c])))
		
		
#		print_vec(coefficients)
#		print([self.sigma, self.omega, self.eta])
		
		res = [0]+[U[i]+PDE.scalar(self.coefficients, U[i-1:i+2])
		                                     for i in range(1, N-1)]+[0]


		(a0, b0, c0, d0) = self.first_eq(self.tau*N)
		(an, bn, cn, dn) = self.last_eq(self.tau*N)
#		print(self.first_eq(self.tau*N))
#		print(self.last_eq(self.tau*N))
		
		res[0]  = (d0 - c0*res[1])  / b0
		res[-1] = (dn - an*res[-2]) / bn
	
		return res

	def implicit_method(self):
		"""Solve with method Progonki"""
		U = self.grid[-1]
		N = len(U)
		t = self.tau*N
	
		M = Tridiagonal_Matrix()
		
		Eq = zip(* [self.first_eq(t)] +
		           [self.middle_eq(i) for i in range(1, N-1)] +
		           [self.last_eq(t)])
		
		[M.a,M.b,M.c,M.d] = list(Eq)
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
