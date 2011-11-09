#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce

from pprint import pprint

eps = 0.0001

def print_mat(mat):
	print("matrix:")
	for row in mat:
		for elem in row:
			print(round(elem, 3), end = "\t")
		print()
	print("--------")

def print_vec(row, k = 4):
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


class Tridiagonal_Matrix:
	def solve(self):
		"""Method progonki"""
		a = self.a
		b = self.b
		c = self.c
		d = self.d
		n = self.n
		
		P = []
		Q = []
		P.append(-c[0]/b[0])
		Q.append( d[0]/b[0])
		
		for i in range(1, n):
			print(P[-1], Q[-1])
			P.append( -c[i] / (b[i]+a[i]*P[i-1]) )
			Q.append( (d[i] - a[i]*Q[i-1]) / (b[i] + a[i]*P[i-1]) )

		x = [0]*n
		x[n-1] = Q[n-1]
		for i in range(n-2, -1, -1):
			x[i] = P[i]*x[i+1] + Q[i]
			
		return x


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
		
		coefficients = PDE.vec_mat([self.sigma, self.omega, self.eta],
		                       [self.coef_a, self.coef_b, self.coef_c])
		# TODO f(x,t) != 0
#		print(list(zip(* [self.coef_a, self.coef_b, self.coef_c])))
		
		print_mat([self.coef_a, self.coef_b, self.coef_c])
		
		print_vec(coefficients)
		print([self.sigma, self.omega, self.eta])
		res = [0] + [U[i] + PDE.scalar(coefficients, U[i-1:i+2])
		             for i in range(1, N-1)] + [0]
		i = 2
		print("Ui-1..Ui+1",  U[i-1:i+2])
		
		
		(a0, b0, c0, d0) = self.first_eq(self.tau*N)
		(an, bn, cn, dn) = self.last_eq(self.tau*N)
		print(self.first_eq(self.tau*N))
		print(self.last_eq(self.tau*N))
		
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
		           [self.middle_eq(k) for k in range(1, N-1)] +
		           [self.last_eq(t)])
		
		Eq = list(Eq)
		
		[M.a,M.b,M.c,M.d] = Eq
		M.n = N
	
		print_mat(Eq)

		x = M.solve()
		return x

	## Description: solve with method 'Progonki'
	def Crank_Nicolson_method(self):
		pass
		# a_i /= Q 
		# b_i /= Q 
		# c_i /= Q 
		# d_i += explicit_method * (1-Q) 
	


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

	def middle_eq(self, i):
		"""Find coefficients of middle equation"""
		U = self.grid[-1]
		N = len(U)
		
		ai = sum(self.coef_a)
		bi = sum(self.coef_b)
		ci = sum(self.coef_c)
		di = -U[i]

		return (ai, bi, ci, di)
	
	def solve(self):
		Us = self.grid
		
#		for t in frange(0, self.t, self.tau):
#			U.append(self.explicit_method())
	
		for t in frange(0, self.t, self.tau):
			Us.append(self.implicit_method())
		
		
		print("ololo")
		return U
	

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

	## Description: find coefficients of last equation
	def middle_eq(u):
		"""Find coefficients of middle equation"""
		pass
		
	def solve(self):
		U = self.grid
		print("ololo")


##====================================================================
## Description: parse all input
def parse_file(f):
	pde = PDE_parser()
	pde.parse_pde(f.readline())
	if pde.type == 'parabolic':
		pde.parse_boundary_condition(f.readline())
		pde.parse_boundary_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		for s in f.readlines():
			pde.parse_stuff(s)
		return Parabolic_PDE(pde)

	elif pde.type == 'hyperbolic':
		pde.parse_boundary_condition(f.readline())
		pde.parse_boundary_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		for s in f.readlines():
			pde.parse_stuff(s)

	elif pde.type == 'elliptic':
		pass

	return pde


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
