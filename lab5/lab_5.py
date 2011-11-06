#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
eps = 0.0001

from pde import *

def frange(x0, x1, d):
	return [ x0 + d*k for k in range(int((x1-x0+eps/10)/d)+1)]

def print_mat(mat):
	print("matrix:")
	for row in mat:
		for elem in row:
			print(round(elem, 3), end = "\t")
		print()
	print("--------")

def print_vec(row):
	print("vector:")
	for elem in row:
		print(round(elem, 4), end = "\t")
	print()
	print("--------")

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
			P.append( -c[i] / (b[i]+a[i]*P[i-1]) ) 
			Q.append( (d[i] - a[i]*Q[i-1]) / (b[i] + a[i]*P[i-1]) )

		x = [0]*n
		x[n-1] = Q[n-1]
		for i in range(n-2, -1, -1):
			x[i] = P[i]*x[i+1] + Q[i]
			
		return x



start_f = lambda x: sin(x)
border1 = lambda x: 0
border2 = lambda x: 0

x0 = 0
x1 = 1

h = 0.1
tau = 0.5

T = 7.0

u0 = [start_f(x) for x in frange(x0, x1, h)]
U1 = [u0]
U2 = [u0]



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
		if left in ['l', 'lx', 'ly', 't', 'w', 'eps']:
			self.__setattr__(left, float(right))

	## Description: for solve init eq
	def parse_initial_condition(self, line):
		[left, right] = line.lower().replace(' ','').split('=')
		if left in ['u(x,0)', 'u(x,y,0)', 'u(x,y)']:
			try: self.init = eval("lambda x,y=0: " + right)
			except: self.init = lambda x,y=0: 0.0
		if left in ['u_t(x,0)', 'u_t(x,y,0)', 'u_t(x,y)']:
			try: self.init1 = eval("lambda x,y=0: " + right)
			except: self.init1 = lambda x,y=0: 0.0

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




s = "15*u_t = 0.1*u_xx + 0.2*u_x + -3.0*u + x*t"



##===========================================================================
## Description: parse all input
def parse_all(s):
	print("s:", s)
	pde = PDE_parser()
	pde.parse_pde(s)
	print("pde =")
	print(pde.type)
	if pde.type == 'parabolic':
		print("f(1,1) =", pde.fun(1,1))
		eq = PDE_parser()
		
		
	elif pde.type == 'hyperbolic':
		pass
	elif pde.type == 'elliptic':
		pass
	
	
	return None







#============================================================================

def main():

	s = open("input").readline()
	
	print(parse_all(s))


#=====================================================================

if __name__ == "__main__":
	main()

