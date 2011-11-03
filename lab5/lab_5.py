#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
eps = 0.0001

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




##==========================================================================

def solve_parabolic(pde):
	
	
	
	







def main2():

	for (step,t) in enumerate(frange(0, T, tau)):
		U1.append( implicid_method(U[-1], step))
		U2.append( explicid_method(U[-1], step))




s = "15*u_t = 0.1*u_xx + 0.2*u_x + -3.0*u + x*t"



##===========================================================================
## Description: parse all input
def parse_all(s):
	print("s:", s)
	pde = parse_pde(s)
	print("pde =")
	print(pde)
	if   pde['type'] == 'parabolic':
		print("f(1,1) =", pde['fun'](1,1))
	elif pde['type'] == 'hyperbolic':
		pass
	elif pde['type'] == 'elliptic':
		pass
	
	
	
## Description: for parsing onr element of summ
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

## Description: for parsing partial differential equation
def parse_pde(s):
	pde = dict()
	[left, right] = s.lower().rstrip().replace(' ','').split('=')
	pde['type'] = None
	if ('u_xx' in right) and ('u_t' in left):
		pde['type'] = 'parabolic'
	elif ('u_xx' in right) and ('u_tt' in left):
		pde['type'] = 'hyperbolic'
	elif ('u_xx' in right) and ('u_yy' in right):
		pde['type'] = 'elliptic'
	tokens = [ x for x in left.split('+') + right.split('+')]
	rest = []
	for x in tokens:
		t = get_token(x)
		if t:
			pde[t[0]] = t[1]
		else:
			rest.append(x)
	try: pde['fun'] = eval("lambda x,t: " + "+".join(rest))
	except: pde['fun'] = lambda x,t: 0.0
	return pde

## Description: for parsing other coefficients
def parse_stuff(s, stuff = dict()):
	[left, right] = s.lower().replace(' ','').split('=')
	if left in ['l', 'lx', 'ly', 't', 'w', 'eps']:
		stuff[left] = float(right)
	return stuff

## Description: for solve start
def parse_initial_condition(s, i_cond = dict()):
	[left, right] = s.lower().replace(' ','').split('=')
	if left in ['u(x,0)', 'u(x,y,0)', 'u(x,y)']:
		try: i_cond['init'] = eval("lambda x,y=0: " + right)
		except: i_cond['init'] = lambda x,y=0: 0.0
	if left in ['u_t(x,0)', 'u_t(x,y,0)', 'u_t(x,y)']:
		try: i_cond['init1'] = eval("lambda x,y=0: " + right)
		except: i_cond['init1'] = lambda x,y=0: 0.0
	return i_cond

## Description: for solve limits
def parse_boundary_condition(s, b_cond = dict()):
	[left, right] = s.lower().replace(' ','').split('=')
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
		b_cond['north'] = result
	elif ('(0,x' in left):
		b_cond['south'] = result
	elif ('(l' in left): # (lx,y,0)
		b_cond['east'] = result
	elif ('(0' in left): # (0,y,0)
		b_cond['west'] = result
	
	try: b_cond['left'] = b_cond['west']
	except: pass
	try: b_cond['right'] = b_cond['east']
	except: pass
	return b_cond














#============================================================================

def main():

	s = open("input").readline()
	
	print(parse_all(s))


#=====================================================================

if __name__ == "__main__":
	main()

