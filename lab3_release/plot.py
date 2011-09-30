#! /usr/bin/python2.6
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from math import sqrt,sin,cos,atan,pi,log

# func = lambda x: cos(x)
eps = 0.001

def func(x,el):
	z = -abs(x - el)/0.33 + 1
	if z<0: z = 0
	return z
	
	
seg = (0.0, 1.0)
n = 3

grid = [ seg[0] + k*(seg[1]-seg[0])/n for k in range(n+1) ] 

def el_func(x, el, grid):
	if not( 0 <= el <= len(grid)): return 0
	if x<grid[0] or x>grid[-1]: return 0
	delta = grid[1] - grid[0]
	
	z = -abs(x - grid[el])/delta + 1.0
	if z<0: z = 0
	return z

def get_func(el, grid):
	return lambda x: el_func(x, el, grid)

Funcs = [ get_func(el, grid) for el in range(n+1)]

def df(func):
	return lambda x: (func(x+eps) - func(x-eps))/(2*eps)

def integral(func, seg):
	d = eps
	res = sum([func(x)*d for x in frange(seg[0], seg[1], d)])
	return res
	
	
Xi = [0.1*pi, 0.2*pi, 0.3*pi, 0.4*pi]
Yi = [sin(xx) for xx in Xi]

# X = np.arange(Xi[0], Xi[-1], 0.01)
# Y = [ func(xx) for xx in X]

X = np.arange(-0.5, 1.5, 0.01)
Y = [ df(Funcs[1])(xx) for xx in X]
# Z = [ func(xx,0.66) for xx in X]

A = [[ f(xx) for xx in X] for f in Funcs]


# [plt.plot(X,A[i],color="grey") for i in range(n+1)]
plt.plot(X,Y,color="red")
plt.show()

	
