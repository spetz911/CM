#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy

segment = [-1.0, 1.0]

def func(x):
	return x/(3*x + 4)**2

def func2(x):
	return x/(2*x + 5)

def roundx(v):
	return [round(xx, 4) for xx in v]

def P_n(A):
	return lambda x: sum([ A[i]*x**i for i in range(0,len(A)) ] )

def MethodOfRectangles(f, seg, h):
	n = round((seg[1] - seg[0])/h )
	return sum( [ h * f( seg[0] + h*i + h/2) for i in range(n)] )
	
def MethodOfTrapezoids(f, seg, h):
	n = round((seg[1] - seg[0])/h )
	y = [ f(seg[0] + h*i) for i in range(n+1)]
	return sum( [ (y[i] + y[i+1])*h/2  for i in range(n)] )
	
def SimpsonMethod(f, seg, h):
	n = round((seg[1] - seg[0])/h )
	n2 = round((seg[1] - seg[0])/h)//2
	y = [ f(seg[0] + h*i) for i in range(n+1)]
	return sum( [( (y[2*i] + 4*y[2*i+1]+y[2*(i+1)] )*h/3 )
											for i in range(n2)] )

def RRR(f1, f2, k, p):
	""" Runge-Romberga-Richardsona """
	# k = h1/h2, p - interpolation level
	return f1 + (f1 - f2) / (k**p - 1)


def main():

	r1 = MethodOfRectangles(func, segment, 0.5)
	print("Rect 0.5:\t", r1)
	t1 = MethodOfTrapezoids(func, segment, 0.5)
	print("Trap 0.5:\t", t1)
	s1 = SimpsonMethod (func, segment, 0.5)
	print("Simspson 0.5:\t", s1)

	r2 = MethodOfRectangles(func, segment, 0.25)
	print("Rect 0.25:\t", r2)
	t2 = MethodOfTrapezoids(func, segment, 0.25)
	print("Trap 0.25:\t", t2)
	s2 = SimpsonMethod (func, segment, 0.25)
	print("Simspson 0.25:\t", s2)

	print("RRR of Rect:\t", RRR(r2, r1, 2, 2) )
	print("RRR of Trap:\t", RRR(t2, t1, 2, 2) )
	print("RRR of Simpson:\t", RRR(s2, s1, 2, 2) )

	print("answer\t", SimpsonMethod(func, segment, 0.0001) )
	
	return 0

main()
