#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce
from matrix_2 import *


#Xi = [-1.0, 0.0, 1.0, 2.0    , 3.0, 4.0    ]
#Yi = [-0.5, 0.0, 0.5, 0.86603, 1.0, 0.86603]

Xi = [0.0, 0.1   , 0.2   , 0.3   , 0.4   ]
Yi = [1.0, 1.1052, 1.2214, 1.3499, 1.4918]

X0 = 0.2
Y0 = 2.4969

def roundx(v):
	return [round(xx, 4) for xx in v]

def P_n(A):
	return lambda x: sum([ A[i]*x**i for i in range(0,len(A)) ] )

def Differetiation( Y, X, n, x):
	if (x<X[0])  : return None
	if (x>=X[-2]-0.001) : return None
	i = 0
	while (x > X[i]+0.0001) : i+=1
	i-=1
	print("segment", [X[i], X[i+1] ], ", x =", x )

	if   (n==1):
		return (Y[i+1]-Y[i])/(X[i+1]-X[i]) + ( ((Y[i+2]-Y[i+1])/(X[i+2]-X[i+1]) - (Y[i+1]-Y[i])/(X[i+1]-X[i]))/(X[i+2]-X[i]) ) * (2*x-X[i]-X[i+1])
	elif (n == 2):
		return 2*(((Y[i+2] - Y[i+1]) / (X[i+2] - X[i+1]) - (Y[i+1] - Y[i]) / (X[i+1] - X[i])) / (X[i+2] - X[i]));
	else:
		return None
		


def diff1(Y, X, n):
	return lambda x: Differetiation( Y, X, n, x)

def main():
	print("table function: ")
	print(Xi)
	print(Yi)

	d1 = diff1(Yi, Xi, 1)
	d2 = diff1(Yi, Xi, 2)
	
	print("diff1", d1(X0) )
	print("diff1", d2(X0) )
	
	return 0

main()
