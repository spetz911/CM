#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

#~ import matplotlib.pyplot as plt
#~ import numpy as np

from math import sqrt,sin,cos,atan,pi,log
from copy import copy, deepcopy
from functools import reduce

Xi = [0.1*pi, 0.2*pi, 0.3*pi, 0.4*pi]
Yi = [sin(xx) for xx in Xi]



Y1_global = []

#Xi = [0.0, 1.0, 2.0, 3.0]
#Yi = [round(sin(xx*pi/6), 5) for xx in Xi]

X0 = pi/4
Y0 = sin(X0)

def product(L):
	mul = lambda x, y: x*y
	return reduce(mul, L)


def l_i(x, i, X):
	""" OK """
	n = len(X)
	S = 1
	for j in range(n):
		if (i!=j):
			S *= (x - X[j])/(X[i] - X[j])
	return S

def L_n(x, Y, X):
	return sum( [(Y[i]*l_i(x, i, X)) for i in range(len(Y))] )


def L(Y, X):
	return lambda x: L_n(x, Y, X)

def Y_m(Y, X):
	m = len(Y)
	if   (m == 1):
		return Y[0]
	elif (m == 2):
		return(Y[0] - Y[1])/(X[0] - X[1])
	else:
		return(Y_m(Y[0:m-1],X[0:m-1])-Y_m(Y[1:m],X[1:m]) )/(X[0]-X[m-1])

def P_m(X):
	if (len(X) == 0): return lambda x: 1
	return lambda x:product( [(x - x_i) for x_i in X])

def N_1(x, Y1, P1):
	return sum( [(Y*P(x)) for (Y,P) in zip(Y1,P1)] )

def N(Y, X):
	m = len(X)
	Y1 = [Y_m(Y[0:i],X[0:i]) for i in range(1,m+1)]
	P1 = [P_m(X[0:i])        for i in range(0,m)]
	
	print("polynom koef:")
	for yy in Y1:  print(round(yy, 5))

	return lambda x: N_1(x,Y1,P1)

def main():
	print("table function: ")
	print(Xi)
	print(Yi)
	
	lagranje = L(Yi, Xi)
	newton   = N(Yi, Xi)
	
	print("Test lagranje", Yi[2], "=", lagranje(Xi[2]) )
	print("Test newton  ", Yi[2], "=", newton  (Xi[2]) )
	
	#~ X = np.arange(Xi[0], Xi[-1], 0.01)
	#~ Y = [ lagranje(xx) for xx in X]
	
	#~ plt.plot(X,Y,color="red")
	#~ plt.plot(Xi,Yi,color="grey")
	#~ plt.show()

	
	
	
	eps_n = abs( Y0 - lagranje(X0) )
	print( "|sin(x) - P_n(x)| = ", eps_n)
	
	return 0

main()
