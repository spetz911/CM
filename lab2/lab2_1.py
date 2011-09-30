#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi,log

def iterations(phi,segment,q,eps):
		xk = 0.0
		xkm = (segment[0]+segment[1])/2.0
		counter = 1
		while 42:
			counter+= 1
			if counter > 9000 : break
			xk = phi(xkm)
			if q/(1-q)*abs(xk-xkm) < eps:
				#print ("finish with %d iters" % counter)
				return xk, counter
			xkm = xk
	
def Newton(func, func1, segment, eps):
		xk = 0.0
		xkm = (segment[0]+segment[1])/2.0
		counter = 1
		while 42:
			counter+= 1
			if counter > 9000 : break
			xk = xkm - func(xkm)/func1(xkm)
			if abs(xk-xkm) < eps:
				#print ("finish with %d iters" % counter)
				return xk, counter
			xkm = xk
		

def func(x):
	return 2**x - x*x - 0.5
	
def phi(x):
	return log(x*x+1/2, 2)
	
def func1(x):
	return (2**x)*log(2) - 2.0*x
	
def main():
	eps = float( input("задайте eps: ") )
	s1 = iterations(phi,[4.0,4.5], 0.9, eps)
	print ("решение итерациями: %f за %d шагов" % s1)
	
	s2 = Newton(func,func1,[4.0,4.5], eps)
	print ("решение Ньютоном: %f за %d итераций" % s2)
	
	print("Зависимость сходимости от числа итераций:")
	print("eps\titer\tNewton")
	for eps in [0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001]:
		s1, c1 = iterations(phi,[4.0,4.5], 0.9, eps)
		s2, c2 = Newton(func,func1,[4.0,4.5], eps)
		print(eps, "\t%d \t%d" %(c1,c2) )

	return

if (__name__ == "__main__"):
	main()
