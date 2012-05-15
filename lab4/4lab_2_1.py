#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import e
from lib import *
from my_matrix_2 import *
#xy"+2y'-xy=0
#y'(1)=0
#1.5y(2)+y'(2)=e**2

#y'=z
#z'=y-2/x*z
#z(1) = 0
#z(2)-1.5y(2) = e**2

x0 = [1., 2.]
#~ koefs = [[0., 1., 0.], [1.5, 1, exp(2)]]
x0 = 1.
a = 1.
b = 2.
eps = .01
h = 0.1
eta0 = 2.
y01 = [eta0,0.]
eta1 = 3.
y02 = [eta1,0.]


def solve(x):
	return e**(x)/x


def F(x,(y,z)):
	return [z,y-2/x*z]

def G(slv):
	return slv[-1][1][1]+1.5*slv[-1][1][0]-e**2

def p(x):
	return 2/x

def q(x):
	return -1

def shooting(h):
	global eta0,eta1
	s = ODY(F,h,x0,y01,a,b)
	s1 = ODY(F,h,x0,y02,a,b)
	slv = s.Rynge_Kytte()
	slv1 = s1.Rynge_Kytte()
	while abs(G(slv1)) > eps:
		eta2 = eta1 - G(slv1)*(eta1-eta0)/(G(slv1)-G(slv))
		s2 = ODY(F,h,x0,[eta2,0.],a,b)
		slv2 = s2.Rynge_Kytte()
		eta0 = eta1
		eta1 = eta2
		slv = slv1
		slv1 = slv2
	return slv2
	
def KRM(h):
	x = frange(1.,2.+h,h)
	f = open("in","w")
	f.write("-100 100 0\n")
	for k in xrange(1,len(x)-1):
		f.write("%f %f %f %f\n" %((1-p(x[k])*h/2)*100, (-2+h**2*q(x[k]))*100,(1+h*p(x[k])/2)*100,0))
	f.write("-100 %f %f\n" %((1.5*h+1)*100,h*e**2*100))
	f.close()
	
	m = Tridiagonal_Matrix("in")
	
	return m.solve()

def main():
	
	global h
	print "Метод стрельбы."
	slv = shooting(h)
	slv2 = shooting(h/2.)
	
	
	for i in xrange(len(slv)):
		print "(%f, %f),%g %g" %(slv[i][0],slv[i][1][0],abs(slv[i][1][0]-solve(slv[i][0])),abs(slv[i][1][0]-slv2[2*i][1][0]))
	
	print "Конечно-разностный метод."
	
	y = KRM(0.1)
	y2 = KRM(0.1/3)
	
	for i in xrange(len(y)):
		print "(%f, %f),%g %g" % (1+i*h,y[i],abs(y[i]-solve(1+i*h)),abs(y[i]-y2[2*i]))

	return 0

main()

