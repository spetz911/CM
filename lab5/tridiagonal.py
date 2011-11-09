#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce

from pprint import pprint

def IMPORT(self, *args):
	return [self.__getattr__(a) for a in args]
	

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


def main():
	pass


#=====================================================================

if __name__ == "__main__":
	main()
