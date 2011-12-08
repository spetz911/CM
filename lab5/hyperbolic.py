#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : hyperbolic.py
### Author  : Oleg Baskakov
### Description : hyperbolic pde solver
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp,sqrt
from copy import copy, deepcopy
from functools import reduce
from tridiagonal import *

from pprint import pprint

from pde import *


class Hyperbolic_PDE(PDE):
	approximate_init = '1lvl'
	approximate_boundary = '1lvl'
	
	u_x = 0.0
	u = 0.0

	coef_a = [+1, -2, +1]
	coef_b = [-3, +4, -1]
	coef_c = [+0, +1, +0]
	coef_t = [+1, -2, +1]
	
	def __init__(self, pde = None):
		super(Hyperbolic_PDE, self).__init__(pde)
		MetaClass.print(self)


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
