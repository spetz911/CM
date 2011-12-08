#! /usr/bin/env python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : lab_5.py
### Author  : Oleg Baskakov
### Description : test cases for 5th lab
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

import sys
from tridiagonal import *
from pde import *
from parabolic import *
from hyperbolic import *
from parser import *


##====================================================================
## Description: parse all input
def parse_file(f):
	pde = PDE_parser()
	pde.parse_pde(f.readline())
	if pde.type == 'parabolic':
		pde.parse_boundary_condition(f.readline())
		pde.parse_boundary_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		for s in f.readlines():
			pde.parse_stuff(s)
		return Parabolic_PDE(pde)

	elif pde.type == 'hyperbolic':
		pde.parse_boundary_condition(f.readline())
		pde.parse_boundary_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		pde.parse_initial_condition(f.readline())
		for s in f.readlines():
			pde.parse_stuff(s)
		return Hyperbolic_PDE(pde)
	
	elif pde.type == 'elliptic':
		pass

	return pde



#=====================================================================
def main():
	if len(sys.argv) < 2: print("please enter test")
	path = sys.argv[1]
	
	f = open(path)
	pde = parse_file(f)
	f.close()

#	print(pde.left[2](0))
#	print(pde.right[2](0))
	
	errs, b1, b2 = pde.check_scheme()
	print("err full analitic solve ->", max(errs))
	print("err bound1 ->")
	print_vec(b1)
	print("err bound2 ->")
	print_vec(b2)

#	print_vec(errs)
	
#	res = pde.solve('explicit')
	res = pde.solve(pde.method)
	grid = list(frange(0, pde.l, pde.h))
	origin = [[pde.res_fun(x,k*pde.tau) for x in grid]
									for k in range(6)]
	
	for k in range(6):
		print("iter ", k)
		print_vec(res[k])
		print_vec(origin[k])
	
	f = open("result", 'w')
	f.write(str({'result':res, 'origin':origin, 'grid':grid}))
#	f.write("\n")
	f.close()
	
#	print(pde.u_x)


#=====================================================================
if __name__ == "__main__":
	main()

