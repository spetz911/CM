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
	print_vec(errs)
	print("err bound1 ->")
	print_vec(b1)
	print("err bound2 ->")
	print_vec(b2)

#	print_vec(errs)
	
#	res = pde.solve('explicit')
	res = pde.solve(pde.method)
	count = len(res)
	grid = list(frange(0, pde.l, pde.h))
	origin = [[pde.res_fun(x,k*pde.tau) for x in grid]
	                                    for k in range(count)]
	
	for k in range(5):
		print("iter ", k)
		print_vec(res[k])
		print_vec(origin[k])

	print("h =", pde.h)
	sigma = pde.a**2 * pde.tau / pde.h**2
	print("SIGMA =", round(sigma, 5))
#	err1 = max([abs(u-v)/abs(v+0.0000001) for u,v in zip(res[-1], origin[-1])])
	err1 = max([abs(u-v) for u,v in zip(res[-1], origin[-1])])
	print("ERR_1 =", err1)
#	err2 = (sum([(abs(u-v)/abs(v+0.0000001))**2 for u,v in zip(res[-1], origin[-1])]) / len(res[-1]))**0.5
	err2 = (sum([abs(u-v)**2 for u,v in zip(res[-1], origin[-1])]) / len(res[-1]))**0.5
	print("ERR_2 =", err2)
	print("ITER =", int(pde.t/pde.tau + 0.00001))
	
	tmp = []
	for us,vs in zip(res, origin):
#		tmp.append( (sum([(abs(u-v+0.000001)/abs(v+0.00000001))**2 for u,v in zip(us, vs)]) / len(us))**0.5)
#		tmp.append(max([abs(u-v)/abs(v+0.0000001) for u,v in zip(us,vs)]))
		tmp.append(max([abs(u-v) for u,v in zip(us,vs)]))
#		tmp.append( min([abs(v) for u,v in zip(us, vs)]))
	
	f = open("result", 'w')
	f.write(str({'result':res, 'origin':origin, 'grid':grid, 'tmp' : tmp, 't' : pde.t, 'tau':pde.tau}))
#	f.write("\n")
	f.close()
	
#	print(pde.u_x)


#=====================================================================
if __name__ == "__main__":
	main()

