#! /usr/bin/env python3
# -*- coding: utf-8 -*-

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

	elif pde.type == 'elliptic':
		pass

	return pde



#=====================================================================

def main():

	f = open("input")
	
	pde = parse_file(f)

#	print(pde.left[2](0))
#	print(pde.right[2](0))
	
	res = pde.solve('implicit')
	
	for k in [0,1,2,3]:
		print("iter ", k)
		print_vec(res[k])
		print_vec([pde.fun(x, k*pde.tau) for x in frange(0, pde.l, pde.h)])
	
#	print(pde.u_x)


#=====================================================================

if __name__ == "__main__":
	main()


