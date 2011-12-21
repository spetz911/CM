#! /usr/bin/env python
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : postprocessor.py
### Author  : Oleg Baskakov
### Description : postprocessor show result
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import sys



#=====================================================================

def main():
	try:
		k = int(sys.argv[1])
	except:
		k = 0
	
	if len(sys.argv) <= 1:
		f = open("result")
		dd = eval(f.read())
		Y = dd['tmp']
		X = list(range(len(Y)))
		plt.plot(X, Y, color="black")
		plt.show()
		return
	
	
	
	
	print("k =", k)

	f = open("result")
	dd = eval(f.read())
	res = dd['result']
	if dd.has_key('grid'):
		grid = dd['grid']
	else:
		grid = range(len(res))
	
	if dd.has_key('origin'):
		origin = dd['origin']
		plt.plot(grid, origin[k], color="blue")

	plt.plot(grid, res[k], color="red")

	plt.show()

	
	
	
#=====================================================================

if __name__ == "__main__":
	main()


