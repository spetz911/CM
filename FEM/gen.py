#! /usr/bin/python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : pde.py
### Author  : Oleg Baskakov
### Description : Partial Differencial Equation
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
from itertools import product as dekart_product
import matrix_2

eps = 0.0001

def frange(x0, x1, d):
	return [ x0 + d*k for k in range(int((x1-x0+eps/10)/d)+1)]

def print_mat(mat):
	print("matrix:")
	for row in mat:
		for elem in row:
			print(" %.2f "%elem, end = '\t')
		print()
	print("--------")

def print_vec(row):
	print("vector:")
	for elem in row:
		print(round(elem, 4), end = "\t")
	print()
	print("--------")

def el_func(x, el, grid):
	if not( 0 <= el <= len(grid)): return 0
	
	z = -abs(x - grid[el])/(grid[1] - grid[0]) + 1.0
	if z<0: z = 0
	return z

points = [(0.5, 1.5),
		  (0.0, 2.0),
		  (1.0, 2.0),
		  (0.0, 1.0),
		  (1.0, 1.0),
		  (0.5, 0.5),
		  (0.0, 0.0),
		  (1.0, 0.0),
		  (1.5, 0.5),
		  (2.0, 1.0),
		  (2.0, 0.0)]

lim_point_tmp = [1,2,4,9,10,7,6,3,1]

lim_points = list(zip(lim_point_tmp, lim_point_tmp[1:]))
lim_points += [(b,a) for a,b in lim_points]





elems_tmp = [(1,2,3),
		     (1,2,4),
		     (1,3,5),
		     (1,4,5),
		     (6,4,5),
		     (6,7,4),
		     (6,5,8),
		     (6,8,7),
		     (9,5,10),
		     (9,8,5),
		     (9,10,11),
		     (9,11,8)]
		     
elems = [(a-1,b-1,c-1) for (a,b,c) in elems_tmp]


def elem_fun(koef):
	return lambda x,y: koef[0] + koef[1]*x + koef[2]*y

# F = elem_plane(elems[0],0)

def elem_fun_all(el, i):
	A = [k for k in [0,1,2] if elems[el][k] == i]
	if A:
		return elem_fun(el, A[0])
	else:
		return lambda x,y: 0.0


def gen_elems():
	"""undefined"""
	return [(a-1,b-1,c-1) for (a,b,c) in elems_tmp]

def gen_lim_points():
	"""undefined"""
	lim_point_tmp = [1,2,4,9,10,7,6,3,1]
	lim_points = list(zip(lim_point_tmp, lim_point_tmp[1:]))
	lim_points += [(b,a) for a,b in lim_points]
	return lim_points



def el_func(x, el, grid):
	if not( 0 <= el <= len(grid)): return 0
	
	z = -abs(x - grid[el])/(grid[1] - grid[0]) + 1.0
	if z<0: z = 0
	return z

def gen_points():
	return   [(0.5, 1.5),
			  (0.0, 2.0),
			  (1.0, 2.0),
			  (0.0, 1.0),
			  (1.0, 1.0),
			  (0.5, 0.5),
			  (0.0, 0.0),
			  (1.0, 0.0),
			  (1.5, 0.5),
			  (2.0, 1.0),
			  (2.0, 0.0)]



alpha = 1
phi = lambda x,y: 1.0
func_start = lambda x,y: 1.0


elems_tmp = [(1,2,3),
		     (1,2,4),
		     (1,3,5),
		     (1,4,5),
		     (6,4,5),
		     (6,7,4),
		     (6,5,8),
		     (6,8,7),
		     (9,5,10),
		     (9,8,5),
		     (9,10,11),
		     (9,11,8)]
		     

def elem_plane(elem, p):
# elem - num of element, p - point of element
	if not(p in elem): return None
	system = [ [1.0, points[p][0], points[p][1]]  for p in elem ]
	b = [0.0] * 3
	b[elem.index(p)] = 1.0
	
	m = matrix_2.Matrix(m = 3, n = 3)
	m.M = system
	m.build_LU()
#	print ("Solve:")
	koef = m.solve(m.shift_b(b))
	koef = [round(xx, 2) for xx in koef]
#	print("koef = ", koef)
	return koef

def elem_fun(koef):
	return lambda x,y: koef[0] + koef[1]*x + koef[2]*y

F = elem_plane(elems[0],0)

def elem_fun_all(el, i):
	A = [k for k in [0,1,2] if elems[el][k] == i]
	if A:
		return elem_fun(el, A[0])
	else:
		return lambda x,y: 0.0


def generate():
	"""Generate params for FEM solver"""
	points = []
	border = []
	elems = []
	
	lx = 1.0
	ly = 1.0
	hx = 0.5
	hy = 0.5
	
	field = None
	
	field = [[(y,x)  for x in frange(0, lx, hx)]
	                 for y in frange(0, ly, hy)]
	
	m = len(field)
	n = len(field[0])
	
	print("n =", n)
	print("m =", m)
	
	#correct last string
	num = lambda i,j: i*(2*m-1) + ((i!=n-1) + 1)*j

	print ("test", num(0,0), num(1,0), num(0,1), num(2,1))

	border.extend(field[0])
	border.extend(list(zip(*field))[-1])
	tmp = list(field[-1])
	tmp.reverse()
	border.extend(tmp)
	tmp = list(list(zip(*field))[0])
	tmp.reverse()
	border.extend(tmp)
	
	points = [(0.0,0.0)] * (n*m + (n-1)*(m-1))
	
	print("len", len(points))
	
	for i in range(n-1):
		for j in range(m-1):
			# add point's coord
			points[num(i+0,j+0)] = field[i+0][j+0]
			points[num(i+0,j+1)] = field[i+0][j+1]
			points[num(i+1,j+0)] = field[i+1][j+0]
			points[num(i+1,j+1)] = field[i+1][j+1]
			
			#add central point with shift = 1
			zx,zy = field[i][j]
			points[num(i, j) + 1] = (zx + hx/2, zy + hx/2)
			
			#add elem index			
			elems.append((num(i,j), num(i,j+1), num(i,j)+1))
			elems.append((num(i,j), num(i+1,j), num(i,j)+1))
			elems.append((num(i,j+1), num(i+1,j+1), num(i,j)+1))
			elems.append((num(i+1,j), num(i+1,j+1), num(i,j)+1))
	
	from pprint import pprint
	pprint(field)
	print(list(enumerate(points))[:9])
	print(border)
	#convert coord to index
	border = list({points.index(p) for p in border})
	
	border.append(border[0])
	
	lim_points = list(zip(border, border[1:]))
#	lim_points += [(b,a) for a,b in lim_points]

	
	return elems, lim_points, points

def main():
	e, p, b = generate()

	print(b)
	
#=====================================================================
if __name__ == "__main__":
	main()



