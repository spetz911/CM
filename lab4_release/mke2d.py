#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,tan,pi,log,exp
from copy import copy, deepcopy
from functools import reduce
from itertools import product as dekart_product
eps = 0.0001

import matrix_2

def square(elem):
	p1 = points[elem[0]]
	p2 = points[elem[1]]
	p3 = points[elem[2]]
	v1 = (p1[0]-p3[0], p1[1]-p3[1])
	v2 = (p2[0]-p3[0], p2[1]-p3[1])
	res = 0.5 * abs(v1[0]*v2[1] - v2[0]*v1[1])
#	print("sqr =", res, elem)
	return res

def median(*v):
	return (sum([x[0] for x in v])/len(v), sum([x[1] for x in v])/len(v))

def dist(A, B):
	return sqrt(sum([ (xi-yi)**2 for (xi,yi) in zip(A,B)]))

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
		     
elems = [(a-1,b-1,c-1) for (a,b,c) in elems_tmp]

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

def matrix_el(elem):
	res = matrix_2.Matrix(m = 3, n = 3)
	k3 = [elem_plane(elem, x) for x in elem]
	beta = [ a[1] for a in k3]
	gama = [ a[2] for a in k3]
	res.M = [[(beta[m]*beta[n] + gama[m]*gama[n]) for n in [0,1,2]]
	                                              for m in [0,1,2]]
	border = [ x for x in dekart_product(elem,elem) if x in lim_points]
	mat = matrix_2.Matrix(m = 3, n = 3)
	if border:
		edge = border[0]
		print("border")
		p1 = points[edge[0]]
		p2 = points[edge[1]]
		tmp_k = alpha/3 * dist(p1, p2)
		mat.M = [[1.0, 0.0, 0.5],
		         [0.0, 0.0, 0.0],
		         [0.5, 0.0, 1.0]]
		G = mat * tmp_k
		f_bord = phi(*median(p1, p2)) * 1/2*dist(p1, p2)

	else:
		print("area", elem)
		mat.M = [[0.0]*3] * 3
		G = mat
		f_bord = 0
	vec = matrix_2.Matrix(m = 3, n = 1)
	center = median(points[elem[0]],points[elem[1]],points[elem[2]])
	f_start = func_start(*center) * square(elem) / 3
	vec.M[0][0] = f_bord - f_start
	vec.M[1][0] = -f_start
	vec.M[2][0] = f_bord - f_start
	
	return (res + mat), vec

def asamble(elems):
	count = len(elems)
	mat = matrix_2.Matrix(m = count, n = count)
	vec = matrix_2.Matrix(m = count, n = 1)
	
	for elem in elems:
		zm, zv = matrix_el(elem)
		
		for i,k in enumerate(elem):
			for j,l in enumerate(elem):
				mat[k][l] = mat[k][l] + zm[i][j]	
		
		for i,k in enumerate(elem):
			vec[k][0] = vec[k][0] + zv[i][0]
	return mat,vec
	
m, v = asamble(elems)
m.pr()
v.pr()




