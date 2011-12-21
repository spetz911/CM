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
from gen import gen_elems, gen_lim_points, gen_points
from gen import generate

eps = 0.0001

import matrix_2

def square(elem, points):
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
	return [ x0 + d*k for k in range(int((x1-x0+1e-13)/d)+1)]

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

alpha = 1
phi = lambda x,y: 1.0
func_start = lambda x,y: 100.0


def elem_plane(i, points):
	""" i - index of point in element"""
	# elem - num of element, p - point of element
#	if not(p in elem): return None
#	print (points)

	system = [ [1.0, p[0], p[1]]  for p in points ]
	b = [0.0] * 3
	b[i] = 1.0
	
	m = matrix_2.Matrix(m = 3, n = 3)
	m.M = system
	m.build_LU()
#	print ("Solve:")
	koef = m.solve(m.shift_b(b))
	koef = [round(xx, 2) for xx in koef]
#	print("koef = ", koef)
	return koef



def matrix_el(elem, lim_points, points):
	res = matrix_2.Matrix(m = 3, n = 3)

	values = [points[x] for x in elem]
	k3 = [elem_plane(i, values) for i in range(len(values))]
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
#		print("area", elem)
		mat.M = [[0.0]*3] * 3
		G = mat
		f_bord = 0
	vec = matrix_2.Matrix(m = 3, n = 1)
	center = median(points[elem[0]],points[elem[1]],points[elem[2]])
	f_start = func_start(*center) * square(elem, points) / 3
	vec.M[0][0] = f_bord - f_start
	vec.M[1][0] = -f_start
	vec.M[2][0] = f_bord - f_start
	
	return (res + mat), vec

def asamble(elems, lim_points, points):
	count = len(elems)
	mat = matrix_2.Matrix(m = count, n = count)
	vec = matrix_2.Matrix(m = count, n = 1)
	
	for elem in elems:
		zm, zv = matrix_el(elem, lim_points, points)
		
		for i,k in enumerate(elem):
			for j,l in enumerate(elem):
				mat[k][l] = mat[k][l] + zm[i][j]	
		
		for i,k in enumerate(elem):
			vec[k][0] = vec[k][0] + zv[i][0]
	return mat,vec
	
# mat, vec = asamble(gen_elems(), gen_lim_points(), gen_points())

mat, vec = asamble(*generate())

mat.pr()
vec.pr()



mat.build_LU()

print ("Solve:")
answ = mat.solve(mat.shift_b(vec.transponate()[0]))
answ = [round(xx, 2) for xx in answ]
print("koef = ", answ)


