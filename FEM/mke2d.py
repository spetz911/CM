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

def vec_prod(p1, p2):
	return p1[0]*p2[1] - p1[1]*p2[0]

def vec_diff(p1,p2):
	return (p1[0]-p2[0], p1[1]-p2[1])

def sign(x):
	if x<0:
		return -1
	else:
		return 1

def in_triangle(p0, p1, p2, p3):
	vp1 = vec_prod(vec_diff(p0,p1), vec_diff(p0,p2))
	vp2 = vec_prod(vec_diff(p0,p2), vec_diff(p0,p3))
	vp3 = vec_prod(vec_diff(p0,p3), vec_diff(p0,p1))
	
	if abs(vp1*vp2*vp3) < 1e-13:
		return False
	elif sign(vp1) == sign(vp2) == sign(vp3):
		return True
	else:
		return False

def elem_fun(x,y, koef3, val3):
	if in_triangle((x,y), *val3):
		return sum(k[0] + k[1]*x + k[2]*y for k in koef3)
	else:
		return 0.0

def test_answ(answ, elems, points):
	funcs3 = []
	values = []
	for el in elems:
		values.append([points[x] for x in el])
		koeffs = [elem_plane(i, values[-1]) for i in range(len(values[-1]))]
		funcs3.append(koeffs)
	
	return lambda x,y: sum(u_k*elem_fun(x,y,koef, val) for u_k, koef, val in zip(answ, funcs3, values))
	
	
	

def square(points):
	p1, p2, p3 = points
	v1 = (p1[0]-p3[0], p1[1]-p3[1])
	v2 = (p2[0]-p3[0], p2[1]-p3[1])
	res = 0.5 * abs(v1[0]*v2[1] - v2[0]*v1[1])
#	print("sqr =", res, points)
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
func_start = lambda x,y: 0.0


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
#	koef2 = [round(xx, 2) for xx in koef]
#	print("koef = ", koef2)
	return koef



def matrix_el(elem, lim_points, points):
	print('matrix_el', elem)
	global alpha
	res = matrix_2.Matrix(m = 3, n = 3)

	values = [points[x] for x in elem]
	print('values', values)
	k3 = [elem_plane(i, values) for i in range(len(values))]
	print('k3 =', k3)
	beta = [ a[1] for a in k3]
	gama = [ a[2] for a in k3]

	# local matrix
	res.M = [[(beta[m]*beta[n] + gama[m]*gama[n]) for n in [0,1,2]]
	                                              for m in [0,1,2]]
	el_square = square(values)
	res = res * el_square

	border = [ x for x in dekart_product(elem,elem) if x in lim_points]

#	print('border =', border)
#	print('dekart_product(elem,elem)', list(dekart_product(elem,elem)))
#	print(lim_points)

	mat = matrix_2.Matrix(m = 3, n = 3)
	if border:
		edge = border[0]
#		print("border")
		p1 = points[edge[0]]
		p2 = points[edge[1]]
		tmp_k = alpha/3 * dist(p1, p2)
		mat.M = [[1.0, 0.0, 0.5],
		         [0.0, 0.0, 0.0],
		         [0.5, 0.0, 1.0]]
		G = mat * tmp_k
#		G.pr()
		f_bord = 0.5 * dist(p1, p2) * phi(*median(p1, p2))

	else:
#		print("area", elem)
		mat.M = [[0.0]*3 for i in range(3)]
		G = mat
		f_bord = 0
		
	vec = matrix_2.Matrix(m = 3, n = 1)
	
	f_start = func_start(*median(*values)) * el_square / 3
	
	vec.M[0][0] = f_bord - f_start
	vec.M[1][0] = -f_start
	vec.M[2][0] = f_bord - f_start
	
	return (res + mat), vec

def asamble(elems, lim_points, points):
	count = len(points)
	mat = matrix_2.Matrix(m = count, n = count)
	vec = matrix_2.Matrix(m = count, n = 1)
	mat.M = [[0.0]*count for i in range(count)]
	vec.M = [[0.0] for i in range(count)]
	
	
	print('elems =', elems)
	
	for elem in elems:
		zm, zv = matrix_el(elem, lim_points, points)
		zm.pr()
		zv.pr()
#		exit(0)
		
		for i,k in enumerate(elem):
			for j,l in enumerate(elem):
				mat[k][l] = mat[k][l] + zm[i][j]	
		
		for i,k in enumerate(elem):
			vec[k][0] = vec[k][0] + zv[i][0]
	return mat,vec




# print(gen_points())

# elems, lim_points, points = gen_elems(), gen_lim_points(), gen_points()
# mat, vec = asamble(elems, lim_points, points)



elems, lim_points, points = generate()

print('old', gen_lim_points())
print('new', lim_points)



mat, vec = asamble(elems, lim_points, points)

mat.pr()
vec.pr()


mat.build_LU()

print ("Solve:")
answ = mat.solve(mat.shift_b(vec.transponate()[0]))
answ = [round(xx, 2) for xx in answ]
print("answ = ", answ)

# answ = [1.0] *len(elems)
# answ[1] = 1.0

# answ = [i/len(elems) for i in range(len(elems))]


Fun = test_answ(answ, elems, points)
Mat = [[Fun(x,y) for y in frange(0.0, 1.0, 0.005)]
                 for x in frange(0.0, 1.0, 0.005)]

# print_mat(Mat)

f = open("result", 'w')
f.write(str(Mat))
f.write("\n")
f.close()



