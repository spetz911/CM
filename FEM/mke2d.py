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
	u_k = []
	for el in elems:
		u_k.append([answ[x] for x in el])
		values.append([points[x] for x in el])
		
		koeffs = [[ answ[p]*x for x in elem_plane(i, values[-1])]
		                      for i,p in enumerate(el)]
		funcs3.append(koeffs)
	
	return lambda x,y: sum(elem_fun(x,y, koef, val) for koef, val in zip(funcs3, values))




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
phi = lambda x,y: x+y
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

#	if abs(k[0] + k[1]*points[i][0] + k[2]*points[i][1] - 1.0) > 0.001 : exit(0)

	return koef # (1.0, koef[1]/koef[0], koef[2]/koef[0])


def mke_solve(elems, lim_points, points):
	global alpha
	mat = matrix_2.Matrix(m = len(points), n = len(points))
	mat = mat * 0
	
	vec = matrix_2.Matrix(m = len(points), n = 1)
	vec = vec * 0

	al =   [0.0]*len(points)
	beta = [0.0]*len(points)
	gama = [0.0]*len(points)
	
	for num,el in enumerate(elems):
		print('matrix_el', el)
		values = [points[x] for x in el]
		el_square = square(values)
		f_start = func_start(*median(*values)) * el_square / 3
		
		for i,e in enumerate(el):
			plane = elem_plane(i, values)
#			print(plane)
			al[e] = plane[0]
			beta[e] = plane[1]
			gama[e] = plane[2]

#		print('alph =', al)
#		print('beta =', beta)
#		print('gama =', gama)

		for p1, p2 in dekart_product(el,el):
			mat[p1][p2] += (beta[p1]*beta[p2] + gama[p1]*gama[p2]) * el_square
			
		border = [ x for x in dekart_product(el,el) if x in lim_points]
		
		print('border', border)
		
		if border:
			i,k = border[0]
			[j] = [x for x in el if x!=i and x!=k]
			bord_dist = dist(points[i], points[k])
			bord_fun = 0.5 * bord_dist * phi(*median(points[i], points[k]))
			
			mat[i][i] += 1.0 * bord_dist * alpha/3
			mat[k][k] += 1.0 * bord_dist * alpha/3
			mat[i][k] += 0.5 * bord_dist * alpha/3
			mat[k][i] += 0.5 * bord_dist * alpha/3
			
			vec[i][0] += bord_fun
			vec[k][0] += bord_fun
			
			print ('f_start =', bord_fun)
	
		else:
			i,j,k = el
			bord_dist = 0.0
			bord_fun = 0.0
			
		vec[i][0] -= f_start
		vec[j][0] -= f_start
		vec[k][0] -= f_start
		
		

	return mat, vec

	


def matrix_el(elem, lim_points, points):
	print('matrix_el', elem)
	global alpha
	res = matrix_2.Matrix(m = 3, n = 3)

	values = [points[x] for x in elem]
#	print('values', values)
	k3 = [elem_plane(i, values) for i in range(len(values))]
#	print('k3 =', k3)
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

	if border:
		edge = border[0]
#		print("border")
		p1 = points[edge[0]]
		p2 = points[edge[1]]
		i = elem.index(edge[0])
		k = elem.index(edge[1])
		bord_dist = dist(p1, p2)
		bord_fun = 0.5 * bord_dist * phi(*median(p1, p2))
	else:
		bord_dist = 0.0
		bord_fun = 0.0
		i = k = 0

	mat = matrix_2.Matrix(m = 3, n = 3)
	mat.M = [[0.0] * 3 for i in range(3)]
	mat[i][i] = 1.0
	mat[k][k] = 1.0
	mat[i][k] = 0.5
	mat[k][i] = 0.5

	mat = mat * (alpha/3 * bord_dist)

	vec = matrix_2.Matrix(m = 3, n = 1)
	
	f_start = func_start(*median(*values)) * el_square / 3
	
	print ('f_start =', bord_fun)
	
	vec.M[0][0] = -f_start
	vec.M[1][0] = -f_start
	vec.M[2][0] = -f_start
	
	vec.M[i][0] += bord_fun
	vec.M[k][0] += bord_fun
	
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

#points = [(0,0), (0,2), (1,1), (2,1)]
# points = [(-1,-1), (0,0), (1,0), (0,1)]

#elems = [(0,1,2), (1,2,3), (2,3,0)]
#lim_points = [(0,1), (1,3), (3,0)]


# print('old', gen_lim_points())
print('new', lim_points)
print('elems', elems)
from pprint import pprint
pprint(list(enumerate(points)))


mat, vec = mke_solve(elems, lim_points, points)
#mat, vec = asamble(elems, lim_points, points)

mat.pr()
vec.transponate().pr()


mat.build_LU()

print ("Solve:")
answ = mat.solve(mat.shift_b(vec.transponate()[0]))

#print ("Test:")
#vec2 = matrix_2.Matrix(m = vec.n, n = vec.m)
#vec2.M[0] = answ
#vec2  = vec2.transponate()
#(mat*vec2).pr()

answ = [round(xx, 5) for xx in answ]
print("answ = ", answ)

etalon = [round(phi(x,y), 2) for x,y in points]
print("etalon = ", etalon)

# answ = [1.0] *len(elems)
# answ[1] = 1.0

# answ = [i/len(elems) for i in range(len(elems))]


Fun = test_answ(answ, elems, points)
Mat = [[Fun(x,y) for y in frange(0.0, 1.0, 0.05)]
                 for x in frange(0.0, 1.0, 0.05)]

# print_mat(Mat)

f = open("result", 'w')
f.write(str(Mat))
f.write("\n")
f.close()

