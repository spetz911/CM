#! /usr/bin/python3
# -*- coding: utf-8 -*-
###-------------------------------------------------------------------
### File    : pde.py
### Author  : Oleg Baskakov
### Description : Partial Differencial Equation
###
### 2011. Written for Moscow Aviation Institute.
###-------------------------------------------------------------------

import numpy
import numpy.linalg as la

import pyopencl as cl

import Tkinter as tk
import Image          # PIL
import ImageTk        # PIL



f = open("result")
Mat = eval(f.read())
f.close()

n = len(Mat)
m = len(Mat[0])

# n = 100
# m = 100

# Mat = [[i+j for j in range(m)]
#             for i in range(n)]


coeff = max(max(el for el in row)
                   for row in Mat)

buffer = numpy.zeros(n*m*4, numpy.uint8)

for i in range(n):
	for j in range(m):
		num = int(Mat[i][j] * 255 / coeff)
		rest = int((Mat[i][j] * 255 / coeff  - num) * 3)
		kk = [0] * 3
		for tt in range(rest):
			kk[tt] = 1
		buffer[4*(i*m + j) + 0] = num + kk[0]
		buffer[4*(i*m + j) + 1] = num + kk[1]
		buffer[4*(i*m + j) + 2] = num + kk[2]
		buffer[4*(i*m + j) + 3] = 255

im = Image.fromstring("RGBA", (n,m), buffer.tostring())
im.save("out.jpg")

