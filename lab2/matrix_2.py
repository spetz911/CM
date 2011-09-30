#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from copy import copy, deepcopy

class Matrix:
	def __init__(self,f = "",m = 0,n = 0):
		if f:
			input1 = open(f, 'r')
			self.M = []
			
			self.b = []
			self.p = 0
			self.p_count = 0
			
			
			for str1 in input1.readlines():
				try:
					z = [float(x) for x in str1.split(' ')]
					
					self.b.append(z.pop()) #deleting b
					self.M.append(z)
				except:
					continue
			#print("LOLA", self.M)
			
			self.m = len(self.M)	
			self.n = len(self.M[0])
		elif m and n:
			l = [0] *  n
			self.M = [copy(l) for i in range(m)]
			self.m = m
			self.n = n
			for i in range(min(n, m)):
				self.M[i][i] = 1
		
		else:
			self.M = []
			self.m = 0
			self.n = 0

	def __getitem__(self,i):
		return self.M[i]
	
	def __setitem__(self,i,y):
		self.M[i] = y

	def __add__(self,y):
		res = Matrix("",self.m,self.n)
		for i in range(self.m):
			for j in range(self.n):
				res[i][j] = self.M[i][j] + y[i][j]
		return res

	def __mul__(self, M2):
		M1 = self
		result = Matrix(m = M1.m, n = M2.n)
		
		for i in range(result.m):
			for j in range(result.n):
				result[i][j] = 0
				for k in range(M1.n):
					result[i][j] += M1[i][k] * M2[k][j]
		return result
					
					
	def __sub__(self,y):
		res = Matrix("",self.m,self.n)
		
		for i in range(self.m):
			for j in range(self.n):
				res[i][j] = self.M[i][j] - y.M[i][j]
		return res


	def norm(self):
		tp = 0.
		mx = (0.0,0)
		for i in range(self.m):
			s = 0.
			for j in range(self.n):
				s+=self.M[i][j]
				#~ print self.M[i][j],
			if abs(mx[0]) < abs(s):
				mx = (s,i)
			#~ print s
		return mx[0]


	def transponate(self):
		res = Matrix("", m = self.n, n = self.m)
		for i in range(self.n):
			for j in range(self.m):
				res[i][j] = self[j][i]
		return res


	def __len__(self):
		return len(self.M)
	
	def pr(self):
		print('Matrix ' + str(self.m)+ 'x' + str(self.n) + " :")
		for xx in self.M :
			for x in xx:
				print("\t%.5f"%x, end = ' ')
			print()
		print("-----")
		return
	
	def det(self):
		"""Calculating det"""
		det = 1
		for i in range(self.n):
			det*=self.U[i][i]
		print("DEDEFED", self.p, self.p_count)
		if self.p_count %2:
			det = -det
		return det

	
#	def pr(self):
#		for i in range(self.m):
#			for j in range(self.n):
#				print("%.2f" % (self.M[i][j]))
#			print("")
	
		
	def swap_rows(self, k, s):
		"""Swap rows k and s"""
		z = self[k]
		self[k] = self[s]
		self[s] = z
		
	def swap_cols(self, k, s):
		"""Swap cols k and s"""
		for j in range(self.n):
			z = self[j][k]
			self[j][k] = self[j][s]
			self[j][s] = z
			

	def build_LU(self):
		n = self.n
		A = Matrix()
		A.M = deepcopy(self.M)
		A.n = self.n
		self.p = [x for x in range(n)]
		self.p_count = 0
#		print ("OLOLO", self.p)
		
		U = Matrix("", n, n)
		L = Matrix("", n, n)
		
		for k in range(n-1):
			
			# Find max
			s = k
			for j in range(k+1, n):
				if (abs(A[s][k])< abs(A[j][k])) : s = j
			
			#print("MAX = ", A[s][k], "k=", k, "s = ", s )
			
			A.swap_rows(k, s)
			L.swap_rows(k, s)
			L.swap_cols(k, s)
			
			# construct vector p
			z = self.p[k]   #equal k
			self.p[k] = self.p[s]
			self.p[s] = z
			if (k!=s): self.p_count += 1
			
			#s = max( [xx for xx in range(k, n)])
			
			
			for i in range(k+1, n): #col number
				koef = A[i][k]/A[k][k]
				A.M[i] = [(xx-xo*koef)for (xx,xo) in zip(A.M[i],A.M[k])]
				L[i][k] = koef
#				print(A.M[i])

		self.L = L
		self.U = A
		
	
	def solve(self,b):
		n = len(b)
		z = [0]*n
		z[0] = b[0]
		for i in range(1,n):
			z[i] = b[i]
			for j in range(i):
				z[i] -= self.L[i][j]*z[j]
		
		x = [0]*n
		x[n-1] = z[n-1]
		i = n-1
		while i > -1 :
			x[i] = z[i]
			for j in range(i+1,n):
				x[i] -= self.U[i][j]*x[j]
			x[i] /= self.U[i][i]
			i-=1
		
		return x
	
	def inverse(self):
		E = Matrix("",self.n,self.n)
		A = Matrix("",self.n, self.n)

		for i in range(self.n):
			e = [0]*self.n
			e[i] = 1.0
			A[self.p[i]] = self.solve(e)
		A.transponate()
		return A

	def shift_b(self, b):
		v = [0]*self.n
		for i in range(self.n):
			v[i] = b[self.p[i]]
		return v
		
