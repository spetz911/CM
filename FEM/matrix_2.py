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
			self.b = []
			for i in range(min(n, m)):
				self.M[i][i] = 1
		
		else:
			self.M = []
			self.m = 0
			self.n = 0

	def __getitem__(self,i):
		return self.M[i]
		
	def __setitem__(self,i,z):
		self.M[i] = z

	def __len__(self):
		return len(self.M)

	def __neg__(self):
		M1 = self
		res = Matrix("",self.m,self.n)
		res.M = [ [(-M1_ij) for M1_ij in M1_i]
							for M1_i  in M1.M]
		return res

	def __add__(self, M2):
		M1 = self
		res = Matrix("",self.m,self.n)
		res.M = [ [(M1_ij + M2_ij)  for (M1_ij, M2_ij) in zip(M1_i,M2_i)]
									for (M1_i, M2_i)   in zip(M1.M,M2.M)]
		return res
		
	def __sub__(self, M2):
		M1 = self
		res = Matrix("",self.m,self.n)
		res.M = [ [(M1_ij - M2_ij)  for (M1_ij, M2_ij) in zip(M1_i,M2_i)]
									for (M1_i, M2_i)   in zip(M1.M,M2.M)]
		return res

	def __mul__(self, M2):
		M1 = self
		if type(M2) in [int,float]:
			res = Matrix("",M1.m, M1.n)
			res.M = [ [(M2*xx) for xx in M1_i] for M1_i in M1]
			return res

		res = Matrix(m = M1.m, n = M2.n)
		res.M = [[ sum([ (M1[i][k] * M2[k][j])  for k in range(M1.n)])
												for j in range(res.n)]
												for i in range(res.m)]
		return res

	def transponate(self):
		res = Matrix("", m = self.n, n = self.m)
		res.M = [ [(self[j][i]) for j in range(res.n)]
								for i in range(res.m)]
		return res

	def pr(self):
		print('Matrix ' + str(self.m)+ 'x' + str(self.n) + " :")
		for xx in self.M :
			for x in xx:
				print(" %.2f "%x, end = '\t')
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
		A = Matrix("", n, n)
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

		return A.transponate()

	def shift_b(self, b):
		v = [0]*self.n
		for i in range(self.n):
			
			v[i] = b[self.p[i]]
		return v
		
		
def main():
	m = Matrix("lab1_1.in")
	print ("Input matrix:")
	m.pr()
	print ("U:")
	m.build_LU()
	m.U.pr()
	print ("L:")
	m.L.pr()
	print ("L*U:")
	(m.L*m.U).pr()
	m1 = m.inverse()
	print ("Solve:")	
	x = m.solve(m.shift_b(m.b))
	x = [round(xx, 2) for xx in x]
	print("x = ", x)
		
	print ("Determinant:")
	print (m.det())
	print ("A^-1:")
	m1.pr()
	print ("A*A^-1:")
	(m*m1).pr()
