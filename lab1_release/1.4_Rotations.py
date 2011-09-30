#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi
from copy import copy, deepcopy

class Matrix:
	def __init__(self,f = "",m = 0,n = 0):
		self.nrm = 0
		if m and n:
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

	def __neg__(self):
		res = Matrix("",self.m,self.n)
		for i in xrange(self.m):
			for j in xrange(self.n):
				res[i][j] = -self.M[i][j]

	def __len__(self):
		return len(self.M)

	def __add__(self,y):
		res = Matrix("",self.m,self.n)
		for i in range(self.m):
			for j in range(self.n):
				res[i][j] = self.M[i][j] + y[i][j]
		return res

	def __sub__(self,y):
		res = Matrix("",self.m,self.n)
		
		for i in range(self.m):
			for j in range(self.n):
				res[i][j] = self.M[i][j] - y.M[i][j]
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
				
	def norm(self):
		tp = 0.
		mx = [0.,0.]
		for i in range(self.m):
			s = 0.
			for j in range(self.n):
				s+=self.M[i][j]
				#~ print self.M[i][j],
			if mx[0] < s:
				mx = [s,i]
			#~ print s
		
		self.nrm = mx[0]
		return self.nrm	

	def transponate(self):
		res = Matrix("", m = self.n, n = self.m)
		for i in range(self.n):
			for j in range(self.m):
				res[i][j] = self[j][i]
		return res

	def pr(self):
		print('Matrix ' + str(self.m)+ 'x' + str(self.n) + " :")
		for xx in self.M :
			for x in xx:
				print("\t%.2f"%x, end = ' ')
			print()
		print("-----")
		return
	
		
class Rotations(Matrix):
	def __init__(self,f="",m=0,n=0):
		if f:
			input1 = open(f, 'r')
			self.M = []
			
			self.b = []
			self.p = 0

			self.eps = float( input1.readline().split(' ').pop() )
			
			for str1 in input1.readlines():
				try:
					z = [float(x) for x in str1.split(' ')]
					self.M.append(z)
				except:
					continue
			
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
	
	def max_nd(self):
		a_m = 0
		i_m = 0
		j_m = 0
		for i in range(self.n):
			for j in range(i):
				if i!=j:
					if (abs(self[i][j])>abs(a_m)):
						a_m = abs(self[i][j])
						i_m = i
						j_m = j
		return (a_m, j_m, i_m)


	def solve(self):
		E  = Matrix("",self.m,self.n)
		Uk = Matrix("",self.m,self.n)
		U_es = []
		
		A = Matrix("",self.m,self.n)
		A.M = deepcopy(self.M)
		Akp = A
		eps_k = 1
		counter = 0
		
		while abs(eps_k) > self.eps:
			counter+=1
			if counter > 9000: break 
			Ak = Akp
			(a, i, j) = Rotations.max_nd(Ak)
			
			if (Ak[i][i] - Ak[j][j]):
				phi = 0.5 * atan(2*Ak[i][j]/(Ak[i][i] - Ak[j][j]))
			else:
				phi = pi/4.
			Uk = deepcopy(E)
			Uk[i][i] = cos(phi)
			Uk[j][j] = cos(phi)
			Uk[i][j] = -sin(phi)
			Uk[j][i] = sin(phi)
			U_es.append(Uk)
			Akp = Uk.transponate()*Ak*Uk
			eps_k = 0.
			for m in range(self.n):
				for l in range(m):
					eps_k += Akp[l][m]*Akp[l][m]
			eps_k = sqrt(eps_k)
		#reduce manual
		for i in range(1, len(U_es)):
			U_es[0] = U_es[0]*U_es[i]
		tp = U_es[0].transponate()
		
		print("finish in %d iters."%counter)
		
		lambdas = []
		vectors = []
		for i in range(Akp.n):
			lambdas.append(Akp[i][i])
			vectors.append(tp[i])
		
		return [lambdas,vectors]


def main():            
	r = Rotations("lab1_4.in")
	
	l, v = r.solve()

	A = Matrix("", r.m, r.n)
	A.M = deepcopy(v)

	for i in range(len(v)):
		v[i] = [round(xx,3) for xx in v[i] ]

	for (l_i, v_i) in zip(l, v): print("  l = \t%.2f\t v = "%l_i, v_i)

	print("test A^t*R*A:")
	(A*r*A.transponate()).pr()


if (__name__ == "__main__"):
	main()
