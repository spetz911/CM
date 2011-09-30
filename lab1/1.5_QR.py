#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi
from copy import copy, deepcopy
from functools import reduce

def product(L):
	mul = lambda x, y: x*y
	return reduce(mul, L)

def sign(x):
	if(x<0):
		return -1
	elif(x>0):
		return 1
	else:
		return 0

#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi
from copy import copy, deepcopy
from functools import reduce

def product(L):
	mul = lambda x, y: x*y
	return reduce(mul, L)

def sign(x):
	if  (x<0): return -1
	elif(x>0): return  1
	else:	  return  0

class Matrix:
	def __init__(self,f = "",m = 0,n = 0):
		if f:
			self.M = None

		elif m and n:
			l = [0] * n
			self.M = [copy(l) for i in range(m)]
			self.m = m
			self.n = n
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
				print("\t %.2f"%x, end = ' ')
			print()
		print("-----")
		return

	def QRdecompositon(self):
		n = self.n
		E = Matrix("",self.m,self.n)
		for i in range(E.n):
			E[i][i] = 1
		Ak = deepcopy(self)
		Hes = []
		for k in range(n):
			vec = Matrix("",n,1)
			
			for j in range(k,n):
				vec[k][0] += Ak[j][k]**2

			vec[k][0] = sqrt(vec[k][0]) * sign(Ak[k][k]) + Ak[k][k] 
			
			for i in range(k+1,n):
				vec[i][0] = Ak[i][k]
				
			vec_t = vec.transponate()
			tp = (vec_t * vec)[0][0]
			tp2 = vec * vec_t

			H = E - tp2*(2.0/tp)
			
			Hes.append(H)
			Ak = H*Ak
		return [product(Hes),Ak]

class QR(Matrix):
	def __init__(self,f="",m=0,n=0):
		if f:
			inp = open(f,"r")
			s = inp.readline().split()
			self.eps = float(s[-1])

			
			self.M = []
			self.b = []
			self.p = 0
			for i in inp.readlines():
				s = [float(x) for x in i.split()]
				self.M.append(s)
						
			self.m = len(self.M)
			self.n = len(self.M[0])
		elif m and n:
			l = [0] * n
			self.M = [copy(l) for i in range(m)]
			self.m = m
			self.n = n
			
		else:
			self.M = []
			self.m = 0
			self.n = 0
		
		self.f = f
		
	def solve(self):
		eps = self.eps
		n = self.n
		m = self.m
		Ak = deepcopy(self)
		lambdas = []

		diag = [False]*n
		
		l1 = [0 + 0j]*n
		l2 = [0 + 0j]*n
		l1p = [0 + 0j]*n
		l2p = [0 + 0j]*n
		count = 0
		while 42 and (count < 9000) and (len(lambdas) < n):
			Q,R = Ak.QRdecompositon()
			Ak = R*Q
			
			for j in range(n):
				s = sum( [Ak[i][j]**2 for i in range(j+1,m)])
				
				if sqrt(s) < eps:
					if j < n-1  and not diag[j]:
						lambdas.append(Ak[j][j] + 0j)
						diag[j] = True
					elif j == n-1 and len(lambdas) == n-1 and not diag[j]:
						lambdas.append(Ak[j][j] + 0j)
						diag[j] = True
				else:
					try: # if j+1 != n
						s-= Ak[j+1][j]**2
					except:
						continue
					if sqrt(s) < eps:
						det = (Ak[j][j]-Ak[j+1][j+1])**2+4*Ak[j+1][j]*Ak[j][j+1]
						b = Ak[j][j] + Ak[j+1][j+1]
						if (det < 0):
							det = sqrt(abs(det))
							l1[j] = b/2.0 + 0.5j*det
							l2[j] = b/2.0 - 0.5j*det
						else:
							det = sqrt(det)
							l1[j] = (b + det)/2.0 + 0j
							l2[j] = (b - det)/2.0 + 0j
						if abs(l1[j]-l1p[j]) < eps and abs(l2[j]-l2p[j]) < eps and not diag[j] and not diag[j+1]:
							lambdas.append(l1[j])
							lambdas.append(l2[j])
							diag[j] = True
							diag[j+1] = True
						else:
							l1p = l1
							l2p = l2
				count += 1
		return lambdas


def round_c(c, n = 2):
	im = c.imag
	re = c.real
	return round(re, n) + 1j*round(im, n)


def main():
	m = QR("lab1_5.in")
	print ("eps = ", m.eps)
	lam = m.solve()
	for xx in lam:
		print(round_c(xx))

if (__name__ == "__main__"):
	main()
