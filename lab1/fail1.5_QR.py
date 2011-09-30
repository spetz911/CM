#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from math import sqrt,sin,cos,atan,pi
from copy import copy, deepcopy
from functools import reduce

def sign(x):
	if(x<0):
		return -1
	elif(x>0):
		return 1
	else:
		return 0

class Matrix:
	def __init__(self,f = "",m = 0,n = 0):
		self.nrm = 0
		if f:
			inp = open(f,"r")
			self.M = []
			self.b = []
			self.p = 0
			for i in inp.readlines():
				try:
					s = map(float,i.split())
					self.M.append(s)
				except:
					continue
						
			for i in range(len(self.M)):
				self.b.append(self.M[i].pop())
			
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

	def __neg__(self):
		M1 = self
		res = Matrix("",self.m,self.n)
		res.M = [ [(-M1_ij) for M1_ij in M1_i]
							for M1_i  in M1.M]
		return res
	
	def __len__(self):
		return len(self.M)

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
				print("\t %.2f"%x, end = ' ')
			print()
		print("-----")
		return
	
	def QR_decompositon(self):
		n = self.n
		E = Matrix("",self.m,self.n)
		for i in range(E.n):
			E.M[i][i] = 1
		Ak = deepcopy(self)
		Hes = []
		for k in range(n):
			v = Matrix("",n,1)
			for j in range(k,n):
				v[k][0]+=Ak[j][k]**2
				#~ print v[k][0]
			v[k][0] = sqrt(v[k][0])
			v[k][0] *= sign(Ak[k][k])
			v[k][0] += Ak[k][k] 
			#~ print v[k][0]
			
			for i in range(k+1,n):
				v[i][0] = Ak[i][k]
				
			vt = v.transponate()
			tp = (vt*v)[0][0]
			#~ print tp
			#~ print
			tp2 = v*vt
			
			#~ tp2.pr()
		
			H = E - tp2*(2.0 / tp)
			
			Hes.append(H)
			Ak = H*Ak
		return [reduce(lambda x,y: x*y,Hes),Ak]


def roundx(v, n):
	return [ round(xx, n) for xx in v]

class QR(Matrix):
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
		eps = self.eps
		n = self.n
		m = self.m
		Ak = self
		lambdas = []
		d1 = 0
		d2 = 0
		d2m = 0
		d1m = 0
		diag = [0]*n
		l1 = [[0.,0.]]*n
		l2 = [[0.,0.]]*n
		l1p = [[0.,0.]]*n
		l2p = [[0.,0.]]*n
		while 1:
			Q,R = Ak.QR_decompositon()
			Ak = R*Q
			test = 0
			for j in range(n):
				
				s = 0.
				for i in range(j+1,m):
					s += Ak[i][j]**2
				s = sqrt(s)
				
				if s < eps:
					if j < n-1  and not diag[j]:
						lambdas.append([Ak[j][j],0.])
						diag[j] = 1
					elif j == n-1 and len(lambdas) == n-1  and not diag[j]:
						lambdas.append([Ak[j][j],0.])	
						diag[j] = 1
				else:
					s = s**2
					try:
						s-= Ak[j+1][j]**2
					except:
						continue
					
					s = sqrt(s)
					if s < eps:
						a1 = Ak[j][j]
						a2 = Ak[j+1][j+1]
						a3 = Ak[j+1][j]
						a4 = Ak[j][j+1]
						try:
							l1[j] = [(a1+a2)/2.+sqrt((a1**2+a2**2 - 2*a1*a2+4*a3*a4))/2.,0.]
							l2[j] = [(a1+a2)/2.-sqrt((a1**2+a2**2 - 2*a1*a2+4*a3*a4))/2.,0.]
						except:
							l1[j] = [(a1+a2)/2.,sqrt(-(a1**2+a2**2 - 2*a1*a2+4*a3*a4))/2.]
							l2[j] = [(a1+a2)/2.,-sqrt(-(a1**2+a2**2 - 2*a1*a2+4*a3*a4))/2.]
						
						if abs(l1[j][0] - l1p[j][0]) < eps and abs(l2[j][0] - l2p[j][0]) < eps and not diag[j] and not diag[j+1]:
							lambdas.append(l1[j])
							lambdas.append(l2[j])
							diag[j] = 1
							diag[j+1] = 1
						else:
							l1p = deepcopy(l1)
							l2p = deepcopy(l2)
			if len(lambdas) == n:
				break
		return lambdas


q = QR("lab1_5.in")

l, v = q.solve()

for (l_i, v_i) in zip(l, v):
	print("  l = \t%.2f\t"%l_i, roundx(v_i, 3))

