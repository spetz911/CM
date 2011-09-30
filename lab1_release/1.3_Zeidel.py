#! /usr/bin/python3.1
# -*- coding: utf-8 -*-

from copy import copy, deepcopy

class Matrix:
	def __init__(self,f = "",m = 0,n = 0):
		self.nrm = 0
		if f:
			input1 = open(f, 'r')
			self.M = []
			self.b = []
			self.p = 0
			
			for str1 in input1.readlines():
				try:
					z = [float(x) for x in str1.split(' ')]
					self.b.append(z.pop()) #deleting b
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
		
class Iteration(Matrix):
	def __init__(self, f="", m=0, n=0):
		inp = open(f, 'r')
		line = inp.readline().split()
		self.eps = float(line[-1])
		Matrix.__init__(self, f, m, n)
	
	def solve(self):
		b = self.b
		A = self.M
		
		alpha = Matrix("",self.n,self.n)
		beta = Matrix("",self.n,1)
		
		for i in range(0, beta.m):
			beta[i][0] = b[i]/A[i][i]
		
		for i in range(0, alpha.n):
			for j in range(0, alpha.n):
				if i==j:
					alpha[i][j] = 0
				else:
					alpha[i][j] = - A[i][j] / A[i][i]
		xk  = beta 
		
		eps_k = 1
		norma = alpha.norm()
		counter = 0
		
		while (eps_k > self.eps):
			xk1 = beta + alpha*xk
			if norma < 1:
				eps_k = norma/(1-norma) * (xk1-xk).norm() 
			else:
				eps_k = (xk1-xk).norm()
			counter+=1
			if counter > 100 : break
			xk = xk1	 
			
		print ("Simple calculate %d iterations; " %counter)
		return xk
		
	def solve_zeidel(self):
		b = self.b
		a = self.M
		alpha = Matrix("",self.m,self.n)
		beta = Matrix("",self.n,1)
		for i in range(self.m):
			beta.M[i] = [b[i]/a[i][i]]
			
			for j in range(self.n):
				if i != j:
					alpha.M[i][j] = - a[i][j] / a[i][i]
				else:
					alpha.M[i][j] = 0		
		C = Matrix("",self.m,self.n)
				
		for i in range(self.m):
			for j in range(self.n):
				if i >= j:
					C.M[i][j] = alpha.M[i][j]
		
		E = Matrix("",self.m,self.n)
		for i in range(E.n):
			E.M[i][i] = 1
		
		ek = 1.
		count = 0
		xk = Matrix("",self.n,1)
		xkm = deepcopy(beta)
		while abs(ek) >= self.eps :
			tp = Matrix("",self.n,1)
			for i in range(self.m):
				tp[i][0] = 0
				for j in range(i):
					tp[i][0]+=alpha[i][j]*xk[j][0]
				for j in range(i,self.n):
					tp[i][0]+=alpha[i][j]*xkm[j][0]
				
				xk[i][0] = beta[i][0] + tp[i][0]
			if alpha.norm() < 1:
				ek = C.norm()/(1-alpha.norm()) * (xk-xkm).norm() 
			else:
				ek = (xk-xkm).norm() 
			count+=1
			xkm = deepcopy(xk)		
		print ("Zeidel finish in %d iteratrions;" %count)
		return xk
		
def main():
	m = Iteration("lab1_3.in")
	print("eps =", m.eps)
	#~ m.pr()
	m.solve().transponate().pr()
	z = m.solve_zeidel()
	z.transponate().pr()
	print("test A*x = b")
	(m*z).transponate().pr()
	
	
	print("b =",m.b)
	
	
	return 

if (__name__ == "__main__"):
	main()
