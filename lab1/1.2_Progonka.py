#! /usr/bin/python3.1
# -*- coding: utf-8 -*-
		
class Tridiagonal_Matrix:
	def __init__(self, f = ""):
		if f:
			inp = open(f, 'r')
			self.a = []
			self.b = []
			self.c = []
			self.d = []
			
			for line in inp.readlines():
				(a, b, c, d) = ( float(x) for x in line.split() )
			#	print("DEBUG:", a,b,c,d)
				self.a.append(a)
				self.b.append(b)
				self.c.append(c)
				self.d.append(d)
			self.n = len(self.a)	
			
		else:
			print ("File wasn't specified")
			exit(1)
			
		return None
	
	def solve(self):
		"""Method progonki"""
		a = self.a
		b = self.b
		c = self.c
		d = self.d
		n = self.n
		
		P = []
		Q = []
		P.append(-c[0]/b[0])
		Q.append( d[0]/b[0])
		
		for i in range(1, n):
			P.append( -c[i] / (b[i]+a[i]*P[i-1]) ) 
			Q.append( (d[i] - a[i]*Q[i-1]) / (b[i] + a[i]*P[i-1]) )

		x = [0]*n
		x[n-1] = Q[n-1]
		for i in range(n-2, -1, -1):
			x[i] = P[i]*x[i+1] + Q[i]
			
		return x
		
def main():
	e = Tridiagonal_Matrix("lab1_2.in")
	print("a =",e.a)
	print("b =",e.b)
	print("c =",e.c)
	print("d =",e.d)
	print()
	x = e.solve()
	x = [round(z,2) for z in x]
	print("x = ", x)
	return 

if (__name__ == "__main__"):
	main()
