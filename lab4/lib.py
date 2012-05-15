#!/usr/bin/env python
# -*- coding: utf-8 -*-
from copy import deepcopy

def frange(start,end,step):
	return list(map(lambda x: x*step, range(int(start*1.0/step),int(end*1.0/step))))


class ODY:
	def __init__(self,f,h,x0,y0,a,b):
		self.f = f
		self.h = h
		self.y0 = y0
		self.x0 = x0
		self.a = a
		self.b = b
		self.tabl = []
		
		#~ self.eps = eps
	
	def euler(self):
		h = self.h
		x0 = self.x0
		y0 = self.y0
		f = self.f
		a = self.a
		b = self.b
		self.tabl = [[x0,y0]]
		N = abs(b-a)/h
		
		for k in xrange(int(N)):
			app = []
			for m in xrange(len(self.tabl[k][1])):
				app.append(self.tabl[k][1][m] + h*(f(self.tabl[k][0],self.tabl[k][1]))[m])
			self.tabl.append([self.tabl[k][0]+h,app])
		return self.tabl
	
	def	Rynge_Kytte(self,itr = -1):
		def mult(a,x):
			c =deepcopy(x)
			for i in xrange(len(x)):
				c[i] = a*x[i]
			return c
		def summ(x,y):
			c =deepcopy(x)
			for j in xrange(len(x)):
				c[j] = x[j]+y[j]
			return c
		
		
		h = self.h
		x0 = self.x0
		y0 = self.y0
		f = self.f
		a = self.a
		b = self.b
		self.tabl = [[x0,y0]]
		N = [abs(b-a)/h,itr][itr!=-1]
		K = [0,0,0,0]
		dy = 0.
		yk = y0
		xk = x0
		self.tabl = [[x0,y0]]
		#~ self.fs = [0]*(int(abs(b-a)/h))
		
		
		for k in xrange(int(N)):
			#~ self.fs.append(f(xk,yk))
			#~ self.fs[k] = f(xk,yk)
			K[0] = mult(h,f(xk,yk))
			K[1] = mult(h,f(xk+h/2.,summ(yk,mult(.5,K[0]))))
			K[2] = mult(h,f(xk+h/2.,summ(yk,mult(.5,K[1]))))
			K[3] = mult(h,f(xk+h,summ(yk,K[2])))
			dy = mult(1./6.,summ(summ(K[0],mult(2.,K[1])),summ(mult(2.,K[2]),K[3])))
			yk = deepcopy(summ(yk,dy))
			
			xk += h
			self.tabl.append([xk,yk])
		return self.tabl
	
	def Addams(self):
		self.tabl = []
		if len(self.tabl) < 4:
			self.Rynge_Kytte(3)
		h = self.h
		f = self.f
		a = self.a
		b = self.b
		x0 = self.x0
		y = deepcopy(self.tabl[3][1])
		N = abs(b-a)/h
		#~ print len(self.tabl)
		for k in xrange(3,int(N)):
			summ = 0
			for m in xrange(len(y)):
				#~ print k
				#~ print self.tabl[k-2][1],self.tabl[k-3][1]
				#~ print self.fs[k][m]
				#~ y[m] += h*(55*self.fs[k][m] - 59*self.fs[k-1][m] + \
				#~ 37*self.fs[k-2][m] - 9*self.fs[k-3][m])/24.
				y[m] += h*(55*f(x0+h*k,self.tabl[k][1])[m] - 59*f(x0+h*(k-1),self.tabl[k-1][1])[m] + \
				37*f(x0+h*(k-2),self.tabl[k-2][1])[m] - 9*f(x0+h*(k-3),self.tabl[k-3][1])[m])/24.
			#~ print k
			#~ print len(self.fs)
				#~ summ +=55.*self.tabl[k][1][m]
				#~ print self.tabl[k-1][1]
				#~ summ -= 59.*self.tabl[k-1][1][m]
				#~ summ += 37.*self.tabl[k-2][1][m]
				#~ summ -= 9.*self.tabl[k-3][1][m]
			#~ 
				#~ y[m] += h*summ/24.
			#~ try:
				#~ self.tabl[k+1] = [x0+h*(k+1),y]
			#~ except IndexError:
			self.tabl.append([x0+h*(k+1),deepcopy(y)])
			#~ print self.tabl[k][1]
		return self.tabl
	def Rynge_Romberg(self,s,p):
		res = []
		for i in xrange(len(self.tabl)):
			print self.tabl[i][1][0],s.tabl[i][1][0]
			res.append((self.tabl[i][1][0] - s.tabl[i][1][0])/(2**p-1))
		return res

#~ class Kraev_solve:
	#~ def __init__(self,x0,koef):
		#~ #koefs =   alpha  betta y0^
		#~ #          delta  gamma y1^
		#~ self.koef = koef
		#~ self.x0 = x0
		#~ pass
	#~ def




