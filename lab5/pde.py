#! /usr/bin/env python3
# -*- coding: utf-8 -*-

interface PDE_Interface:
	def __init__():
		"TODO description"
	def first_eq(u):
		"TODO description"
	def middle_eq(u):
		"TODO description"
	def last_eq(u):
		"TODO description"
	

class PDE:
	def scalar(v1, v2):
		return sum(x1*x2 for (x1,x2) in zip(v1,v2))
	
	def mat_vec(m, v):
		return [scalar(vx, v) for vx in m]

	def vec_mat(v, m):
		return mat_vec(zip(*m), v)
	
	## Description: solve equation
	def explicit_method(u, u_k_1 = 0): #TODO add threads
		sigma = tau * alpha**2 / h**2
		omega = tau * b / (2*h)
		eta = c * tau
		# move it to inheritance class
		coef_a = [+1, -2, +1]
		coef_b = [-3, +4, -1]
		coef_c = [+0, +1, +0]
		
		coef_mat = vec_mat([sigma, omega, eta], [coef_a, coef_b, coef_c])
		
		res = [0] + [mat_vec(coef_mat, u[i-1:(i+1)+1]) for i in range(1, N-1)] + [0]
		
		(a0, b0, c0, d0) = first_eq(u)
		(an, bn, cn, dn) = last_eq(u)
		res[0]  = (d0 - c0*res[1])  / b0
		res[-1] = (dn - bn*res[-2]) / an
	
		return res

	## Description: solve with method 'Progonki'
	def implicit_method():
		N = len(u)
	
		M = Tridiagonal_Matrix();
	
		M = zip(* first_eq(pde) + middle_eq(pde) * (N-2) + last_eq(pde) )
	

		M.n = N
	
		print_mat([M.a, M.b, M.c, M.d])

		x = M.solve()
		print_vec( [u0] + x + [u1])

	## Description: solve with method 'Progonki'
	def Crank_Nicolson_method():
	
		# a_i /= Q 
		# b_i /= Q 
		# c_i /= Q 
		# d_i += explicit_method * (1-Q) 
	





class Parabolic_PDE(pde) implements PDE_Interface:
	def __init__(pde = dict()):
		pde.__init__(pde)
		
	
	def first_eq(pde):
		"""Find coefficients of first equation"""
#		alpha = pde['init'][0]
#		beta  = pde['init'][1]
#		phi1  = pde['init'][2]
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
	
		a0 = 0
		b0 = self.alpha * (2*a*a/h + h/tau - c*h) - self.beta * (2*a*a - b*h)
		c0 = self.alpha * (2*a*a/h)
		d0 = self.alpha * (u[0][k] * h/tau) - self.phi1(tau*(k+1)) * (2*a*a - b*h)
		return (a0, b0, c0, d0)
	
	def last_eq(u):
		"""Find coefficients of last equation"""
#		alpha = pde['init1'][0]
#		beta  = pde['init1'][1]
#		phi1  = pde['init1'][2]
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']

		an = self.alpha * (-2*a*a/h)
		bn = self.alpha * (2*a*a/h + h/tau - c*h) + self.beta * (2*a*a + b*h)
		cn = 0
		dn = self.alpha * (u[N][k] * h/tau) + self.phi1(tau*(k+1)) * (2*a*a + b*h)
		return (an, bn, cn, dn)

	## Description: find coefficients of last equation
	def middle_eq(u):
		"""Find coefficients of middle equation"""
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
	
		sigma = a*a * tau / (h*h)
		omega = b * tau / (2*h)
		eta   = c*tau #TODO add f(x,t) here
	
		ai = sigma - 3*omega
		bi = 1 - 2*sigma + 4*omega + eta
		ci = sigma - omega
		di = -u[N][k]

		return (ai, bi, ci, di)

class Hyperbolic_PDE(pde) implements PDE_Interface:
	def __init__():
		pde.__init__()
		pass
	
	## Description: find coefficients of last equation
	def middle_eq(u):
		"""Find coefficients of middle equation"""
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
	
		sigma = tau**2 * alpha**2 / h**2
		omega = tau**2 * b / (2*h)
		eta = c * tau**2 + 2 # TODO add f(x,t)

		ai = sigma - 3*omega
		bi = -2*sigma + 4*omega + eta + 1
		ci = sigma - omega
		di = -u[N][k]

		return (ai, bi, ci, di)
	def last_eq(u):
		"""Find coefficients of last equation"""
		pass

	## Description: find coefficients of last equation
	def middle_eq(u):
		"""Find coefficients of middle equation"""
		pass



def main():

	s = open("input").readline()
	
	print(parse_all(s))


#=====================================================================

if __name__ == "__main__":
	main()

