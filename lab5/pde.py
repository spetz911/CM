#! /usr/bin/python3.1
# -*- coding: utf-8 -*-




class pde:
	def __init__():
		pass

	def scalar(v1, v2):
		return sum(x1*x2 for (x1,x2) in zip(v1,v2))
	
	def mat_vec(m, v):
		return [scalar(vx, v) for vx in m]

	def vec_mat(v, m):
		return mat_vec(zip(*m), v)
	
	## Description: solve equation
	def explicit_method(u): #TODO add threads
		sigma = tau * alpha**2 / h**2
		omega = tau * b / (2*h)
		eta = c * tau
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
	
		M = zip(* first_eq(pde) + midle_eq(pde) * (N-2) + last_eq(pde) )
	
		# a_i /= Q 
		# b_i /= Q 
		# c_i /= Q 
		# d_i += explicit_method * (1-Q) 
	
		M.n = N
	
		print_mat([M.a, M.b, M.c, M.d])

		x = M.solve()
		print_vec( [u0] + x + [u1])
	
	





class parabolic_pde(pde):
	def __init__():
		pde.__init__()
		pass
	
	## Description: find coefficients of first equation
	def first_eq(pde):
		alpha = pde['init'][0]
		beta  = pde['init'][1]
		phi1  = pde['init'][2]
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
	
		a0 = 0
		b0 = alpha * (2*a*a/h + h/tau - c*h) - beta * (2*a*a - b*h)
		c0 = alpha * (2*a*a/h)
		d0 = alpha * (u[0][k] * h/tau) - phi1(tau*(k+1)) * (2*a*a - b*h)
		return (a0, b0, c0, d0)
	
	## Description: find coefficients of last equation
	def last_eq(u):
		alpha = pde['init1'][0]
		beta  = pde['init1'][1]
		phi1  = pde['init1'][2]
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']

		an = alpha * (-2*a*a/h)
		bn = alpha * (2*a*a/h + h/tau - c*h) + beta * (2*a*a + b*h)
		cn = 0
		dn = alpha * (u[N][k] * h/tau) + phi2(tau*(k+1)) * (2*a*a + b*h)
		return (an, bn, cn, dn)

	## Description: find coefficients of last equation
	def midle_eq(pde):
		a = pde['u_xx']
		b = pde['u_x']
		c = pde['u']
		u = pde['grid']
	
		sigma = a*a * tau / (h*h)
		omega = b * tau / (2*h)
		eta   = c*tau #TODO add f(x,t) here
	
		ai = sigma - 3*omega
		bi = 1 - 2*sigma + 4*omega + eta
		ci = sigma - omega
		di = -u[N][k]

		return (ai, bi, ci, di)

	




def main():

	s = open("input").readline()
	
	print(parse_all(s))


#=====================================================================

if __name__ == "__main__":
	main()

