#! /usr/bin/python3.1
# -*- coding: utf-8 -*-


## 

class pde:
	## Description: solve equation
	def explicit_method( u, step): #TODO add threads
		N = len(u)
		sigma = tau * alpha**2 / h**2
		omega = tau * b / (2*h)
		eta = c * tau
		(a0, b0, c0, d0) = first_eq(pde)
		(an, bn, cn, dn) = last_eq(pde)
	
		res = [0] + [sigma*(u[i-1] - 2*u[i] + u[i+1]) +
			         omega*(-3*u[i-1] + 4*u[i] - u[i+1]) +
			         eta * u[i]   for i in range(1, N-1)] + [0]
		res[0] = (d0 - c0*res[1]) / b0
		res[N-1] = (dn - bn*res[N-2]) / an
	
		return res

	## Description: solve with method 'Progonki'
	def implicit_method():
		N = len(u)
		sigma = tau * alpha**2 / h**2
		omega = tau * b / (2*h)
		eta = c * tau
	
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
	__init__():
		pass
	
	




def main():

	s = open("input").readline()
	
	print(parse_all(s))


#=====================================================================

if __name__ == "__main__":
	main()

