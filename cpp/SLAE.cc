/**
 * Solver for System of linear equations.
 * 
 * @author Oleg Baskakov
 * 2012. Written by NatSys Lab. (info@natsys-lab.com).
 */
#ifndef __SLAE_CC__
#define __SLAE_CC__

#include <assert.h>
#include <vector>

#include "vector.cc"
#include "matrix.cc"

#define INVARIANT(cond)	assert(cond)
#define RANGE(i, a, b) for (int i = a; i < b; ++i)
#define DOWNRANGE(i, a, b) for (int i = a; i > b; --i)


#define ABS(A) ((A)<0?(A):-(A))


class SLAE {
private:
	Matrix U;
	Matrix L;
	std::vector<int> p;
	bool flag;
public:
	SLAE(Matrix &A) : U(Matrix(A)) {}
	~SLAE() {}

	void
	LU_decomposition()
	{
		INVARIANT(U.size_n() == U.size_m());
		L = Matrix(U.size(), U.size());
		p = std::vector<int>(U.size());
		LU_decomposition(U, L, p);
		flag = true;
	}
	
	static int
	LU_decomposition(Matrix &U, Matrix &L, std::vector<int> &p)
	{
		int p_count = 0;
		int n = U.size_n();
		int m = U.size_m();
		// set permutation vector
		RANGE(i, 0, n) {
			p[i] = i;
		}
		// for_each row
		RANGE(k, 0, m-1) {
			// arg_max(column)
			int s = k;
			RANGE(j, k+1, n) {
				if (ABS(U[s][k]) < ABS(U[j][k]))
					s = j;
			}
			// swap when s != k
			if (s != k) {
				U.swap_rows(k, s);
				L.swap_rows(k, s);
				L.swap_cols(k, s);
			
				std::swap(p[s], p[k]);
				
				p_count += 1;
			}
			// Gauss method one pass
			RANGE(i, k+1, n) {
				L[i][k] = U[i][k] / U[k][k];
				U[i] += U[k] * (-L[i][k]);
				//FIXME don't create new vector
			}
		}
		return p_count;
	}

	Vector<double> &
	solve(Vector<double> &b)
	{
		INVARIANT(U.size() == b.size());
		INVARIANT(flag);
		int n = b.size();
		Vector<double> &x = *(new Vector<double>(b));
		
		// forward Gauss pass
		RANGE(i, 0, n) {
			RANGE(j, 0, i)
				x[i] -= L[i][j] * x[j];
		}
		
		// reverse pass
		DOWNRANGE(i, n-1, -1) {
			RANGE(j, i+1, n) {
				x[i] -= U[i][j]*x[j];
			}
			x[i] /= U[i][i];
		}
		return x;
	}
	
	Matrix &
	inverse()
	{
		Matrix &res = *(new Matrix(U.size()));
		Matrix &E = Matrix::identity(U.size());
		int n = U.size();
		
		RANGE(i, 0, n)
			res[p[i]] = solve(E[i]);
		
		return res.transpose();
	}


};
#endif // __SLAE_CC__

