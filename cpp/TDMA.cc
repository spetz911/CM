/**
 *  Tridiagonal matrix algorithm (TDMA) for SLAE.
 * 
 * @author Oleg Baskakov
 * 2012. Written by NatSys Lab. (info@natsys-lab.com).
 */
#ifndef __TDMA_CC__
#define __TDMA_CC__

#include <assert.h>

#include "vector.cc"

#define INVARIANT(cond)	assert(cond)
#define RANGE(i, a, b) for (int i = a; i < b; ++i)
#define DOWNRANGE(i, a, b) for (int i = a; i > b; --i)


#define ABS(A) ((A)<0?(A):-(A))




class TDM {
private:
	Vector<double> a;
	Vector<double> b;
	Vector<double> c;
	Vector<double> d;
public:
	TDM(size_t size)
		:a(Vector<double>(size)),
		 b(Vector<double>(size)),
		 c(Vector<double>(size)),
		 d(Vector<double>(size))
	{}
	
	TDM(Vector<double> &a, Vector<double> &b, Vector<double> &c, Vector<double> &d)
		:a(Vector<double>(a)),
		 b(Vector<double>(b)),
		 c(Vector<double>(c)),
		 d(Vector<double>(d))
	{}
	
	~TDM() {}

	double & get_a(size_t i) { return a[i];}
	double & get_b(size_t i) { return b[i];}
	double & get_c(size_t i) { return c[i];}
	double & get_d(size_t i) { return d[i];}

	inline Vector<double> &
	solve()
	{
		TDM::solve(a,b,c,d);
	}

	// Thomas algorithm for Tridiagonal matrix
	static Vector<double> &
	solve(Vector<double> &a, Vector<double> &b, Vector<double> &c, Vector<double> &d)
	{
		size_t n = d.size();
		Vector<double> &x = *(new Vector<double>(n));
	
		for (int i = 1; i < n; i++)
		{
			b[i] -= a[i]/b[i-1] * c[i-1];
			d[i] -= a[i]/b[i-1] * d[i-1];
		}

		x[n-1] = d[n-1]/b[n-1];
	
		for (int i = n - 2; i >= 0; i--)
			x[i] = (d[i] - c[i]*x[i+1]) / b[i];
		return x;
	}


};
#endif // __TDMA_CC__

