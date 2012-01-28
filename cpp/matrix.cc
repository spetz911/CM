/**
 * Solver for System of linear equations.
 * 
 * @author Oleg Baskakov
 * 2012. Written by NatSys Lab. (info@natsys-lab.com).
 */
#ifndef __MATRIX_CC__
#define __SLAE_CC__


#include <iostream>
#include <fstream>
#include <vector>

#include "vector.cc"




class Matrix : public Vector<Vector<double> > {
private:
public:
	Matrix()
		:Vector<Vector<double> >() {}
	Matrix(size_t n)
		:Vector<Vector<double> >(n) {}
	Matrix(size_t n, size_t m)
		:Vector<Vector<double> >(n)
	{
		for (int i=0; i<n; ++i)
			(*this)[i] = *(new Vector<double> (m));
	}
	inline
	Matrix(Vector<Vector<double> > mat)
		:Vector<Vector<double> >(mat) {}
	
	static Matrix &
	identity(size_t n)
	{
		Matrix &mat = *(new Matrix(n, n));
		for (int i=0; i<n; ++i)
			mat[i][i] = 1.0;
		return mat;
	}
	
	Matrix(const Vector<double> &vec)
		:Vector<Vector<double> >(vec.size())
	{
		for (int i=0; i<vec.size(); ++i) {
			(*this)[i] = *(new Vector<double> (1));
			(*this)[i][0] = vec[i];
		}
//		Matrix m(1);
//		m[0] = vec;
//		*this = transpose(m);
	}

	size_t inline
	size_n() const
	{
		return size();
	}
	
	size_t inline
	size_m() const
	{
		if (empty())
			return 0;
		else
			return (*this)[0].size();
	}
	
	void
	resize(size_t n, size_t m)
	{
		std::vector<Vector<double> >::resize(n);
		for (int i=0; i<n; ++i)
			(*this)[i].resize(m);
	}
	
	void
	swap_rows(int i, int j)
	{
		std::swap((*this)[i], (*this)[j]);
	}

	void
	swap_cols(int i, int j)
	{
		for (int k=0; k < size(); ++k)
			std::swap((*this)[k][i], (*this)[k][j]);
	}
	
	void
	operator*=(double a)
	{
		for (int i=0; i<size(); ++i)
			(*this)[i] *= a;
	}

	Matrix &
	operator*(double a) const
	{
		Matrix &res = *(new Matrix(size_n(), size_m()));
		for (int i=0; i < size_n(); ++i)
			for (int j=0; j < size_m(); ++j)
				res[i][j] = (*this)[i][j] * a;
		return res;
	}

	Vector<double>  &
	operator*(const Vector<double> &vec) const
	{
		const Matrix &mat = *this;
		Vector<double> &res = *(new Vector<double> (vec.size()));
		for (int i=0; i < size_n(); ++i)
			res[i] = mat[i] * vec;
		return res;
	}
	
	
	Matrix &
	operator*(const Matrix &mat0) const
	{
		const Matrix &mat1 = *this;
		const Matrix &mat2 = transpose(mat0);
		Matrix &res = *(new Matrix(mat1.size_n(), mat2.size_n()));
		
		for (int i=0; i < res.size_n(); ++i)
			for (int j=0; i < res.size_m(); ++i)
			res[i][j] = Vector<double>::dot(mat1[i], mat2[i]);
		return res;
	}

	Matrix &
	transpose()
	{
		Matrix &res = *this;
		for (int i=0; i < res.size_n(); ++i)
			for (int j=0; j < res.size_m(); ++j)
				std::swap(res[i][j], res[j][i]);
		return res;
	}

	static Matrix &
	transpose(const Matrix &mat)
	{
		Matrix &res = *(new Matrix(mat.size_m(), mat.size_n()));
		for (int i=0; i < res.size_n(); ++i)
			for (int j=0; j < res.size_m(); ++j)
				res[i][j] = mat[j][i];
		return res;
	}

	std::string &
	str()
	{
		std::string *res;
		for (int i=0; i<size_n(); ++i)
			*res += (*this)[i].str() + std::endl;
		return *res;
	}

	friend std::ostream& operator<<(std::ostream& cout_, const Matrix& mat)
	{
		int n = mat.size_n();
		int m = mat.size_m();
		
		for (int i=0; i<n; ++i) {
			for (int j=0; j<m; ++j)
				cout_ << mat[i][j] << "  ";
			cout_ << std::endl;
		}
		return cout_;
	}
	
	friend std::istream& operator>>(std::istream& cin_, Matrix& mat)
	{
		int n, m;
		cin_ >> n >> m;
		mat.resize(n, m);
		
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				cin_ >> mat[i][j];
		return cin_;
	}

};
#endif // __MATRIX_CC__


