#include <stdlib.h> // for exit(1)

#include "vector.cc"

using namespace std;


class Matrix {
private:
	Vector *data;
	int n, m;
public:
	Matrix() : data(NULL), n(0) {}
	Matrix(int size) : data(new Vector[size]), n(size) {}
	Matrix(int _n, int _m) : data(new Vector[_n]), n(_n), m(_m)
	{
		for (int i=0; i<n; ++i)
			data[i] = *(new Vector(m));
	}
	~Matrix()
	{
		if (data) {
			delete[] data;
			data = NULL;
			n = 0;
			m = 0;
		}
	}
	
	int
	size_n() const
	{
		return n;
	}
	int
	size_m() const
	{
		return m;
	}
	
	Vector &
	resize(int _n, int _m)
	{
		n = _n;
		m = _m;
		if (data) delete[] data;
		data = new Vector[_n];
		
		for (int i=0; i<n; ++i)
			data[i].resize(m);
		return *this;
	}
	
	Vector &
	operator[](int i)
	{
		return data[i];
	}
	const Vector &
	operator[](int i) const
	{
		return data[i];
	}
	
	Matrix &
	operator=(const Matrix &old)
	{
		n = old.size_n();
		m = old.size_m();
		resize(n, m);
		
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				data[i][j] = old[i][j]; //FIXME operator= doesn't work
		return *this;
	}
	
	Matrix &
	operator+(const Matrix &v2) const
	{
		Matrix *res = new Matrix(n, m);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
			(*res)[i][j] = data[i][j] + v2[i][j];
		return *res;
	}

	Vector &
	operator+=(const Matrix &v2) const
	{
		for (int i=0; i<n; ++i)
			data[i] += v2[i];
		return *this;
	}

	Matrix &
	operator*(double a) const
	{
		Matrix *res = new Matrix(n, m);
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				res[i][j] = data[i][j] * a;
		return *res;
	}

	Vector &
	operator*(const Vector &vec) const
	{
		Vector *res = new Vector(vec.size());
		for (int i=0; i<n; ++i)
			res[i] = data[i] * vec;
		return *res;
	}
	
	void
	operator*=(double a)
	{
		for (int i=0; i<n; ++i)
			data[i] *= a;
		return *this;
	}

	double
	operator*(const Matrix &mat2) const
	{
		double res;
		for (int i=0; i<n; ++i)
			res += data[i] * mat2[i];
		return res;
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
		mat.resize(n,m);
		
		for (int i=0; i<n; ++i)
			for (int j=0; j<m; ++j)
				cin_ >> mat[i][j];
		return cin_;
	}
};


