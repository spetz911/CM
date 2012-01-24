#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h> // for exit(1)

// using namespace std;

class Vector {
private:
	double *data;
	int n;
public:
	Vector() : data(NULL){}
	Vector(double *ptr, int size) : data(ptr), n(size) {}
	Vector(int size) : data(new double[size]), n(size) {}
	~Vector()
	{
		if (data) {
			delete[] data;
			data = NULL;
			n = 0;
		}
	}
	
	int
	size() const
	{
		return n;
	}
	
	void
	resize(int size)
	{
		n = size;
		if (data) delete[] data;
		data = new double[size];
	}
	
	Vector &
	operator=(const Vector &old)
	{
		resize(old.size());
		for (int i=0; i<old.n; ++i)
			data[i] = old[i];
	}
	
	inline double &
	operator[](int i)
	{
		return data[i];
	}
	inline const double &
	operator[](int i) const
	{
		return data[i];
	}

	void
	operator+=(const Vector &v2)
	{
		for (int i=0; i<n; ++i)
			data[i] += v2[i];
	}
	
	Vector &
	operator+(const Vector &v2) const
	{
		Vector *res = new Vector(n);
		for (int i=0; i<n; ++i)
			(*res)[i] = data[i] + v2[i];
		return *res;
	}

	void
	operator*=(double a)
	{
		for (int i=0; i<n; ++i)
			data[i] *= a;
	}

	Vector &
	operator*(double a) const
	{
		Vector *res = new Vector(n);
		for (int i=0; i<n; ++i)
			res[i] = a * data[i];
		return *res;
	}


	double
	operator*(const Vector &v2) const
	{
		double res;
		for (int i=0; i<n; ++i)
			res += data[i] * v2[i];
		return res;
	}
	
	friend std::ostream& operator<<(std::ostream& cout_, const Vector& vec)
	{
		cout_ << "vector " << vec.size() << ":" << std::endl;
		for (int i=0; i<vec.size(); ++i)
			cout_ << vec[i] << "  ";
		cout_ << std::endl;
		return cout_;
	}
	
	friend std::istream& operator>>(std::istream& cin_, Vector& vec)
	{
		int n;
		cin_ >> n;
		vec.resize(n);

		for (int i=0; i<n; ++i)
			cin_ >> vec[i];
		return cin_;
	}
};

/*
int
main()
{
	double a1[7] = {1,2,3,4,5,6,7};
	double a2[7] = {1,2,3,4,5,6,7};
	
	double *v1 = new double[7];
	double *v2 = new double[7];
	
	for (int i=0; i<7; ++i) {
		v1[i] = a1[i];
		v2[i] = a2[i];
	}
	
	Vector qq(v1, 7);
	Vector ww(v1, 7);
	
	
	
	return 0;
}

*/
