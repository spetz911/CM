/**
 * Solver for System of linear equations.
 * 
 * @author Oleg Baskakov
 * 2012. Written by NatSys Lab. (info@natsys-lab.com).
 */
#ifndef __VECTOR_CC__
#define __VECTOR_CC__

#include <iostream>
#include <fstream>
#include <vector>

// using namespace std;


template <typename Type>

class Vector : public std::vector<Type> {
public:
	Vector()
		:std::vector<Type>() {}

	Vector(int size)
		:std::vector<Type>(size) {}
	
	size_t size() const
	{
		std::vector<Type>::size();
	}

	void
	operator+=(const Vector &v2)
	{
		for (int i=0; i<size(); ++i)
			(*this)[i] += v2[i];
	}
	
	Vector &
	operator+(const Vector &v2) const
	{
		const Vector &v1 = *this;
		Vector &res = *(new Vector(size()));
		for (int i=0; i<size(); ++i)
			res[i] = v1[i] + v2[i];
		return res;
	}

	void
	operator*=(Type a)
	{
		Vector &vec = *this;
		for (int i=0; i<size(); ++i)
			vec[i] *= a;
	}

	Vector &
	operator*(Type a) const
	{
		const Vector &vec = *this;
		Vector &res = *(new Vector(size()));
		for (int i=0; i<size(); ++i)
			res[i] = vec[i] * a;
		return res;
	}

	Type
	operator*(const Vector<Type> &v2) const
	{
		const Vector &v1 = *this;
		Type res;
		for (int i=0; i<v1.size(); ++i)
			res += v1[i] * v2[i];
		return res;
	}

	static inline Type
	dot(const Vector<Type> &v1, const Vector<Type> &v2)
	{
		return v1 * v2;
	}

	std::string &
	str()
	{
		std::string *res;
		char s[100];
		for (int i=0; i<size(); ++i) {
			int n = snprintf(s, 100, "%f", (*this)[i]);
			*res += std::string(s, n) + "\t";
		}
		return *res;
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
#endif // __VECTOR_CC__

