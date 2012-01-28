#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h> // for exit(1)

#include "vector.cc"
#include "matrix.cc"
#include "SLAE.cc"

using namespace std;

int
test()
{
	Matrix A(2,2);
	Matrix B;
	Matrix C;
	
	A[0][0] = 1.0;
	A[1][0] = 2.0;
	
	cin >> B;
	
//	std::cout << A.size() << std::endl;
//	std::cout << B.size() << std::endl;
	
	C = A + B;
	
	cout << C;

}



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
	
	Vector<double> qq(7);
	Vector<double> ww(7);
	
	for (int i=0; i<7; ++i) {
		qq[i] = a1[i];
		ww[i] = a2[i];
	}
	
	Matrix A;
	cin >> A;
	
	Matrix B = A;
	
	cout << "A:" << endl << A;

	cout << "B:" << endl << B;
	
	Matrix C = A * B;
	
	cout << "A*B:" << endl << C;
	
	
	
	return 0;
}

