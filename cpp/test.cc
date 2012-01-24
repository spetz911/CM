#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h> // for exit(1)


#include "matrix.cc"

using namespace std;

int
main()
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
