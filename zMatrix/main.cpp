#include <iostream>
#include "config_default.h"
#include "zmatrix.h"

using namespace std;

void test()
{
	//Matrix mat1(3, 2);
	////mat1.eye();
	////cout << mat1;
	////mat1.ones(5, 4);
	////cout << mat1;
	////Matrix mat2(mat1);
	////Matrix mat3;
	////mat3 = mat1;
	//mat1.eye();
	//cout << mat1;
	//Matrix mat2;
	//mat2 = mat1.clone();
}

void myTest()
{
	//cout << "----------A---------------------" << endl;
	//Matrix A(3, 3);
	//A.zeros();
	//cout << A;
	//A.ones();
	//cout << A;
	//A.eye();
	//cout << A;
	//
	//cout << "----------CA---------------------" << endl;
	//Matrix CA(A);
	//cout << CA;

	//cout << "----------B---------------------" << endl;
	//Matrix B = A.clone();	//构造tempM,拷贝构造B，析构tempM
	//cout << B;

	//cout << "----------C---------------------" << endl;
	//Matrix C;
	//A.copyTo(C);
	//cout << C;

	//cout << "----------D---------------------" << endl;
	//Matrix D;	//构造D
	//D = A.clone();	//构造tempM,拷贝构造临时对象，析构tempM，赋值，析构临时对象。
	//cout << D;
	
	Matrix matA = { 1,2,3 };
	cout << matA;
	double dA[2] = { 1,2 };
	Matrix matB;
	matB(dA, 2);
	cout << matB;
	Matrix matC = matA;
	cout << matC[0][0] << endl;
	cout << ((matA == matC) ? "matA == matC" : "matA != matC" )<< endl;
	cout << ((matB != matC) ? "matB != matC" : "matB == matC" )<< endl;
	Matrix sumM;
//	sumM = matA + matB;
	sumM = matA + matC;
	cout << sumM;
	Matrix resM;
	resM = matA - matC;
	cout << resM;
	Matrix matD(3, 1);
	matD = { 3,
		2,
		1 };
	cout << matD;
	Matrix dotM = matA * matD;
	cout << dotM;
	cout << matA.at(0, 0) << endl;
}

int main()
{
//	_log_("log test");
//	__log__("wfefw");
	myTest();
	
	system("pause");

	return 0;
}
