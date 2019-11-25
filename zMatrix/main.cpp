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
	//---------------------------------------------		列举各种构造方法	--------------------------------------------//
	//cout << "----------A---------------------" << endl;
	//Matrix A(3, 3);
	//A.zeros();
	//cout << A;
	//A.ones();
	//cout << A;
	//A.eye();
	//cout << A;
	//
	////拷贝（浅）构造
	//cout << "----------CA---------------------" << endl;
	//Matrix CA(A);
	//cout << CA;

	////复制1
	//cout << "----------B---------------------" << endl;
	//Matrix B = A.clone();	//构造tempM,拷贝构造B，析构tempM
	//cout << B;

	////复制2
	//cout << "----------C---------------------" << endl;
	//Matrix C;
	//A.copyTo(C);
	//cout << C;

	////赋值
	//cout << "----------D---------------------" << endl;
	//Matrix D;	//构造D
	//D = A.clone();	//构造tempM,拷贝构造临时对象，析构tempM，赋值，析构临时对象。
	//cout << D;

	////初始化列表
	//cout << "-------------------	E	---------------------------" << endl;
	//Matrix E{ 1,2,3 };
	//cout << E;

	//test operation--------------------------------------
//	Matrix matA = { 1,2,3 };
//	cout << matA;
//	double dA[2] = { 1,2 };
//	Matrix matB;
//	matB(dA, 2);
//	cout << matB;
//	Matrix matC = matA;
//	cout << matC[0][0] << endl;
//	cout << ((matA == matC) ? "matA == matC" : "matA != matC" )<< endl;
//	cout << ((matB != matC) ? "matB != matC" : "matB == matC" )<< endl;
//	Matrix sumM;
////	sumM = matA + matB;
//	sumM = matA + matC;
//	cout << sumM;
//	Matrix resM;
//	resM = matA - matC;
//	cout << resM;
//	Matrix matD(3, 1);
//	matD = { 3,
//		2,
//		1 };
//	cout << matD;
//	Matrix dotM = matA * matD;
//	cout << dotM;
//	cout << matA.at(0, 0) << endl;

	//test rank
	Matrix swapMat(3, 2);
	swapMat = { 1,2,
		3,4,
		5,6
	};
	cout << swapMat;
	swapMat.swap_rows(0,2);
	cout << "after swap:" << endl << swapMat;
	//cout << swapMat.argMax(0, 1, 2) << endl;

	Matrix rrefMat_src(4, 5);
	rrefMat_src =
	{
		1,1,0,-3,-1,
		1,-1,2,-1,0,
		4,-2,6,3,-4,
		2,4,-2,4,-7
	};
	cout << rrefMat_src;
	Matrix rrefMat_dst = rrefMat_src.rref_wiki();
	cout << rrefMat_dst;
	cout << "rank = " <<rrefMat_src.rank() << endl;

	Matrix invMat_src(3, 3);
	invMat_src =
	{
		0,2,-1,
		1,1,2,
		-1,-1,-1
	};
	Matrix invMat_dst = invMat_src.inv();
	cout << "inv:" << endl << invMat_dst;

	Matrix transMat = rrefMat_src.transpose();
	cout << "transMat" << endl << transMat;

	cout << "-----------------------------------------	test end	-----------------------------------------------" << endl;
}

int main()
{
	_log_("log test");
	//OutputDebugString("输出信息到输出窗口");
	myTest();
	
	system("pause");

	return 0;
}
