#ifndef _ZMATRIX_H
#define _ZMATRIX_H

#include <iostream>
#include <iomanip>
#include "config_default.h"
using namespace std;

class Matrix
{
public:
	//2019.5.30
	//构造函数
	Matrix();	//default constructor
	Matrix(int rows, int cols);

	Matrix(const Matrix& m);
	//6.1
	Matrix(initializer_list<double> il);//使用数组初始化，形如a[2] = {0,1};mat(a,2);	//zmz

	~Matrix();
	
	//资源控制
	void create(int _rows, int _cols);
	int refAdd(int *addr, int delta);
	void release();


	//初始化
	void zeros();
	void zeros(int _rows, int _cols);
	void ones();
	void ones(int _rows, int _cols);
	void eye();
	void eye(int _rows, int _cols);
	
	//2019.5.31
	//深拷贝函数	
	void copyTo(Matrix& outputMatrix) const;
	Matrix clone() const;

	//2019.6.1
	//运算符
	Matrix& operator = (const Matrix& m);
	Matrix& operator = (initializer_list<double> il);	//使用初始化列表形如mat = {1,2,3};
	Matrix& operator () (double *inputArr, size_t n);//使用数组赋值，形如a[2] = {0,1};mat(a,2);
	Matrix& operator () (double *inputArr, int rows, int cols);
	inline double* operator [] (size_t n) { return &data[n * cols]; }	//返回第n行首地址//无越界处理
	inline const double* operator [] (size_t n) const  { return &data[n * cols]; }
	Matrix& operator += (const Matrix& m);

	//空判断和大小
	inline bool empty() const { return data == nullptr; }
	inline size_t size() const { return _size; }
	
	//访问数据
	double at(int _rows, int _cols);	//越界返回0.0

	//其他矩阵操作
	//待完善
	//20191125------------------------------------
	size_t rank();	//矩阵的秩
	Matrix inv();	//逆	//bug
	Matrix transpose();	//转置
	Matrix dot(const Matrix& m);	//点乘
	Matrix cross(const Matrix &m);	//叉乘
	//------------------------------------

	//20191126------------------------------------
	Matrix conv(const Matrix &m);	//卷积,
	Matrix converse() const;	//对矩阵（卷积核）进行翻转
	//----------------------------------------

	//20191121-------------------------------
	Matrix rref_bad(const Matrix& M);	//化矩阵为行阶梯矩阵
	Matrix rref_wiki();	//化矩阵为行阶梯矩阵
	Matrix rref_wiki_2();
	int argMax(int col, int row_start, int row_end);	//寻找列主元
	void swap_rows(int row1, int row2);		//交换两行数据	
	//-------------------------------20191122

	//成员变量
	int rows, cols;
	double *data;

private:
	size_t _size;	//size of matrix
	int *refcount;	//pointer to the reference counter
	void initEmpty();	
	
};

//重载矩阵运算符
//2019.5.31
ostream &operator << (ostream& os, const Matrix& m);

//2019.6.1
bool operator == (const Matrix& m1,const Matrix &m2);
bool operator != (const Matrix& m1,const Matrix &m2);
Matrix operator +(const Matrix& m1, const Matrix &m2);
Matrix operator -(const Matrix &m1, const Matrix &m2);
Matrix operator *(const Matrix &m1, const Matrix &m2);




#endif
