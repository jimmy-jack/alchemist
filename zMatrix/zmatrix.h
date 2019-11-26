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
	//���캯��
	Matrix();	//default constructor
	Matrix(int rows, int cols);

	Matrix(const Matrix& m);
	//6.1
	Matrix(initializer_list<double> il);//ʹ�������ʼ��������a[2] = {0,1};mat(a,2);	//zmz

	~Matrix();
	
	//��Դ����
	void create(int _rows, int _cols);
	int refAdd(int *addr, int delta);
	void release();


	//��ʼ��
	void zeros();
	void zeros(int _rows, int _cols);
	void ones();
	void ones(int _rows, int _cols);
	void eye();
	void eye(int _rows, int _cols);
	
	//2019.5.31
	//�������	
	void copyTo(Matrix& outputMatrix) const;
	Matrix clone() const;

	//2019.6.1
	//�����
	Matrix& operator = (const Matrix& m);
	Matrix& operator = (initializer_list<double> il);	//ʹ�ó�ʼ���б�����mat = {1,2,3};
	Matrix& operator () (double *inputArr, size_t n);//ʹ�����鸳ֵ������a[2] = {0,1};mat(a,2);
	Matrix& operator () (double *inputArr, int rows, int cols);
	inline double* operator [] (size_t n) { return &data[n * cols]; }	//���ص�n���׵�ַ//��Խ�紦��
	inline const double* operator [] (size_t n) const  { return &data[n * cols]; }
	Matrix& operator += (const Matrix& m);

	//���жϺʹ�С
	inline bool empty() const { return data == nullptr; }
	inline size_t size() const { return _size; }
	
	//��������
	double at(int _rows, int _cols);	//Խ�緵��0.0

	//�����������
	//������
	//20191125------------------------------------
	size_t rank();	//�������
	Matrix inv();	//��	//bug
	Matrix transpose();	//ת��
	Matrix dot(const Matrix& m);	//���
	Matrix cross(const Matrix &m);	//���
	//------------------------------------

	//20191126------------------------------------
	Matrix conv(const Matrix &m);	//���,
	Matrix converse() const;	//�Ծ��󣨾���ˣ����з�ת
	//----------------------------------------

	//20191121-------------------------------
	Matrix rref_bad(const Matrix& M);	//������Ϊ�н��ݾ���
	Matrix rref_wiki();	//������Ϊ�н��ݾ���
	Matrix rref_wiki_2();
	int argMax(int col, int row_start, int row_end);	//Ѱ������Ԫ
	void swap_rows(int row1, int row2);		//������������	
	//-------------------------------20191122

	//��Ա����
	int rows, cols;
	double *data;

private:
	size_t _size;	//size of matrix
	int *refcount;	//pointer to the reference counter
	void initEmpty();	
	
};

//���ؾ��������
//2019.5.31
ostream &operator << (ostream& os, const Matrix& m);

//2019.6.1
bool operator == (const Matrix& m1,const Matrix &m2);
bool operator != (const Matrix& m1,const Matrix &m2);
Matrix operator +(const Matrix& m1, const Matrix &m2);
Matrix operator -(const Matrix &m1, const Matrix &m2);
Matrix operator *(const Matrix &m1, const Matrix &m2);




#endif
