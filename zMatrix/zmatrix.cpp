#include "zmatrix.h"

double Matrix::at(int _rows, int _cols)
{
	if (_rows < 0 || _rows > rows || _cols < 0 || _cols > cols)
		return 0.0;
	else
		return (*this)[_rows][_cols];
}

/**
* @berif 将矩阵初始化为空矩阵 
*/
void Matrix::initEmpty()
{
	rows = cols = _size = 0;
	data = nullptr;
	refcount = nullptr;
}

//寻找列主元，对每一列消元时首先要找到哪一行元素最大，将其作为主元
int Matrix::argMax(int col, int row_start, int row_end)
{
	int i_max = row_start;
	double data_max = 0;
	for (int i = row_start; i <= row_end; ++i)
	{
		if (fabs(data[i * cols + col]) >= data_max)
		{
			i_max = i;
			data_max = fabs(data[i * cols + col]);
		}
	}
	return i_max;
}

//交换两行数据，inPlace
void Matrix::swap_rows(int row1, int row2)
{
	double *temp = new double[cols];
	memcpy(temp, data + row1 * cols, cols * sizeof(double));
	memcpy(data + row1 * cols, data + row2 * cols, cols * sizeof(double));
	memcpy(data + row2 * cols, temp, cols * sizeof(double));
	delete[] temp;
}

/**
* @berif 真正的创建矩阵，分配内存
* @attention 所有矩阵数据的分配都应该通过调用该函数实现（调用该函数一般意味着重新创建函数）
* @param[in] _rows，行数
* @param[in] _cols，列数
*/
void Matrix::create(int _rows, int _cols)
{
	release();

	rows = _rows;
	cols = _cols;
	_size = rows * cols;
	refcount = new int(1);
	data = new double[rows * cols];

	_log_("Matrix create.");
}
/**
* @berif 控制引用计数的值
*/
int Matrix::refAdd(int *addr, int delta)
{
	int temp = *addr;
	*addr += delta;
	return temp;
}

/**
* @berif 释放资源
* @attention  矩阵的资源由该函数控制并释放
*/
void Matrix::release()
{
	if (refcount && refAdd(refcount, -1) == 1)	//引用计数为1，释放
	{
		delete[] data;
		data = nullptr;
		delete refcount;
		refcount = nullptr;
		rows = cols = _size = 0;
		_log_("Matrix release.");
	}
}

//将矩阵初始化为零矩阵
void Matrix::zeros()
{
	for (auto i = 0; i < _size;++i)
	{
		data[i] = 0;
	}
}
//重新分配内存并初始化为零矩阵
void Matrix::zeros(int _rows, int _cols)
{
	create(_rows, _cols);
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 0;
	}
}
//初始化为1
void Matrix::ones()
{
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 1;
	}
}
//重新分配内存并初始化为1
void Matrix::ones(int _rows, int _cols)
{
	create(_rows, _cols);
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 1;
	}
}
//初始化为单位矩阵
void Matrix::eye()
{
	size_t pos = 0;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if (i == j)
				data[pos] = 1;
			else
				data[pos] = 0;
			++pos;
		}
	}
}
//重新分配内存并初始化为单位矩阵
void Matrix::eye(int _rows, int _cols)
{
	create(_rows, cols);
	size_t pos = 0;
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			if (i == j)
				data[pos] = 1;
			else
				data[pos] = 0;
			++pos;
		}
	}
}

//两个深拷贝函数
//为目标矩阵分配内存，并拷贝数据到目标矩阵
void Matrix::copyTo(Matrix & outputMatrix) const
{
	outputMatrix.create(rows, cols);
	memcpy(outputMatrix.data, data, sizeof(double) * _size);
}
//返回临时矩阵的拷贝
Matrix Matrix::clone() const
{
	Matrix tempM;
	copyTo(tempM);
	return tempM;
}

//赋值运算符，浅拷贝
Matrix & Matrix::operator=(const Matrix & m)
{
	// TODO: 在此处插入 return 语句
	if (this != &m)	//如果右值不是自己
	{
		if (m.refcount)
			refAdd(m.refcount, 1);	//引用计数加1
		release();	//释放左值
	
		//重新赋值
		rows = m.rows;
		cols = m.cols;
		_size = m._size;
		data = m.data;
		refcount = m.refcount;
	}//如果右值是自己，则什么都不做
	_log_("Matrix assignment function.");
	return *this;
}

Matrix & Matrix::operator=(initializer_list<double> il)
{
	// TODO: 在此处插入 return 语句
	size_t num = sizeof(il) / sizeof(double);
	if (rows != 0 && cols != 0)
	{
		create(rows, cols);
	}
	else
	{
		create(1, num);
	}

	auto index = il.begin();
	auto end = il.end();
	for (int i = 0; i < _size; ++i, ++index)
	{
		if (index < end)
		{
			data[i] = *index;
		}
		else
		{
			data[i] = 0.0;
		}
	}

	return *this;
}

Matrix & Matrix::operator()(double * inputArr, size_t n)
{
	// TODO: 在此处插入 return 语句
	create(1, n);
	for (int i = 0; i < n; ++i)
	{
		data[i] = inputArr[i];
	}
	return *this;
}

Matrix & Matrix::operator()(double * inputArr, int rows, int cols)
{
	// TODO: 在此处插入 return 语句
	create(rows, cols);
	for (int i = 0; i < rows * cols; ++i)
	{
		data[i] = inputArr[i];
	}
	return *this;
}

Matrix & Matrix::operator+=(const Matrix & m)
{
	// TODO: 在此处插入 return 语句
	try
	{
		if (rows != m.rows || cols != m.cols)
		{
			//_log_("Matrix's size are not the same.");
			throw runtime_error("Matrix's size are not the same.");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return *this;
	}
	for (auto i = 0; i < _size; ++i)
	{
		data[i] += m.data[i];
	}
	return *this;
}

//无参构造函数
Matrix::Matrix()
{
	initEmpty();
	_log_("Matrix construct without params.");
}
//带参构造函数
Matrix::Matrix(int rows, int cols)
{
	initEmpty();
	create(rows, cols);
	_log_("Matrix construct with params.");
}

//浅拷贝构造函数
Matrix::Matrix(const Matrix & m) : rows(m.rows),cols(m.cols),_size(m._size),data(m.data),refcount(m.refcount)
{
	if (refcount)
	{
		refAdd(refcount, 1);
	}
	_log_("Matrix shallow copy construct function.");
}

//初始化列表 构造函数
Matrix::Matrix(initializer_list<double> il)
{
	//size_t num = sizeof(il) / sizeof(double);
	if (rows != 0 && cols != 0)
	{
		create(rows, cols);
	}
	else
	{
		create(1, il.size());
	}

	auto index = il.begin();
	auto end = il.end();
	for (int i = 0; i < _size; ++i, ++index)
	{
		if (index < end)
		{
			data[i] = *index;
		}
		else
		{
			data[i] = 0.0;
		}
	}
	
	_log_("Matrix construct with initializeList.");
}
//析构函数
Matrix::~Matrix()
{
	release();
	_log_("Matrix deconstruct.");
}

//重载输出运算符
ostream &operator << (ostream& os, const Matrix& m)
{
	os.flags(ios::left); //左对齐
	size_t pos = 0;
	os  << '[';
	for (int i = 0; i < m.rows; ++i) {
		for (int j = 0; j < m.cols; ++j) {
			os << setw(10) << m.data[pos];
			//if (m.cols != j + 1)
			//	os << ',';
			++pos;
		}
		if (m.rows != i + 1) 
			os  << ';' << endl << ' ';
		else
			os << ']' << endl;
	}
	return os;
}

bool operator == (const Matrix& m1, const Matrix &m2)
{
	//1.没有分配内存的情况
	if (m1.data == nullptr && m1.data == m2.data)
	{
		return true;
	}
	//2.有分配内存，判断行列和数据
	else if(m1.data != nullptr)
	{
		//地址相同，直接相等
		if (m1.data == m2.data)
			return true;
		//地址不相同
		if (m1.rows == m2.rows && m1.cols == m2.cols)
		{
			int i = 0;
			for (; i < m1.size(); ++i)
			{
				if (m1.data[i] != m2.data[i])
					break;
			}
			if (i == m1.size())
				return true;
		}
	}
	return false;
}
bool operator != (const Matrix& m1,const Matrix &m2)
{
	return !(m1 == m2);
}
Matrix operator +(const Matrix& m1, const Matrix &m2) 
{
	try
	{
		if (m1.rows != m2.rows || m1.cols != m2.cols)
			throw runtime_error("m1.rows != m2.rows || m1.cols != m2.cols");
	}
	catch (runtime_error err)
	{
		//_log_("m1.rows != m2.rows || m1.cols != m2.cols");
		cout << err.what() << endl;
		return Matrix();
	}
	Matrix temp(m1.rows,m1.cols);
	for (int i = 0; i < m1.size(); ++i)
	{
		temp.data[i] = m1.data[i] + m2.data[i];
	}
	
	return temp;
}
Matrix operator -(const Matrix &m1, const Matrix &m2) 
{
	try
	{
		if (m1.rows != m2.rows || m1.cols != m2.cols)
			throw runtime_error("m1.rows != m2.rows || m1.cols != m2.cols");
	}
	catch (runtime_error err)
	{
		//_log_("m1.rows != m2.rows || m1.cols != m2.cols");
		cout << err.what() << endl;
		return Matrix();
	}
	Matrix temp(m1.rows,m1.cols);
	for (int i = 0; i < m1.size(); ++i)
	{
		temp.data[i] = m1.data[i] - m2.data[i];
	}
	
	return temp;
}
Matrix operator *(const Matrix &m1, const Matrix &m2) 
{
	try
	{
		if (m1.cols != m2.rows)
		{
			throw runtime_error("m1.cols != m2.rows");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}
	Matrix temp(m1.rows,m2.cols);
	temp.zeros();
	for (int i = 0; i < m1.rows; ++i)
	{
		for (int j = 0; j < m2.cols; ++j)
		{
			for (int k = 0; k < m1.cols; ++k)
			{
				temp[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
	return temp;
}

//将矩阵化为行阶梯矩阵
Matrix Matrix::rref_bad(const Matrix& M)
{
	Matrix R(M);
	int m = R.rows;
	int n = R.cols;
	//判断奇异矩阵？

	double *temp = new double[n];
	try
	{
		//找主元
		double divisor = R[0][0];
		if (divisor != 0)
		{
			for (int j = 0; j < n; ++j)
			{
				R[0][j] /= divisor;
			}
		
			//消元
			for (int j = 0; j < n; ++j)
			{
				for (int i = 1; i < m; ++i)
				{
												
				}
			}

		}
		else
		{
			throw runtime_error("divisor is 0");	
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		delete[] temp;
		return R;
	}

	delete[] temp;
	return R;
}



/* --------------------	from wiki   ----------------------------------
// h := 1 /* Initialization of the pivot row */
// k := 1 /* Initialization of the pivot column */
// while h ≤ m and k ≤ n
//   /* Find the k-th pivot: */
//   i_max := argmax (i = h ... m, abs(A[i, k]))
//   if A[i_max, k] = 0
//     /* No pivot in this column, pass to next column */
//     k := k+1
//   else
//      swap rows(h, i_max)
//      /* Do for all rows below pivot: */
//      for i = h + 1 ... m:
//         f := A[i, k] / A[h, k]
//         /* Fill with zeros the lower part of pivot column: */
//         A[i, k]  := 0
//         /* Do for all remaining elements in current row: */
//         for j = k + 1 ... n:
//            A[i, j] := A[i, j] - A[h, j] * f
//      /* Increase pivot row and column */
//      h := h+1 
//      k := k+1	*/
Matrix Matrix::rref_wiki()
{
	Matrix R = this->clone();
	int m = R.rows;
	int n = R.cols;
	
	int h = 0, k = 0;
	while (h < m && k < n)
	{
		//用第h行消第k列
		//1,寻找第k列主元
		int i_max = R.argMax(k, h, m - 1);
		//2,将主元行置顶，消去主元行以下元素
		if (R[i_max][k] == 0)
		{
			k++;
		}
		else
		{
			R.swap_rows(h, i_max);
			for (int i = h + 1; i < m; ++i)
			{
				double f = R[i][k] / R[h][k];
				R[i][k] = 0;
				for (int j = k + 1; j < n; ++j)
				{
					R[i][j] -= R[h][j] * f;
				}
			}
			h++;
			k++;
		}
	}

	return R;
}



