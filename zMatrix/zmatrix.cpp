#include "zmatrix.h"

double Matrix::at(int _rows, int _cols)
{
	if (_rows < 0 || _rows > rows || _cols < 0 || _cols > cols)
		return 0.0;
	else
		return (*this)[_rows][_cols];
}

//�������ȣ��н��ݾ���ȫΪ0�еĸ���
size_t Matrix::rank()
{
	size_t count = 0;
	Matrix rrefMat = this->rref_wiki();
	double eps = 1e-10;

	for (int i = 0; i < rrefMat.rows; ++i)
	{
		for (int j = 0; j < rrefMat.cols; ++j)
		{
			if (rrefMat[i][j] > eps)
			{
				count++;
				break;
			}
		}
	}

	return count;
}

//��������
Matrix Matrix::inv()
{	
	try
	{
		if (this->rows != this->cols)
			throw runtime_error("This is not a square Matrix");
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}
	
	Matrix zgMat(this->rows, this->cols + this->cols);
	Matrix eyeMat(this->rows, this->cols);
	eyeMat.eye();
	for (int i = 0; i < zgMat.rows; ++i)
	{
		for (int j = 0; j < zgMat.cols; ++j)
		{
			if (j < this->cols)
			{
				zgMat[i][j] = data[i * cols + j];
			}
			else
			{
				zgMat[i][j] = eyeMat[i][j - this->cols];
			}
		}	
	}

	Matrix invMat_zg = zgMat.clone();
	int m = invMat_zg.rows;
	int n = invMat_zg.cols / 2;
	
	int h = 0, k = 0;
	while (h < m && k < n)
	{
		//�õ�h������k��
		//1,Ѱ�ҵ�k����Ԫ
		int i_max = invMat_zg.argMax(k, h, m - 1);
		//2,����Ԫ���ö�����ȥ��Ԫ������Ԫ��
		if (invMat_zg[i_max][k] == 0)
		{
			k++;
		}
		else
		{
			invMat_zg.swap_rows(h, i_max);
			for (int i = h + 1; i < m; ++i)
			{
				double f = invMat_zg[i][k] / invMat_zg[h][k];
				invMat_zg[i][k] = 0;
				for (int j = k + 1; j < n * 2; ++j)
				{
					invMat_zg[i][j] -= invMat_zg[h][j] * f;
				}
			}
			h++;
			k++;
		}
	}

	h = m - 1;
	k = n - 1;
	while (h >= 0 && k >= 0)
	{
		int i_max = invMat_zg.argMax(k, 0, h);
		if (invMat_zg[i_max][k] == 0)
		{
			k--;
		}
		else
		{
			invMat_zg.swap_rows(h, i_max);
			for (int i = h - 1; i > 0; --i)
			{
				double f = invMat_zg[i][k] / invMat_zg[h][k];
				invMat_zg[i][k] = 0;
				for (int j = k - 1; j > 0; --j)
				{
					invMat_zg[i][j] -= invMat_zg[h][j] * f;
				}

				for (int j = k - 1; j < invMat_zg.cols; ++j)
				{
					invMat_zg[i][j] -= invMat_zg[h][j] * f;
				}
			}
			h--;
			k--;
		}
	}

	return invMat_zg;
}

Matrix Matrix::transpose()
{
	Matrix T(cols, rows);
	for (int i = 0; i < T.rows; ++i)
	{
		for (int j = 0; j < T.cols; ++j)
		{
			T[i][j] = data[j * cols + i];
		}
	}

	return T;
}

Matrix Matrix::dot(const Matrix & m)
{
	try
	{
		if (cols != m.cols || rows != m.rows)
		{
			throw runtime_error("two different dimension matrixs");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}
	Matrix temp(rows, cols);
	temp.zeros();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < m.cols; ++j)
		{
			temp[i][j] = data[i * cols + j] * m[i][j];
		}
	}
	return temp;
}

Matrix Matrix::cross(const Matrix & m)
{
	try
	{
		if (cols != m.rows)
		{
			throw runtime_error("m1.cols != m2.rows");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}
	Matrix temp(rows, m.cols);
	temp.zeros();
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < m.cols; ++j)
		{
			for (int k = 0; k < cols; ++k)
			{
				temp[i][j] += data[i * cols + k] * m[k][j];
			}
		}
	}
	return temp;
}

/**
* @berif �������ʼ��Ϊ�վ��� 
*/
void Matrix::initEmpty()
{
	rows = cols = _size = 0;
	data = nullptr;
	refcount = nullptr;
}

//Ѱ������Ԫ����ÿһ����Ԫʱ����Ҫ�ҵ���һ��Ԫ����󣬽�����Ϊ��Ԫ
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

//�����������ݣ�inPlace
void Matrix::swap_rows(int row1, int row2)
{
	double *temp = new double[cols];
	memcpy(temp, data + row1 * cols, cols * sizeof(double));
	memcpy(data + row1 * cols, data + row2 * cols, cols * sizeof(double));
	memcpy(data + row2 * cols, temp, cols * sizeof(double));
	delete[] temp;
}

/**
* @berif �����Ĵ������󣬷����ڴ�
* @attention ���о������ݵķ��䶼Ӧ��ͨ�����øú���ʵ�֣����øú���һ����ζ�����´���������
* @param[in] _rows������
* @param[in] _cols������
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
* @berif �������ü�����ֵ
*/
int Matrix::refAdd(int *addr, int delta)
{
	int temp = *addr;
	*addr += delta;
	return temp;
}

/**
* @berif �ͷ���Դ
* @attention  �������Դ�ɸú������Ʋ��ͷ�
*/
void Matrix::release()
{
	if (refcount && refAdd(refcount, -1) == 1)	//���ü���Ϊ1���ͷ�
	{
		delete[] data;
		data = nullptr;
		delete refcount;
		refcount = nullptr;
		rows = cols = _size = 0;
		_log_("Matrix release.");
	}
}

//�������ʼ��Ϊ�����
void Matrix::zeros()
{
	for (auto i = 0; i < _size;++i)
	{
		data[i] = 0;
	}
}
//���·����ڴ沢��ʼ��Ϊ�����
void Matrix::zeros(int _rows, int _cols)
{
	create(_rows, _cols);
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 0;
	}
}
//��ʼ��Ϊ1
void Matrix::ones()
{
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 1;
	}
}
//���·����ڴ沢��ʼ��Ϊ1
void Matrix::ones(int _rows, int _cols)
{
	create(_rows, _cols);
	for (auto i = 0; i < _size; ++i)
	{
		data[i] = 1;
	}
}
//��ʼ��Ϊ��λ����
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
//���·����ڴ沢��ʼ��Ϊ��λ����
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

//�����������
//ΪĿ���������ڴ棬���������ݵ�Ŀ�����
void Matrix::copyTo(Matrix & outputMatrix) const
{
	outputMatrix.create(rows, cols);
	memcpy(outputMatrix.data, data, sizeof(double) * _size);
}
//������ʱ����Ŀ���
Matrix Matrix::clone() const
{
	Matrix tempM;
	copyTo(tempM);
	return tempM;
}

//��ֵ�������ǳ����
Matrix & Matrix::operator=(const Matrix & m)
{
	// TODO: �ڴ˴����� return ���
	if (this != &m)	//�����ֵ�����Լ�
	{
		if (m.refcount)
			refAdd(m.refcount, 1);	//���ü�����1
		release();	//�ͷ���ֵ
	
		//���¸�ֵ
		rows = m.rows;
		cols = m.cols;
		_size = m._size;
		data = m.data;
		refcount = m.refcount;
	}//�����ֵ���Լ�����ʲô������
	_log_("Matrix assignment function.");
	return *this;
}

Matrix & Matrix::operator=(initializer_list<double> il)
{
	// TODO: �ڴ˴����� return ���
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
	// TODO: �ڴ˴����� return ���
	create(1, n);
	for (int i = 0; i < n; ++i)
	{
		data[i] = inputArr[i];
	}
	return *this;
}

Matrix & Matrix::operator()(double * inputArr, int rows, int cols)
{
	// TODO: �ڴ˴����� return ���
	create(rows, cols);
	for (int i = 0; i < rows * cols; ++i)
	{
		data[i] = inputArr[i];
	}
	return *this;
}

Matrix & Matrix::operator+=(const Matrix & m)
{
	// TODO: �ڴ˴����� return ���
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

//�޲ι��캯��
Matrix::Matrix()
{
	initEmpty();
	_log_("Matrix construct without params.");
}
//���ι��캯��
Matrix::Matrix(int rows, int cols)
{
	initEmpty();
	create(rows, cols);
	_log_("Matrix construct with params.");
}

//ǳ�������캯��
Matrix::Matrix(const Matrix & m) : rows(m.rows),cols(m.cols),_size(m._size),data(m.data),refcount(m.refcount)
{
	if (refcount)
	{
		refAdd(refcount, 1);
	}
	_log_("Matrix shallow copy construct function.");
}

//��ʼ���б� ���캯��
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
//��������
Matrix::~Matrix()
{
	release();
	_log_("Matrix deconstruct.");
}

//������������
ostream &operator << (ostream& os, const Matrix& m)
{
	os.flags(ios::left); //�����
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
	//1.û�з����ڴ�����
	if (m1.data == nullptr && m1.data == m2.data)
	{
		return true;
	}
	//2.�з����ڴ棬�ж����к�����
	else if(m1.data != nullptr)
	{
		//��ַ��ͬ��ֱ�����
		if (m1.data == m2.data)
			return true;
		//��ַ����ͬ
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

//valid�;����same�ͺ�full����ʱûʵ�֣����0��
Matrix Matrix::conv(const Matrix & m)
{
	try
	{
		if (m.cols != m.rows)
		{
			throw runtime_error("�˾���Ϊ����");
		}
		if (this->cols < m.cols || this->rows < m.rows)
		{
			throw runtime_error("�˾���̫��");
		}
		if (this->cols % 2 == 0 || m.cols % 2 == 0)
		{
			throw runtime_error("����ά�Ȳ�Ϊ����");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}
	
	Matrix converMat = m.converse();
	int src_halfLen = this->rows / 2;
	int ker_halfLen = m.rows / 2;
	int res_halfLen = src_halfLen - ker_halfLen;
	int resDim = res_halfLen * 2 + 1;
	Matrix resMat(resDim, resDim);
	Matrix temp(rows, cols);
	temp.zeros();
	
	int u_start = ker_halfLen;
	int u_end = this->rows - 1 - ker_halfLen;
	int v_start = ker_halfLen;
	int v_end = this->cols - 1 - ker_halfLen;
	for (int u = u_start; u <= u_end; ++u)
	{
		for (int v = v_start; v <= v_end; ++v)
		{
			for (int i = -ker_halfLen; i <= ker_halfLen; i++)
			{
				for (int j = -ker_halfLen; j <= ker_halfLen; j++)
				{
					temp[u][v] += (*this)[u + i][v + j] * converMat[i + ker_halfLen][j + ker_halfLen];
				}
			}
		}
	}
	
	for (int i = 0; i < resDim; i++)
	{
		for (int j = 0; j < resDim; j++)
		{
			resMat[i][j] = temp[i + ker_halfLen][j + ker_halfLen];
		}
	}

	return resMat;
}

Matrix Matrix::converse() const
{
	try
	{
		if (this->cols != this->rows)
		{
			throw runtime_error("�˾���Ϊ����");
		}
	}
	catch (runtime_error err)
	{
		cout << err.what() << endl;
		return Matrix();
	}

	int r = this->rows;
	int c = this->cols;
	Matrix temp(this->rows,this->cols);
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			temp[i][j] = (*this)[r - 1 - i][c - 1 - j];
		}
	}
	return temp;
}

//������Ϊ�н��ݾ���
Matrix Matrix::rref_bad(const Matrix& M)
{
	Matrix R(M);
	int m = R.rows;
	int n = R.cols;
	//�ж��������

	double *temp = new double[n];
	try
	{
		//����Ԫ
		double divisor = R[0][0];
		if (divisor != 0)
		{
			for (int j = 0; j < n; ++j)
			{
				R[0][j] /= divisor;
			}
		
			//��Ԫ
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
// while h �� m and k �� n
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
		//�õ�h������k��
		//1,Ѱ�ҵ�k����Ԫ
		int i_max = R.argMax(k, h, m - 1);
		//2,����Ԫ���ö�����ȥ��Ԫ������Ԫ��
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

Matrix Matrix::rref_wiki_2()
{
	Matrix R = this->clone();
	int m = R.rows;
	int n = R.cols;

	int h = m - 1;
	int k = n - 1;
	while (h >= 0 && k >= 0)
	{
		int i_max = R.argMax(k, 0, h);
		if (R[i_max][k] == 0)
		{
			k--;
		}
		else
		{
			R.swap_rows(h, i_max);
			for (int i = h - 1; i > 0; --i)
			{
				double f = R[i][k] / R[h][k];
				R[i][k] = 0;
				for (int j = k - 1; j > 0; --j)
				{
					R[i][j] -= R[h][j] * f;
				}
			}
			h--;
			k--;
		}
	}
	return R;
}



