#include "Matrix.h"

Matrix::Matrix()
{

}

Matrix::Matrix(int rows, int cols)
{
	coefs.resize(rows, vector<float>(cols));
}

vector<vector<float>> Matrix::getMatrix()
{
	return coefs;
}

float Matrix::getData(int row, int col)
{
	return coefs[row][col];
}

void Matrix::setData(int row, int col, float value)
{
	coefs[row][col] = value;
}

float Matrix::determinant_2x2(vector<vector<float>> mat)
{
	return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

float Matrix::determinant(vector<vector<float>> mat)
{
	if (mat.empty())
	{
		throw exception();
	}

	int n = mat.size();

	if (n == 1)
	{
		return mat[0][0];
	}

	if (n == 2)
	{
		return determinant_2x2(mat);
	}

	float result = 0;
	for (int i = 0; i < n; i++)
	{
		vector<vector<float>> submat;
		for (int j = 1; j < n; j++)
		{
			vector<float> row;
			for (int k = 0; k < n; k++)
			{
				if (k != i)
				{
					row.push_back(mat[j][k]);
				}
			}
			if (!row.empty())
			{
				submat.push_back(row);
			}
		}
		int sign = (i % 2 == 0) ? 1 : -1;
		result += sign * mat[0][i] * determinant(submat);
	}

	return result;
}

float Matrix::getDetMinor(int row, int col)
{
	int n = coefs.size();
	vector<vector<float>> minor(n - 1, vector<float>(n - 1));

	int minor_row = 0, minor_col = 0;

	for (int i = 0; i < n; i++)
	{
		if (i != row)
		{
			minor_col = 0;
			for (int j = 0; j < n; j++)
			{
				if (j != col)
				{
					minor[minor_row][minor_col] = coefs[i][j];
					minor_col++;
				}
			}
			minor_row++;
		}
	}

	return determinant(minor);
}

vector<vector<float>> Matrix::getInverseMatrix()
{
	float det = determinant(coefs);

	if (det == 0)
	{
		throw exception();
	}
	else
	{
		int n = coefs.size();
		vector<vector<float>> result(n, vector<float>(n));

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				int sign = pow(-1, i + j);
				result[j][i] = sign * getDetMinor(i, j) / det;
			}
		}

		return result;
	}
}

Matrix::~Matrix()
{
	coefs.clear();
	coefs.~vector();
}