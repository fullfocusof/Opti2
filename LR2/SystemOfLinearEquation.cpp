#include "SystemOfLinearEquation.h"

SystemOfLinearEquation::SystemOfLinearEquation()
{
	
}

SystemOfLinearEquation::SystemOfLinearEquation(vector<vector<float>>& coefs, vector<float>& consts)
{
	int rows = coefs.size();
	int cols = coefs[0].size();

	coefsSLE = Matrix(rows, cols);
	constants.resize(rows);
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			coefsSLE.setData(i,j, coefs[i][j]);
		}
		constants[i] = consts[i];
	}
}

Matrix SystemOfLinearEquation::getMatrix()
{
	return coefsSLE;
}

//void SystemOfLinearEquation::readFromFile(string filename)
//{
//	ifstream ifs(filename);
//
//	if (ifs.is_open())
//	{
//		int rows, cols;
//		ifs >> rows;
//		cols = rows;
//		coefsSLE = Matrix(rows, cols);
//		constants.resize(rows);
//
//		float value;
//		for (int i = 0; i < rows; i++)
//		{
//			for (int j = 0; j < cols; j++)
//			{
//				ifs >> value;
//				coefsSLE.setData(i, j, value);
//			}
//			ifs >> constants[i];
//		}
//
//		ifs.close();
//
//		cout << "�������� ������ ������";
//	}
//	else
//	{
//		cerr << ">>>������ �������� �����<<<";
//	}
//}

void SystemOfLinearEquation::printData()
{
	if (coefsSLE.getMatrix().empty() || constants.empty())
	{
		cout << "������ �����������";
	}
	else
	{
		int size = coefsSLE.getMatrix().size();
		cout << "���������� ��������� - " << size << endl << "���������� ���������� - " << size << endl;
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size - 1; j++)
			{
				cout << coefsSLE.getMatrix()[i][j] << "*x" << j + 1 << " + " << " ";
			}
			cout << coefsSLE.getMatrix()[i][size - 1] << "*x" << size << " = " << constants[i];
			cout << endl;
		}
	}
}

float SystemOfLinearEquation::getDet()
{
	return coefsSLE.determinant(coefsSLE.getMatrix());
}

vector<float> SystemOfLinearEquation::gaussian_solve()
{
	if (coefsSLE.getMatrix().empty())
	{
		cout << "������ �����������";
		return vector<float>();
	}

	int n = coefsSLE.getMatrix().size();
	vector<vector<float>> coefsTemp = coefsSLE.getMatrix();
	vector<float> constTemp = constants;

	for (int k = 0; k < n; k++) // ������ ���
	{
		for (int i = k + 1; i < n; i++)
		{
			float factor = coefsTemp[i][k] / coefsTemp[k][k];
			for (int j = k; j < n; j++)
			{
				coefsTemp[i][j] -= factor * coefsTemp[k][j];
			}
			constTemp[i] -= factor * constTemp[k];
		}
	}

	bool hasUniqueSolution = true; // ��������
	for (int i = 0; i < n; i++)
	{
		if (coefsTemp[i][i] == 0 && constTemp[i] != 0)
		{
			hasUniqueSolution = false;
			break;
		}
	}

	if (!hasUniqueSolution)
	{
		cout << "������� ��������� �� ����� �������" << endl << endl;
		return vector<float>();
	}
	else
	{
		int rank = 0;
		for (int i = 0; i < n; i++)
		{
			if (coefsTemp[i][i] != 0)
			{
				rank++;
			}
		}

		if (rank < n)
		{
			cout << "������� ��������� ����� ����������� ��������� �������" << endl << endl;
			return vector<float>();
		}
		else
		{
			cout << "������� ��������� ����� ���� �������:" << endl;

			vector<float> solution(n); // �������� ���
			for (int i = n - 1; i >= 0; i--)
			{
				float sum = 0;
				for (int j = i + 1; j < n; ++j) 
				{
					sum += coefsTemp[i][j] * solution[j];
				}
				solution[i] = (constTemp[i] - sum) / coefsTemp[i][i];	
			}

			for (int i = 0; i < solution.size(); i++)
			{
				cout << "x" << i + 1 << " = " << solution[i] << endl;
			}

			cout << endl;

			return solution;
		}
	}
}

vector<float> SystemOfLinearEquation::gaussian_solveLeadElem()
{
	if (coefsSLE.getMatrix().empty())
	{
		cout << "������ �����������";
		return vector<float>();
	}

	int n = coefsSLE.getMatrix().size();
	vector<vector<float>> coefsTemp = coefsSLE.getMatrix();
	vector<float> constTemp = constants;

	for (int i = 0; i < n; i++) // ������ ���
	{
		int max_row = i;
		for (int j = i + 1; j < n; j++)
		{
			if (abs(coefsTemp[j][i]) > abs(coefsTemp[max_row][i]))
			{
				max_row = j;
			}
		}

		swap(coefsTemp[i], coefsTemp[max_row]);
		swap(constTemp[i], constTemp[max_row]);

		for (int j = i + 1; j < n; j++)
		{
			float factor = coefsTemp[j][i] / coefsTemp[i][i];
			constTemp[j] -= factor * constTemp[i];
			for (int k = i; k < n; k++)
			{
				coefsTemp[j][k] -= factor * coefsTemp[i][k];
			}
		}
	}

	bool hasUniqueSolution = true; // ��������
	for (int i = 0; i < n; i++) 
	{
		if (coefsTemp[i][i] == 0 && constTemp[i] != 0)
		{
			hasUniqueSolution = false;
			break;
		}
	}

	if (!hasUniqueSolution) 
	{
		cout << "������� ��������� �� ����� �������" << endl << endl;
		return vector<float>();
	}
	else
	{
		int rank = 0;
		for (int i = 0; i < n; i++) 
		{
			if (coefsTemp[i][i] != 0) 
			{
				rank++;
			}
		}

		if (rank < n)
		{
			cout << "������� ��������� ����� ����������� ��������� �������" << endl << endl;
			return vector<float>();
		}
		else
		{
			cout << "������� ��������� ����� ���� �������" << endl;

			vector<float> solution(n); // �������� ���
			for (int i = n - 1; i >= 0; i--)
			{
				solution[i] = constants[i];
				for (int j = i + 1; j < n; j++)
				{
					solution[i] -= coefsTemp[i][j] * solution[j];
				}
				solution[i] /= coefsTemp[i][i];
			}

			for (int i = 0; i < n; i++)
			{
				if (solution[i] == -0.0f) solution[i] = 0.0f;
			}

			for (int i = 0; i < solution.size(); i++)
			{
				cout << "x" << i + 1 << " = " << solution[i] << endl;
			}

			return solution;
		}
	}
}

vector<float> SystemOfLinearEquation::inverseMatrix_Solve()
{
	try
	{
		Matrix invMatrix = coefsSLE.getInverseMatrix();
		vector<float> constTemp = constants;
		vector<float> solution;

		int matColsSize = invMatrix.getMatrix()[0].size();
		int constSize = constTemp.size();
		if (matColsSize != constSize) throw exception("������ ������� ������������� ������");

		int n = invMatrix.getMatrix().size();
		for (int i = 0; i < n; i++)
		{
			float value = 0.0f;
			int k = invMatrix.getMatrix()[i].size();
			for (int j = 0; j < k; j++)
			{
				value += invMatrix.getData(i, j) * constTemp[j];
			}
			solution.push_back(value);
		}

		return solution;
	}
	catch (const exception&)
	{
		cout << "������������ ������� ������� ����� ����, ��� �������� �������" << endl << endl;
	}
	
}

void SystemOfLinearEquation::simpleIteration(int choice)
{
	vector<vector<float>> coefsTemp = coefsSLE.getMatrix();
	vector<float> constTemp = constants;
	int n = coefsTemp.size();
	vector<int> idX(n, 0);

	for (int i = 0; i < n; i++) // ���������� � ����
	{
		float max_in_row = FLT_MIN;
		int colID = 0;
		for (int j = 0; j < n; j++)
		{
			if (coefsTemp[i][j] > max_in_row)
			{
				max_in_row = coefsTemp[i][j];
				colID = j;
			}
		}

		idX[i] = colID;

		for (int j = 0; j < n; j++)
		{
			if (j == colID)
			{
				coefsTemp[i][j] /= max_in_row;
			}
			else
			{
				coefsTemp[i][j] /= -max_in_row;
			}			
		}
		constTemp[i] /= max_in_row;
	}

	float p1 = 0, p2 = 0, p3 = 0, square_sum = 0;

	for (int i = 0; i < n; i++) // �������
	{
		float sum_in_row = 0, sum_in_col = 0;
		for (int j = 0; j < n; j++)
		{
			if (coefsTemp[i][j] != 1.f)
			{
				sum_in_row += abs(coefsTemp[i][j]);
				square_sum += pow(coefsTemp[i][j], 2);				
			}
			if (coefsTemp[j][i] != 1.f)
			{
				sum_in_col += abs(coefsTemp[j][i]);
			}
		}

		if (sum_in_row > p1)
		{
			p1 = sum_in_row;
		}
		if (sum_in_col > p2)
		{
			p2 = sum_in_col;
		}
	}

	p3 = sqrt(square_sum);

	if (p1 > 1 && p2 > 1 && p3 > 1)
	{
		cout << "���������� ����� �������";
		return;
	}

	float eps = 0.00001, favMetrics = min(p1, p2, p3), compressionCoef = eps * (1 - favMetrics) / favMetrics; // �� ������

	vector<vector<float>> values(2, vector<float>(n, 0));
	float dValues = FLT_MAX;

	if (choice == 0) // ������� ��������
	{
		do
		{
			for (int i = 0; i < n; i++)
			{
				values[0][i] = values[1][i];
				values[1][i] = 0;
			}

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (j != idX[i])
					{
						values[1][idX[i]] += coefsTemp[i][j] * values[0][j];
					}
				}
				values[1][idX[i]] += constTemp[i];
			}

			dValues = max(values[1][0] - values[0][0], values[1][1] - values[0][1], values[1][2] - values[0][2]);

		} while (dValues > compressionCoef);
	}
	else // ����� �������
	{
		do
		{
			for (int i = 0; i < n; i++)
			{
				values[0][i] = values[1][i];
				values[1][i] = 0;
			}

			//values[1][0] += coefsTemp[1][1] * values[0][1] + coefsTemp[1][2] * values[0][2] + constTemp[1]; // x1  1 ������   2
			//values[1][1] += coefsTemp[2][0] * values[1][0] + coefsTemp[2][2] * values[0][2] + constTemp[2]; // x2  2 ������   0
			//values[1][2] += coefsTemp[0][0] * values[1][0] + coefsTemp[0][1] * values[1][1] + constTemp[0]; // x3  0 ������   1

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (i == idX[j]) // i = 0 j = 1 idX[j] = 0
					{
						for (int k = 0; k < n; k++)
						{
							if (k != idX[j])
							{
								values[1][i] += coefsTemp[j][k] * values[(i - 1) < 0 ? i : (i - 1)][k];
							}				
						}		
						values[1][i] += constTemp[j];
					}	
				}
				
			}

			dValues = max(values[1][0] - values[0][0], values[1][1] - values[0][1], values[1][2] - values[0][2]);

		} while (dValues > compressionCoef);
	}

	cout << "������� ��������� ����� ���� �������" << endl;
	for (int i = 0; i < values[0].size(); i++)
	{
		cout << "x" << i + 1 << " = " << values[0][i] << endl;
	}
}

SystemOfLinearEquation::~SystemOfLinearEquation()
{
	coefsSLE.~Matrix();
	constants.~vector();
}