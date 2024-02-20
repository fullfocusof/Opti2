#include "SystemOfLinearEquation.h"

SystemOfLinearEquation::SystemOfLinearEquation()
{
	cntEquation = cntVar = 0;
}

SystemOfLinearEquation::SystemOfLinearEquation(vector<vector<float>> matrix, int n)
{
	cntEquation = n;
	cntVar = n + 1;
	SLE = matrix;
}

int SystemOfLinearEquation::getCntEq()
{
	return cntEquation;
}

int SystemOfLinearEquation::getCntVar()
{
	return cntVar;
}

vector<vector<float>> SystemOfLinearEquation::getMatrix()
{
	return SLE;
}

void SystemOfLinearEquation::getFromFile(string filename)
{
	ifstream ifs(filename);

	if (ifs.is_open())
	{
		ifs >> cntEquation;
		cntVar = cntEquation + 1;
		SLE.resize(cntEquation, vector<float>(cntVar + 1));

		while (!ifs.eof())
		{
			for (int i = 0; i < cntEquation; i++)
			{
				for (int j = 0; j < cntVar; j++)
				{
					ifs >> SLE[i][j];
				}
			}

			ifs.close();
		}

		cout << "Успешный импорт данных";
	}
	else
	{
		cerr << ">>>Ошибка открытия файла<<<";
	}
}

void SystemOfLinearEquation::printData()
{
	if (SLE.empty())
	{
		cout << "Данные отсутствуют";
	}
	else
	{
		cout << "Количество уравнений - " << cntEquation << endl << "Количество переменных - " << cntVar - 1 << endl << endl;
		for (int i = 0; i < cntEquation; i++)
		{
			for (int j = 0; j < cntVar; j++)
			{
				cout << SLE[i][j] << "\t";
			}
			cout << endl;
		}
	}
}

pair<vector<float>, int> SystemOfLinearEquation::gaussian_solve()
{
	vector<float> solution(cntVar - 1);

	for (int i = 0; i < cntEquation; i++) // приведение к треугольной
	{
		for (int j = i + 1; j < cntEquation; j++)
		{
			float factor = SLE[j][i] / SLE[i][i];
			for (int k = i; k < cntVar; k++)
			{
				SLE[j][k] -= factor * SLE[i][k];
			}
		}
	}

	for (int i = cntEquation - 1; i >= 0; i--) // обратный ход
	{
		solution[i] = SLE[i][cntEquation];
		for (int j = i + 1; j < cntEquation; j++)
		{
			solution[i] -= SLE[i][j] * solution[j];
		}
		solution[i] /= SLE[i][i];
	}

	bool hasUniqueSolution = true;  // проверка
	bool hasInfiniteSolutions = false;
	bool noSolution = false;

	for (int i = 0; i < SLE.size(); i++) 
	{
		bool allZero = true;
		bool hasNonZero = false;
		for (int j = 0; j < SLE[i].size() - 1; j++) 
		{
			if (SLE[i][j] != 0) {
				allZero = false;
				hasNonZero = true;
				break;
			}
		}
		if (allZero && SLE[i].back() != 0) 
		{
			noSolution = true;
			break;
		}
		else if (allZero && SLE[i].back() == 0) 
		{
			hasInfiniteSolutions = true;
		}
		else if (!allZero && !hasNonZero) 
		{
			hasUniqueSolution = false;
		}
	}

	return make_pair(solution, hasUniqueSolution ? 0 : (hasInfiniteSolutions ? 1 : 2));
}

SystemOfLinearEquation::~SystemOfLinearEquation()
{
	SLE.clear();
	SLE.~vector();
}

void SystemOfLinearEquation::printQuit()
{
	cout << endl << endl << "Backspace - возврат в меню";

	int answ = _getch();
	bool press = false;

	while (!press)
	{
		if (answ == 8)
		{
			press = true;
			break;
		}
		else
		{
			answ = _getch();
		}
	}

	system("cls");
}