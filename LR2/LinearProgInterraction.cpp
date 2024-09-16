#include "LinearProgInterraction.h"

void LinearProgInterraction::readFromFile(string filename)
{
	ifstream ifs(filename);

	if (ifs.is_open())
	{
		int vars, limits;
		ifs >> vars >> limits;

		float value;
		for (int i = 0; i < vars; i++)
		{
			ifs >> value;
			objFunc.push_back(value);
		}

		limitsCoefs.resize(limits);

		string sign;
		for (int i = 0; i < limits; i++)
		{
			for (int j = 0; j < vars; j++)
			{
				ifs >> value;
				limitsCoefs[i].first.push_back(value);
			}
			ifs >> sign;
			limitsCoefs[i].second = sign == "<=" ? false : true;
			ifs >> value;
			limitsConstants.push_back(value);
		}

		ifs.close();

		cout << "Успешный импорт данных";
	}
	else
	{
		cerr << ">>>Ошибка открытия файла<<<";
	}
}

void LinearProgInterraction::printData()
{
	if (objFunc.empty() || limitsCoefs.empty() || limitsConstants.empty())
	{
		cout << "Данные неполные или отсутствуют";
	}
	else
	{
		int limits = limitsCoefs.size();
		int vars = objFunc.size();
		cout << "Количество ограничений - " << limits << endl << "Количество переменных - " << vars << endl << endl;

		cout << "Целевая функция: F(X) = ";
		for (int i = 0; i < vars; i++)
		{
			if (i == vars - 1) cout << objFunc[i] << "*x" << i + 1;
			else cout << objFunc[i] << "*x" << i + 1 << " + ";
		}
		cout << endl << endl;

		cout << "Ограничения:" << endl;
		for (int i = 0; i < limits; i++)
		{
			for (int j = 0; j < vars; j++)
			{
				cout << limitsCoefs[i].first[j] << "*x" << j + 1 << " + ";
			}
			string sign = (limitsCoefs[i].second == false) ? "<=" : ">=";
			cout << sign << " " << limitsConstants[i] << endl;
		}
	}
}

void LinearProgInterraction::printQuit()
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

void LinearProgInterraction::genComb(vector<vector<int>>& result, int n, int k)
{
	if (k == n)
	{
		vector<int> temp(n);
		for (int j = 0; j < temp.size(); j++)
		{
			temp[j] = j;
		}
		//reverse(temp.begin(), temp.end());
		result.push_back(temp);

		return;
	}

	if (k == 0)
	{
		result.push_back(vector<int>());

		return;
	}

	vector<int> comb(n + 2);

	for (int j = 1; j <= k; j++)
	{
		comb[j] = j - 1;
	}
	comb[k + 1] = n;
	comb[k + 2] = 0;

	while (true)
	{
		vector<int> temp;
		copy(comb.begin() + 1, comb.begin() + k + 1, back_inserter(temp));
		//reverse(temp.begin(), temp.end());
		result.push_back(temp);

		int j = 1;
		while (comb[j + 1] == comb[j] + 1)
		{
			comb[j] = j - 1;
			j++;
		}

		if (j > k)
		{
			break;
		}

		comb[j] += 1;
	}
}

float LinearProgInterraction::getObjFuncValue(vector<float>& solution)
{
	float result = 0.0f;
	int sizeFunc = objFunc.size();
	int sizeSolution = solution.size();
	if (sizeFunc != sizeSolution)
	{
		cout << "Размеры векторов целевой функции и решения не совпадают" << endl;
		return 0.0f;
	}

	for (int i = 0; i < sizeFunc; i++)
	{
		result += objFunc[i] * solution[i];
	}

	return result;
}

bool LinearProgInterraction::compareBySecond(const pair<vector<float>, float>& a, const pair<vector<float>, float>& b)
{
	return a.second < b.second;
}

void LinearProgInterraction::bruteForceMethod(vector<pair<vector<float>, float>>& solutions)
{
	if (objFunc.empty() || limitsCoefs.empty() || limitsConstants.empty())
	{
		cerr << "Данные неполные или отсутствуют";
	}

	vector<vector<int>> combs;
	genComb(combs, limitsCoefs.size(), objFunc.size());

	int iterations = combs.size();
	int limits = limitsCoefs.size();
	int vars = objFunc.size();
	for (int i = 0; i < iterations; i++)
	{
		vector<vector<float>> curCoefs;
		vector<float> curConsts;
		vector<int> curComb = combs[i];
		bool notSolution = false;
		int combSize = combs[i].size();
		for (int j = 0; j < combSize; j++)
		{
			curCoefs.push_back(limitsCoefs[curComb[j]].first);
			curConsts.push_back(limitsConstants[curComb[j]]);
		}

		cout << "Текущая комбинация: ";
		for (auto& el : curComb)
		{
			cout << el + 1 << " ";
		}
		cout << endl;

		SystemOfLinearEquation curSLE(curCoefs, curConsts);
		curSLE.printData();

		vector<float> solution = curSLE.inverseMatrix_Solve();
		int solutionSize = solution.size();
		
		if (!solution.empty())
		{
			for (int j = 0; j < limits; j++)
			{
				if (find(curComb.begin(), curComb.end(), j) == curComb.end())
				{
					float left = 0.0f;
					float right = limitsConstants[j];
					for (int k = 0; k < vars; k++)
					{
						left += limitsCoefs[j].first[k] * solution[k];
					}

					bool sign = limitsCoefs[j].second;
					if ((sign && left < right) || (!sign && left > right))
					{
						notSolution = true;
						break;
					}
				}
			}

			if (!notSolution)
			{
				float objFuncValue = getObjFuncValue(solution);
				cout << "Значение целевой функции: F(";
				for (int i = 0; i < solutionSize; i++)
				{
					if (i == solutionSize - 1) cout << solution[i] << ") = " << objFuncValue;
					else cout << solution[i] << ", ";
				}
				cout << endl << endl;

				pair<vector<float>, float> tempAnswer = make_pair(solution, objFuncValue);
				if (find(solutions.begin(), solutions.end(), tempAnswer) == solutions.end()) solutions.push_back(tempAnswer);				
			}
			else
			{
				cout << "Не проходит по всем ограничениям" << endl << endl;
			}
		}
	}

	sort(solutions.begin(), solutions.end(), compareBySecond);

	int answerSize = solutions.size();
	cout << "Ответ: ОДР имеет " << answerSize << " вершин(ы):" << endl;
	for (int i = 0; i < answerSize; i++)
	{
		int curAnswer = solutions[i].first.size();
		cout << "F(";
		for (int j = 0; j < curAnswer; j++)
		{
			if (j == curAnswer - 1) cout << solutions[i].first[j] << ") = " << solutions[i].second;
			else cout << solutions[i].first[j] << ", ";
		}
		if (i == 0) cout << " - максимум";
		else if (i == answerSize - 1) cout << " - минимум";
		cout << endl;
	}

}