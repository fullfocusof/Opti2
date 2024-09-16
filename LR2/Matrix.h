#pragma once
#include <iostream>
#include <fstream>

#include <Windows.h>
#include <conio.h>
#include <cmath>

#include <vector>
#include <string>
#include <algorithm>
#include <utility>

using namespace std;

class Matrix
{
	vector<vector<float>> coefs;

public:

	Matrix();
	Matrix(int rows, int cols);

	vector<vector<float>> getMatrix();
	float getData(int row, int col);
	void setData(int row, int col, float value);

	float determinant_2x2(vector<vector<float>> mat);
	float determinant(vector<vector<float>> mat);

	float getDetMinor(int row, int col);
	Matrix getInverseMatrix();

	~Matrix();
};