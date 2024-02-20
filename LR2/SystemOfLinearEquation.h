#pragma once
#include <iostream>
#include <fstream>

#include <Windows.h>
#include <conio.h>
#include <cmath>

#include <vector>
#include <string>

using namespace std;

class SystemOfLinearEquation
{
	int cntEquation;
	int cntVar;
	vector<vector<float>> SLE;

public:

	SystemOfLinearEquation();
	SystemOfLinearEquation(vector<vector<float>> matrix, int n);

	int getCntEq();
	int getCntVar();
	vector<vector<float>> getMatrix();

	void getFromFile(string filename);
	void printData();

	pair<vector<float>, int> gaussian_solve();

	~SystemOfLinearEquation();

	void printQuit();
};