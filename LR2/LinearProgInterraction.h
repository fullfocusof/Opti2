#pragma once
#include "SystemOfLinearEquation.h"

class LinearProgInterraction
{
	vector<float> objFunc;
	vector<pair<vector<float>, bool>> limitsCoefs; // 0 - <=; 1 - >=
	vector<float> limitsConstants;

public:

	void readFromFile(string filename);
	void printData();
	void printQuit();

	void genComb(vector<vector<int>>& result, int n, int k);

	float getObjFuncValue(vector<float>& solution);
	static bool compareBySecond(const pair<vector<float>, float>& a, const pair<vector<float>, float>& b);

	void bruteForceMethod(vector<pair<vector<float>, float>>& solutions);
};