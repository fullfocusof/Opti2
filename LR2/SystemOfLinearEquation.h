#pragma once
#include "Matrix.h"

class SystemOfLinearEquation
{
	Matrix coefsSLE;
	vector<float> constants;

public:

	SystemOfLinearEquation();
	SystemOfLinearEquation(vector<vector<float>>& coefs, vector<float>& consts);

	Matrix getMatrix();

	//void readFromFile(string filename);
	void printData();

	float getDet();

	vector<float> gaussian_solve();
	vector<float> gaussian_solveLeadElem();

	vector<float> inverseMatrix_Solve();

	void simpleIteration(int choice);

	//void printQuit();

	~SystemOfLinearEquation();
};