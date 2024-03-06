#pragma once
#include "Matrix.h"

class SystemOfLinearEquation
{
	Matrix coefsSLE;
	vector<float> constants;

public:

	SystemOfLinearEquation();

	Matrix getMatrix();

	void readFromFile(string filename);
	void printData();

	float getDet();

	void gaussian_solve();
	void gaussian_solveLeadElem();

	void simpleIteration(int choice);

	void printQuit();

	~SystemOfLinearEquation();
};