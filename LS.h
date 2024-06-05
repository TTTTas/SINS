#pragma once
#include <Eigen/Dense>
#include <iostream>
using namespace std;
using namespace Eigen;
#include"Configure.h"

class Least_Squares
{
public:
	MatrixXd X;
	MatrixXd x;
	MatrixXd B;
	MatrixXd l;
	MatrixXd P;
	MatrixXd Qxx;
	double sigma;

	Least_Squares();
	Least_Squares(Configure cfg);
	void LS();
	void ELS();
	int set_B_Pos(MatrixXd B_);
	int set_B_Vel(MatrixXd B_);
	int set_l(MatrixXd l_);
	int set_P(MatrixXd P_);
	void reset();
};