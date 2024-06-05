#include"LS.h"
#include"cal.h"
#include"GNSS_data.h"

Least_Squares::Least_Squares()
{
	X = MatrixXd::Zero(4, 1);
}

Least_Squares::Least_Squares(Configure cfg)
{
	if (cfg.SYS_num == 1)
		X = MatrixXd::Zero(4, 1);
	else if (cfg.SYS_num == 2)
		X = MatrixXd::Zero(5, 1);
}


void Least_Squares::ELS()
{
	Qxx = (B.transpose() * P * B).inverse();
	x = Qxx * B.transpose() * P * l;
	X = X + x;
	MatrixXd v = B * x - l;
	int m = B.rows() - B.cols();
	sigma = sqrt(((v.transpose() * P * v) / m)(0, 0));
}

void Least_Squares::LS()
{
	Qxx = (B.transpose() * P * B).inverse();
	X = Qxx * B.transpose() * P * l;
	MatrixXd v = B * X - l;
	int m = B.rows() - B.cols();
	sigma = sqrt(((v.transpose() * P * v) / m)(0, 0));
}

int Least_Squares::set_B_Pos(MatrixXd B_)
{
	if (B.rows() == 0)
	{
		B = B_;
		return B.rows();
	}
	int Last_Row = B.rows();
	int Last_Col = B.cols();
	int temp_Row = B_.rows();
	int temp_Col = B_.cols();
	B.conservativeResize(Last_Row + temp_Row, Last_Col + temp_Col - 3);
	B.block(Last_Row, 0, temp_Row, 3) = B_.leftCols(3);
	B.block(0, Last_Col, Last_Row, temp_Col - 3) = MatrixXd::Zero(Last_Row, temp_Col - 3);
	B.block(Last_Row, Last_Col, temp_Row, temp_Col - 3) = B_.rightCols(temp_Col - 3);
	B.block(Last_Row, 3, temp_Row, Last_Col - 3) = MatrixXd::Zero(temp_Row, Last_Col - 3);
	return B.rows();
}

int Least_Squares::set_B_Vel(MatrixXd B_)
{
	if (B.rows() == 0)
	{
		B = B_;
		return B.rows();
	}
	B.conservativeResize(B.rows() + B_.rows(), 4);
	B.bottomRows(B_.rows()) = B_;
	return B.rows();
}
int Least_Squares::set_l(MatrixXd l_)
{
	if (l.rows() == 0)
	{
		l = l_;
		return l.rows();
	}
	l.conservativeResize(l.rows() + l_.rows(), 1);
	l.bottomRows(l_.rows()) = l_;
	return l.rows();
}
int Least_Squares::set_P(MatrixXd P_)
{
	if (P.rows() == 0)
	{
		P = P_;
		return P.rows();
	}
	int Last_Row = P.rows();
	int Last_Col = P.cols();
	int temp_Row = P_.rows();
	int temp_Col = P_.cols();
	P.conservativeResize(P.rows() + P_.rows(), P.cols() + P_.cols());
	P.block(Last_Row, Last_Col, temp_Row, temp_Col) = P_;
	P.block(0, Last_Col, Last_Row, temp_Col) = MatrixXd::Zero(Last_Row, temp_Col);
	P.block(Last_Row, 0, temp_Row, Last_Col) = MatrixXd::Zero(temp_Row, Last_Col);
	return P.rows();
}

void Least_Squares::reset()
{
	B = MatrixXd::Zero(0, 0);
	l = MatrixXd::Zero(0, 0);
	P = MatrixXd::Zero(0, 0);
}