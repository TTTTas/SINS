#include "KF.h"
#include "cal.h"
#include "GNSS_data.h"

KalmanFilter::KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance)
{
	A_ = MatrixXd::Identity(initial_state.rows(), initial_state.rows());
	H_ = MatrixXd::Zero(0, 0);
	z_ = MatrixXd::Zero(0, 0);
	R_ = MatrixXd::Zero(0, 0);
	x_hat_ = initial_state;
	P_ = initial_covariance;
}

void KalmanFilter::predict()
{
	x_hat_minus_ = A_ * x_hat_;
	P_minus_ = A_ * P_ * A_.transpose() + Q_;
}

void KalmanFilter::update()
{
	K_ = P_minus_ * H_.transpose() * (H_ * P_minus_ * H_.transpose() + R_).inverse();
	x_hat_ = x_hat_minus_ + K_ * (z_);
	P_ = (MatrixXd::Identity(x_hat_.rows(), x_hat_.rows()) - K_ * H_) * P_minus_;
}

void KalmanFilter::set_A(MatrixXd A)
{
	A_ = A;
}

void KalmanFilter::set_H(MatrixXd H, int mode)
{
	if(H_.rows()==0)
		H_ = H;
	else
	{
		int Last_Row, Last_Col, temp_Row, temp_Col, sys_num;
		MatrixXd B_old_one, B_pos_add, B_vel_add, B_add_one, B_vel, B_pos;
		switch (mode)
		{
		case MODE_SPP:
			Last_Row = H_.rows() / 2;
			Last_Col = H.cols();
			temp_Row = H.rows() / 2;
			temp_Col = H.cols();
			B_pos = H_.block(0, 0, Last_Row, 3);
			B_vel = H_.block(Last_Row, 3, Last_Row, 3);
			sys_num = Last_Col - 7;
			B_old_one = H.block(0, 6, Last_Row, sys_num);
			B_pos_add = H.block(0, 0, temp_Row, 3);
			B_vel_add = H.block(temp_Row, 3, temp_Row, 3);
			B_add_one = H.block(0, 6, temp_Row, 1);

			B_pos.conservativeResize(Last_Row + temp_Row, 3);
			B_pos.bottomRows(temp_Row) = B_pos_add;

			B_vel.conservativeResize(Last_Row + temp_Row, 3);
			B_vel.bottomRows(B_vel_add.rows()) = B_vel_add;

			H_ = MatrixXd::Zero(B_pos.rows() + B_vel.rows(), 6 + sys_num + 2);
			H_.block(0, 0, B_pos.rows(), 3) = B_pos;
			H_.block(B_pos.rows(), 3, B_vel.rows(), 3) = B_vel;
			H_.block(0, 6, Last_Row, sys_num) = B_old_one;
			H_.block(Last_Row, 6 + sys_num, temp_Row, 1) = B_add_one;
			H_.block(B_pos.rows(), 7 + sys_num, temp_Row, 1) = B_add_one;
			break;
		case MODE_RTK:
			Last_Row = H_.rows();
			Last_Col = H_.cols();
			temp_Row = H.rows();
			temp_Col = H.cols();
			H_.conservativeResize(Last_Row + temp_Row, Last_Col + temp_Col - 3);
			H_.block(Last_Row, 0, temp_Row, 3) = H.leftCols(3);
			H_.block(0, Last_Col, Last_Row, temp_Col - 3) = MatrixXd::Zero(Last_Row, temp_Col - 3);
			H_.block(Last_Row, Last_Col, temp_Row, temp_Col - 3) = H.rightCols(temp_Col - 3);
			H_.block(Last_Row, 3, temp_Row, Last_Col - 3) = MatrixXd::Zero(temp_Row, Last_Col - 3);
		}
		
	}

}

void KalmanFilter::set_Z(MatrixXd z)
{
	if (z_.rows() == 0)
		z_ = z;
	else
	{
		z_.conservativeResize(z_.rows() + z.rows(), 1);
		z_.bottomRows(z.rows()) = z;
	}
}

void KalmanFilter::set_Q(MatrixXd Q)
{
	Q_ = Q;
}

void KalmanFilter::set_R(MatrixXd R, int mode)
{
	if (R_.rows() == 0)
		R_ = R;
	else
	{
		int Last_Row, Last_Col, temp_Row, temp_Col;
		MatrixXd P_pos, P_vel, P_pos_add, P_vel_add;
		switch (mode)
		{
		case MODE_SPP:
			Last_Row = R_.rows() / 2;
			temp_Row = R.rows() / 2;
			P_pos = R_.block(0, 0, Last_Row, Last_Row);
			P_vel = R_.block(Last_Row, Last_Row, Last_Row, Last_Row);
			P_pos_add = R.block(0, 0, temp_Row, temp_Row);
			P_vel_add = R.block(temp_Row, temp_Row, temp_Row, temp_Row);

			P_pos.conservativeResize(Last_Row + temp_Row, Last_Row + temp_Row);
			P_pos.block(Last_Row, Last_Row, temp_Row, temp_Row) = P_pos_add;
			P_pos.block(0, Last_Row, Last_Row, temp_Row) = MatrixXd::Zero(Last_Row, temp_Row);
			P_pos.block(Last_Row, 0, temp_Row, Last_Row) = MatrixXd::Zero(temp_Row, Last_Row);

			P_vel.conservativeResize(Last_Row + temp_Row, Last_Row + temp_Row);
			P_vel.block(Last_Row, Last_Row, temp_Row, temp_Row) = P_vel_add;
			P_vel.block(0, Last_Row, Last_Row, temp_Row) = MatrixXd::Zero(Last_Row, temp_Row);
			P_vel.block(Last_Row, 0, temp_Row, Last_Row) = MatrixXd::Zero(temp_Row, Last_Row);

			R_.conservativeResize(P_pos.rows() + P_vel.rows(), P_pos.cols() + P_vel.cols());
			R_.block(0, 0, P_pos.rows(), P_pos.cols()) = P_pos;
			R_.block(P_pos.rows(), P_pos.cols(), P_vel.rows(), P_vel.cols()) = P_vel;
			R_.block(0, P_pos.cols(), P_pos.rows(), P_vel.cols()) = MatrixXd::Zero(P_pos.rows(), P_vel.cols());
			R_.block(P_pos.rows(), 0, P_vel.rows(), P_pos.cols()) = MatrixXd::Zero(P_vel.rows(), P_pos.cols());
			break;
		case MODE_RTK:
			Last_Row = R_.rows();
			Last_Col = R_.cols();
			temp_Row = R.rows();
			temp_Col = R.cols();
			R_.conservativeResize(R_.rows() + R.rows(), R_.cols() + R.cols());
			R_.block(Last_Row, Last_Col, temp_Row, temp_Col) = R;
			R_.block(0, Last_Col, Last_Row, temp_Col) = MatrixXd::Zero(Last_Row, temp_Col);
			R_.block(Last_Row, 0, temp_Row, Last_Col) = MatrixXd::Zero(temp_Row, Last_Col);
			break;
		default:break;
		}

	}
}

void KalmanFilter::setState(MatrixXd state)
{
	x_hat_ = state;
}

void KalmanFilter::reset()
{
	H_ = MatrixXd::Zero(0, 0);
	z_ = MatrixXd::Zero(0, 0);
	R_ = MatrixXd::Zero(0, 0);
}

MatrixXd KalmanFilter::getState() const
{
	return x_hat_;
}

MatrixXd KalmanFilter::getState_minus()
{
	return x_hat_minus_;
}

MatrixXd getA(double T, GNSS_Configure cfg)
{
	int row = cfg.SYS_num + 7;
	MatrixXd A = MatrixXd::Identity(row, row);
	A.block(0, 3, 3, 3) = MatrixXd::Identity(3, 3) * T;
	if (cfg.SYS_num == 1)
		A(6, 7) = T;
	else if (cfg.SYS_num == 2)
	{
		A(6, 8) = T;
		A(7, 8) = T;
	}
	return A;
}

MatrixXd getH(MatrixXd B)
{
	MatrixXd H = MatrixXd::Zero(B.rows() + B.rows(), B.cols() + B.cols());
	H.block(0, 0, B.rows(), B.cols() - 1) = B.leftCols(3);
	H.block(0, 6, B.rows(), 1) = B.rightCols(1);
	H.block(B.rows(), 3, B.rows(), 3) = B.leftCols(3);
	H.block(B.rows(), 7, B.rows(), 1) = B.rightCols(1);
	return H;
}

MatrixXd getQ(double T, GNSS_Configure cfg)
{
	MatrixXd Sv = MatrixXd::Identity(3, 3) * 0.05;
	double St = 0.05;
	double Sf = 0.05;
	int row = cfg.SYS_num + 7;
	MatrixXd Q = MatrixXd::Zero(row, row);
	Q.block(0, 0, 3, 3) = Sv * pow(T, 3) / 3;
	Q.block(0, 3, 3, 3) = Sv * pow(T, 2) / 2;
	Q.block(3, 0, 3, 3) = Sv * pow(T, 2) / 2;
	Q.block(3, 3, 3, 3) = Sv * T;
	if (cfg.SYS_num == 1)
	{
		Q(6, 6) = St * T + Sf * pow(T, 3) / 3;
		Q(6, 7) = Sf * pow(T, 2) / 2;
		Q(7, 6) = Sf * pow(T, 2) / 2;
		Q(7, 7) = Sf * T;
	}
	else if (cfg.SYS_num == 2)
	{
		Q(6, 6) = St * T + Sf * pow(T, 3) / 3;
		Q(6, 8) = Sf * pow(T, 2) / 2;
		Q(8, 6) = Sf * pow(T, 2) / 2;
		Q(7, 7) = St * T + Sf * pow(T, 3) / 3;
		Q(7, 8) = Sf * pow(T, 2) / 2;
		Q(8, 7) = Sf * pow(T, 2) / 2;
		Q(8, 8) = Sf * T;
	}
	return Q;
}

MatrixXd getR(double ROW)
{
	MatrixXd R = MatrixXd::Identity(ROW + ROW, ROW + ROW);
	R.block(0, 0, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 4.5;
	R.block(ROW, ROW, ROW, ROW) = MatrixXd::Identity(ROW, ROW) * 0.18;
	return R;
}

MatrixXd getz(MatrixXd l_P, MatrixXd l_V)
{
	MatrixXd z = MatrixXd::Zero(l_P.rows() + l_V.rows(), 1);
	z.block(0, 0, l_P.rows(), 1) = l_P;
	z.block(l_P.rows(), 0, l_V.rows(), 1) = l_V;
	return z;
}

int RTK_getA(vector<int> old_prn, map<int, int>& new_prn, MatrixXd& A, int new_ref)
{
	A = MatrixXd::Zero(new_prn.size(), old_prn.size());
	map<int, int>::iterator old2new_it;
	int val = -1;
	for(int i=0;i<old_prn.size();i++)
	{
		old2new_it = new_prn.find(old_prn[i]);
		if(old2new_it==new_prn.end())
		{
			if (old_prn[i] == new_ref)
				val = i;
			continue;
		}		
		A(old2new_it->second, i) = 1;
		new_prn.erase(old2new_it);
	}
	for(old2new_it=new_prn.begin();old2new_it!=new_prn.end();old2new_it++)
	{
		if (find(old_prn.begin(), old_prn.end(), old2new_it->first) != old_prn.end())
			new_prn.erase(old2new_it);
	}
	return val;
}


MatrixXd Cal_V(RTK_DATA* rtk, GNSS_Configure cfg, XYZ* rove_pos)
{
	MatrixXd V = MatrixXd::Zero(0, 1);
	int ROWS = 0;
	vector<Satellate*> Rove_Sates;
	vector<Satellate*> Base_Sates;
	if (cfg.GPS_Cfg.used)
	{
		Rove_Sates = rtk->Rove_GPS;
		Base_Sates = rtk->Base_GPS;
		MatrixXd B_Rove = MatrixXd::Zero(0, 3);
		MatrixXd l_Rove = MatrixXd::Zero(0, 1);
		MatrixXd Q_Rove = MatrixXd::Zero(1, 1);
		MatrixXd B_Base = MatrixXd::Zero(0, 3);
		MatrixXd l_Base = MatrixXd::Zero(0, 1);
		MatrixXd Q_Base = MatrixXd::Zero(1, 1);
		vector<int> prn_used;
		get_mat(B_Base, l_Base, Q_Base, rtk->OBSTIME, rtk->Base_appro_pos, Base_Sates, SYS_GPS, cfg, rtk, prn_used);
		prn_used.clear();
		int elev_row = get_mat(B_Rove, l_Rove, Q_Rove, rtk->OBSTIME, rove_pos, Rove_Sates, SYS_GPS, cfg, rtk, prn_used);
		vector<MatrixXd> l_Rove_max_elev_row;
		vector<MatrixXd> l_Base_max_elev_row;
		ROWS = B_Rove.rows(); 
		for (int j = 0; j < 2 * cfg.phase_num; j++)
		{
			l_Rove_max_elev_row.push_back(l_Rove.row(elev_row + j * ROWS));
			l_Base_max_elev_row.push_back(l_Base.row(elev_row + j * ROWS));
		}
		for (int i = 0; i < ROWS; i++)
		{
			for (int j = 0; j < 2 * cfg.phase_num; j++)
			{
				l_Rove.row(i + j * ROWS) = -l_Rove_max_elev_row[j] + l_Rove.row(i + j * ROWS);
				l_Base.row(i + j * ROWS) = -l_Base_max_elev_row[j] + l_Base.row(i + j * ROWS);
			}
		}
		for (int i = 0; i < 2 * cfg.phase_num; i++)
		{
			RemoveRow(l_Rove, elev_row + i * ROWS - i);
			RemoveRow(l_Base, elev_row + i * ROWS - i);
		}
		MatrixXd l_dd = l_Rove - l_Base;
		if (V.rows() == 0)
			V = l_dd;
		else
		{
			V.conservativeResize(V.rows() + l_dd.rows(), 1);
			V.bottomRows(l_dd.rows()) = l_dd;
		}
	}
	if (cfg.BDS_Cfg.used)
	{
		Rove_Sates = rtk->Rove_BDS;
		Base_Sates = rtk->Base_BDS;
		MatrixXd B_Rove = MatrixXd::Zero(0, 3);
		MatrixXd l_Rove = MatrixXd::Zero(0, 1);
		MatrixXd Q_Rove = MatrixXd::Zero(1, 1);
		MatrixXd B_Base = MatrixXd::Zero(0, 3);
		MatrixXd l_Base = MatrixXd::Zero(0, 1);
		MatrixXd Q_Base = MatrixXd::Zero(1, 1);
		vector<int> prn_used;
		get_mat(B_Base, l_Base, Q_Base, rtk->OBSTIME, rtk->Base_appro_pos, Base_Sates, SYS_BDS, cfg, rtk, prn_used);
		prn_used.clear();
		int elev_row = get_mat(B_Rove, l_Rove, Q_Rove, rtk->OBSTIME, rove_pos, Rove_Sates, SYS_BDS, cfg, rtk, prn_used);
		vector<MatrixXd> l_Rove_max_elev_row;
		vector<MatrixXd> l_Base_max_elev_row;
		ROWS = B_Rove.rows();
		for (int j = 0; j < 2 * cfg.phase_num; j++)
		{
			l_Rove_max_elev_row.push_back(l_Rove.row(elev_row + j * ROWS));
			l_Base_max_elev_row.push_back(l_Base.row(elev_row + j * ROWS));
		}
		for (int i = 0; i < ROWS; i++)
		{
			for (int j = 0; j < 2 * cfg.phase_num; j++)
			{
				l_Rove.row(i + j * ROWS) = -l_Rove_max_elev_row[j] + l_Rove.row(i + j * ROWS);
				l_Base.row(i + j * ROWS) = -l_Base_max_elev_row[j] + l_Base.row(i + j * ROWS);
			}
		}
		for (int i = 0; i < 2 * cfg.phase_num; i++)
		{
			RemoveRow(l_Rove, elev_row + i * ROWS - i);
			RemoveRow(l_Base, elev_row + i * ROWS - i);
		}
		MatrixXd l_dd = l_Rove - l_Base;
		if (V.rows() == 0)
			V = l_dd;
		else
		{
			V.conservativeResize(V.rows() + l_dd.rows(), 1);
			V.bottomRows(l_dd.rows()) = l_dd;
		}
	}
	return V;
}