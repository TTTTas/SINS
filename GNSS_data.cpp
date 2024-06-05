#include "GNSS_data.h"
#include "cal.h"
#include <io.h>
#include <direct.h>

EPHEMERIS* DATA_SET::GPS_eph[GPS_SAT_QUAN];
EPHEMERIS* DATA_SET::BDS_eph[BDS_SAT_QUAN];

void DATA_SET::inital_eph()
{
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_eph[i] = new EPHEMERIS();
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_eph[i] = new EPHEMERIS();
	}
}

DATA_SET::DATA_SET(Configure cfg)
{
	OBSTIME = new GPSTIME();

	LS_SATES = new string();
	KF_SATES = new string();

	LS_GPS_num = 0;
	LS_BDS_num = 0;
	KF_GPS_num = 0;
	KF_BDS_num = 0;
	*LS_SATES = "";
	*KF_SATES = "";

	range = new OBS_DATA();
	range->OBS_TIME = new GPSTIME();
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_GF[i] = 0;
		GPS_MW[i] = 0;
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_GF[i] = 0;
		BDS_MW[i] = 0;
	}
	inital_eph();
	LS_result = UN_Solve;
	KF_result = UN_Solve;

	LS_first = true;
	KF_first = true;
	temp_ref = MatrixXd::Zero(4, 1);
	Pos = MatrixXd::Zero(3, 1);

	LS_Pos = new Least_Squares(cfg);
	LS_Vel = new Least_Squares();

	int row = 7 + cfg.SYS_num;
	MatrixXd P_ = MatrixXd::Zero(row, row);
	P_.block(0, 0, 3, 3) = MatrixXd::Identity(3, 3) * 0.05;
	P_.block(3, 3, 3, 3) = MatrixXd::Identity(3, 3) * 0.05;
	if (cfg.SYS_num == 1)
	{
		P_(6, 6) = 0.05;
		P_(7, 7) = 0.05;
	}
	else if (cfg.SYS_num == 2)
	{
		P_(6, 6) = 0.05;
		P_(7, 7) = 0.05;
		P_(8, 8) = 0.05;
	}
	KF = new KalmanFilter(MatrixXd::Zero(7 + cfg.SYS_num, 1), P_);
	Real_Pos = new XYZ(-2267807.853, 5009320.431, 3221020.875);
}

void DATA_SET::Sate_pos_pre(double t, Configure cfg)
{
	for(Satellate* sate:range->GPS_SATE)
	{
		int prn = sate->PRN;

		// 检查星历数据是否匹配
		if (GPS_eph[prn - 1]->PRN != prn)
			continue;

		// 计算卫星位置、钟差
		double ts = OBSTIME->SecOfWeek - sate->PSERA[0] / velocity_c;
		int week = OBSTIME->Week;
		double dt = abs(ts - GPS_eph[prn - 1]->toe_tow + (week - GPS_eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double clk = 0;
		double clk1 = 0;
		XYZ sate_pos, sate_pos1;
		// 计算卫星位置和速度
		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, GPS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk, GPS_eph[prn - 1], &sate_pos, clk, sate->PSERA[0] / velocity_c, sate->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, GPS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + 1e-3, GPS_eph[prn - 1], &sate_pos1, clk1, sate->PSERA[0] / velocity_c,
				sate->SYS);
		}

		// 计算卫星速度
		sate->pos = sate_pos;
		sate->vel = (sate_pos1 - sate_pos) / 1e-3;
		sate->clk = clk;
		sate->clk_vel = (clk1 - clk) / 1e-3;
	}
	for (Satellate* sate : range->BDS_SATE)
	{
		int prn = sate->PRN;

		// 检查星历数据是否匹配
		if (BDS_eph[prn - 1]->PRN != prn)
			continue;

		// 计算卫星位置、钟差
		double ts = OBSTIME->SecOfWeek - sate->PSERA[0] / velocity_c;
		int week = OBSTIME->Week;
		ts -= 14;
		week -= 1356;
		double dt = abs(ts - BDS_eph[prn - 1]->toe_tow + (week - BDS_eph[prn - 1]->toe_wn) * 604800);
		if (dt > 14400)
			continue;

		double clk = 0;
		double clk1 = 0;
		XYZ sate_pos, sate_pos1;
		// 计算卫星位置和速度
		for (int j = 0; j < 3; j++)
		{
			clk = CORRECT_CLK(ts - clk, BDS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk, BDS_eph[prn - 1], &sate_pos, clk, sate->PSERA[0] / velocity_c, sate->SYS);
			clk1 = CORRECT_CLK(ts - clk1 + 1e-3, BDS_eph[prn - 1]);
			SAT_POS_CAL(ts - clk1 + 1e-3, BDS_eph[prn - 1], &sate_pos1, clk1, sate->PSERA[0] / velocity_c,
				sate->SYS);
		}

		// 计算卫星速度
		sate->pos = sate_pos;
		sate->vel = (sate_pos1 - sate_pos) / 1e-3;
		sate->clk = clk;
		sate->clk_vel = (clk1 - clk) / 1e-3;
	}
}

void DATA_SET::reset()
{
	if (!range->GPS_SATE.empty() && !range->BDS_SATE.empty())
	{
		for (vector<Satellate*>::iterator it = range->GPS_SATE.begin(); it != range->GPS_SATE.end(); it++)
		{
			if (NULL != *it)
			{
				delete*it;
				*it = NULL;
			}
		}
		range->GPS_SATE.clear();
		for (vector<Satellate*>::iterator it = range->BDS_SATE.begin(); it != range->BDS_SATE.end(); it++)
		{
			if (NULL != *it)
			{
				delete*it;
				*it = NULL;
			}
		}
		range->BDS_SATE.clear();
	}
}

int DATA_SET::LS_print(Configure cfg)
{
	XYZ xyz = get_XYZ(Pos);
	double dt_G = 0;
	double dt_C = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = LS_Pos->X(3, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = LS_Pos->X(3, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = LS_Pos->X(3, 0);
		dt_C = LS_Pos->X(4, 0);
	}
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V;
	switch (LS_result)
	{
	case UN_Solve:
		printf("GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Success_Solve:
		Q_local = R * ((LS_Pos->Qxx).block(0, 0, 3, 3)) * R.transpose();
		m_H = (LS_Pos->sigma) * sqrt(Q_local(2, 2));
		m_V = (LS_Pos->sigma) * sqrt(Q_local(0, 0) + Q_local(1, 1));
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %8.4f\t%8.4f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
		       OBSTIME->Week, OBSTIME->SecOfWeek,
		       xyz.X, xyz.Y, xyz.Z,
		       blh.Lat, blh.Lon, blh.Height,
		       enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			printf("GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			printf("BDS Clk : % 7.4f\t", dt_C);
		printf(
			"m_H : % 7.4f\tm_V : % 7.4f\tVelocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tsigma_P: % 6.4f\tsigma_V : % 6.4f\tPDOP : % 6.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			m_H, m_V,
			LS_Vel->X(0, 0), LS_Vel->X(1, 0), LS_Vel->X(2, 0), LS_Vel->X(3, 0),
			LS_Pos->sigma, LS_Vel->sigma, Cal_PDOP(LS_Pos->Qxx),
			LS_GPS_num, LS_BDS_num,
			LS_SATES->c_str());
		return 1;
		break;
	case OBS_DATA_Loss:
		printf("GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Epoch_Loss:
		printf("GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Set_UP_B_fail:
		printf("GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	default:
		break;
	}
	return 0;
}

int DATA_SET::LS_Filewrite(FILE* fpr, Configure cfg)
{
	XYZ xyz = get_XYZ(Pos);
	double dt_G = 0;
	double dt_C = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = LS_Pos->X(3, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = LS_Pos->X(3, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = LS_Pos->X(3, 0);
		dt_C = LS_Pos->X(4, 0);
	}
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	MatrixXd R = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	MatrixXd Q_local;
	double m_H, m_V;
	switch (LS_result)
	{
	case UN_Solve:
		//fprintf(fpr, "GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Success_Solve:
		Q_local = R * ((LS_Pos->Qxx).block(0, 0, 3, 3)) * R.transpose();
		m_H = (LS_Pos->sigma) * sqrt(Q_local(2, 2));
		m_V = (LS_Pos->sigma) * sqrt(Q_local(0, 0) + Q_local(1, 1));
		fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %10.6f\t%10.6f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
		        OBSTIME->Week, OBSTIME->SecOfWeek,
		        xyz.X, xyz.Y, xyz.Z,
		        blh.Lat, blh.Lon, blh.Height,
		        enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			fprintf(fpr, "GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			fprintf(fpr, "BDS Clk : % 7.4f\t", dt_C);
		fprintf(
			fpr,
			"m_H : % 7.4f\tm_V : % 7.4f\tVelocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tsigma_P: % 6.4f\tsigma_V : % 6.4f\tPDOP : % 6.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
			m_H, m_V,
			LS_Vel->X(0, 0), LS_Vel->X(1, 0), LS_Vel->X(2, 0), LS_Vel->X(3, 0),
			LS_Pos->sigma, LS_Vel->sigma, Cal_PDOP(LS_Pos->Qxx),
			LS_GPS_num, LS_BDS_num,
			LS_SATES->c_str());
		return 1;
		break;
	case OBS_DATA_Loss:
		//fprintf(fpr, "GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Epoch_Loss:
		//fprintf(fpr, "GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	case Set_UP_B_fail:
		//fprintf(fpr, "GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		return 0;
	default:
		break;
	}
	return 0;
}

void DATA_SET::KF_Print(FILE* fpr, Configure cfg)
{
	XYZ xyz = get_XYZ(KF->getState().block(0, 0, 3, 1));
	BLH blh = XYZ2BLH(xyz, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Real_Pos, xyz, SYS_GPS);
	double dt_G = 0;
	double dt_C = 0;
	XYZ Vel = get_XYZ(KF->getState().block(3, 0, 3, 1));
	double Rcv_t_v = 0;
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = KF->getState()(6, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = KF->getState()(6, 0);
		Rcv_t_v = KF->getState()(7, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = KF->getState()(6, 0);
		dt_C = KF->getState()(7, 0);
		Rcv_t_v = KF->getState()(8, 0);
	}
	switch (KF_result)
	{
	case UN_Solve:
		printf("GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	//fprintf(fpr, "GPSTIME: %d\t%.3f\tUN_Solve\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		break;
	case Success_Solve:
		printf("GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %10.6f\t%10.6f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
		       OBSTIME->Week, OBSTIME->SecOfWeek,
		       xyz.X, xyz.Y, xyz.Z,
		       blh.Lat, blh.Lon, blh.Height,
		       enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			printf("GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			printf("BDS Clk : % 7.4f\t", dt_C);
		printf("Velocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
		       Vel.X, Vel.Y, Vel.Z, Rcv_t_v,
		       KF_GPS_num, KF_BDS_num,
		       KF_SATES->c_str());

		fprintf(fpr, "GPSTIME: %d\t%.3f\tXYZ: %.4f\t%.4f\t%.4f\tBLH: %10.6f\t%10.6f\t%7.4f\tENU: %7.4f\t%7.4f\t%7.4f\t",
		        OBSTIME->Week, OBSTIME->SecOfWeek,
		        xyz.X, xyz.Y, xyz.Z,
		        blh.Lat, blh.Lon, blh.Height,
		        enu.X, enu.Y, enu.Z);
		if (cfg.GPS_Cfg.used)
			fprintf(fpr, "GPS Clk : % 7.4f\t", dt_G);
		if (cfg.BDS_Cfg.used)
			fprintf(fpr, "BDS Clk : % 7.4f\t", dt_C);
		fprintf(fpr, "Velocity : % 8.4f\t % 8.4f\t % 8.4f\t % 8.4f\tGPS : % 2d\tBDS : % 2d\t % s\n",
		        Vel.X, Vel.Y, Vel.Z, Rcv_t_v,
		        KF_GPS_num, KF_BDS_num,
		        KF_SATES->c_str());
		break;
	case OBS_DATA_Loss:
		printf("GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	//fprintf(fpr, "GPSTIME: %d\t%.3f\tOBS_DATA_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		break;
	case Epoch_Loss:
		printf("GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	//fprintf(fpr, "GPSTIME: %d\t%.3f\tEpoch_Loss\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		break;
	case Set_UP_B_fail:
		printf("GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
	//fprintf(fpr, "GPSTIME: %d\t%.3f\tSet_UP_B_fail\n", OBSTIME->Week, OBSTIME->SecOfWeek);
		break;
	default:
		break;
	}
}


int createDirectory(string path)
{
	int len = path.length();
	char tmpDirPath[256] = {0};
	for (int i = 0; i < len; i++)
	{
		tmpDirPath[i] = path[i];
		if (tmpDirPath[i] == '\\' || tmpDirPath[i] == '/')
		{
			if (_access(tmpDirPath, 0) == -1)
			{
				int ret = _mkdir(tmpDirPath);
				if (ret == -1)
					return ret;
			}
		}
	}
	return 0;
}

int DATA_SET::CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag)
{
	//函数可以优化
	int prn = sate->PRN;
	int Index = decode_SYN(sys, sate->SYG_TYPE[index]);
	double f = CODE2FREQ(Index);
	double lamda = 1e-6 * velocity_c / f;
	double C_D = 0;
	double mean_D = 0;
	double C_P = 0;
	double C_L = 0;
	if (sate->DOPPLER[0] == 0)
	{
		cout << "DOPPLER DATA LOSS!" << endl;
		return 0;
	}
	switch (sys)
	{
	case SYS_GPS:
		Index -= 1;
		if (GPS_PHA[Index][prn - 1] == 0 && GPS_PSE[Index][prn - 1] == 0 && GPS_DOP[Index][prn - 1] == 0)
		{
			GPS_PHA[Index][prn - 1] = sate->PHASE[index];
			GPS_PSE[Index][prn - 1] = sate->PSERA[index];
			GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(GPS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (GPS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - GPS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - GPS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		GPS_PSE[Index][prn - 1] = sate->PSERA[index];
		GPS_PHA[Index][prn - 1] = sate->PHASE[index];
		GPS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	case SYS_BDS:
		Index -= 7;
		if (BDS_PHA[Index][prn - 1] == 0 && BDS_PSE[Index][prn - 1] == 0 && BDS_DOP[Index][prn - 1] == 0)
		{
			BDS_PHA[Index][prn - 1] = sate->PHASE[index];
			BDS_PSE[Index][prn - 1] = sate->PSERA[index];
			BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];
			return 1;
		}
		C_D = lamda * abs(BDS_DOP[Index][prn - 1] - sate->DOPPLER[index]);
		if (C_D > 20)
			return -1;
		mean_D = (BDS_DOP[Index][prn - 1] + sate->DOPPLER[index]) / 2;
		C_P = abs((sate->PSERA[index] - BDS_PSE[Index][prn - 1]) + lamda * mean_D * t);
		C_L = abs(lamda * (sate->PHASE[index] - BDS_PHA[Index][prn - 1]) + lamda * mean_D * t);
		if (C_P > 8)
			PSE_flag = false;
		if (C_L > 0.5)
			PHA_flag = false;
		BDS_PSE[Index][prn - 1] = sate->PSERA[index];
		BDS_PHA[Index][prn - 1] = sate->PHASE[index];
		BDS_DOP[Index][prn - 1] = sate->DOPPLER[index];

		return 1;
	default:
		return 0;
		break;
	}
}

int DATA_SET::DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2)
{
	int prn = sate->PRN;
	double f1 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index1]));
	double f2 = CODE2FREQ(decode_SYN(sys, sate->SYG_TYPE[index2]));
	double lamda1 = 1e-6 * velocity_c / f1;
	double lamda2 = 1e-6 * velocity_c / f2;
	bool PSE_flag1 = true;
	bool PHA_flag1 = true;
	bool PSE_flag2 = true;
	bool PHA_flag2 = true;
	double GF = 0;
	double MW = 0;
	double dGF = 0;
	double dMW = 0;
	sate->Cycle_jump = false;
	switch (sys)
	{
	case SYS_GPS:
		GF = f1 * sate->PHASE[index1] - f2 * sate->PHASE[index2];
		MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[
			index1] - sate->PHASE[index2]);
		if (GPS_GF[prn - 1] == 0 && GPS_MW[prn - 1] == 0)
		{
			GPS_GF[prn - 1] = GF;
			GPS_MW[prn - 1] = MW;
			GPS_COUNT[prn - 1]++;
			return 1;
		}
		dGF = abs(GF - GPS_GF[prn - 1]);
		dMW = abs(MW - GPS_MW[prn - 1]);
		GPS_GF[prn - 1] = GF;
		GPS_MW[prn - 1] = (GPS_MW[prn - 1] * (GPS_COUNT[prn - 1]++) + MW);
		GPS_MW[prn - 1] /= GPS_COUNT[prn - 1];
		if (dGF > GF_THRESH)
			sate->Cycle_jump = true;
		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(
			sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{

			return 1;
		}
		else
		{
			return 0;
		}
	case SYS_BDS:
		GF = f1 * sate->PHASE[index1] - f2 * sate->PHASE[index2];
		MW = (f1 - f2) * (sate->PSERA[index1] / lamda1 + sate->PSERA[index2] / lamda2) / (f1 + f2) - (sate->PHASE[
			index1] - sate->PHASE[index2]);
		if (BDS_GF[prn - 1] == 0 && BDS_MW[prn - 1] == 0)
		{
			BDS_GF[prn - 1] = GF;
			BDS_MW[prn - 1] = MW;
			BDS_COUNT[prn - 1]++;
			return 1;
		}
		dGF = abs(GF - BDS_GF[prn - 1]);
		dMW = abs(MW - BDS_MW[prn - 1]);
		BDS_GF[prn - 1] = GF;
		BDS_MW[prn - 1] = (BDS_MW[prn - 1] * (BDS_COUNT[prn - 1]++) + MW);
		BDS_MW[prn - 1] /= BDS_COUNT[prn - 1];
		if (dGF > GF_THRESH)
			sate->Cycle_jump = true;

		if (CheckOBSConsist(sate, sys, t, index1, PSE_flag1, PHA_flag1) && CheckOBSConsist(
			sate, sys, t, index2, PSE_flag2, PHA_flag2) && PSE_flag1 && PHA_flag1 && PSE_flag2 && PHA_flag2)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	default:
		return 0;
		break;
	}
}

void DATA_SET::DetectOut(Configure cfg, double dt_e)
{
	for (int i = 0; i < range->GPS_SATE.size(); i++)
	{
		if (range->GPS_SATE[i]->Phase_NUM < 2) //需修改
			continue;
		int prn = range->GPS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(range->GPS_SATE[i]->SYS, range->GPS_SATE[i]->SYG_TYPE[Index1])) != cfg.GPS_Cfg.f1 &&
			Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(range->GPS_SATE[i]->SYS, range->GPS_SATE[i]->SYG_TYPE[Index2])) != cfg.GPS_Cfg.f2 &&
			Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(range->GPS_SATE[i], range->GPS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM ||
			Index2 == MAXNUM)
		{
			range->GPS_SATE[i]->Outlier = true;
			continue;
		}
	}
	for (int i = 0; i < range->BDS_SATE.size(); i++)
	{
		if (range->BDS_SATE[i]->Phase_NUM < 2)
			continue;
		int prn = range->BDS_SATE[i]->PRN;
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(range->BDS_SATE[i]->SYS, range->BDS_SATE[i]->SYG_TYPE[Index1])) != cfg.BDS_Cfg.f1 &&
			Index1 < MAXNUM)
			Index1++;
		while (CODE2FREQ(decode_SYN(range->BDS_SATE[i]->SYS, range->BDS_SATE[i]->SYG_TYPE[Index2])) != cfg.BDS_Cfg.f2 &&
			Index2 < MAXNUM)
			Index2++;

		if (!DetectOutlier(range->BDS_SATE[i], range->BDS_SATE[i]->SYS, dt_e, Index1, Index2) || Index1 == MAXNUM ||
			Index2 == MAXNUM)
		{
			range->BDS_SATE[i]->Outlier = true;
			continue;
		}
	}
}

RTK_DATA::RTK_DATA(Configure cfg)
{
	OBSTIME = new GPSTIME();
	LS_SATES = new string();
	KF_SATES = new string();
	LS_GPS_num = 0;
	LS_BDS_num = 0;

	LS_IS_FIXED = false;
	KF_IS_FIXED = false;

	Ratio_LS = 0;
	Ratio_KF = 0;

	Ref_PRN_GPS_NEW = -1;
	Ref_PRN_BDS_NEW = -1;
	Ref_PRN_BDS_OLD = 0;
	Ref_PRN_GPS_OLD = 0;

	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		GPS_GF[i] = 0;
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		BDS_GF[i] = 0;
	}

	*LS_SATES = "";
	*KF_SATES = "";

	LS_result = 0;
	KF_result = 0;

	Base_appro_pos = new XYZ();
	Rove_appro_pos = new XYZ();
	KF_pos = new XYZ();
	Q_baseLine = MatrixXd::Zero(3, 3);
	Fix_N = MatrixXd::Zero(0, 1);
	Stable_Num = 0;

	LS = new Least_Squares();
	KF = new KalmanFilter(MatrixXd::Zero(3, 1), MatrixXd::Identity(3, 3) * SQR(cfg.initial_covariance));
}

void RTK_DATA::reset()
{
	Base_GPS.clear();
	Base_BDS.clear();
	Rove_GPS.clear();
	Rove_BDS.clear();

	LS_GPS_num = 0;
	LS_BDS_num = 0;

	LS_IS_FIXED = false;
	KF_IS_FIXED = false;
}

void RTK_DATA::LS_Output(FILE* fpr, Configure cfg)
{
	XYZ base_line = *Rove_appro_pos - *Base_appro_pos;
	BLH blh = XYZ2BLH(*Rove_appro_pos, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Base_appro_pos, *Rove_appro_pos, SYS_GPS);
	MatrixXd Qxx = LS->Qxx.block(0, 0, 3, 3);
	double Base_Len = Len(&base_line);
	MatrixXd R = get_XYZ(base_line / Base_Len).transpose();
	MatrixXd K = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	double Q_s = (R * Qxx * R.transpose())(0, 0);
	MatrixXd Qnn = K * Qxx * K.transpose();
	MatrixXd V = Cal_V(this, cfg, Rove_appro_pos);
	V -= LS->B.rightCols(LS->B.cols() - 3) * LS->X.bottomRows(LS->X.rows() - 3);
	double sigma = sqrt((V.transpose() * LS->P * V / (LS->B.rows() - LS->B.cols()))(0, 0));
	MatrixXd m_xx = Qxx * sigma;
	MatrixXd m_nn = Qnn * sigma;
	int gps = cfg.GPS_Cfg.used;
	int bds = cfg.BDS_Cfg.used;
	MatrixXd V_L = MatrixXd::Zero((gps * LS_GPS_num + bds * LS_BDS_num) * cfg.phase_num, 1);	
	MatrixXd V_P = MatrixXd::Zero((gps * LS_GPS_num + bds * LS_BDS_num) * cfg.phase_num, 1);
	if (gps)
	{
		//GPS系统
		V_L.topRows(KF_GPS_num * cfg.phase_num) = V.topRows(KF_GPS_num * cfg.phase_num);
		V_P.topRows(KF_GPS_num * cfg.phase_num) = V.block(KF_GPS_num * cfg.phase_num, 0, KF_GPS_num * cfg.phase_num, 1);
	}
	if (bds)
	{
		//BDS系统
		V_L.bottomRows(KF_BDS_num * cfg.phase_num) = V.block(2 * gps * KF_GPS_num * cfg.phase_num, 0, KF_BDS_num * cfg.phase_num, 1);
		V_P.bottomRows(KF_BDS_num * cfg.phase_num) = V.block(2 * gps * KF_GPS_num * cfg.phase_num + LS_BDS_num * cfg.phase_num, 0, KF_BDS_num * cfg.phase_num, 1);
	}
	double RDOP = sqrt(Qxx.trace());
	double RMS = sqrt((V.transpose() * V)(0, 0) / V.rows());
	double RMS_L = sqrt((V_L.transpose() * V_L)(0, 0) / V_L.rows());
	double RMS_P = sqrt((V_P.transpose() * V_P)(0, 0) / V_P.rows());
	double PDOP = Cal_PDOP(LS->Qxx.topLeftCorner(3, 3));
	string is_fixed;
	if (LS_IS_FIXED)
		is_fixed = "Fix";
	else
		is_fixed = "Float";
	printf("Mode: RTK\tMethod: LS\tGPSTIME: %d\t%.3f\t  %s\n", OBSTIME->Week, OBSTIME->SecOfWeek, is_fixed.c_str());
	fprintf(fpr, "Mode: RTK\tMethod: LS\tGPSTIME: %d\t%.3f\t  %s\n", OBSTIME->Week, OBSTIME->SecOfWeek,
	        is_fixed.c_str());
	switch (LS_result)
	{
	case UN_Solve:
		printf("UN_Solve\n");
		break;
	case Success_Solve:
		printf("Rove XYZ: %.4f\t%.4f\t%.4f\tBase_Line XYZ: %.6f\t%.6f\t%.6f\n"
		       "Rove BLH: %10.6f\t%10.6f\t%7.4f\t\tBase_Line ENU: %.6f\t%.6f\t%.6f\n"
		       "m_x: % 1.6f\tm_y: % 1.6f\t\m_z: % 1.6f\tm_e: % 1.6f\tm_n: % 1.6f\t\m_u: % 1.6f\n"
		       "Ratio: %.4f\tQ_s: %1.6f\tRMS: %.4f\tRMS_L: %.4f\tRMS_P: %.4f\tRDOP: %.6f\tPDOP: %.6f\tSigma: %.6f\n"
			"GPS_REF: %2d\tBDS_REF: %2d\tGPS: % 2d\tBDS: % 2d\t % s\n",
		       Rove_appro_pos->X, Rove_appro_pos->Y, Rove_appro_pos->Z,
		       base_line.X, base_line.Y, base_line.Z,
		       blh.Lat, blh.Lon, blh.Height,
		       enu.X, enu.Y, enu.Z,
		       sqrt(m_xx(0, 0)), sqrt(m_xx(1, 1)), sqrt(m_xx(2, 2)),
		       sqrt(m_nn(0, 0)), sqrt(m_nn(1, 1)), sqrt(m_nn(2, 2)),
		       Ratio_LS, Q_s, RMS, RMS_L, RMS_P, RDOP, PDOP, sigma, 
		       Ref_PRN_GPS_NEW, Ref_PRN_BDS_NEW,
		       LS_GPS_num, LS_BDS_num,
		       LS_SATES->c_str());

		fprintf(fpr,"Rove XYZ: %.4f\t%.4f\t%.4f\tBase_Line XYZ: %.6f\t%.6f\t%.6f\n"
		       "Rove BLH: %10.6f\t%10.6f\t%7.4f\t\tBase_Line ENU: %.6f\t%.6f\t%.6f\n"
		       "m_x: % 1.6f\tm_y: % 1.6f\t\m_z: % 1.6f\tm_e: % 1.6f\tm_n: % 1.6f\t\m_u: % 1.6f\n"
		       "Ratio: %.4f\tQ_s: %1.6f\tRMS: %.4f\tRMS_L: %.4f\tRMS_P: %.4f\tRDOP: %.6f\tPDOP: %.6f\tSigma: %.6f\n"
			"GPS_REF: %2d\tBDS_REF: %2d\tGPS: % 2d\tBDS: % 2d\t % s\n",
			Rove_appro_pos->X, Rove_appro_pos->Y, Rove_appro_pos->Z,
			base_line.X, base_line.Y, base_line.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z,
			sqrt(m_xx(0, 0)), sqrt(m_xx(1, 1)), sqrt(m_xx(2, 2)),
			sqrt(m_nn(0, 0)), sqrt(m_nn(1, 1)), sqrt(m_nn(2, 2)),
			Ratio_LS, Q_s, RMS, RMS_L, RMS_P, RDOP, PDOP, sigma,
			Ref_PRN_GPS_NEW, Ref_PRN_BDS_NEW,
			LS_GPS_num, LS_BDS_num,
			LS_SATES->c_str());
		break;
	case OBS_DATA_Loss:
		printf("OBS_DATA_Loss\n");
		break;
	case Epoch_Loss:
		printf("Epoch_Loss\n");
		break;
	case Set_UP_B_fail:
		printf("Set_UP_B_fail\n");
		break;
	default:
		break;
	}
	printf("------------------------------------End of Epoch------------------------------------\n");
	fprintf(fpr, "------------------------------------End of Epoch------------------------------------\n");
}

void RTK_DATA::KF_Output(FILE* fpr, Configure cfg)
{
	XYZ rove_pos = *KF_pos;
	XYZ base_line = rove_pos - *Base_appro_pos;
	BLH blh = XYZ2BLH(rove_pos, WGS84_e2, WGS84_a);
	XYZ enu = XYZ2ENU(*Base_appro_pos, rove_pos, SYS_GPS);
	MatrixXd Qxx = Q_baseLine;
	double Base_Len = Len(&base_line);
	MatrixXd R = get_XYZ(base_line / Base_Len).transpose();
	MatrixXd K = get_Rot(degree2rad(blh.Lat), degree2rad(blh.Lon));
	double Q_s = (R * Qxx * R.transpose())(0, 0);
	MatrixXd Qnn = K * Qxx * K.transpose();
	MatrixXd V = Cal_V(this, cfg, KF_pos);
	MatrixXd N;
	if (KF_IS_FIXED)N = Fix_N;
	else N = KF->x_hat_.bottomRows(KF->x_hat_.rows() - 3);
	V -= KF->H_.rightCols(KF->H_.cols() - 3) * N;
	int gps = cfg.GPS_Cfg.used;
	int bds = cfg.BDS_Cfg.used;
	MatrixXd V_L = MatrixXd::Zero((gps * KF_GPS_num + bds * KF_BDS_num) * cfg.phase_num, 1);
	MatrixXd V_P = MatrixXd::Zero((gps * KF_GPS_num + bds * KF_BDS_num) * cfg.phase_num, 1);
	if (gps)
	{
		//GPS系统
		V_L.topRows(KF_GPS_num * cfg.phase_num) = V.topRows(KF_GPS_num * cfg.phase_num);
		V_P.topRows(KF_GPS_num * cfg.phase_num) = V.block(KF_GPS_num * cfg.phase_num, 0, KF_GPS_num * cfg.phase_num, 1);
	}
	if (bds)
	{
		//BDS系统
		V_L.bottomRows(KF_BDS_num * cfg.phase_num) = V.block(2 * gps * KF_GPS_num * cfg.phase_num, 0, KF_BDS_num * cfg.phase_num, 1);
		V_P.bottomRows(KF_BDS_num * cfg.phase_num) = V.block(2 * gps * KF_GPS_num * cfg.phase_num + KF_BDS_num * cfg.phase_num, 0, KF_BDS_num * cfg.phase_num, 1);
	}
	//std :: cout << V << endl << endl << V_L << endl << endl << V_P << endl;
	double sigma= sqrt((V.transpose() * KF->R_.inverse() * V)(0, 0) / (KF->H_.rows()-KF->H_.cols()));
	MatrixXd m_xx = Qxx * sigma;
	MatrixXd m_nn = Qnn * sigma;	
	double RMS = sqrt((V.transpose() * V)(0, 0) / V.rows());
	double RMS_L = sqrt((V_L.transpose() * V_L)(0, 0) / V_L.rows());
	double RMS_P = sqrt((V_P.transpose() * V_P)(0, 0) / V_P.rows());
	double RDOP = sqrt(Qxx.trace());
	double PDOP = Cal_PDOP(Q_baseLine);
	string is_fixed;
	if (KF_IS_FIXED)
		is_fixed = "Fix";
	else
		is_fixed = "Float";
	printf("Mode: RTK\tMethod: KF\tGPSTIME: %d\t%.3f\t  %s\n", OBSTIME->Week, OBSTIME->SecOfWeek, is_fixed.c_str());
	fprintf(fpr, "Mode: RTK\tMethod: KF\tGPSTIME: %d\t%.3f\t  %s\n", OBSTIME->Week, OBSTIME->SecOfWeek,
	        is_fixed.c_str());
	switch (KF_result)
	{
	case UN_Solve:
		printf("UN_Solve\n");
		break;
	case Success_Solve:
		printf("Rove XYZ: %.4f\t%.4f\t%.4f\tBase_Line XYZ: %.6f\t%.6f\t%.6f\n"
		       "Rove BLH: %10.6f\t%10.6f\t%7.4f\t\tBase_Line ENU: %.6f\t%.6f\t%.6f\n"
		       "m_x: % 1.6f\tm_y: % 1.6f\t\m_z: % 1.6f\tm_e: % 1.6f\tm_n: % 1.6f\t\m_u: % 1.6f\n"
		       "Ratio: %.4f\tQ_s: %1.6f\tRMS: %.4f\tRMS_L: %.4f\tRMS_P: %.4f\tRDOP: %.6f\tPDOP: %.6f\tSigma: %.6f\n"
		       "GPS_REF: %2d\tBDS_REF: %2d\tGPS: % 2d\tBDS: % 2d\t % s\n",
		       KF->x_hat_(0, 0), KF->x_hat_(1, 0), KF->x_hat_(2, 0),
		       base_line.X, base_line.Y, base_line.Z,
		       blh.Lat, blh.Lon, blh.Height,
		       enu.X, enu.Y, enu.Z,
		       sqrt(m_xx(0, 0)), sqrt(m_xx(1, 1)), sqrt(m_xx(2, 2)),
		       sqrt(m_nn(0, 0)), sqrt(m_nn(1, 1)), sqrt(m_nn(2, 2)),
		       Ratio_KF, Q_s, RMS, RMS_L, RMS_P, RDOP, PDOP, sigma, 
		       Ref_PRN_GPS_NEW, Ref_PRN_BDS_NEW,
		       KF_GPS_num, KF_BDS_num,
		       KF_SATES->c_str());

		fprintf(fpr,"Rove XYZ: %.4f\t%.4f\t%.4f\tBase_Line XYZ: %.6f\t%.6f\t%.6f\n"
		       "Rove BLH: %10.6f\t%10.6f\t%7.4f\t\tBase_Line ENU: %.6f\t%.6f\t%.6f\n"
		       "m_x: % 1.6f\tm_y: % 1.6f\t\m_z: % 1.6f\tm_e: % 1.6f\tm_n: % 1.6f\t\m_u: % 1.6f\n"
		       "Ratio: %.4f\tQ_s: %1.6f\tRMS: %.4f\tRMS_L: %.4f\tRMS_P: %.4f\tRDOP: %.6f\tPDOP: %.6f\tSigma: %.6f\n"
			"GPS_REF: %2d\tBDS_REF: %2d\tGPS: % 2d\tBDS: % 2d\t % s\n",
			KF->x_hat_(0, 0), KF->x_hat_(1, 0), KF->x_hat_(2, 0),
			base_line.X, base_line.Y, base_line.Z,
			blh.Lat, blh.Lon, blh.Height,
			enu.X, enu.Y, enu.Z,
			sqrt(m_xx(0, 0)), sqrt(m_xx(1, 1)), sqrt(m_xx(2, 2)),
			sqrt(m_nn(0, 0)), sqrt(m_nn(1, 1)), sqrt(m_nn(2, 2)),
			Ratio_KF, Q_s, RMS, RMS_L, RMS_P, RDOP, PDOP, sigma,
			Ref_PRN_GPS_NEW, Ref_PRN_BDS_NEW,
			KF_GPS_num, KF_BDS_num,
			KF_SATES->c_str());
		break;
	case OBS_DATA_Loss:
		printf("OBS_DATA_Loss\n");
		break;
	case Epoch_Loss:
		printf("Epoch_Loss\n");
		break;
	case Set_UP_B_fail:
		printf("Set_UP_B_fail\n");
		break;
	default:
		break;
	}
	printf("------------------------------------End of Epoch------------------------------------\n");
	fprintf(fpr, "------------------------------------End of Epoch------------------------------------\n");
}
