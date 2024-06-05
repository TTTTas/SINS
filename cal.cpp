#include "cal.h"
#include "GNSS_data.h"
#include "LAMBDA.h"
#include "comm.h"
#include <map>

/* 函数名: SQR
	输入参数:
	 - x: 要计算平方的值
	返回值: 输入值的平方
*/
double SQR(double x)
{
	return x * x;
}

// 函数名: Len
// 输入参数:
//   - pos: 包含三维坐标(X, Y, Z)的结构体
// 返回值: 坐标点的长度（模）
double Len(XYZ* pos)
{
	return sqrt(SQR(pos->X) + SQR(pos->Y) + SQR(pos->Z));
}

double Len(XYZ* pos1, XYZ* pos2)
{
	return sqrt(SQR(pos1->X - pos2->X) + SQR(pos1->Y - pos2->Y) + SQR(pos1->Z - pos2->Z));
}

// 函数名: Cal_PDOP
// 输入参数:
//   - Qxx: 协方差矩阵
// 返回值: PDOP值
double Cal_PDOP(MatrixXd Qxx)
{
	if (Qxx.rows() < 3 || Qxx.cols() < 3)
		return 0;
	return sqrt(Qxx(0, 0) + Qxx(1, 1) + Qxx(2, 2));
}

// 函数名: Cal_LEAST_SQR
// 输入参数:
//   - B: 设计矩阵
//   - l: 观测残差
//   - P: 观测权重矩阵
//   - Qxx: 输出参数，协方差矩阵
//   - x: 输出参数，最小二乘解
//   - thegma: 输出参数，单位权中误差
//   - DOP: 输出参数，PDOP值
// 返回值: 操作结果，成功返回1，否则返回0
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd* Qxx, MatrixXd& x, double* thegma, double* DOP)
{
	// 检查输入矩阵的维度是否合法
	if (B.rows() != P.rows() || B.rows() != l.rows())
		return 0;
	if (B.rows() < B.cols())
		return 0;

	// 计算协方差矩阵
	*Qxx = (B.transpose() * P * B).inverse();

	// 计算最小二乘解
	x = *Qxx * B.transpose() * P * l;

	// 计算观测残差
	MatrixXd v = B * x - l;

	// 计算单位权中误差
	int m = B.rows() - B.cols();
	*thegma = sqrt(((v.transpose() * P * v) / m)(0, 0));

	// 计算PDOP值
	*DOP = Cal_PDOP(*Qxx);

	return 1; // 成功
}

// 函数名: decode_SYN
// 输入参数:
//   - sys: 系统编号
//   - signal: 信号编号
// 返回值: 解码后的信号编号，未知信号返回UNKOWN
unsigned int decode_SYN(int sys, int signal)
{
	switch (sys)
	{
	case SYS_GPS:
		switch (signal)
		{
		case 0:
			return CODE_L1C; /* L1C/A */
		case 5:
			return CODE_L2P; /* L2P    (OEM7) */
		case 9:
			return CODE_L2W; /* L2P(Y),semi-codeless */
		case 14:
			return CODE_L5Q; /* L5Q    (OEM6) */
		case 16:
			return CODE_L1L; /* L1C(P) (OEM7) */
		case 17:
			return CODE_L2S; /* L2C(M) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	case SYS_BDS:
		switch (signal)
		{
		case 0:
			return CODE_L2I; /* B1I with D1 (OEM6) */
		case 1:
			return CODE_L7I; /* B2I with D1 (OEM6) */
		case 2:
			return CODE_L6I; /* B3I with D1 (OEM7) */
		case 4:
			return CODE_L2I; /* B1I with D2 (OEM6) */
		case 5:
			return CODE_L7I; /* B2I with D2 (OEM6) */
		case 6:
			return CODE_L6I; /* B3I with D2 (OEM7) */
		case 7:
			return CODE_L1P; /* B1C(P) (OEM7) */
		case 9:
			return CODE_L5P; /* B2a(P) (OEM7) */
		default:
			return UNKOWN;
			break;
		}
	default:
		return UNKOWN;
		break;
	}
}

// 函数名: CODE2FREQ
// 输入参数:
//   - code: 信号编号
// 返回值: 对应信号编号的频率值
double CODE2FREQ(int code)
{
	switch (code)
	{
	case 0:
		return 0;
	case 1:
		return Freq_L1;
	case 2:
	case 3:
		return Freq_L2;
	case 4:
		return Freq_L5;
	case 5:
	case 6:
		return Freq_L1;
	case 7:
		return Freq_B1;
	case 8:
		return Freq_B2;
	case 9:
		return Freq_B3;
	case 10:
		return Freq_B1_C;
	case 11:
		return Freq_B2_a;
	default:
		return 0;
		break;
	}
}

// 函数名: CORRECT_CLK
// 输入参数:
//   - t: 当前时刻
//   - eph: 指向EPHEMERIS结构的指针
// 返回值: 改正后的钟差
double CORRECT_CLK(double t, EPHEMERIS* eph)
{
	double correct_clk = t;
	for (int i = 0; i < 10; i++)
	{
		correct_clk = t - (eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (
			correct_clk - eph->toc));
	}
	return eph->a_f0 + eph->a_f1 * (correct_clk - eph->toc) + eph->a_f2 * (correct_clk - eph->toc) * (correct_clk - eph
		->toc);
}

// 函数名: TGD
// 输入参数:
//   - e: 指向EPHEMERIS结构的指针
//   - f: 频率值
//   - sys: 卫星系统
// 返回值: TGD值
double TGD(EPHEMERIS* e, double f, int sys)
{
	switch (sys)
	{
	case SYS_GPS:
		return SQR(f / Freq_L1) * e->T_GD1;
		break;
	case SYS_BDS:
		if (f == Freq_B1)
			return e->T_GD1;
		if (f == Freq_B2)
			return e->T_GD2;
		if (f == Freq_B3)
			return 0.0;
	default:
		return -1;
		break;
	}
}

// 函数名: SAT_POS_CAL
// 输入参数:
//   - t: 观测时刻
//   - eph: 卫星星历数据结构指针
//   - Sat_Pos: 卫星位置结果结构指针
//   - clk: 卫星钟差
//   - dt: 传播时间修正值
//   - SYS: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: 1 表示计算成功

unsigned int SAT_POS_CAL(double t, EPHEMERIS* eph, XYZ* Sat_Pos, double& clk, double dt, int SYS)
{
	double n, delt_t, M, E, E0, V, u_, u, r, i, dt0, F;
	switch (SYS)
	{
	case SYS_GPS:
		n = sqrt(WGS84_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(WGS84_GM) / (velocity_c * velocity_c);
		break;
	case SYS_BDS:
		n = sqrt(CGCS2000_GM) / (eph->sqrt_A * eph->sqrt_A * eph->sqrt_A) + eph->delt_n;
		dt0 = delt_t = t - eph->toe_tow;
		F = -2 * sqrt(CGCS2000_GM) / (velocity_c * velocity_c);
		break;
	default:
		break;
	}

	double dtr = 0;
	while (abs(delt_t) > 302400)
	{
		if (delt_t > 302400)
		{
			delt_t -= 604800;
		}
		else if (delt_t < -302400)
		{
			delt_t += 604800;
		}
	}

	for (int i = 0; i < 10; i++)
	{
		M = eph->M0 + n * delt_t;
		E = M;
		E0 = M;
		int count = 0;
		do
		{
			E0 = E;
			E = M + eph->e * sin(E0);
			count++;
		}
		while (abs(E - E0) > 0.000000000001 && count < 10);
		dtr = F * eph->e * eph->sqrt_A * sin(E);
		delt_t = dt0 - dtr;
	}
	clk += dtr;
	dt += clk;
	V = atan2((sqrt(1 - eph->e * eph->e) * sin(E)), (cos(E) - eph->e));
	u_ = eph->omiga + V;
	u = u_ + eph->Cuc * cos(2 * u_) + eph->Cus * sin(2 * u_);
	r = eph->sqrt_A * eph->sqrt_A * (1 - eph->e * cos(E)) + eph->Crc * cos(2 * u_) + eph->Crs * sin(2 * u_);
	i = eph->i0 + eph->Cic * cos(2 * u_) + eph->Cis * sin(2 * u_) + eph->dot_i * delt_t;
	double x = r * cos(u);
	double y = r * sin(u);
	double z = 0;
	double L;
	double x0, y0, z0;
	switch (SYS)
	{
	case SYS_GPS:
		L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
		Sat_Pos->X = (x * cos(L) - y * cos(i) * sin(L)) + omiga_earth * dt * (x * sin(L) + y * cos(i) * cos(L));
		Sat_Pos->Y = x * sin(L) + y * cos(i) * cos(L) - omiga_earth * dt * (x * cos(L) - y * cos(i) * sin(L));
		Sat_Pos->Z = y * sin(i);
		break;
	case SYS_BDS:
		if (fabs(eph->i0 - 0.0873) < 0.1 && fabs(eph->sqrt_A - 6493) < 1)
		// 通过轨道倾角和轨道根数判断是否为GEO卫星  i: 5/deg sqrt_A: 6493/sqrt_meter
		{
			L = eph->Omiga0 + eph->dot_Omiga * delt_t - omiga_earth * eph->toe_tow;
			x0 = x * cos(L) - y * cos(i) * sin(L);
			y0 = x * sin(L) + y * cos(i) * cos(L);
			z0 = y * sin(i);
			MatrixXd P_GK(3, 1);
			MatrixXd R_Z(3, 3);
			MatrixXd R_X(3, 3);
			MatrixXd P(3, 1);
			P_GK << x0,
				y0,
				z0;
			R_X << 1, 0, 0,
				0, cos(-5 * Pi / 180), sin(-5 * Pi / 180),
				0, -sin(-5 * Pi / 180), cos(-5 * Pi / 180);
			R_Z << cos(omiga_earth * delt_t), sin(omiga_earth * delt_t), 0,
				-sin(omiga_earth * delt_t), cos(omiga_earth * delt_t), 0,
				0, 0, 1;
			P = R_Z * R_X * P_GK;
			Sat_Pos->X = cos(omiga_earth * dt) * P(0, 0) + sin(omiga_earth * dt) * P(1, 0);
			Sat_Pos->Y = cos(omiga_earth * dt) * P(1, 0) - sin(omiga_earth * dt) * P(0, 0);
			Sat_Pos->Z = P(2, 0);
		}
		else
		{
			L = eph->Omiga0 + (eph->dot_Omiga - omiga_earth) * t - eph->dot_Omiga * eph->toe_tow;
			Sat_Pos->X = cos(omiga_earth * dt) * (x * cos(L) - y * cos(i) * sin(L)) + sin(omiga_earth * dt) * (x *
				sin(L) + y * cos(i) * cos(L));
			Sat_Pos->Y = cos(omiga_earth * dt) * (x * sin(L) + y * cos(i) * cos(L)) - sin(omiga_earth * dt) * (x *
				cos(L) - y * cos(i) * sin(L));
			Sat_Pos->Z = y * sin(i);
		}
		break;
	default:
		break;
	}

	return 1;
}

// 函数名: Ele_Angle
// 输入参数:
//   - SatPos: 卫星直角坐标系下的位置
//   - RcvPos: 接收机直角坐标系下的位置
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: 卫星高度角

double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys)
{
	// 将卫星坐标转换到ENU坐标系
	XYZ Satenu = XYZ2ENU(RcvPos, SatPos, sys);

	// 计算卫星高度角
	return asin(Satenu.Z / Len(&Satenu));
}

// 函数名: Hopefield
// 输入参数:
//   - E: 卫星高度角
//   - H: 接收机海拔高度
// 返回值: Hopefield对流层改正值

double Hopefield(double E, double H)
{
	// 计算标准大气压力和温度
	double Ts = T0 - 0.0065 * (H - H0);
	double hd = 40136 + 148.72 * (Ts - 273.16);

	// 计算标准大气压和相对湿度
	double Ps = P0 * pow(1 - 0.0000226 * (H - H0), 5.225);
	double RH = RH0 * exp(-0.0006396 * (H - H0));

	// 计算饱和水汽压力
	double es = RH * exp(-37.2465 + 0.213166 * Ts - 0.000256908 * Ts * Ts);

	// 计算方位角对应的方向余弦
	double md = sin(deg2rad(sqrt(E * E + 6.25)));
	double mw = sin(deg2rad(sqrt(E * E + 2.25)));

	// 计算水汽分压和水汽刻度高
	double hw = 11000;
	double ZHD = 155.2 * 1e-7 * Ps * (hd - H) / Ts;
	double ZWD = 155.2 * 1e-7 * 4810 * es * (hw - H) / (Ts * Ts);

	// 返回Hopefield对流层改正值
	return ZHD / md + ZWD / mw;
}

// 函数名: Hopefield
// 输入参数:
//   - SatPos: 卫星直角坐标系下的位置
//   - RcvPos: 接收机直角坐标系下的位置
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: Hopefield对流层改正值

double Hopefield(XYZ SatPos, XYZ RcvPos, int sys)
{
	// 若卫星坐标为零向量，返回0
	if (SatPos.X == 0 && SatPos.Y == 0 && SatPos.Z == 0)
	{
		return 0;
	}

	// 计算卫星高度角
	double E = rad2deg(Ele_Angle(SatPos, RcvPos, sys));

	// 获取接收机海拔高度
	double H = 0;
	switch (sys)
	{
	case SYS_GPS:
		H = XYZ2BLH(RcvPos, WGS84_e2, WGS84_a).Height;
		break;
	case SYS_BDS:
		H = XYZ2BLH(RcvPos, CGCS2000_e2, CGCS2000_a).Height;
		break;
	default:
		break;
	}

	// 判断海拔高度是否在有效范围内，计算Hopefield对流层改正值
	if (H < 20e3 && H > -100)
	{
		return Hopefield(E, H);
	}
	else
	{
		return 0;
	}
}


double Klobuchar(XYZ RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys)
{
	if (!(alpha[0] * alpha[1] * alpha[2] * alpha[3] * beta[0] * beta[1] * beta[2] * beta[3]))
		return -1;
	BLH RcvBLH;
	double T_g = 0;
	double EA = 0;
	double B_IPP = 0;
	double L_IPP = 0;
	double B_m = 0;
	double t = 0;
	double A_I = 0;
	double P_I = 0;
	double Phase_I = 0;
	double F = 0;
	switch (sys)
	{
	case SYS_GPS:
		RcvBLH = XYZ2BLH(RcvPos, WGS84_e2, WGS84_a);
		EA = 0.0137 / (E + 0.11) - 0.022;
		B_IPP = RcvBLH.Lat + EA * cos(A);
		if (B_IPP < -0.416)
			B_IPP = -0.416;
		if (B_IPP > 0.416)
			B_IPP = 0.416;
		L_IPP = RcvBLH.Lon + EA * sin(A) / cos(B_IPP);
		B_m = B_IPP + 0.064 * cos(L_IPP - 1.617);
		t = 43200 * L_IPP + UT;
		while (t > 86400 || t < 0)
		{
			if (t > 86400)
				t -= 86400;
			if (t < 0)
				t += 86400;
		}
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_m, i);
			P_I += beta[i] * pow(B_m, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		Phase_I = 2 * Pi * (t - 50400) / P_I;
		F = 1 + 16 * pow((0.53 - E), 3);
		if (abs(Phase_I) < 1.57)
		{
			T_g = F * (5e-9 + A_I * (1 - pow(Phase_I, 2) / 2 + pow(Phase_I, 4) / 24));
		}
		else
		{
			T_g = F * 5e-9;
		}
		return pow(Freq_L1 / code, 2) * T_g;
		break;
	case SYS_BDS:
		RcvBLH = XYZ2BLH(RcvPos, CGCS2000_e2, CGCS2000_a);
		EA = Pi / 2 - E - asin(cos(E) * 6378 / (6378 + 375));
		B_IPP = asin(sin(RcvBLH.Lat) * cos(EA) + cos(RcvBLH.Lat) * sin(EA) * cos(A));
		L_IPP = RcvBLH.Lon + asin(sin(EA) * sin(A) / cos(B_IPP));
		t = UT + L_IPP * 43200 / Pi;
		A_I = 0;
		P_I = 0;
		for (int i = 0; i < 4; i++)
		{
			A_I += alpha[i] * pow(B_IPP / Pi, i);
			P_I += beta[i] * pow(B_IPP / Pi, i);
		}
		if (A_I < 0)
			A_I = 0;
		if (P_I < 72000)
			P_I = 72000;
		if (P_I > 172800)
			P_I = 172800;
		if (abs(t - 50400) < P_I / 4)
		{
			T_g = 5e-9 + A_I * cos(2 * Pi * (t - 50400) / P_I);
		}
		else
		{
			5e-9;
		}
		return T_g / sqrt(1 - pow(cos(E) * 6378 / (6378 + 375), 2));
		break;
	default:
		return 0;
		break;
	}
}

// 函数名: setup_LS
// 输入参数:
//   - data: 数据集结构体指针
//   - cfg: 配置信息结构体
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: LS定位矩阵行数

unsigned int setup_LS(DATA_SET* data, Configure cfg, int sys)
{
	// 声明变量
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;

	// 根据卫星系统类型选择相应的卫星数据和星历数据
	switch (sys)
	{
	case SYS_GPS:
		Sates = data->range->GPS_SATE;
		eph = data->GPS_eph;
		f = cfg.GPS_Cfg.f1;
		break;
	case SYS_BDS:
		Sates = data->range->BDS_SATE;
		eph = data->BDS_eph;
		f = cfg.BDS_Cfg.f1;
		break;
	default:
		return 0;
		break;
	}
	if (Sates.size() < cfg.SYS_num + 3)
	{
		data->LS_result = OBS_DATA_Loss;
		return 0;
	}
	// 获取接收机位置和接收机钟差
	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	double dt_Rcv = data->temp_ref(3, 0);

	// 初始化变量
	MatrixXd B, l_Pos, P_Pos;
	MatrixXd l_Vel, P_Vel;
	B = MatrixXd::Zero(0, 4);
	l_Vel = l_Pos = MatrixXd::Zero(0, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);

	double p_pos = 0;

	// 循环遍历卫星
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;

		// 检查星历数据是否匹配
		if (eph[prn - 1]->PRN != prn)
			continue;

		// 查找卫星对应的频率索引
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;

		// 获取测量值
		double measure = get_measure(Sates[i], cfg, eph[prn - 1], p_pos);
		if (!measure)
			continue;

		if (Sates[i]->pos.X == 0)
			continue;
		double elev = Ele_Angle(Sates[i]->pos, RcvPos, Sates[i]->SYS);
		Sates[i]->ELEV = elev;
		// 检查是否为首次定位
		if (!data->LS_first)
		{
			if (elev < deg2rad(cfg.Ele_Mask))
				continue;
		}

		// 计算方位角对应的波长和长度
		double lamda = (1e-6 * velocity_c / f);
		double len = sqrt(
			SQR(RcvPos.X - Sates[i]->pos.X) + SQR(RcvPos.Y - Sates[i]->pos.Y) + SQR(RcvPos.Z - Sates[i]->pos.Z));

		// 计算预测值
		double w_pos = measure - (len + dt_Rcv - velocity_c * Sates[i]->clk + Hopefield(
			Sates[i]->pos, RcvPos, Sates[i]->SYS));
		double v0 = ((Sates[i]->pos.X - RcvPos.X) * Sates[i]->vel.X + (Sates[i]->pos.Y - RcvPos.Y) * Sates[i]->vel.Y + (
			Sates[i]->pos.Z -
			RcvPos.Z) * Sates[i]->vel.Z) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - Sates[i]->clk_vel);

		// 检查是否为首次定位
		//if (!data->LS_first)
		//{
		//	if (abs(w_pos) > cfg.w_thresh)
		//		continue;
		//}

		// 构造最小二乘问题的设计矩阵B和观测矩阵l
		double l = (RcvPos.X - Sates[i]->pos.X) / len;
		double m = (RcvPos.Y - Sates[i]->pos.Y) / len;
		double n = (RcvPos.Z - Sates[i]->pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		ROWS++;

		// 更新B、l、P矩阵
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		if (ROWS == 1)
		{
			P_Pos = MatrixXd::Identity(1, 1);
			P_Pos(0, 0) = p_pos;
		}
		else
		{
			MatrixXd temp_P = MatrixXd::Identity(ROWS, ROWS);
			temp_P.block(0, 0, P_Pos.rows(), P_Pos.cols()) = P_Pos;
			temp_P(ROWS - 1, ROWS - 1) = p_pos;
			P_Pos = temp_P;
		}

		// 根据卫星系统类型记录卫星PRN
		switch (sys)
		{
		case SYS_GPS:
			*(data->LS_SATES) += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*(data->LS_SATES) += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
	}

	// 初始化权阵矩阵
	P_Pos = MatrixXd::Identity(ROWS, ROWS);
	P_Vel = MatrixXd::Identity(ROWS, ROWS);


	// 若有卫星数据，则更新LS_Pos和LS_Vel矩阵
	if (ROWS >= cfg.SYS_num + 3)
	{
		//std::cout << B << endl;
		//std::cout << l_Pos << endl;
		data->LS_Pos->set_B_Pos(B);
		data->LS_Pos->set_l(l_Pos);
		data->LS_Pos->set_P(P_Pos);
		data->LS_Vel->set_B_Vel(B);
		data->LS_Vel->set_l(l_Vel);
		data->LS_Vel->set_P(P_Vel);
	}
	else
	{
		data->LS_result = OBS_DATA_Loss;
		return 0;
	}

	// 根据卫星系统类型记录卫星数量
	switch (sys)
	{
	case SYS_GPS:
		data->LS_GPS_num = ROWS;
		break;
	case SYS_BDS:
		data->LS_BDS_num = ROWS;
		break;
	default:
		break;
	}

	// 返回LS定位矩阵行数
	return ROWS;
}


// 函数名: setup_KF
// 输入参数:
//   - data: 数据集结构体指针
//   - cfg: 配置信息结构体
//   - sys: 卫星系统类型 (SYS_GPS 或 SYS_BDS)
// 返回值: KF定位矩阵行数

unsigned int setup_KF(DATA_SET* data, Configure cfg, int sys)
{
	// 声明变量
	int ROWS = 0;
	vector<Satellate*> Sates;
	EPHEMERIS** eph;
	double f = 0;

	// 根据卫星系统类型选择相应的卫星数据和星历数据
	switch (sys)
	{
	case SYS_GPS:
		Sates = data->range->GPS_SATE;
		eph = data->GPS_eph;
		f = cfg.GPS_Cfg.f1;
		break;
	case SYS_BDS:
		Sates = data->range->BDS_SATE;
		eph = data->BDS_eph;
		f = cfg.BDS_Cfg.f1;
		break;
	default:
		return 0;
		break;
	}

	// 获取接收机位置和接收机钟差
	XYZ RcvPos = get_XYZ(data->temp_ref.block(0, 0, 3, 1));
	MatrixXd RcvVel(4, 4);
	RcvVel.block(0, 0, 3, 1) = data->KF->getState_minus().block(3, 0, 3, 1);
	RcvVel(3, 0) = data->KF->getState_minus()(cfg.SYS_num + 6, 0);
	double dt_Rcv = data->temp_ref(3, 0);
	MatrixXd B, l_Pos, l_Vel;
	B = MatrixXd::Zero(0, 4);
	l_Vel = l_Pos = MatrixXd::Zero(0, 1);
	MatrixXd B_pos_new = MatrixXd::Zero(1, 4);
	MatrixXd l_pos_new = MatrixXd::Zero(1, 1);
	MatrixXd l_vel_new = MatrixXd::Zero(1, 1);
	MatrixXd P_Pos;
	double p_pos = 0;
	// 循环遍历卫星
	for (int i = 0; i < Sates.size(); i++)
	{
		int prn = Sates[i]->PRN;

		// 检查星历数据是否匹配
		if (eph[prn - 1]->PRN != prn)
			continue;

		// 查找卫星对应的频率索引
		int Index = 0;
		while (CODE2FREQ(decode_SYN(Sates[i]->SYS, Sates[i]->SYG_TYPE[Index])) != f && Index < MAXNUM)
			Index++;
		if (Index == MAXNUM)
			continue;

		// 获取测量值
		double measure = get_measure(Sates[i], cfg, eph[prn - 1], p_pos);
		if (!measure)
			continue;

		if (Sates[i]->pos.X == 0)
			continue;
		double elev = Ele_Angle(Sates[i]->pos, RcvPos, Sates[i]->SYS);
		Sates[i]->ELEV = elev;
		// 检查仰角是否小于阈值
		if (elev < deg2rad(cfg.Ele_Mask))
			continue;

		// 计算方位角对应的波长和长度
		double lamda = (1e-6 * velocity_c / f);
		double len = Len(&RcvPos, &Sates[i]->pos);
		// 构造卡尔曼滤波问题的设计矩阵B和观测矩阵l
		double l = (RcvPos.X - Sates[i]->pos.X) / len;
		double m = (RcvPos.Y - Sates[i]->pos.Y) / len;
		double n = (RcvPos.Z - Sates[i]->pos.Z) / len;
		B_pos_new(0, 0) = l;
		B_pos_new(0, 1) = m;
		B_pos_new(0, 2) = n;
		B_pos_new(0, 3) = 1;
		// 计算预测值
		double w_pos = measure - (len + dt_Rcv - velocity_c * Sates[i]->clk + Hopefield(
			Sates[i]->pos, RcvPos, Sates[i]->SYS));
		double v0 = ((Sates[i]->pos.X - RcvPos.X) * Sates[i]->vel.X + (Sates[i]->pos.Y - RcvPos.Y) * Sates[i]->vel.Y + (
			Sates[i]->pos.Z -
			RcvPos.Z) * Sates[i]->vel.Z) / len;
		double w_vel = -lamda * Sates[i]->DOPPLER[Index] - (v0 - Sates[i]->clk_vel) - (B_pos_new * RcvVel)(0, 0);

		// 检查预测值是否大于阈值
		if (abs(w_pos) > cfg.w_thresh)
			continue;

		// 检查预测速度是否大于阈值
		if (abs(w_vel) > cfg.w_thresh)
			continue;

		l_pos_new(0, 0) = w_pos;
		l_vel_new(0, 0) = w_vel;
		ROWS++;

		// 更新B、l矩阵
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_pos_new;
		l_Pos.conservativeResize(l_Pos.rows() + 1, l_Pos.cols());
		l_Pos.bottomRows(1) = l_pos_new;
		l_Vel.conservativeResize(l_Vel.rows() + 1, l_Vel.cols());
		l_Vel.bottomRows(1) = l_vel_new;
		if (ROWS == 1)
		{
			P_Pos = MatrixXd::Identity(1, 1);
			P_Pos(0, 0) = p_pos;
		}
		else
		{
			MatrixXd temp_P = MatrixXd::Identity(ROWS, ROWS);
			temp_P.block(0, 0, P_Pos.rows(), P_Pos.cols()) = P_Pos;
			temp_P(ROWS - 1, ROWS - 1) = p_pos;
			P_Pos = temp_P;
		}

		// 根据卫星系统类型记录卫星PRN
		switch (sys)
		{
		case SYS_GPS:
			*(data->KF_SATES) += "G" + to_string(Sates[i]->PRN);
			break;
		case SYS_BDS:
			*(data->KF_SATES) += "C" + to_string(Sates[i]->PRN);
			break;
		default:
			break;
		}
	}

	// 如果卫星数为0，返回0
	if (ROWS == 0)
		return 0;

	MatrixXd R = MatrixXd::Identity(ROWS + ROWS, ROWS + ROWS);
	R.block(0, 0, ROWS, ROWS) = P_Pos;
	R.block(ROWS, ROWS, ROWS, ROWS) = MatrixXd::Identity(ROWS, ROWS) * 0.18;

	// 更新Kalman Filter对象的H、R、Z矩阵
	data->KF->set_H(getH(B), MODE_SPP);
	data->KF->set_R(getR(ROWS), MODE_SPP);
	data->KF->set_Z(getz(l_Pos, l_Vel));

	// 根据卫星系统类型记录卫星数量
	switch (sys)
	{
	case SYS_GPS:
		data->KF_GPS_num = ROWS;
		break;
	case SYS_BDS:
		data->KF_BDS_num = ROWS;
		break;
	default:
		break;
	}

	// 返回KF定位矩阵行数
	return ROWS;
}


// 函数名: get_measure
// 输入参数:
//   - Sate: 卫星信息结构体指针
//   - cfg: 配置信息结构体
//   - eph: 星历数据结构体指针
//   - p: 权重
// 返回值: 测量值

double get_measure(Satellate* Sate, Configure cfg, EPHEMERIS* eph, double& p)
{
	// 声明变量
	double measure = 0;
	double f1 = 0;
	double f2 = 0;

	// 根据卫星系统类型选择相应的频率
	switch (Sate->SYS)
	{
	case SYS_GPS:
		f1 = cfg.GPS_Cfg.f1;
		f2 = cfg.GPS_Cfg.f2;
		break;
	case SYS_BDS:
		f1 = cfg.BDS_Cfg.f1;
		f2 = cfg.BDS_Cfg.f2;
		break;
	default:
		break;
	}

	// 如果卫星标记为异常值，返回测量值为0
	if (Sate->Outlier)
		return 0;

	int prn = Sate->PRN;
	int Index1 = 0;
	int Index2 = 0;

	// 查找卫星对应的频率索引
	while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
		Index1++;
	while (CODE2FREQ(decode_SYN(Sate->SYS, Sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
		Index2++;

	// 如果未找到对应频率索引，返回测量值为0
	if (Index1 == MAXNUM)
		return 0;

	if (Sate->SNR[Index1] < cfg.SNR_thresh)
		return 0;

	if (cfg.phase_num == 2)
		if (Sate->SNR[Index2] < cfg.SNR_thresh) return 0;

	// 如果指定的频率1没有相应的相位和伪距观测值，返回测量值为0
	if (!(Sate->LOCK_PSE[Index1] && Sate->LOCK_PHA[Index1]))
		return 0;

	// 根据配置的系统数目计算测量值
	if (cfg.phase_num == 1)
	{
		measure = Sate->PSERA[Index1];
		p = SQR(Sate->PSE_PREC[Index1]);
	}
	else if (cfg.phase_num == 2)
	{
		// 如果相位观测值数量小于2，返回测量值为0
		if (Sate->Phase_NUM < 2)
			return 0;

		// 如果未找到对应频率索引，返回测量值为0
		if (Index2 == MAXNUM)
			return 0;

		// 如果指定的频率2没有相应的相位和伪距观测值，返回测量值为0
		if (!(Sate->LOCK_PSE[Index2] && Sate->LOCK_PHA[Index2]))
			return 0;

		// 计算测量值（无电离层组合）
		measure = Sate->PSERA[Index1] * SQR(f1) / (SQR(f1) - SQR(f2)) - SQR(f2) * Sate->PSERA[Index2] / (SQR(f1) -
			SQR(f2));
		p = SQR(SQR(f1) * Sate->PSE_PREC[Index1] / (SQR(f1) - SQR(f2))) + SQR(
			SQR(f2) * Sate->PSE_PREC[Index2] / (SQR(f1) - SQR(f2)));

		// 对于BDS系统，加上TGD1对应的钟差改正
		if (Sate->SYS == SYS_BDS)
			measure += velocity_c * SQR(f2 / f1) * eph->T_GD1 / (1 - SQR(f2 / f1));
	}

	// 返回计算得到的测量值
	return measure;
}

double P_height1(double ELEV)
{
	return SQR(4 * 1e-3) + SQR(3 * 1e-3 / sin(ELEV));
}

double P_height2(double ELEV)
{
	return SQR(0.003) * (1 + 1.5 * SQR(cos(ELEV)));
}


// 函数名: LS_SPV
// 输入参数:
//   - data: 数据集，包含卡尔曼滤波器、最小二乘解等数据
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
unsigned int LS_SPV(DATA_SET* data, Configure cfg)
{
	int val = 0;

	int i = 0;
	// 迭代4次
	MatrixXd x_last = data->LS_Pos->X.topRows(3);
	for (i = 0; i < 10; i++)
	{
		// 清空卫星信息
		*data->LS_SATES = "";
		// 重置最小二乘解和速度
		data->LS_Pos->reset();
		data->LS_Vel->reset();
		// 设置最小二乘
		val = data->Set_LS(cfg);

		// 如果成功设置最小二乘，则执行最小二乘解算
		if (val)
		{
			data->LS_Pos->LS();
			data->Pos += data->LS_Pos->X.topRows(3);
			data->LS_Vel->LS();
			data->LS_result = Success_Solve;
			double dx = (MatrixXd::Ones(1, 3) * (data->LS_Pos->X.topRows(3) - x_last))(0, 0);
			if (abs(dx) < 1e-6)
				break;
			x_last = data->LS_Pos->X.topRows(3);
		}
	}
	std::cout << "LS迭代" << i << "次" << endl;
	return val;
}

// 函数名: Set_LS
// 输入参数:
//   - cfg: 配置信息，包括系统数和传感器配置
// 返回值: 最小二乘问题的行数
int DATA_SET::Set_LS(Configure cfg)
{
	// 将最小二乘解的前三行赋值给临时参考数据
	temp_ref.topRows(3) = Pos;

	// 初始化GPS和BDS的接收机钟差
	double dt_G = 0;
	double dt_C = 0;

	// 根据系统数量设置接收机钟差
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

	int val1 = 1, val2 = 1;
	// 设置GPS的接收机钟差
	if (cfg.GPS_Cfg.used)
	{
		val1 = 0;
		temp_ref(3, 0) = dt_G;
		val1 = setup_LS(this, cfg, SYS_GPS);
	}

	// 设置BDS的接收机钟差
	if (cfg.BDS_Cfg.used)
	{
		val2 = 0;
		temp_ref(3, 0) = dt_C;
		val2 = setup_LS(this, cfg, SYS_BDS);
	}

	// 返回最小二乘的行数
	return val1 * val2 ? val1 + val2 : 0;
}


// 函数名: KF_SPV
// 输入参数:
//   - data: 数据集，包含卡尔曼滤波器、最小二乘解等数据
//   - dt_e: 时间步长
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
unsigned int KF_SPV(DATA_SET* data, double dt_e, Configure cfg)
{
	int val = 0;

	// 如果是第一次滤波且最小二乘解成功
	if (data->KF_first && data->LS_result == Success_Solve)
	{
		// 从最小二乘解中提取状态向量
		MatrixXd state(8, 1);
		state.block(0, 0, 3, 1) = data->LS_Pos->X.block(0, 0, 3, 1);
		state.block(3, 0, 3, 1) = data->LS_Vel->X.block(0, 0, 3, 1);
		// 根据系统编号配置状态向量
		if (cfg.SYS_num == 1)
		{
			state(6, 0) = data->LS_Pos->X(3, 0);
			state(7, 0) = data->LS_Vel->X(3, 0);
		}
		else if (cfg.SYS_num == 2)
		{
			state.conservativeResize(9, 1);
			state(6, 0) = data->LS_Pos->X(3, 0);
			state(7, 0) = data->LS_Pos->X(4, 0);
			state(8, 0) = data->LS_Vel->X(3, 0);
		}

		// 设置卡尔曼滤波器状态
		data->KF->setState(state);
		data->KF_GPS_num = data->LS_GPS_num;
		data->KF_BDS_num = data->LS_BDS_num;
		*data->KF_SATES = *(data->LS_SATES);
		data->KF_result = Success_Solve;
		return 1;
	}

	// 重置卡尔曼滤波器并进行预测
	*data->KF_SATES = "";
	data->KF->reset();
	data->KF->set_A(getA(dt_e, cfg));
	data->KF->set_Q(getQ(dt_e, cfg));
	data->KF->predict();

	// 设置卡尔曼滤波器的状态
	val = data->Set_KF(cfg);

	// 如果设置成功，执行更新步骤
	if (val)
	{
		data->KF->update();
		data->KF_result = Success_Solve;
	}

	return val;
}

// 函数名: DATA_SET::Set_KF
// 输入参数:
//   - cfg: 配置信息，包括系统编号和传感器配置
// 返回值: 操作结果，成功返回0，否则返回非零值
int DATA_SET::Set_KF(Configure cfg)
{
	int val = 0;
	temp_ref.topRows(3) = KF->getState_minus().block(0, 0, 3, 1);
	double dt_G = 0;
	double dt_C = 0;

	// 根据系统编号配置GPS和BDS的接收机钟差
	if (cfg.SYS_num == 1)
	{
		if (cfg.GPS_Cfg.used)
			dt_G = KF->getState_minus()(6, 0);
		if (cfg.BDS_Cfg.used)
			dt_C = KF->getState_minus()(6, 0);
	}
	else if (cfg.SYS_num == 2)
	{
		dt_G = KF->getState_minus()(6, 0);
		dt_C = KF->getState_minus()(7, 0);
	}

	// 根据配置信息设置GPS的状态向量
	if (cfg.GPS_Cfg.used)
	{
		temp_ref(3, 0) = dt_G;
		val = setup_KF(this, cfg, SYS_GPS);
	}

	// 根据配置信息设置BDS的状态向量
	if (cfg.BDS_Cfg.used)
	{
		temp_ref(3, 0) = dt_C;
		val = setup_KF(this, cfg, SYS_BDS);
	}

	return val;
}

void processString(std::string* str, int* gps, int* bds)
{
	char prevChar = '\0'; // 上一个字符初始化为非数字
	// 遍历字符串
	for (size_t i = 0; i < str->length(); ++i)
	{
		// 如果当前字符是数字
		if (isdigit((*str)[i]))
		{
			int num = 0;
			// 解析数字
			while (i < str->length() && isdigit((*str)[i]))
			{
				num = num * 10 + ((*str)[i] - '0');
				++i;
			}
			// 根据上一个字符将对应的数组的值设置为1
			if (prevChar != '\0')
			{
				if (prevChar == 'G' && num > 0 && num <= GPS_SAT_QUAN)
				{
					gps[num - 1] = 1;
				}
				else if (prevChar == 'C' && num > 0 && num <= BDS_SAT_QUAN)
				{
					bds[num - 1] = 1;
				}
			}
		}
		// 更新上一个字符
		prevChar = (*str)[i];
	}
}


unsigned int Select_Common_Sates(DATA_SET* rove, DATA_SET* base, RTK_DATA* rtk, Configure cfg)
{
	int Rove_GPS[GPS_SAT_QUAN];
	int Rove_BDS[BDS_SAT_QUAN];
	int Base_GPS[GPS_SAT_QUAN];
	int Base_BDS[BDS_SAT_QUAN];
	for (int i = 0; i < GPS_SAT_QUAN; i++)
	{
		Rove_GPS[i] = 0;
		Base_GPS[i] = 0;
	}
	for (int i = 0; i < BDS_SAT_QUAN; i++)
	{
		Rove_BDS[i] = 0;
		Base_BDS[i] = 0;
	}

	processString(rove->LS_SATES, Rove_GPS, Rove_BDS);
	processString(base->LS_SATES, Base_GPS, Base_BDS);

	if (cfg.GPS_Cfg.used)
	{
		for (Satellate* sate : rove->range->GPS_SATE)
		{
			if (Rove_GPS[sate->PRN - 1] * Base_GPS[sate->PRN - 1] == 1)
			{
				rtk->Rove_GPS.push_back(sate);
			}
		}
		for (Satellate* sate : base->range->GPS_SATE)
		{
			if (Rove_GPS[sate->PRN - 1] * Base_GPS[sate->PRN - 1] == 1)
			{
				rtk->Base_GPS.push_back(sate);
			}
		}
	}

	if(cfg.BDS_Cfg.used)
	{
		for (Satellate* sate : rove->range->BDS_SATE)
		{
			if (Rove_BDS[sate->PRN - 1] * Base_BDS[sate->PRN - 1] == 1)
			{
				if (sate->PRN < 6 || (sate->PRN > 58 && sate->PRN < 63))
					continue;
				rtk->Rove_BDS.push_back(sate);
			}
		}
		for (Satellate* sate : base->range->BDS_SATE)
		{
			if (Rove_BDS[sate->PRN - 1] * Base_BDS[sate->PRN - 1] == 1)
			{
				if (sate->PRN < 6 || (sate->PRN > 58 && sate->PRN < 63))
					continue;
				rtk->Base_BDS.push_back(sate);
			}
		}
	}
	
	XYZ xyz = get_XYZ(base->Pos);
	rtk->Base_appro_pos->X = xyz.X;
	rtk->Base_appro_pos->Y = xyz.Y;
	rtk->Base_appro_pos->Z = xyz.Z;
	xyz = get_XYZ(rove->Pos);
	rtk->Rove_appro_pos->X = xyz.X;
	rtk->Rove_appro_pos->Y = xyz.Y;
	rtk->Rove_appro_pos->Z = xyz.Z;
	if (rtk->KF->x_hat_.rows() == 3 && abs(rtk->KF->x_hat_(0, 0)) < 1e-4)
	{
		rtk->KF->x_hat_(0, 0) = xyz.X;
		rtk->KF->x_hat_(1, 0) = xyz.Y;
		rtk->KF->x_hat_(2, 0) = xyz.Z;
		*rtk->KF_pos = xyz;
	}
	if (cfg.GPS_Cfg.used)sort(rtk->Rove_GPS.begin(), rtk->Rove_GPS.end(), Satellate::ELEV_comp);
	if (cfg.BDS_Cfg.used)sort(rtk->Rove_BDS.begin(), rtk->Rove_BDS.end(), Satellate::ELEV_comp);
	if (cfg.GPS_Cfg.used)rtk->Ref_PRN_GPS_NEW = rtk->Rove_GPS[0]->PRN;
	if (cfg.BDS_Cfg.used)rtk->Ref_PRN_BDS_NEW = rtk->Rove_BDS[0]->PRN;
	if (cfg.GPS_Cfg.used)sort(rtk->Rove_GPS.begin(), rtk->Rove_GPS.end(), Satellate::compare);
	if (cfg.BDS_Cfg.used)sort(rtk->Rove_BDS.begin(), rtk->Rove_BDS.end(), Satellate::compare);
	if (cfg.GPS_Cfg.used)sort(rtk->Base_GPS.begin(), rtk->Base_GPS.end(), Satellate::compare);
	if (cfg.BDS_Cfg.used)sort(rtk->Base_BDS.begin(), rtk->Base_BDS.end(), Satellate::compare);
	rtk->OBSTIME->SecOfWeek = base->OBSTIME->SecOfWeek;
	rtk->OBSTIME->Week = base->OBSTIME->Week;

	return 1;
}

int get_mat(MatrixXd& B, MatrixXd& L, MatrixXd& Q, GPSTIME* time, XYZ* RcvPos, vector<Satellate*> sates,
            int sys, Configure cfg, RTK_DATA* rtk, vector<int>& PRN_used)
{
	MatrixXd B_new = MatrixXd::Zero(1, 3);
	MatrixXd l_new = MatrixXd::Zero(1, 1);
	MatrixXd l_p = MatrixXd::Zero(0, 1);
	MatrixXd L_2 = MatrixXd::Zero(0, 1);
	MatrixXd L_p_2 = MatrixXd::Zero(0, 1);
	double f1 = 0;
	double f2 = 0;
	double p_pos = 0;
	int Max_elev_row = 0;
	int ref_prn;
	int ROWS = 0;
	switch (sys)
	{
	case SYS_GPS:
		f1 = cfg.GPS_Cfg.f1;
		f2 = cfg.GPS_Cfg.f2;
		ref_prn = rtk->Ref_PRN_GPS_NEW;
		break;
	case SYS_BDS:
		f1 = cfg.BDS_Cfg.f1;
		f2 = cfg.BDS_Cfg.f2;
		ref_prn = rtk->Ref_PRN_BDS_NEW;
		break;
	default:
		break;
	}
	for (Satellate* sate : sates)
	{
		int prn = sate->PRN;

		// 查找卫星对应的频率索引
		int Index1 = 0;
		int Index2 = 0;
		while (CODE2FREQ(decode_SYN(sate->SYS, sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
			Index1++;
		if (Index1 == MAXNUM)
			continue;
		while (CODE2FREQ(decode_SYN(sate->SYS, sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
			Index2++;
		if (Index2 == MAXNUM)
			continue;
		//筛选SNR
		if (sate->SNR[Index1] < 30)
			continue;
		if (cfg.phase_num == 2)
		{
			if (sate->SNR[Index2] < 30)
				continue;
		}
		if (sate->pos.X == 0)
			continue;
		// 获取测量值
		double measure_p1 = sate->PSERA[Index1];
		double measure_l1 = sate->PHASE[Index1] * 1e-6 * velocity_c / f1;
		double measure_p2 = sate->PSERA[Index2];
		double measure_l2 = sate->PHASE[Index2] * 1e-6 * velocity_c / f2;
		if (!measure_p1)
			continue;
		if (!measure_l1)
			continue;
		if (cfg.phase_num == 2)
		{
			if (!measure_p2)
				continue;
			if (!measure_l2)
				continue;
		}

		// 计算方位角对应的波长和长度
		double len = Len(RcvPos, &sate->pos);
		double w_pos = measure_l1 - len;
		// 构造最小二乘问题的设计矩阵B和观测矩阵l
		double l = (RcvPos->X - sate->pos.X) / len;
		double m = (RcvPos->Y - sate->pos.Y) / len;
		double n = (RcvPos->Z - sate->pos.Z) / len;
		B_new(0, 0) = l;
		B_new(0, 1) = m;
		B_new(0, 2) = n;
		l_new(0, 0) = w_pos;
		double temp_elev = sate->ELEV;
		if (sate->PRN==ref_prn)
		{
			Max_elev_row = ROWS;
		}
		ROWS++;

		// 更新B、L、P矩阵
		B.conservativeResize(B.rows() + 1, B.cols());
		B.bottomRows(1) = B_new;
		L.conservativeResize(L.rows() + 1, L.cols());
		L.bottomRows(1) = l_new;
		w_pos = measure_p1 - len;
		l_new(0, 0) = w_pos;
		l_p.conservativeResize(l_p.rows() + 1, l_p.cols());
		l_p.bottomRows(1) = l_new;
		if (cfg.phase_num == 2)
		{
			w_pos = measure_l2 - len;
			l_new(0, 0) = w_pos;
			L_2.conservativeResize(L_2.rows() + 1, L_2.cols());
			L_2.bottomRows(1) = l_new;
			w_pos = measure_p2 - len;
			l_new(0, 0) = w_pos;
			L_p_2.conservativeResize(L_p_2.rows() + 1, L_p_2.cols());
			L_p_2.bottomRows(1) = l_new;
		}
		if (ROWS == 1)
		{
			Q = MatrixXd::Identity(1, 1);
			Q(0, 0) = P_height1(temp_elev);
		}
		else
		{
			MatrixXd temp_Q = MatrixXd::Identity(ROWS, ROWS);
			temp_Q.block(0, 0, Q.rows(), Q.cols()) = Q;
			temp_Q(ROWS - 1, ROWS - 1) = P_height1(temp_elev);
			Q = temp_Q;
		}

		PRN_used.push_back(sate->PRN);
	}
	L.conservativeResize(L.rows() * 2 * cfg.phase_num, 1);
	if (cfg.phase_num == 1)
	{
		L.bottomRows(ROWS) = l_p;
	}
	else if (cfg.phase_num == 2)
	{
		L.block(ROWS, 0, ROWS, 1) = L_2;
		L.block(ROWS * 2, 0, ROWS, 1) = l_p;
		L.bottomRows(ROWS) = L_p_2;
	}
	return Max_elev_row;
}

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows)
	{
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
			matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
	}

	matrix.conservativeResize(numRows, numCols);
}

void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;

	if (colToRemove < numCols)
	{
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
			matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
	}

	matrix.conservativeResize(numRows, numCols);
}


unsigned int RTK_Set(DATA_SET* base, RTK_DATA* rtk, Configure cfg, int sys, int method)
{
	int ROWS = 0;
	vector<Satellate*> Rove_Sates;
	vector<Satellate*> Base_Sates;
	double f1 = 0;
	double f2 = 0;
	// 根据卫星系统类型选择相应的卫星数据和星历数据
	switch (sys)
	{
	case SYS_GPS:
		Rove_Sates = rtk->Rove_GPS;
		Base_Sates = rtk->Base_GPS;
		f1 = cfg.GPS_Cfg.f1;
		f2 = cfg.GPS_Cfg.f2;
		break;
	case SYS_BDS:
		Rove_Sates = rtk->Rove_BDS;
		Base_Sates = rtk->Base_BDS;
		f1 = cfg.BDS_Cfg.f1;
		f2 = cfg.BDS_Cfg.f2;
		break;
	default:
		return 0;
		break;
	}

	XYZ* rove_pos = new XYZ();
	switch (method)
	{
	case Method_LS:
		rove_pos->X = rtk->Rove_appro_pos->X;
		rove_pos->Y = rtk->Rove_appro_pos->Y;
		rove_pos->Z = rtk->Rove_appro_pos->Z;
		break;
	case Method_KF:
		rove_pos->X = rtk->KF->x_hat_minus_(0, 0);
		rove_pos->Y = rtk->KF->x_hat_minus_(1, 0);
		rove_pos->Z = rtk->KF->x_hat_minus_(2, 0);
		break;
	default: break;
	}

	MatrixXd B_Rove = MatrixXd::Zero(0, 3);
	MatrixXd l_Rove = MatrixXd::Zero(0, 1);
	MatrixXd Q_Rove = MatrixXd::Zero(1, 1);
	MatrixXd B_Base = MatrixXd::Zero(0, 3);
	MatrixXd l_Base = MatrixXd::Zero(0, 1);
	MatrixXd Q_Base = MatrixXd::Zero(1, 1);
	vector<int> prn_used;
	get_mat(B_Base, l_Base, Q_Base, rtk->OBSTIME, rtk->Base_appro_pos, Base_Sates, sys, cfg, rtk, prn_used);
	prn_used.clear();
	int elev_row = get_mat(B_Rove, l_Rove, Q_Rove, rtk->OBSTIME, rove_pos, Rove_Sates, sys, cfg, rtk, prn_used);
	for (int prn : prn_used)
	{
		switch (method)
		{
		case Method_LS:
			switch (sys)
			{
			case SYS_GPS:
				if (prn == rtk->Ref_PRN_GPS_NEW)
					continue;
				*(rtk->LS_SATES) += "G" + std::to_string(prn);
				break;
			case SYS_BDS:
				if (prn == rtk->Ref_PRN_BDS_NEW)
					continue;
				*(rtk->LS_SATES) += "C" + std::to_string(prn);
				break;
			default:
				break;
			}
			break;
		case Method_KF:
			switch (sys)
			{
			case SYS_GPS:
				if (prn == rtk->Ref_PRN_GPS_NEW)
					continue;
				*(rtk->KF_SATES) += "G" + std::to_string(prn);
				rtk->GPS_PRN.push_back(prn);
				break;
			case SYS_BDS:
				if (prn == rtk->Ref_PRN_BDS_NEW)
					continue;
				*(rtk->KF_SATES) += "C" + std::to_string(prn);
				rtk->BDS_PRN.push_back(prn);
				break;
			default: break;
			}
			break;
		default: break;
		}
	}
	switch (method)
	{
	case Method_LS:
		switch (sys)
		{
		case SYS_GPS:
			rtk->LS_GPS_num = prn_used.size() - 1;
			break;
		case SYS_BDS:
			rtk->LS_BDS_num = prn_used.size() - 1;
			break;
		default:
			break;
		}
		break;
	case Method_KF:
		switch (sys)
		{
		case SYS_GPS:
			rtk->KF_GPS_num = prn_used.size() - 1;
			break;
		case SYS_BDS:
			rtk->KF_BDS_num = prn_used.size() - 1;
			break;
		default: break;
		}
		break;
	default: break;
	}
	//参考星单差
	ROWS = B_Rove.rows();
	MatrixXd B_Max_elev_row;
	B_Max_elev_row = B_Rove.row(elev_row);
	vector<MatrixXd> l_Rove_max_elev_row;
	vector<MatrixXd> l_Base_max_elev_row;
	for (int j = 0; j < 2 * cfg.phase_num; j++)
	{
		l_Rove_max_elev_row.push_back(l_Rove.row(elev_row + j * ROWS));
		l_Base_max_elev_row.push_back(l_Base.row(elev_row + j * ROWS));
	}
	for (int i = 0; i < ROWS; i++)
	{
		B_Rove.row(i) = B_Rove.row(i) - B_Max_elev_row;
		for (int j = 0; j < 2 * cfg.phase_num; j++)
		{
			l_Rove.row(i + j * ROWS) = -l_Rove_max_elev_row[j] + l_Rove.row(i + j * ROWS);
			l_Base.row(i + j * ROWS) = -l_Base_max_elev_row[j] + l_Base.row(i + j * ROWS);
		}
	}

	RemoveRow(B_Rove, elev_row);
	for (int i = 0; i < 2 * cfg.phase_num; i++)
	{
		RemoveRow(l_Rove, elev_row + i * ROWS - i);
		RemoveRow(l_Base, elev_row + i * ROWS - i);
	}

	ROWS--;
	B_Rove.conservativeResize(ROWS * 2 * cfg.phase_num, 3);
	for (int i = 1; i < cfg.phase_num * 2; i++)
	{
		B_Rove.block(i * ROWS, 0, ROWS, 3) = B_Rove.topRows(ROWS);
	}

	//函数模型
	MatrixXd B_dd = B_Rove;
	B_dd.conservativeResize(B_dd.rows(), 3 + ROWS * cfg.phase_num);
	B_dd.block(0, 3, B_dd.rows(), ROWS * cfg.phase_num) = MatrixXd::Zero(B_dd.rows(), ROWS * cfg.phase_num);
	B_dd.block(0, 3, ROWS, ROWS) = 1e-6 * MatrixXd::Identity(ROWS, ROWS) * velocity_c / f1;
	B_dd.block(ROWS * cfg.phase_num, 3, ROWS * cfg.phase_num, ROWS * cfg.phase_num) = MatrixXd::Zero(
		ROWS * cfg.phase_num, ROWS * cfg.phase_num);
	if (cfg.phase_num == 2)
		B_dd.block(ROWS, 3 + ROWS, ROWS, ROWS) = MatrixXd::Identity(ROWS, ROWS) * 1e-6 * velocity_c / f2;
	MatrixXd l_dd = l_Rove - l_Base;

	//随机模型
	MatrixXd Q_d = Q_Base + Q_Rove;
	MatrixXd Q_dd = MatrixXd::Ones(ROWS, ROWS) * Q_d(elev_row, elev_row);
	RemoveColumn(Q_d, elev_row);
	RemoveRow(Q_d, elev_row);
	Q_dd = Q_d + Q_dd;
	//std::cout << Q_dd << endl;
	MatrixXd temp_P = MatrixXd::Zero(ROWS * 2 * cfg.phase_num, ROWS * 2 * cfg.phase_num);
	MatrixXd temp_Q = MatrixXd::Zero(ROWS * 2 * cfg.phase_num, ROWS * 2 * cfg.phase_num);
	MatrixXd P_dd = Q_dd.inverse();
	//std::cout << P_dd << endl;
	if (cfg.phase_num == 1)
	{
		temp_P.block(0, 0, ROWS, ROWS) = P_dd;
		temp_P.block(ROWS, ROWS, ROWS, ROWS) = P_dd / cfg.LrateP;
		temp_Q.block(0, 0, ROWS, ROWS) = Q_dd;
		temp_Q.block(ROWS, ROWS, ROWS, ROWS) = Q_dd * cfg.LrateP;
	}
	else if (cfg.phase_num == 2)
	{
		temp_P.block(0, 0, ROWS, ROWS) = P_dd;
		temp_P.block(ROWS, ROWS, ROWS, ROWS) = P_dd;
		temp_P.block(ROWS * 2, ROWS * 2, ROWS, ROWS) = P_dd / cfg.LrateP;
		temp_P.block(ROWS * 3, ROWS * 3, ROWS, ROWS) = P_dd / cfg.LrateP;
		temp_Q.block(0, 0, ROWS, ROWS) = Q_dd;
		temp_Q.block(ROWS, ROWS, ROWS, ROWS) = Q_dd;
		temp_Q.block(ROWS * 2, ROWS * 2, ROWS, ROWS) = Q_dd * cfg.LrateP;
		temp_Q.block(ROWS * 3, ROWS * 3, ROWS, ROWS) = Q_dd * cfg.LrateP;
	}
	P_dd = temp_P;
	Q_dd = temp_Q;

	switch (method)
	{
	case Method_LS:
		rtk->LS->set_B_Pos(B_dd);
		rtk->LS->set_l(l_dd);
		rtk->LS->set_P(P_dd);
		break;
	case Method_KF:
		rtk->KF->set_H(B_dd, MODE_RTK);
		rtk->KF->set_Z(l_dd);
		rtk->KF->set_R(Q_dd, MODE_RTK);
		break;
	default:
		break;
	}


	delete rove_pos;
	return 1;
}

unsigned int RTK_Solve(DATA_SET* base, RTK_DATA* rtk, Configure cfg)
{
	if (cfg.RTK_LS_used)
	{
		for (int i = 0; i < 4; i++)
		{
			rtk->LS->reset();
			*rtk->LS_SATES = "";

			if (cfg.GPS_Cfg.used)
				RTK_Set(base, rtk, cfg, SYS_GPS, Method_LS);
			/*std::cout << rtk->LS->B << endl;
			std::cout << rtk->LS->l << endl;
			std::cout << rtk->LS->P << endl;*/
			if (cfg.BDS_Cfg.used)
				RTK_Set(base, rtk, cfg, SYS_BDS, Method_LS);
			/*std::cout << rtk->LS->B << endl;
			std::cout << rtk->LS->l << endl;
			std::cout << rtk->LS->P << endl;*/
			rtk->LS->LS();
			rtk->Rove_appro_pos->X += rtk->LS->X(0, 0);
			rtk->Rove_appro_pos->Y += rtk->LS->X(1, 0);
			rtk->Rove_appro_pos->Z += rtk->LS->X(2, 0);
			rtk->LS_result = Success_Solve;
		}
		int n = rtk->LS->X.rows() - 3;
		MatrixXd Qaa = rtk->LS->Qxx.bottomRightCorner(n, n);
		double* dQaa = matrixXdToDoublePtr(Qaa);
		double* dF = new double[n * 2]; // 分配 n*2 个连续的 double 内存空间
		MatrixXd F = MatrixXd::Zero(n, 2);
		double* ds = new double[2]; // 分配 2 个连续的 double 内存空间
		VectorXd s = VectorXd::Zero(2);
		VectorXd a = rtk->LS->X.bottomRows(n);
		double* da = vectorXdToDoublePtr(a);
		int val = lambda(n, 2, da, dQaa, dF, ds);
		F = doublePtrToMatrixXd1(dF, n, 2);
		s = doublePtrToVectorXd(ds, 2);
		Qaa = doublePtrToMatrixXd1(dQaa, n, n);
		free(dF);
		free(da);
		free(ds);
		free(dQaa);
		if (!val)
		{
			rtk->Ratio_LS = s[0] > s[1] ? s[0] / s[1] : s[1] / s[0];
			if (rtk->Ratio_LS < cfg.Ratio_thresh)
			{
				rtk->LS_IS_FIXED = false;
			}
			else
			{
				rtk->LS_IS_FIXED = true;
				MatrixXd Fix_pos = get_XYZ(*rtk->Rove_appro_pos) - rtk->LS->Qxx.topRightCorner(3, n) * Qaa.inverse() * (
					a -
					F.col(0));
				rtk->LS->Qxx.topLeftCorner(3, 3) += -rtk->LS->Qxx.topRightCorner(3, n) * Qaa.inverse() * rtk->LS->Qxx.
					bottomLeftCorner(n, 3);
				*rtk->Rove_appro_pos = get_XYZ(Fix_pos);
				rtk->LS->X.bottomRows(n) = F.col(0);
			}
		}
	}

	if (cfg.RTK_KF_used)
	{
		*rtk->KF_SATES = "";

		if (rtk->KF->x_hat_.rows() == 3)
		{
			rtk->Ref_PRN_GPS_OLD = rtk->Ref_PRN_GPS_NEW;
			rtk->Ref_PRN_BDS_OLD = rtk->Ref_PRN_BDS_NEW;
			if (cfg.GPS_Cfg.used)
			{
				rtk->GPS_PRN.clear();
				int last_row = rtk->KF->x_hat_.rows();
				rtk->KF->x_hat_.conservativeResize(last_row + (rtk->Rove_GPS.size() - 1) * cfg.phase_num, 1);
				double f1 = cfg.GPS_Cfg.f1;
				double f2 = cfg.GPS_Cfg.f2;
				double lambda1 = 1e-6 * velocity_c / f1;
				double lambda2 = 1e-6 * velocity_c / f2;
				int num = rtk->Rove_GPS.size() - 1;
				int i = 0;
				for (Satellate* sate : rtk->Rove_GPS)
				{
					int Index1 = 0;
					int Index2 = 0;

					// 查找卫星对应的频率索引
					while (CODE2FREQ(decode_SYN(SYS_GPS, sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
						Index1++;
					while (CODE2FREQ(decode_SYN(SYS_GPS, sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
						Index2++;
					if (sate->PRN == rtk->Ref_PRN_GPS_NEW)
						continue;
					
					rtk->KF->x_hat_(i + last_row, 0) = sate->PHASE[Index1] - sate->PSERA[Index1] / lambda1;
					if (cfg.phase_num == 2)
						rtk->KF->x_hat_(last_row + num + i, 0) = sate->PHASE[Index2] - sate->PSERA[Index2] / lambda2;
					rtk->GPS_PRN.push_back(sate->PRN);
					i++;
				}
				MatrixXd P_new = MatrixXd::Zero(rtk->KF->x_hat_.rows(), rtk->KF->x_hat_.rows());
				P_new.block(0, 0, last_row, last_row) = rtk->KF->P_;
				P_new.block(last_row, last_row, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda1);
				if (cfg.phase_num == 2)
					P_new.block(last_row + num, last_row + num, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda2);
				rtk->KF->P_ = P_new;
			}
			if (cfg.BDS_Cfg.used)
			{
				rtk->BDS_PRN.clear();
				int last_row = rtk->KF->x_hat_.rows();
				rtk->KF->x_hat_.conservativeResize(last_row + (rtk->Rove_BDS.size() - 1) * cfg.phase_num, 1);
				double f1 = cfg.BDS_Cfg.f1;
				double f2 = cfg.BDS_Cfg.f2;
				double lambda1 = 1e-6 * velocity_c / f1;
				double lambda2 = 1e-6 * velocity_c / f2;
				int num = rtk->Rove_BDS.size() - 1;
				int i = 0;
				for (Satellate* sate : rtk->Rove_BDS)
				{
					int Index1 = 0;
					int Index2 = 0;

					// 查找卫星对应的频率索引
					while (CODE2FREQ(decode_SYN(SYS_BDS, sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
						Index1++;
					while (CODE2FREQ(decode_SYN(SYS_BDS, sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
						Index2++;
					if (sate->PRN == rtk->Ref_PRN_BDS_NEW)
						continue;

					rtk->KF->x_hat_(i + last_row, 0) = sate->PHASE[Index1] - sate->PSERA[Index1] / lambda1;
					if (cfg.phase_num == 2)
						rtk->KF->x_hat_(last_row + num + i, 0) = sate->PHASE[Index2] - sate->PSERA[Index2] / lambda2;
					rtk->BDS_PRN.push_back(sate->PRN);
					i++;
				}
				MatrixXd P_new = MatrixXd::Zero(rtk->KF->x_hat_.rows(), rtk->KF->x_hat_.rows());
				P_new.block(0, 0, last_row, last_row) = rtk->KF->P_;
				P_new.block(last_row, last_row, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda1);
				if (cfg.phase_num == 2)
					P_new.block(last_row + num, last_row + num, num, num) = MatrixXd::Identity(num, num) *
						SQR(cfg.initial_covariance / lambda2);
				rtk->KF->P_ = P_new;
			}
			rtk->KF->Q_ = MatrixXd::Zero(rtk->KF->x_hat_.rows(), rtk->KF->x_hat_.rows());
		}
		map<int, int> new_prn_GPS, new_prn_BDS;
		map<int, int> new_sate_GPS, new_sate_BDS;
		MatrixXd GPS_A, BDS_A, RTK_A;
		MatrixXd GPS_D, BDS_D, RTK_D;
		RTK_A = MatrixXd::Identity(3, 3);
		RTK_D = RTK_A;
		bool GPS_reset = false;
		bool BDS_reset = false;
		if (cfg.GPS_Cfg.used)
		{
			int i = 0;
			for (Satellate* sate : rtk->Rove_GPS)
			{
				if (sate->PRN == rtk->Ref_PRN_GPS_NEW)
				{
					if (sate->Cycle_jump)
						GPS_reset = true;
					continue;
				}
				new_prn_GPS.insert({sate->PRN, i});
				i++;
			}
			new_sate_GPS = new_prn_GPS;
			int new_ref_col = RTK_getA(rtk->GPS_PRN, new_prn_GPS, GPS_A, rtk->Ref_PRN_GPS_NEW);
			if (rtk->Ref_PRN_GPS_OLD != rtk->Ref_PRN_GPS_NEW && new_ref_col == -1)
			{
				rtk->Fix_N = MatrixXd::Zero(0, 1);
				int last_row = 3;
				int num = rtk->Rove_GPS.size() - 1;
				MatrixXd trans = MatrixXd::Zero(rtk->KF->x_hat_.rows() + (num - rtk->KF_GPS_num) * cfg.phase_num, rtk->KF->x_hat_.rows());
				trans.topLeftCorner(3, 3) = MatrixXd::Identity(3, 3);
				trans.bottomRightCorner(rtk->KF_BDS_num* cfg.phase_num, rtk->KF_BDS_num* cfg.phase_num) = MatrixXd::Identity(rtk->KF_BDS_num * cfg.phase_num, rtk->KF_BDS_num * cfg.phase_num);
				rtk->KF->x_hat_ = trans * rtk->KF->x_hat_;
				rtk->KF->P_ = trans * rtk->KF->P_ * trans.transpose();
				double f1 = cfg.GPS_Cfg.f1;
				double f2 = cfg.GPS_Cfg.f2;
				double lambda1 = 1e-6 * velocity_c / f1;
				double lambda2 = 1e-6 * velocity_c / f2;
				int i = 0;
				for (Satellate* sate : rtk->Rove_GPS)
				{
					int Index1 = 0;
					int Index2 = 0;

					// 查找卫星对应的频率索引
					while (CODE2FREQ(decode_SYN(SYS_GPS, sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
						Index1++;
					while (CODE2FREQ(decode_SYN(SYS_GPS, sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
						Index2++;
					if (sate->PRN == rtk->Ref_PRN_GPS_NEW)
						continue;

					rtk->KF->x_hat_(i, 0) = sate->PHASE[Index1] - sate->PSERA[Index1] / lambda1;
					if (cfg.phase_num == 2)
						rtk->KF->x_hat_(num + i, 0) = sate->PHASE[Index2] - sate->PSERA[Index2] / lambda2;
					i++;
				}
				rtk->KF->P_.block(last_row, last_row, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda1);
				if (cfg.phase_num == 2)
					rtk->KF->P_.block(last_row + num, last_row + num, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda2);
				GPS_A = MatrixXd::Identity(num, num);
			}
			GPS_D = GPS_A;
			GPS_A = MatrixXd::Identity(GPS_A.rows(), GPS_A.rows());
			if (rtk->Ref_PRN_GPS_OLD != rtk->Ref_PRN_GPS_NEW && new_ref_col!=-1)
				GPS_D.col(new_ref_col) = MatrixXd::Ones(GPS_D.rows(), 1) * -1;
			rtk->Ref_PRN_GPS_OLD = rtk->Ref_PRN_GPS_NEW;
			MatrixXd temp_A = MatrixXd::Zero(RTK_A.rows() + GPS_A.rows() * cfg.phase_num,
			                                 RTK_A.cols() + GPS_A.cols() * cfg.phase_num);
			MatrixXd temp_D = MatrixXd::Zero(RTK_D.rows() + GPS_D.rows() * cfg.phase_num,
			                                 RTK_D.cols() + GPS_D.cols() * cfg.phase_num);
			temp_A.block(0, 0, RTK_A.rows(), RTK_A.cols()) = RTK_A;
			temp_D.block(0, 0, RTK_D.rows(), RTK_D.cols()) = RTK_D;
			for (int i = 0; i < cfg.phase_num; i++)
			{
				temp_A.block(RTK_A.rows() + GPS_A.rows() * i, RTK_A.cols() + GPS_A.cols() * i, GPS_A.rows(),
				             GPS_A.cols()) = GPS_A;
				temp_D.block(RTK_D.rows() + GPS_D.rows() * i, RTK_D.cols() + GPS_D.cols() * i, GPS_D.rows(),
				             GPS_D.cols()) = GPS_D;
			}
			RTK_A = temp_A;
			RTK_D = temp_D;
		}
		if (cfg.BDS_Cfg.used)
		{
			int i = 0;
			for (Satellate* sate : rtk->Rove_BDS)
			{
				if (sate->PRN == rtk->Ref_PRN_BDS_NEW)
				{
					if (sate->Cycle_jump)
						BDS_reset = true;
					continue;
				}
				new_prn_BDS.insert({sate->PRN, i});
				i++;	
			}
			new_sate_BDS = new_prn_BDS;
			int new_ref_col = RTK_getA(rtk->BDS_PRN, new_prn_BDS, BDS_A, rtk->Ref_PRN_BDS_NEW);
			if (rtk->Ref_PRN_BDS_OLD != rtk->Ref_PRN_BDS_NEW && new_ref_col == -1)
			{
				rtk->Fix_N = MatrixXd::Zero(0, 1);
				int last_row = 3;
				int num = rtk->Rove_BDS.size() - 1;
				if (cfg.GPS_Cfg.used)
					last_row += (rtk->Rove_GPS.size() - 1) * cfg.phase_num;
				MatrixXd trans = MatrixXd::Zero(rtk->KF->x_hat_.rows() + (num - rtk->KF_BDS_num) * cfg.phase_num, rtk->KF->x_hat_.rows());
				trans.topLeftCorner(last_row, last_row) = MatrixXd::Identity(last_row, last_row);
				rtk->KF->x_hat_ = trans * rtk->KF->x_hat_;
				rtk->KF->P_ = trans * rtk->KF->P_ * trans.transpose();
				double f1 = cfg.BDS_Cfg.f1;
				double f2 = cfg.BDS_Cfg.f2;
				double lambda1 = 1e-6 * velocity_c / f1;
				double lambda2 = 1e-6 * velocity_c / f2;
				int i = 0;
				for (Satellate* sate : rtk->Rove_BDS)
				{
					int Index1 = 0;
					int Index2 = 0;

					// 查找卫星对应的频率索引
					while (CODE2FREQ(decode_SYN(SYS_BDS, sate->SYG_TYPE[Index1])) != f1 && Index1 < MAXNUM)
						Index1++;
					while (CODE2FREQ(decode_SYN(SYS_BDS, sate->SYG_TYPE[Index2])) != f2 && Index2 < MAXNUM)
						Index2++;
					if (sate->PRN == rtk->Ref_PRN_BDS_NEW)
						continue;

					rtk->KF->x_hat_(i + last_row, 0) = sate->PHASE[Index1] - sate->PSERA[Index1] / lambda1;
					if (cfg.phase_num == 2)
						rtk->KF->x_hat_(last_row + num + i, 0) = sate->PHASE[Index2] - sate->PSERA[Index2] / lambda2;
					i++;
				}
				rtk->KF->P_.block(last_row, last_row, num, num) = MatrixXd::Identity(num, num) * SQR(cfg.initial_covariance / lambda1);
				if (cfg.phase_num == 2)
					rtk->KF->P_.block(last_row + num, last_row + num, num, num) = MatrixXd::Identity(num, num) *
					SQR(cfg.initial_covariance / lambda2);
				BDS_A = MatrixXd::Identity(num, num);
			}
			BDS_D = BDS_A;
			BDS_A = MatrixXd::Identity(BDS_A.rows(), BDS_A.rows());
			if (rtk->Ref_PRN_BDS_OLD != rtk->Ref_PRN_BDS_NEW && new_ref_col!=-1)
				BDS_D.col(new_ref_col) = MatrixXd::Ones(BDS_D.rows(), 1) * -1;
			rtk->Ref_PRN_BDS_OLD = rtk->Ref_PRN_BDS_NEW;
			MatrixXd temp_A = MatrixXd::Zero(RTK_A.rows() + BDS_A.rows() * cfg.phase_num,
			                                 RTK_A.cols() + BDS_A.cols() * cfg.phase_num);
			MatrixXd temp_D = MatrixXd::Zero(RTK_D.rows() + BDS_D.rows() * cfg.phase_num,
			                                 RTK_D.cols() + BDS_D.cols() * cfg.phase_num);
			temp_A.block(0, 0, RTK_A.rows(), RTK_A.cols()) = RTK_A;
			temp_D.block(0, 0, RTK_D.rows(), RTK_D.cols()) = RTK_D;
			for (int i = 0; i < cfg.phase_num; i++)
			{
				temp_A.block(RTK_A.rows() + BDS_A.rows() * i, RTK_A.cols() + BDS_A.cols() * i, BDS_A.rows(),
				             BDS_A.cols()) = BDS_A;
				temp_D.block(RTK_D.rows() + BDS_D.rows() * i, RTK_D.cols() + BDS_D.cols() * i, BDS_D.rows(),
				             BDS_D.cols()) = BDS_D;
			}
			RTK_A = temp_A;
			RTK_D = temp_D;
		}
		rtk->KF->P_ = RTK_D * rtk->KF->P_ * RTK_D.transpose();
		rtk->KF->x_hat_ = RTK_D * rtk->KF->x_hat_;

		if (cfg.GPS_Cfg.used)
		{
			double f[2] = {cfg.GPS_Cfg.f1, cfg.GPS_Cfg.f2};
			double lambda[2] = {1e-6 * velocity_c / f[0], 1e-6 * velocity_c / f[1]};
			int Index[2] = {0, 0};
			if (!new_prn_GPS.empty())
			{
				for (auto prn : new_prn_GPS)
				{
					Index[0] = 0;
					Index[1] = 0;
					while (CODE2FREQ(decode_SYN(SYS_GPS, rtk->Rove_GPS[prn.second]->SYG_TYPE[Index[0]])) != f[0] &&
						Index[0]
						< MAXNUM)
						Index[0]++;
					while (CODE2FREQ(decode_SYN(SYS_GPS, rtk->Rove_GPS[prn.second]->SYG_TYPE[Index[1]])) != f[1] &&
						Index[1]
						< MAXNUM)
						Index[1]++;
					for (int i = 0; i < cfg.phase_num; i++)
					{
						rtk->KF->x_hat_(3 + (rtk->Rove_GPS.size() - 1) * i + prn.second, 0) =
							rtk->Rove_GPS[prn.second]->PHASE[Index[i]] - rtk->Rove_GPS[prn.second]->PSERA[Index[i]] / lambda[i];
						rtk->KF->P_(3 + (rtk->Rove_GPS.size() - 1) * i + prn.second,
						            3 + (rtk->Rove_GPS.size() - 1) * i + prn.second) = SQR(cfg.initial_covariance / lambda[i]);
					}
				}
			}
			int num = new_sate_GPS.size();
			for (auto prn : new_sate_GPS)
			{
				Index[0] = 0;
				Index[1] = 0;
				while (CODE2FREQ(decode_SYN(SYS_GPS, rtk->Rove_GPS[prn.second]->SYG_TYPE[Index[0]])) != f[0] &&
					Index[0]
					< MAXNUM)
					Index[0]++;
				while (CODE2FREQ(decode_SYN(SYS_GPS, rtk->Rove_GPS[prn.second]->SYG_TYPE[Index[1]])) != f[1] &&
					Index[1]
					< MAXNUM)
					Index[1]++;
				if (GPS_reset || rtk->Rove_GPS[prn.second]->Cycle_jump || rtk->Base_GPS[prn.second]->Cycle_jump)
				{
					for (int i = 0; i < cfg.phase_num; i++)
					{
						rtk->KF->x_hat_(3 + num * i + prn.second, 0) = rtk->Rove_GPS[prn.second]->PHASE[Index[i]] - rtk
							->Rove_GPS[prn.second]->PSERA[Index[i]] / lambda[i];
						rtk->KF->P_.row(3 + num * i + prn.second).setZero();
						rtk->KF->P_.col(3 + num * i + prn.second).setZero();
						rtk->KF->P_(3 + num * i + prn.second,
						            3 + num * i + prn.second) = SQR(cfg.initial_covariance / lambda[i]);
					}
				}
			}
		}
		if (cfg.BDS_Cfg.used)
		{
			double f[2] = {cfg.BDS_Cfg.f1, cfg.BDS_Cfg.f2};
			double lambda[2] = {1e-6 * velocity_c / f[0], 1e-6 * velocity_c / f[1]};
			int Index[2] = {0, 0};
			int pre_index = 0;
			if (cfg.GPS_Cfg.used)
				pre_index = rtk->Rove_GPS.size() - 1;
			if (!new_prn_BDS.empty())
			{
				for (auto prn : new_prn_BDS)
				{
					Index[0] = 0;
					Index[1] = 0;
					while (CODE2FREQ(decode_SYN(SYS_BDS, rtk->Rove_BDS[prn.second]->SYG_TYPE[Index[0]])) != f[0] &&
						Index[0]
						< MAXNUM)
						Index[0]++;
					while (CODE2FREQ(decode_SYN(SYS_BDS, rtk->Rove_BDS[prn.second]->SYG_TYPE[Index[1]])) != f[1] &&
						Index[1]
						< MAXNUM)
						Index[1]++;
					for (int i = 0; i < cfg.phase_num; i++)
					{
						rtk->KF->x_hat_(3 + pre_index * cfg.phase_num + (rtk->Rove_BDS.size() - 1) * i + prn.second, 0)
							= rtk->Rove_BDS[prn.second]->PHASE[Index[i]] - rtk->Rove_BDS[prn.second]->PSERA[Index[i]] / lambda[i];
						rtk->KF->P_(3 + pre_index * cfg.phase_num + (rtk->Rove_BDS.size() - 1) * i + prn.second,
						            3 + pre_index * cfg.phase_num + (rtk->Rove_BDS.size() - 1) * i + prn.second) = SQR(cfg.initial_covariance / lambda[i]);
					}
				}
			}
			int num = new_sate_BDS.size();
			for (auto prn : new_sate_BDS)
			{
				Index[0] = 0;
				Index[1] = 0;
				while (CODE2FREQ(decode_SYN(SYS_BDS, rtk->Rove_BDS[prn.second]->SYG_TYPE[Index[0]])) != f[0] &&
					Index[0]
					< MAXNUM)
					Index[0]++;
				while (CODE2FREQ(decode_SYN(SYS_BDS, rtk->Rove_BDS[prn.second]->SYG_TYPE[Index[1]])) != f[1] &&
					Index[1]
					< MAXNUM)
					Index[1]++;
				if (BDS_reset || rtk->Rove_BDS[prn.second]->Cycle_jump || rtk->Base_BDS[prn.second]->Cycle_jump)
				{
					for (int i = 0; i < cfg.phase_num; i++)
					{
						rtk->KF->x_hat_(3 + pre_index * cfg.phase_num + num * i + prn.second, 0)
							= rtk->Rove_BDS[prn.second]->PHASE[Index[i]] - rtk->Rove_BDS[prn.second]->PSERA[Index[i]] / lambda[i];
						rtk->KF->P_.row(3 + pre_index * cfg.phase_num + num * i + prn.second).
						     setZero();
						rtk->KF->P_.col(3 + pre_index * cfg.phase_num + num * i + prn.second).
						     setZero();
						rtk->KF->P_(3 + pre_index * cfg.phase_num + num * i + prn.second,
						            3 + pre_index * cfg.phase_num + num * i + prn.second) = SQR(cfg.initial_covariance / lambda[i]);
					}
				}
			}
		}
		rtk->KF->set_A(RTK_A);
		MatrixXd Q = MatrixXd::Identity(rtk->KF->A_.rows(), rtk->KF->A_.rows()) * SQR(cfg.process_noise_abg);
		Q.block(0, 0, 3, 3) = MatrixXd::Identity(3, 3) * SQR(cfg.process_noise_pos);
		rtk->KF->set_Q(Q);
		rtk->KF->predict();

		rtk->GPS_PRN.clear();
		rtk->BDS_PRN.clear();

		if (cfg.GPS_Cfg.used)
			RTK_Set(base, rtk, cfg, SYS_GPS, Method_KF);
		if (cfg.BDS_Cfg.used)
			RTK_Set(base, rtk, cfg, SYS_BDS, Method_KF);

		rtk->KF->z_ -= rtk->KF->H_.rightCols(rtk->KF->H_.cols() - 3) * rtk->KF->x_hat_minus_.bottomRows(
			rtk->KF->x_hat_minus_.rows() - 3);

		rtk->KF->update();
		rtk->KF_result = Success_Solve;

		*rtk->KF_pos = get_XYZ(rtk->KF->x_hat_.topRows(3));
		rtk->Q_baseLine = rtk->KF->P_.topLeftCorner(3, 3);

		

		int n = rtk->KF->x_hat_.rows() - 3;
		MatrixXd Qaa = rtk->KF->P_.bottomRightCorner(n, n);
		double* dQaa = matrixXdToDoublePtr(Qaa);
		double* dF = new double[n * 2]; // 分配 n*2 个连续的 double 内存空间
		MatrixXd F = MatrixXd::Zero(n, 2);
		double* ds = new double[2]; // 分配 2 个连续的 double 内存空间
		VectorXd s = VectorXd::Zero(2);
		VectorXd a = rtk->KF->x_hat_.bottomRows(n);
		double* da = vectorXdToDoublePtr(a);
		int val = lambda(n, 2, da, dQaa, dF, ds);
		F = doublePtrToMatrixXd1(dF, n, 2);
		s = doublePtrToVectorXd(ds, 2);
		Qaa = doublePtrToMatrixXd1(dQaa, n, n);
		free(dF);
		free(da);
		free(ds);
		free(dQaa);
		if (rtk->Fix_N.rows() != 0)
			rtk->Fix_N = RTK_D.bottomRightCorner(RTK_D.rows()-3, RTK_D.cols() - 3) * rtk->Fix_N;
		if (!val)
		{
			rtk->Ratio_KF = s[0] > s[1] ? s[0] / s[1] : s[1] / s[0];
			if (rtk->Ratio_KF < cfg.Ratio_thresh)
			{
				rtk->KF_IS_FIXED = false;
				rtk->Fix_N = MatrixXd::Zero(rtk->Fix_N.rows(), 1);
			}
			else
			{
				int num = 0;
				if (rtk->Fix_N.rows() != 0)
				{
					for (int i = 0; i < F.rows(); i++)
					{
						if (abs(rtk->Fix_N(i, 0) - F(i, 0)) < 1e-3 && rtk->Fix_N(i, 0) * F(i, 0) != 0)num++;
					}
				}
				rtk->KF_IS_FIXED = true;
				rtk->Fix_N = F.col(0);
				*rtk->KF_pos = *rtk->KF_pos -
					get_XYZ(rtk->KF->P_.topRightCorner(3, n) * Qaa.inverse() * (a - F.col(0)));
				rtk->Q_baseLine = rtk->KF->P_.topLeftCorner(3, 3) - rtk->KF->P_.topRightCorner(3, n) * Qaa.inverse() *
					rtk->KF->P_.bottomLeftCorner(n, 3);
				if (rtk->Ratio_KF>5&&(rtk->KF_GPS_num+rtk->KF_BDS_num)>7)
				{
					rtk->KF->x_hat_.topRows(3) = get_XYZ(*rtk->KF_pos);
					rtk->KF->P_.topLeftCorner(3, 3) = rtk->Q_baseLine;
					rtk->KF->P_.topRightCorner(3, n) = MatrixXd::Zero(3, n);
					rtk->KF->P_.bottomLeftCorner(n, 3) = MatrixXd::Zero(n, 3);
				}
				if (rtk->Ratio_KF > 10 && Cal_PDOP(rtk->Q_baseLine) < 2.0 && (rtk->KF_GPS_num + rtk->KF_BDS_num) > 10 && num > 5)
				{
					rtk->Stable_Num++;
					if (rtk->Stable_Num > 3)
					{
						rtk->KF->x_hat_.topRows(3) = get_XYZ(*rtk->KF_pos);
						rtk->KF->x_hat_.bottomRows(n) = F.col(0);
						rtk->KF->P_.topLeftCorner(3, 3) = rtk->Q_baseLine;
						rtk->KF->P_.topRightCorner(3, n) = MatrixXd::Zero(3, n);
						rtk->KF->P_.bottomLeftCorner(n, 3) = MatrixXd::Zero(n, 3);
						rtk->KF->P_.bottomRightCorner(n, n) = MatrixXd::Identity(n, n) * 1e-6;
					}
				}
				else
				{
					rtk->Stable_Num = 0;
				}
			}
		}
	}


	return 1;
}
