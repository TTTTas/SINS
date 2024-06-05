#pragma once
#include <Eigen/Dense>
#include "transform.h"
#include "KF.h"
#include "LS.h"
#include "Configure.h"
#include <string>
#include <vector>
using namespace Eigen;
using namespace std;

#define MAXRAWLEN 40960		 // 最大读取数据长度
#define MAXNUM 8			 // 波段数

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*解算结果*/
#define UN_Solve 0
#define Success_Solve 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

/*卫星数据*/
struct Satellate
{
	unsigned short PRN = 0;	 // 标识码
	unsigned short SYS = -1; // 系统
	int Phase_NUM = 0;		 // 波段数

	/*波段*/
	int SYG_TYPE[MAXNUM] = { -1, -1, -1, -1, -1, -1, -1, -1 };
	/*伪距(m)*/
	double PSERA[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*相位(cycle)*/
	double PHASE[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*多普勒(HZ)*/
	double DOPPLER[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*载噪比*/
	double SNR[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*伪距精度*/
	double PSE_PREC[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*载波相位精度*/
	double PHA_PREC[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*通道状态*/
	int LOCK_PSE[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	int LOCK_PHA[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	int PARITY[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };

	/*是否有粗差*/
	bool Outlier = false;
	bool Cycle_jump = false;

	static bool compare(const Satellate* a1, const Satellate* a2) {
		return a1->PRN < a2->PRN;
	}

	static bool ELEV_comp(const Satellate* a1, const Satellate* a2)
	{
		return a1->ELEV > a2->ELEV;
	}

	double ELEV = 0;
	XYZ pos;
	XYZ vel;
	double clk;
	double clk_vel;
};

struct EPHEMERIS
{
	/* data */
	// 卫星编号
	int PRN = 0;
	// 基本轨道参数
	double toe_tow;
	double toe_wn;
	double toc;
	double sqrt_A;
	double e;
	double M0;
	double omiga;
	double i0;
	double Omiga0;
	// 摄动改正数
	double delt_n;
	double dot_i;
	double dot_Omiga;
	double Cus;
	double Cuc;
	double Crs;
	double Crc;
	double Cis;
	double Cic;
	// 时间误差参数
	/*钟飘*/
	double a_f0;
	double a_f1;
	double a_f2;
	/*群延差*/
	double T_GD1;
	double T_GD2;
	/*数据龄期*/
	double AODC;
	int AODE;
	int IODC;
	// 健康指数
	unsigned int health;
	double URA;
	GPSTIME* Z_T;
	int IODE1;
	int IODE2;
};

struct OBS_DATA
{
	GPSTIME* OBS_TIME = new GPSTIME();
	int Sate_Num = 0;
	vector<Satellate*> GPS_SATE;
	vector<Satellate*> BDS_SATE;
};

/*解算结果存储类*/
class DATA_SET
{
public:
	GPSTIME* OBSTIME;	// 观测时间

	int LS_GPS_num;		// GPS卫星数
	int LS_BDS_num;		// BDS卫星数
	int KF_GPS_num;		// GPS卫星数
	int KF_BDS_num;		// BDS卫星数
	string* LS_SATES;		// 使用的卫星集合
	string* KF_SATES;
	int LS_result;	// 解算结果
	int KF_result;
	XYZ* Real_Pos;		// 参考真值

	/*卫星粗差探测存储变量*/
	double GPS_GF[GPS_SAT_QUAN];
	double GPS_MW[GPS_SAT_QUAN];
	double GPS_PSE[6][GPS_SAT_QUAN];
	double GPS_PHA[6][GPS_SAT_QUAN];
	double GPS_DOP[6][GPS_SAT_QUAN];
	int GPS_COUNT[GPS_SAT_QUAN];

	double BDS_GF[BDS_SAT_QUAN];
	double BDS_MW[BDS_SAT_QUAN];
	double BDS_PSE[5][BDS_SAT_QUAN];
	double BDS_PHA[5][BDS_SAT_QUAN];
	double BDS_DOP[5][BDS_SAT_QUAN];
	int BDS_COUNT[BDS_SAT_QUAN];

	/*观测值数据*/
	OBS_DATA* range;

	/*单历元星历数据存储*/
	static EPHEMERIS* GPS_eph[GPS_SAT_QUAN];
	static EPHEMERIS* BDS_eph[BDS_SAT_QUAN];

	/*求解过程参量*/
	bool LS_first;
	bool KF_first;
	MatrixXd Pos;
	MatrixXd temp_ref;

	/*Least_Square*/
	Least_Squares* LS_Pos;
	Least_Squares* LS_Vel;

	/*KalmanFilter*/
	KalmanFilter* KF;

	DATA_SET(GNSS_Configure cfg);

	static void inital_eph();

	void reset();

	int LS_print(GNSS_Configure cfg);				// 控制台输出
	int LS_Filewrite(FILE* fpr, GNSS_Configure cfg); // 文件输出
	void KF_Print(FILE* fpr, GNSS_Configure cfg);

	int Set_KF(GNSS_Configure cfg);
	int Set_LS(GNSS_Configure cfg);

	/*一致性检验*/
	int CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag);

	/*粗差探测*/
	int DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2);

	void DetectOut(GNSS_Configure cfg, double dt_e);

	void Sate_pos_pre(double t, GNSS_Configure cfg);

};

class RTK_DATA
{
public:
	GPSTIME* OBSTIME;
	int LS_GPS_num;		// GPS卫星数
	int LS_BDS_num;		// BDS卫星数
	int KF_GPS_num;		// GPS卫星数
	int KF_BDS_num;		// BDS卫星数
	string* LS_SATES;
	string* KF_SATES;
	vector<int>GPS_PRN;
	vector<int>BDS_PRN;
	int Ref_PRN_GPS_OLD;
	int Ref_PRN_BDS_OLD;
	int Ref_PRN_GPS_NEW;
	int Ref_PRN_BDS_NEW;
	int LS_result;	// 解算结果
	int KF_result;

	bool LS_IS_FIXED;
	bool KF_IS_FIXED;

	double Ratio_LS;
	double Ratio_KF;

	double GPS_GF[GPS_SAT_QUAN];
	double BDS_GF[BDS_SAT_QUAN];

	vector<Satellate*> Rove_GPS;
	vector<Satellate*> Rove_BDS;
	vector<Satellate*> Base_GPS;
	vector<Satellate*> Base_BDS;
	XYZ* Base_appro_pos;
	XYZ* Rove_appro_pos;

	XYZ* KF_pos;
	MatrixXd Q_baseLine;
	MatrixXd Fix_N;
	int Stable_Num;

	Least_Squares* LS;
	KalmanFilter* KF;

	RTK_DATA(GNSS_Configure cfg);
	void reset();
	void LS_Output(FILE* fpr, GNSS_Configure cfg);
	void KF_Output(FILE* fpr, GNSS_Configure cfg);
};

/*创建文件夹*/
int createDirectory(string path);

