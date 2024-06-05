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

#define MAXRAWLEN 40960		 // ����ȡ���ݳ���
#define MAXNUM 8			 // ������

#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*������*/
#define UN_Solve 0
#define Success_Solve 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

/*��������*/
struct Satellate
{
	unsigned short PRN = 0;	 // ��ʶ��
	unsigned short SYS = -1; // ϵͳ
	int Phase_NUM = 0;		 // ������

	/*����*/
	int SYG_TYPE[MAXNUM] = { -1, -1, -1, -1, -1, -1, -1, -1 };
	/*α��(m)*/
	double PSERA[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*��λ(cycle)*/
	double PHASE[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*������(HZ)*/
	double DOPPLER[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*�����*/
	double SNR[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*α�ྫ��*/
	double PSE_PREC[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*�ز���λ����*/
	double PHA_PREC[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	/*ͨ��״̬*/
	int LOCK_PSE[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	int LOCK_PHA[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };
	int PARITY[MAXNUM] = { 0, 0, 0, 0, 0, 0, 0, 0 };

	/*�Ƿ��дֲ�*/
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
	// ���Ǳ��
	int PRN = 0;
	// �����������
	double toe_tow;
	double toe_wn;
	double toc;
	double sqrt_A;
	double e;
	double M0;
	double omiga;
	double i0;
	double Omiga0;
	// �㶯������
	double delt_n;
	double dot_i;
	double dot_Omiga;
	double Cus;
	double Cuc;
	double Crs;
	double Crc;
	double Cis;
	double Cic;
	// ʱ��������
	/*��Ʈ*/
	double a_f0;
	double a_f1;
	double a_f2;
	/*Ⱥ�Ӳ�*/
	double T_GD1;
	double T_GD2;
	/*��������*/
	double AODC;
	int AODE;
	int IODC;
	// ����ָ��
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

/*�������洢��*/
class DATA_SET
{
public:
	GPSTIME* OBSTIME;	// �۲�ʱ��

	int LS_GPS_num;		// GPS������
	int LS_BDS_num;		// BDS������
	int KF_GPS_num;		// GPS������
	int KF_BDS_num;		// BDS������
	string* LS_SATES;		// ʹ�õ����Ǽ���
	string* KF_SATES;
	int LS_result;	// ������
	int KF_result;
	XYZ* Real_Pos;		// �ο���ֵ

	/*���Ǵֲ�̽��洢����*/
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

	/*�۲�ֵ����*/
	OBS_DATA* range;

	/*����Ԫ�������ݴ洢*/
	static EPHEMERIS* GPS_eph[GPS_SAT_QUAN];
	static EPHEMERIS* BDS_eph[BDS_SAT_QUAN];

	/*�����̲���*/
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

	int LS_print(GNSS_Configure cfg);				// ����̨���
	int LS_Filewrite(FILE* fpr, GNSS_Configure cfg); // �ļ����
	void KF_Print(FILE* fpr, GNSS_Configure cfg);

	int Set_KF(GNSS_Configure cfg);
	int Set_LS(GNSS_Configure cfg);

	/*һ���Լ���*/
	int CheckOBSConsist(Satellate* sate, int sys, double t, int index, bool& PSE_flag, bool& PHA_flag);

	/*�ֲ�̽��*/
	int DetectOutlier(Satellate* sate, int sys, double t, int index1, int index2);

	void DetectOut(GNSS_Configure cfg, double dt_e);

	void Sate_pos_pre(double t, GNSS_Configure cfg);

};

class RTK_DATA
{
public:
	GPSTIME* OBSTIME;
	int LS_GPS_num;		// GPS������
	int LS_BDS_num;		// BDS������
	int KF_GPS_num;		// GPS������
	int KF_BDS_num;		// BDS������
	string* LS_SATES;
	string* KF_SATES;
	vector<int>GPS_PRN;
	vector<int>BDS_PRN;
	int Ref_PRN_GPS_OLD;
	int Ref_PRN_BDS_OLD;
	int Ref_PRN_GPS_NEW;
	int Ref_PRN_BDS_NEW;
	int LS_result;	// ������
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

/*�����ļ���*/
int createDirectory(string path);

