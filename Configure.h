#pragma once

/*FREQUNCY*/ // MHz
/*GPS*/
#define Freq_L1 1575.42
#define Freq_L2 1227.60
#define Freq_L5 1176.45
/*BDS*/
#define Freq_B1 1561.098
#define Freq_B1_C 1575.42
#define Freq_B2 1207.14
#define Freq_B2_a 1176.45
#define Freq_B3 1268.52
#include "INS_types.h"
#include "transform.h"

class Sate_Configure
{
public:
	bool used;

	double f1;
	double f2;

	Sate_Configure();
	Sate_Configure(double f1, double f2);
};

/*���ںʹ洢�ṹ�������*/
class GNSS_Configure
{
public:
	const char* NetIP_1;			 // IP
	unsigned short NetPort_1; // �˿�
	const char* NetIP_2;			 // IP
	unsigned short NetPort_2; // �˿�
	const char* ObsDatFile_1;				 // log�ļ�·��
	const char* ObsDatFile_2;				 // log�ļ�·��
	const char* ResDatFile;				 // pos�ļ�·��
	const char* KFDatFile;
	/*Ƶ������ϵͳ������*/
	int phase_num; // ��Ƶor˫Ƶ
	int SYS_num;   // ����or˫��
	int Hop_used;  // �Ƿ����ö��������

	/*GPSϵͳ����*/
	Sate_Configure GPS_Cfg;

	/*BDSϵͳ����*/
	Sate_Configure BDS_Cfg;

	/*����ģʽ*/
	bool SPP_LS_used;
	bool SPP_KF_used;
	bool RTK_LS_used;
	bool RTK_KF_used;
	/*��ֵ*/
	double w_thresh;
	double Ele_Mask;
	double SNR_thresh;
	double Ratio_thresh;
	double GF_thresh;

	double LrateP;

	/*�������˲�����*/
	double initial_covariance;
	double process_noise_pos;
	double process_noise_abg;
	GNSS_Configure();
	void Load_cfg();
};

class INS_Configure
{
public:
	// �ļ�·��
	const char* Imu_path;
	const char* GNSS_path;
	const char* ODO_path;
	string Out_Folder;
	// ��������
	bool use_GNSS_vel;
	bool use_ZUPT;
	bool use_ODONHC;
	bool use_GNSS_file;
	double odo_update_rate;
	double zupt_update_rate;
	// ��ʼ����Ϣ
	int Init_time;
	int Samp_rate;
	double start_time;
	double end_time;
	GINSOptions gins_options;
	// ��װ����
	Vector3d odo_lever;			//m
	// ����
	Vector3d ODONHC_mean_noise;	//m/s

	INS_Configure();
};