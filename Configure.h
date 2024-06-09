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

/*网口和存储结构相关配置*/
class GNSS_Configure
{
public:
	const char* NetIP_1;			 // IP
	unsigned short NetPort_1; // 端口
	const char* NetIP_2;			 // IP
	unsigned short NetPort_2; // 端口
	const char* ObsDatFile_1;				 // log文件路径
	const char* ObsDatFile_2;				 // log文件路径
	const char* ResDatFile;				 // pos文件路径
	const char* KFDatFile;
	/*频点数、系统数配置*/
	int phase_num; // 单频or双频
	int SYS_num;   // 单星or双星
	int Hop_used;  // 是否启用对流层改正

	/*GPS系统配置*/
	Sate_Configure GPS_Cfg;

	/*BDS系统配置*/
	Sate_Configure BDS_Cfg;

	/*解算模式*/
	bool SPP_LS_used;
	bool SPP_KF_used;
	bool RTK_LS_used;
	bool RTK_KF_used;
	/*阈值*/
	double w_thresh;
	double Ele_Mask;
	double SNR_thresh;
	double Ratio_thresh;
	double GF_thresh;

	double LrateP;

	/*卡尔曼滤波配置*/
	double initial_covariance;
	double process_noise_pos;
	double process_noise_abg;
	GNSS_Configure();
	void Load_cfg();
};

class INS_Configure
{
public:
	// 文件路径
	const char* Imu_path;
	const char* GNSS_path;
	const char* ODO_path;
	string Out_Folder;
	// 解算配置
	bool use_GNSS_vel;
	bool use_ZUPT;
	bool use_ODONHC;
	bool use_GNSS_file;
	double odo_update_rate;
	double zupt_update_rate;
	// 初始化信息
	int Init_time;
	int Samp_rate;
	double start_time;
	double end_time;
	GINSOptions gins_options;
	// 安装参数
	Vector3d odo_lever;			//m
	// 噪声
	Vector3d ODONHC_mean_noise;	//m/s

	INS_Configure();
};