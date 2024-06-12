#include"Configure.h"
#include<iostream>

#include "comm.h"

GNSS_Configure::GNSS_Configure()
{
	Load_cfg();
	int count = 0;
	if (GPS_Cfg.used)
		count++;
	if (BDS_Cfg.used)
		count++;
	SYS_num = count;
	std::cout << "SYS_num: " << SYS_num << "\n"
		<< "Phase_num: " << phase_num << "\n"
		<< std::endl;
}

void GNSS_Configure::Load_cfg()
{
	NetIP_1 = "8.140.46.126";
	NetPort_1 = 4002;

	NetIP_2 = "8.140.46.126";
	NetPort_2 = 4002;

	//NetIP_1 = "123.57.30.145";
	//NetPort_1 = 4003;

	//NetIP_2 = "123.57.30.145";
	//NetPort_2 = 4003;

	phase_num = 1;
	SYS_num = 0;
	Hop_used = 1;

	GPS_Cfg = Sate_Configure(Freq_L1, Freq_L2);
	BDS_Cfg = Sate_Configure(Freq_B3, Freq_B1);
	GPS_Cfg.used = true;
	BDS_Cfg.used = true;

	SPP_LS_used = true;
	SPP_KF_used = false;
	RTK_LS_used = true;
	RTK_KF_used = true;

	w_thresh = 10;
	Ele_Mask = 15;
	SNR_thresh = 30;
	Ratio_thresh = 2;
	GF_thresh = 0.0131;

	LrateP = 10000;

	initial_covariance = 30;
	process_noise_pos = 1;
	process_noise_abg = 1e-4;
}

Sate_Configure::Sate_Configure(double f1_, double f2_)
{
	used = false;
	f1 = f1_;
	f2 = f2_;
}

Sate_Configure::Sate_Configure()
{
	used = false;
	f1 = 0;
	f2 = 0;
}

INS_Configure::INS_Configure()
{
	Imu_path = "./test_data/IMU_data.ASC";
	GNSS_path = "./test_data/RTK.nav";
	ODO_path = "ODO.bin";
	Out_Folder="./test_data/result";

	use_GNSS_vel = false;
	use_ODONHC = false;
	use_ZUPT = false;
	odo_update_rate = 0.5;
	zupt_update_rate = 0.5;
	use_GNSS_file = true;

	Init_time = 5 * 60;
	Samp_rate = 100;

	start_time = 355467.76;
	end_time = -1;

	gins_options.initstate.pos << -2267683.401, 5009281.160, 3221140.845;
	gins_options.initstate.vel << 0.132, -0.318, 0.600;
	gins_options.initstate.euler << 0.971825, 0.148966, 2.762344;

	gins_options.initstate_std.pos << 1.00000000e+2, 1.00000000e+2, 1.00000000e+2;
	gins_options.initstate_std.vel << 1.00000000e+0, 1.00000000e+0, 1.00000000e+0;
	gins_options.initstate_std.euler << 5.00000000e-1, 5.00000000e-1, 5.00000000e-1;

	gins_options.initstate.imuerror.gyrbias << 0, 0, 0;
	gins_options.initstate.imuerror.accbias << 0, 0, 0;
	gins_options.initstate_std.imuerror.gyrbias << 8.33300000e-5, 8.33300000e-5, 8.33300000e-5;
	gins_options.initstate_std.imuerror.accbias << 9.80000000e-4, 9.80000000e-4, 9.80000000e-4;

	gins_options.imunoise.pos_prw << 1.00000000e-4, 1.00000000e-4, 1.00000000e-4;
	gins_options.imunoise.acc_vrw << 1.00000000e-3, 1.00000000e-3, 1.00000000e-3;
	gins_options.imunoise.gyr_arw << 2.75000000e-3, 2.75000000e-3, 2.75000000e-3;
	gins_options.imunoise.gyrbias_std << 7.71600069e-10, 7.71600069e-10, 7.71600069e-10;
	gins_options.imunoise.accbias_std << 6.00000000e-8, 6.00000000e-8, 6.00000000e-8;
	gins_options.imunoise.corr_time = 4;

	gins_options.antlever << -0.0450, 0.0050, 0.8880;
	odo_lever << -0.0450, 0.0050, 0.8880;

	ODONHC_mean_noise << 0.1, 0.1, 0.1;

	gins_options.initstate.euler *= DEG2RAD;
	gins_options.initstate_std.euler *= DEG2RAD;

	gins_options.initstate.imuerror.gyrbias *= DEG2RAD;
	gins_options.initstate.imuerror.accbias *= 1;
	gins_options.imunoise.gyrbias_std *= DEG2RAD;
	gins_options.imunoise.accbias_std *= 1;

	gins_options.imunoise.gyr_arw *= DEG2RAD;

	gins_options.imunoise.corr_time *= 3600;
}
