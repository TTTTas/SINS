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
	Imu_path = "Leador-A15.txt";
	GNSS_path = "GNSS-RTK.txt";
	ODO_path = "ODO.bin";
	Out_Folder="./Ins/result";

	use_GNSS_vel = false;
	use_ODONHC = false;
	use_ZUPT = false;
	odo_update_rate = 0.5;
	zupt_update_rate = 0.5;
	use_GNSS_file = true;

	Init_time = 5 * 60;
	Samp_rate = 200;

	start_time = 456300;
	end_time = 457000;

	gins_options.initstate.pos << deg2rad(30.4447873701), deg2rad(114.4718632047), 20.899;
	gins_options.initstate.vel << 0.0, 0.0, 0.0;
	gins_options.initstate.euler << 0.85421502, -2.03480295, 185.70235133;

	gins_options.initstate_std.pos << 0.005, 0.004, 0.008;
	gins_options.initstate_std.vel << 0.003, 0.004, 0.004;
	gins_options.initstate_std.euler << 0.003, 0.003, 0.023;

	gins_options.initstate.imuerror.gyrbias << 0, 0, 0;
	gins_options.initstate.imuerror.accbias << 0, 0, 0;
	gins_options.initstate.imuerror.gyrscale << 0, 0, 0;
	gins_options.initstate.imuerror.accscale << 0, 0, 0;
	gins_options.imunoise.gyrbias_std << 0.027, 0.027, 0.027;
	gins_options.imunoise.accbias_std << 15.0, 15.0, 15.0;
	gins_options.imunoise.gyrscale_std << 300.0, 300.0, 300.0;
	gins_options.imunoise.accscale_std << 300.0, 300.0, 300.0;

	gins_options.imunoise.gyr_arw << 0.003, 0.003, 0.003;
	gins_options.imunoise.acc_vrw << 0.03, 0.03, 0.03;
	gins_options.imunoise.corr_time = 4;

	gins_options.antlever << 0.136, -0.301, -0.184;
	odo_lever << -0.522, -0.47, 1.797;

	ODONHC_mean_noise << 0.1, 0.1, 0.1;

	gins_options.initstate.euler *= DEG2RAD;
	gins_options.initstate_std.euler *= DEG2RAD;

	gins_options.initstate.imuerror.gyrbias *= DEG2RAD / 3600;
	gins_options.initstate.imuerror.accbias *= 1e-5;
	gins_options.initstate.imuerror.gyrscale *= 1e-6;
	gins_options.initstate.imuerror.accscale *= 1e-6;
	gins_options.imunoise.gyrbias_std *= DEG2RAD / 3600;
	gins_options.imunoise.accbias_std *= 1e-5;
	gins_options.imunoise.gyrscale_std *= 1e-6;
	gins_options.imunoise.accscale_std *= 1e-6;

	gins_options.imunoise.gyr_arw *= DEG2RAD / 60;
	gins_options.imunoise.acc_vrw /= 60;

	gins_options.imunoise.corr_time *= 3600;
}
