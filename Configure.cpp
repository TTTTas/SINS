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
	Imu_path = "IMU_data.ASC";
	GNSS_path = "GNSS_data.pos";
	ODO_path = "ODO.bin";
	Out_Folder="D:\\INS\\result\\";

	use_GNSS_vel = false;
	use_ODONHC = false;
	use_ZUPT = false;
	odo_update_rate = 0.5;
	zupt_update_rate = 0.5;
	use_GNSS_file = true;

	Init_time = 5 * 60;
	init_pos = XYZ(deg2rad(30.53009), deg2rad(114.35609), 26.704);
	init_vel = XYZ(0.0, 0.0, 0.0);
	init_att = XYZ(0.0, 0.0, 0.0);

	init_pos_std = XYZ(0.05, 0.05, 0.1);
	init_vel_std = XYZ(0.05, 0.05, 0.05);
	init_att_std = XYZ(0.1, 0.1, 0.5);

	init_gyr_bias = XYZ(-3700, 3400, 1000);
	init_acc_bias = XYZ(-10000, 3500, -6700);
	init_gyr_scale = XYZ(10000, 1700, -500);
	init_acc_scale = XYZ(0, 2000, 0);
	init_gyr_bias_std = XYZ(50, 50, 50);
	init_acc_bias_std = XYZ(250, 250, 250);
	init_gyr_scale_std = XYZ(1000, 1000, 1000);
	init_acc_scale_std = XYZ(1000, 1000, 1000);

	gyr_ARW = 0.24;
	acc_VRW = 0.24;
	gyr_bias_std = 50;
	acc_bias_std = 250;
	gyr_scale_std = 1000;
	acc_scale_std = 1000;
	corr_time = 1;

	ant_lever = XYZ(0.045, 0.46,-0.238);
	odo_lever = XYZ(-0.522, -0.47, 1.797);
	install_angle = XYZ(0, -0.2, 1.2);

	ODONHC_mean_noise = XYZ(0.1, 0.1, 0.1);

	init_att *= DEG2RAD;
	init_att_std *= DEG2RAD;

	init_gyr_bias *= DEG2RAD / 3600;
	init_acc_bias *= 1e-5;
	init_gyr_scale *= 1e-6;
	init_acc_scale *= 1e-6;
	init_gyr_bias_std *= DEG2RAD / 3600;
	init_acc_bias_std *= 1e-5;
	init_gyr_scale_std *= 1e-6;
	init_acc_scale_std *= 1e-6;

	gyr_ARW *= DEG2RAD / 60;
	acc_VRW /= 60;
	gyr_bias_std *= DEG2RAD / 3600;
	acc_bias_std *= 1e-5;
	gyr_scale_std *= 1e-6;
	acc_scale_std *= 1e-6;
	corr_time *= 3600;

	install_angle *= DEG2RAD;
}
