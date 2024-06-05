#include"Configure.h"
#include<iostream>

Configure::Configure()
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

void Configure::Load_cfg()
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