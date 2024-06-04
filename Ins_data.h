#pragma once
#include <vector>
#include <string>
#include "imu_data.h"
#include <cmath>

#include "comm.h"
using namespace Eigen;

class INS_Eigen
{
private:
	double gravity;
	double omega_e;
	double Samp_rate;
	int Init_time;

	double Time;
	Vector3d Atti;
	Vector3d Vel;
	Vector3d Pos;			//ECEF    m   m   m
	Vector3d BLH;			//CGCS200 rad rad m

public:
	std::vector<IMU_data*> Imu;


	std::string IMU_file_path;
	std::string Pos_path;
	std::string BLH_path;
	INS_Eigen()
	{
		IMU_file_path = "IMU_data.ASC";
		Pos_path = "ECEF_pos.pos";
		BLH_path = "CGCS2000_pos.pos";

		Samp_rate = 100;
		gravity = 9.7936174;				//m/(s^2)
		omega_e = 7.292115 * 1e-5;			//rad/s
		Init_time = 5 * 60;					//s

		Time = 0;
		Atti = Vector3d::Zero();
		Vel = Vector3d::Zero();
		Pos = Vector3d::Zero();
		BLH = Vector3d(deg2rad(30.53009), deg2rad(114.35609), 26.704);
	}

	void Init_Yaw();

	void Pure_IMU();
};
