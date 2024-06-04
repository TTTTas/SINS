#pragma once
#define ACC_SCALE 1.5258789063E-06
#define GYR_SCALE 1.0850694444E-07
#include <Eigen/Dense>

class IMU_data
{
public:
	double time;
	double gyro_x;
	double gyro_y;
	double gyro_z;
	double accel_x;
	double accel_y;
	double accel_z;

	IMU_data()
	{
		time = 0;
		gyro_x = gyro_y = gyro_z = 0;
		accel_x = accel_y = accel_z = 0;
	}

	IMU_data(double t, double g_x,double g_y, double g_z, double a_x, double a_y, double a_z)
	{
		time = t;
		gyro_x = g_x;
		gyro_y = g_y;
		gyro_z = g_z;
		accel_x = a_x;
		accel_y = a_y;
		accel_z = a_z;
	}

public:
	Eigen::Vector3d GYR_Vector() { return Eigen::Vector3d(gyro_x, gyro_y, gyro_z); }
	Eigen::Vector3d ACC_Vector() { return Eigen::Vector3d(accel_x, accel_y, accel_z); }
};
