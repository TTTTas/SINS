#include <iostream>
#include <Eigen/Dense>

#include "comm.h"
#include "file_load.h"
#include "Ins_data.h"


int main()
{
	INS_Eigen ins_eigen;
	read_imu_asc(ins_eigen);
	ins_eigen.Init_Yaw();
	ins_eigen.Pure_IMU();


	system("pause");

	return 0;
}