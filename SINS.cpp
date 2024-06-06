#include <iostream>
#include <Eigen/Dense>

#include "comm.h"
#include "GNSS_data.h"
#include "file_load.h"
#include "Ins_data.h"


int main()
{
	/*INS_Eigen ins_eigen;
	read_imu_asc(ins_eigen);
	ins_eigen.Init_Yaw();
	ins_eigen.Pure_IMU();*/

	GNSS_Configure cfg;
	DATA_SET* data = new DATA_SET(cfg);
	read_EPHEMERIS_Rnx(data->GPS_eph, data->BDS_eph, "nav.24p");


	system("pause");

	return 0;
}