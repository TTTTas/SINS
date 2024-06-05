#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

#include "GNSS_data.h"
#include "Ins_data.h"

IMU_data* read_line_data(const std::string& line);

int read_imu_asc(INS_Eigen& ins);

int read_OBS_Rnx(OBS_DATA* obs);

int read_EPHEMERIS_Rnx(EPHEMERIS* gps_eph, EPHEMERIS* bds_eph);