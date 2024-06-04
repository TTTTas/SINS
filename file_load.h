#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

#include "Ins_data.h"

IMU_data* read_line_data(const std::string& line);

int read_imu_asc(INS_Eigen& ins);