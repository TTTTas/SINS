#pragma once
#include <string>
#include <vector>

#include "GNSS_data.h"
#include "Ins_data.h"
#include "INS_types.h"

#define ACC_SCALE 1.5258789063E-06
#define GYR_SCALE 1.0850694444E-07

IMU read_line_imu_txt(const std::string& line);

IMU read_line_imu(const std::string& line);

int read_imu_asc(vector<IMU*>& Imu, INS_Configure cfg);

vector<double> read_line_gnss(const std::string& line);

int read_GNSS(vector<GNSS*>& Gnss, INS_Configure cfg);

bool isAllSpaces(const std::string& str);

Satellate* read_OBS_line(string line, vector<pair<string, int>> pairs, vector<string>signals);

int read_OBS_Rnx(vector<OBS_DATA*>* obs, const char* path, XYZ* appro_pos);

int read_EPHEMERIS_Rnx(EPHEMERIS** gps_eph, EPHEMERIS** bds_eph, const char* path);