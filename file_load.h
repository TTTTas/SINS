#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <unordered_map>

#include "GNSS_data.h"
#include "Ins_data.h"

IMU_data* read_line_data(const std::string& line);

int read_imu_asc(INS_Eigen& ins);

bool isAllSpaces(const std::string& str);

Satellate* read_OBS_line(string line, vector<pair<string, int>> pairs, vector<string>signals);

int read_OBS_Rnx(vector<OBS_DATA*>* obs, const char* path, XYZ* appro_pos);

int read_EPHEMERIS_Rnx(EPHEMERIS** gps_eph, EPHEMERIS** bds_eph, const char* path);

int read_EPHEMERIS_One(EPHEMERIS* eph, FILE& fpr, int prn);