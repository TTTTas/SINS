#include "file_load.h"

IMU_data* read_line_data(const std::string& line)
{
    double g_x, g_y, g_z, a_x, a_y, a_z;
    g_x = g_y = g_z = a_x = a_y = a_z = 0.0;
    std::stringstream lineStream(line);
    std::string temp;

    std::string prefix;
    char comma, semicolon;
    int unknown;
    double timestamp1, timestamp2;


    // 忽略 "%RAWIMUSA,"
    lineStream.ignore(10);

    // 读取 prefix 和 timestamp1
    lineStream >> unknown >> comma >> timestamp1 >> semicolon;
    // 读取 prefix 和 timestamp2
    lineStream >> unknown >> comma >> timestamp2 >> comma;

    // 读取 unknown 和 a_x, a_y, a_z, g_x, g_y, g_z
    lineStream >> unknown >> comma >> a_z >> comma >> a_y >> comma >> a_x >> comma;
    lineStream >> g_z >> comma >> g_y >> comma >> g_x;

    IMU_data* imu = new IMU_data(timestamp1, -g_y * GYR_SCALE, g_x * GYR_SCALE, -g_z * GYR_SCALE, -a_y * ACC_SCALE, a_x * ACC_SCALE, -a_z * ACC_SCALE);
    return imu;
}

int read_imu_asc(INS_Eigen& ins)
{
    std::ifstream inputFile(ins.IMU_file_path);
    std::string line;
    int val = 0;
    if (!inputFile) {
        std::cerr << "无法打开文件!" << std::endl;
        return 1;
    }
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue;
        }
        else
        {
            val = 1;
            ins.Imu.push_back(read_line_data(line));
        }
    }
    return val;
}