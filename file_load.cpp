#include "file_load.h"
#include "progressbar.hpp"

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


    // ���� "%RAWIMUSA,"
    lineStream.ignore(10);

    // ��ȡ prefix �� timestamp1
    lineStream >> unknown >> comma >> timestamp1 >> semicolon;
    // ��ȡ prefix �� timestamp2
    lineStream >> unknown >> comma >> timestamp2 >> comma;

    // ��ȡ a_x, a_y, a_z, g_x, g_y, g_z
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
        std::cerr << "�޷����ļ�!" << std::endl;
        return 1;
    }

    // ��ȡ�ļ��ܴ�С���ֽڣ�
    inputFile.seekg(0, std::ios::end);
    std::streamsize totalSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // ��ʼ��������
    progressbar bar(totalSize);
    bar.set_todo_char(" ");
    bar.set_done_char("��");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");
    std::streamsize bytesRead = 0;

    std::printf("��ȡ�ļ�:\n");
    std::streamsize bytes_last = 0;
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue;
        }
        else
        {
            val = 1;
            ins.Imu.push_back(read_line_data(line));

            // �����Ѷ�ȡ�ֽ���
            bytesRead = inputFile.tellg();
            // ������ȰٷֱȲ����½�����
            if(bytesRead-bytes_last>0.005*totalSize)
            {
                bar.update(bytesRead);
                bytes_last = bytesRead;
            }
        }

    }
    return val;
}