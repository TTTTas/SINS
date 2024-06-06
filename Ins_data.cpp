#include "Ins_data.h"
#include "comm.h"
#include <fstream>
#include "progressbar.hpp"

#define _CRT_SECURE_NO_WARNINGS


void INS_Eigen::Init_Yaw()
{
	int total_time = Samp_rate * Init_time;
	double mean_x, mean_y, mean_z;
    double meanfx = 0.0, meanfy = 0.0, meanfz = 0.0;
    mean_x = mean_y = mean_z = 0.0;
	for(int i=0;i<total_time;i++)
	{
		mean_x += Imu[i]->gyro_x;
		mean_y += Imu[i]->gyro_y;
		mean_z += Imu[i]->gyro_z;
        meanfx += Imu[i]->accel_x;
        meanfy += Imu[i]->accel_y;
        meanfz += Imu[i]->accel_z;
	}
	mean_x /= total_time;
	mean_y /= total_time;
	mean_z /= total_time;
    meanfx /= total_time;
    meanfy /= total_time;
    meanfz /= total_time;

    Eigen::Vector3d g_n(0, 0, gravity);
    Eigen::Vector3d v_g = g_n.normalized();
    Eigen::Vector3d omiga_n_ie(OMEGA_E * cos(Atti(0)),0,-OMEGA_E * sin(Atti(0)));
    Eigen::Vector3d v_omiga = (g_n.cross(omiga_n_ie)).normalized();
    Eigen::Vector3d v_gomiga = (g_n.cross(omiga_n_ie).cross(g_n)).normalized();

    // 计算omiga_b_ie
    Eigen::Vector3d omiga_b_ie(mean_x, mean_y, mean_z);

    // 计算g_b
    Eigen::Vector3d g_b(meanfx, meanfy, meanfz);
    g_b = -g_b;

    // 计算omiga_g
    Eigen::Vector3d omiga_g = g_b.normalized();

    // 计算omiga_omiga
    Eigen::Vector3d omiga_omiga = (g_b.cross(omiga_b_ie)).normalized();

    // 计算omiga_gomiga
    Eigen::Vector3d omiga_gomiga = (g_b.cross(omiga_b_ie).cross(g_b)).normalized();

    // 构建C_n_b矩阵
    Eigen::Matrix3d V;
    Eigen::Matrix3d O;
    V.col(0) = v_g;
    V.col(1) = v_omiga;
    V.col(2) = v_gomiga;
    O.row(0) = omiga_g.transpose();
    O.row(1) = omiga_omiga.transpose();
    O.row(2) = omiga_gomiga.transpose();

    Eigen::Matrix3d C_n_b = V * O;

    // 计算俯仰角、横滚角和航向角
    Atti(0) = std::atan(-C_n_b(2, 0) / std::sqrt(C_n_b(2, 1) * C_n_b(2, 1) + C_n_b(2, 2) * C_n_b(2, 2)));
    Atti(1) = std::atan2(C_n_b(2, 1), C_n_b(2, 2));
    Atti(2) = std::atan2(C_n_b(1, 0), C_n_b(0, 0));


    std::cout << std::endl << "初始姿态角:\n" << rad2deg(Atti).transpose() << std::endl;
}

void INS_Eigen::INS_Mech()
{
    // Prepare result containers
    int data_size = Imu.size();
    Eigen::Matrix3d E = Euler2C(Atti);
    Eigen::Matrix3d e = Eigen::Matrix3d::Zero();
    Eigen::Vector3d v = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    Eigen::Vector3d v_last_two = Eigen::Vector3d::Zero();
    Eigen::Vector3d pos_last_two = Eigen::Vector3d::Zero();

    double time_last = Imu[Samp_rate * Init_time]->time;

    FILE* BLH_Fobs;
    try
    {
        if ((fopen_s(&BLH_Fobs,this->BLH_path.c_str(), "w")) != 0)printf("Failed to Open File!");
    }
    catch (std::exception exception)
    {
        printf("Failed to Open File!");
    }

    // 初始化进度条
    progressbar bar(data_size - Samp_rate * Init_time);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");
    int last_epoch = Samp_rate * Init_time;

    printf("惯导机械编排解算:\n");
    for (int i = Samp_rate * Init_time + 1; i < data_size; ++i) {
        Time = Imu[i]->time;
        double delt_t = Time - time_last;
        time_last = Time;
        // 解算
        e = Update_Euler_C(E, Imu[i - 1]->GYR_Vector(), Imu[i]->GYR_Vector(), BLH, Vel, delt_t);
        v = Update_velocity(v_last_two, Vel, pos_last_two, BLH, delt_t, Imu[i]->ACC_Vector(), Imu[i]->ACC_Vector(), E, Imu[i - 1]->GYR_Vector(), Imu[i]->GYR_Vector());
        pos = Update_pos(BLH, Vel, v, delt_t);

        // 更新
        v_last_two = Vel;
        pos_last_two = BLH;
        E = e;
        Vel = v;
        BLH = pos;

        std::fprintf(BLH_Fobs, "%.4f\t%10.6f\t%10.6f\t%7.6f\n", Time, rad2deg(BLH(0)), rad2deg(BLH(1)), BLH(2));

        if(i- last_epoch>0.005*(data_size- Samp_rate * Init_time))
        {
            last_epoch = i;
            bar.update(i - Samp_rate * Init_time);
        }
    }
    Atti = C2Euler(E);
}

