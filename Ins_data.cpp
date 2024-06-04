#include "Ins_data.h"
#include "comm.h"
#include <fstream>

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

    Vector3d g_n(0, 0, gravity);
    Vector3d v_g = g_n.normalized();
    Vector3d omiga_n_ie(OMEGA_E * cos(Atti(0)),0,-OMEGA_E * sin(Atti(0)));
    Vector3d v_omiga = (g_n.cross(omiga_n_ie)).normalized();
    Vector3d v_gomiga = (g_n.cross(omiga_n_ie).cross(g_n)).normalized();

    // 计算omiga_b_ie
    Vector3d omiga_b_ie(mean_x, mean_y, mean_z);

    // 计算g_b
    Vector3d g_b(meanfx, meanfy, meanfz);
    g_b = -g_b;

    // 计算omiga_g
    Vector3d omiga_g = g_b.normalized();

    // 计算omiga_omiga
    Vector3d omiga_omiga = (g_b.cross(omiga_b_ie)).normalized();

    // 计算omiga_gomiga
    Vector3d omiga_gomiga = (g_b.cross(omiga_b_ie).cross(g_b)).normalized();

    // 构建C_n_b矩阵
    Matrix3d V;
    Matrix3d O;
    V.col(0) = v_g;
    V.col(1) = v_omiga;
    V.col(2) = v_gomiga;
    O.row(0) = omiga_g.transpose();
    O.row(1) = omiga_omiga.transpose();
    O.row(2) = omiga_gomiga.transpose();

    Matrix3d C_n_b = V * O;

    // 计算俯仰角、横滚角和航向角
    Atti(0) = std::atan(-C_n_b(2, 0) / std::sqrt(C_n_b(2, 1) * C_n_b(2, 1) + C_n_b(2, 2) * C_n_b(2, 2)));
    Atti(1) = std::atan2(C_n_b(2, 1), C_n_b(2, 2));
    Atti(2) = std::atan2(C_n_b(1, 0), C_n_b(0, 0));

    cout << rad2deg(Atti).transpose() << endl;
}

void INS_Eigen::Pure_IMU()
{
    // Prepare result containers
    int data_size = Imu.size();
    Matrix3d E = Euler2C(Atti);
    Matrix3d e = Matrix3d::Zero();
    Vector3d v = Vector3d::Zero();
    Vector3d pos = Vector3d::Zero();
    Vector3d v_last_two = Vector3d::Zero();
    Vector3d pos_last_two = Vector3d::Zero();

    double time_last = Imu[Samp_rate * Init_time]->time;

    fstream sw;
    try
    {
        sw.open(this->BLH_path);
    }
    catch (exception exception)
    {
        printf("Failed to Open File!");
    }

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

        std::cout << Time << "\t" << BLH.transpose() << "\n";
        sw << Time << "\t" << BLH.transpose() << "\n";
    }
    Atti = C2Euler(E);
}

