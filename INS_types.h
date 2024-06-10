#pragma once
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include "comm.h"

typedef struct GNSS {
    double time;

    Eigen::Vector3d blh;
    Eigen::Vector3d std;
    Eigen::Vector3d vel;
    Eigen::Vector3d vel_std;

    bool isvalid;
} GNSS;

typedef struct IMU {
    double time;
    double dt;

    Eigen::Vector3d dtheta;
    Eigen::Vector3d dvel;

    double odovel;
} IMU;

typedef struct Attitude {
    Eigen::Quaterniond qbe;
    Eigen::Matrix3d cbe;
    Eigen::Vector3d euler;
} Attitude;

typedef struct PVA {
    Eigen::Vector3d pos;
    Eigen::Vector3d vel;
    Attitude att;
} PVA;

typedef struct ImuError {
    Eigen::Vector3d gyrbias;
    Eigen::Vector3d accbias;
} ImuError;

typedef struct NavState {
    Eigen::Vector3d pos;
    Eigen::Vector3d vel;
    Eigen::Vector3d euler;

    ImuError imuerror;
} NavState;

typedef struct ImuNoise {
    Eigen::Vector3d gyr_arw;
    Eigen::Vector3d acc_vrw;
    Eigen::Vector3d gyrbias_std;
    Eigen::Vector3d accbias_std;
    double corr_time;
} ImuNoise;

typedef struct GINSOptions {

    // 初始状态和状态标准差
    NavState initstate;
    NavState initstate_std;

    // IMU噪声参数
    ImuNoise imunoise;

    // 安装参数
    Eigen::Vector3d antlever = { 0, 0, 0 };

    void print_options() {
        std::cout << "---------------SINS Options:---------------" << std::endl;

        // 打印初始状态
        std::cout << " - Initial State: " << std::endl;
        std::cout << '\t' << "- initial position: ";
        std::cout << std::setprecision(12) << initstate.pos[0] * RAD2DEG << "  ";
        std::cout << std::setprecision(12) << initstate.pos[1] * RAD2DEG << "  ";
        std::cout << std::setprecision(6) << initstate.pos[2] << " [deg, deg, m] " << std::endl;
        std::cout << '\t' << "- initial velocity: " << initstate.vel.transpose() << " [m/s] " << std::endl;
        std::cout << '\t' << "- initial attitude: " << initstate.euler.transpose() * RAD2DEG << " [deg] " << std::endl;
        std::cout << '\t' << "- initial gyrbias : " << initstate.imuerror.gyrbias.transpose() * RAD2DEG * 3600
            << " [deg/h] " << std::endl;
        std::cout << '\t' << "- initial accbias : " << initstate.imuerror.accbias.transpose() * 1e5 << " [mGal] "
            << std::endl;

        // 打印初始状态标准差
        std::cout << " - Initial State STD: " << std::endl;
        std::cout << '\t' << "- initial position std: " << initstate_std.pos.transpose() << " [m] " << std::endl;
        std::cout << '\t' << "- initial velocity std: " << initstate_std.vel.transpose() << " [m/s] " << std::endl;
        std::cout << '\t' << "- initial attitude std: " << initstate_std.euler.transpose() * RAD2DEG << " [deg] "
            << std::endl;
        std::cout << '\t' << "- initial gyrbias std: " << initstate_std.imuerror.gyrbias.transpose() * RAD2DEG * 3600
            << " [deg/h] " << std::endl;
        std::cout << '\t' << "- initial accbias std: " << initstate_std.imuerror.accbias.transpose() * 1e5 << " [mGal] "
            << std::endl;

        // 打印IMU噪声参数
        std::cout << " - IMU noise: " << std::endl;
        std::cout << '\t' << "- arw: " << imunoise.gyr_arw.transpose() * RAD2DEG * 60 << " [deg/sqrt(h)] " << std::endl;
        std::cout << '\t' << "- vrw: " << imunoise.acc_vrw.transpose() * 60 << " [m/s/sqrt(h)] " << std::endl;
        std::cout << '\t' << "- gyrbias  std: " << imunoise.gyrbias_std.transpose() * RAD2DEG * 3600 << " [deg/h] "
            << std::endl;
        std::cout << '\t' << "- accbias  std: " << imunoise.accbias_std.transpose() * 1e5 << " [mGal] " << std::endl;
        std::cout << '\t' << "- correlation time: " << imunoise.corr_time / 3600.0 << " [h] " << std::endl;

        // 打印GNSS天线杆臂
        std::cout << " - Antenna leverarm: " << antlever.transpose() << " [m] " << std::endl << std::endl;
    }

} GINSOptions;