#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

#define M_PI 3.14159265358979
constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180.0 / M_PI;
constexpr double OMEGA_E = 7.292115e-5;

// Function to convert degrees to radians
 double deg2rad(double degrees);

Eigen::Vector3d deg2rad(Eigen::Vector3d& degrees);

// Function to convert radians to degrees
 double rad2deg(double radians);

Eigen::Vector3d rad2deg(Eigen::Vector3d& radians);

// 线性外推函数
Eigen::Vector3d Extrapol(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1);

// 向量元素平方
Eigen::Vector3d SQR_Mat(Eigen::Vector3d v);

// Function to calculate the local gravity value
 double GRS80_g(const Eigen::Vector3d& pos);

// Function to calculate the meridian radius
 double Cal_RM(double B);

// Function to calculate the prime vertical radius
 double Cal_RN(double B);

// Function to create a skew-symmetric matrix
Eigen::Matrix3d Skew(const Eigen::Vector3d& vector);

// Function to convert rotation matrix to quaternion
Eigen::Quaterniond C2q(const Eigen::Matrix3d& C);

// Function to convert quaternion to rotation matrix
Eigen::Matrix3d q2C(const Eigen::Quaterniond& q);

// Function to convert quaternion to equivalent rotation vector
Eigen::Vector3d q2Phi(const Eigen::Quaterniond& q);

// Function to convert equivalent rotation vector to rotation matrix
Eigen::Matrix3d Phi2C(const Eigen::Vector3d& Phi);

// Function to convert equivalent rotation vector to quaternion
Eigen::Quaterniond Phi2q(const Eigen::Vector3d& Phi);

// Function to convert a rotation matrix to Euler angles
Eigen::Vector3d C2Euler(const Eigen::Matrix3d& C);

// Function to convert Euler angles to a rotation matrix
Eigen::Matrix3d Euler2C(const Eigen::Vector3d& Euler);

// Function to convert Euler angles to quaternion
Eigen::Quaterniond Euler2q(const Eigen::Vector3d& Euler);

/* n系(导航坐标系)到e系(地心地固坐标系)转换四元数 */
Eigen::Quaterniond q_n2e(const Eigen::Vector3d& blh);

/* 从n系到e系转换四元数得到纬度和经度 */
Eigen::Vector3d q_n2e_2_blh(const Eigen::Quaterniond& qne, double height);

/* n系相对位置转大地坐标相对位置 */
Eigen::Matrix3d DRi(const Eigen::Vector3d& blh);

/* 大地坐标相对位置转n系相对位置 */
Eigen::Matrix3d DR(const Eigen::Vector3d& blh);

// Function to update Euler angles to quaternion
Eigen::Quaterniond Update_Euler_q(const Eigen::Vector3d& E, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& pos, const Eigen::Vector3d& v, double dt);

// Function to convert BLH to NE 
std::vector<Eigen::Vector2d> BLH2NE(const std::vector<Eigen::Vector3d>& BLH, const Eigen::Vector3d& BLH0);