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

// Function to calculate the local gravity value
 double GRS80_g(const Eigen::Vector3d& pos);

// Function to calculate the meridian radius
 double Cal_RM(double B);

// Function to calculate the prime vertical radius
 double Cal_RN(double B);

// Function to create a skew-symmetric matrix
Eigen::Matrix3d Cross_vector(const Eigen::Vector3d& vector);

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

// Function to update Euler angles to quaternion
Eigen::Quaterniond Update_Euler_q(const Eigen::Vector3d& E, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& pos, const Eigen::Vector3d& v, double dt);

// Function to update the position
Eigen::Vector3d Update_pos(const Eigen::Vector3d& pos0, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double dt);

// Function to convert BLH to NE coordinates
std::vector<Eigen::Vector2d> BLH2NE(const std::vector<Eigen::Vector3d>& BLH, const Eigen::Vector3d& BLH0);

// Function to update attitude matrix
Eigen::Matrix3d Update_Euler_C(const Eigen::Matrix3d& E, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& pos, const Eigen::Vector3d& v, double dt);

// Function to calculate gravity/Coriolis integral term
 Eigen::Vector3d Cal_deltv_g_cor(const Eigen::Vector3d g, const Eigen::Vector3d& omiga_ie, const Eigen::Vector3d& omiga_en, const Eigen::Vector3d& v0, double dt);

// Function to calculate specific force integral term
 Eigen::Vector3d Cal_deltv_f(const Eigen::Vector3d& omiga_ie, const Eigen::Vector3d& omiga_en, const Eigen::Matrix3d& C, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double dt);

// Function to update velocity
Eigen::Vector3d Update_velocity(const Eigen::Vector3d& v00, const Eigen::Vector3d& v0, const Eigen::Vector3d& pos, const Eigen::Vector3d& pos0, double dt, const Eigen::Vector3d& dv0, const Eigen::Vector3d& dv1, const Eigen::Matrix3d& C, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1);