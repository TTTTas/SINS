#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

#define M_PI 3.14159265358979
constexpr double DEG2RAD = M_PI / 180.0;
constexpr double RAD2DEG = 180.0 / M_PI;
constexpr double OMEGA_E = 7.292115e-5;

// Function to convert degrees to radians
 double deg2rad(double degrees);

 Vector3d deg2rad(Vector3d& degrees);

// Function to convert radians to degrees
 double rad2deg(double radians);

 Vector3d rad2deg(Vector3d& radians);

// 线性外推函数
 Vector3d Extrapol(const Vector3d& v0, const Vector3d& v1);

// Function to calculate the local gravity value
 double GRS80_g(const Vector3d& pos);

// Function to calculate the meridian radius
 double Cal_RM(double B);

// Function to calculate the prime vertical radius
 double Cal_RN(double B);

// Function to create a skew-symmetric matrix
 Matrix3d Cross_vector(const Vector3d& vector);

// Function to convert rotation matrix to quaternion
 Quaterniond C2q(const Matrix3d& C);

// Function to convert quaternion to rotation matrix
 Matrix3d q2C(const Quaterniond& q);

// Function to convert quaternion to equivalent rotation vector
 Vector3d q2Phi(const Quaterniond& q);

// Function to convert equivalent rotation vector to rotation matrix
 Matrix3d Phi2C(const Vector3d& Phi);

// Function to convert equivalent rotation vector to quaternion
 Quaterniond Phi2q(const Vector3d& Phi);

// Function to convert a rotation matrix to Euler angles
 Vector3d C2Euler(const Matrix3d& C);

// Function to convert Euler angles to a rotation matrix
Matrix3d Euler2C(const Vector3d& Euler);

// Function to convert Euler angles to quaternion
Quaterniond Euler2q(const Vector3d& Euler);

// Function to update Euler angles to quaternion
 Quaterniond Update_Euler_q(const Vector3d& E, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& pos, const Vector3d& v, double dt);

// Function to update the position
Vector3d Update_pos(const Vector3d& pos0, const Vector3d& v0, const Vector3d& v1, double dt);

// Function to convert BLH to NE coordinates
 vector<Vector2d> BLH2NE(const vector<Vector3d>& BLH, const Vector3d& BLH0);

// Function to update attitude matrix
Matrix3d Update_Euler_C(const Matrix3d& E, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& pos, const Vector3d& v, double dt);

// Function to calculate gravity/Coriolis integral term
 Vector3d Cal_deltv_g_cor(const Vector3d g, const Vector3d& omiga_ie, const Vector3d& omiga_en, const Vector3d& v0, double dt);

// Function to calculate specific force integral term
 Vector3d Cal_deltv_f(const Vector3d& omiga_ie, const Vector3d& omiga_en, const Matrix3d& C, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& v0, const Vector3d& v1, double dt);

// Function to update velocity
Vector3d Update_velocity(const Vector3d& v00, const Vector3d& v0, const Vector3d& pos, const Vector3d& pos0, double dt, const Vector3d& dv0, const Vector3d& dv1, const Matrix3d& C, const Vector3d& theta0, const Vector3d& theta1);