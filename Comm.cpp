#include "comm.h"

// Function to convert degrees to radians
 double deg2rad(double degrees) {
    return degrees * DEG2RAD;
}

 Eigen::Vector3d deg2rad(Eigen::Vector3d& degrees)
{
    return degrees * DEG2RAD;
}

// Function to convert radians to degrees
 double rad2deg(double radians) {
    return radians * RAD2DEG;
}

 Eigen::Vector3d rad2deg(Eigen::Vector3d& radians)
{
    return radians * RAD2DEG;
}


// 线性外推函数
 Eigen::Vector3d Extrapol(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1) {
    // 输入 v0     历元0数据
    //      v1     历元1数据
    // 输出 result 外推历元2数据
    Eigen::Vector3d result = 1.5 * v1 - 0.5 * v0;
    return result;
}

Eigen::Matrix3d SQR_Mat(Eigen::Vector3d v)
{
    Eigen::Matrix3d sqr;
    sqr = v * Eigen::MatrixXd::Identity(3, 3) * v.transpose();
    return sqr;
}

// Function to calculate the local gravity value
 double GRS80_g(const Eigen::Vector3d& pos) {
    double B = pos(0);
    double H = pos(2);
    double g0 = 9.7803267715 * (1 + 0.0052790414 * sin(B) * sin(B) + 0.0000232718 * sin(B) * sin(B) * sin(B) * sin(B));
    return g0 - (3.087691089e-6 - 4.397731e-9 * sin(B) * sin(B)) * H + 0.721e-12 * H * H;
}

Eigen::Vector3d g_e(Eigen::Vector3d pos)
{
    double g = GRS80_g(pos);
    Eigen::Vector3d g_e;
    g_e << cos(pos[1]) * cos(pos[0]) * (-g), sin(pos[1])* cos(pos[0])* (-g), sin(pos[0])* (-g);
    return g_e;
}

// Function to calculate the meridian radius
 double Cal_RM(double B) {
    const double a = 6378137.0;
    const double e = 0.08181919104;
    return a * (1 - e * e) / pow((1 - e * e * sin(B) * sin(B)), 1.5);
}

// Function to calculate the prime vertical radius
 double Cal_RN(double B) {
    const double a = 6378137.0;
    const double e = 0.08181919104;
    return a / sqrt(1 - e * e * sin(B) * sin(B));
}

// Function to create a skew-symmetric matrix
 Eigen::Matrix3d Skew(const Eigen::Vector3d& vector) {
    Eigen::Matrix3d cross_matrix;
    cross_matrix << 0, -vector(2), vector(1),
        vector(2), 0, -vector(0),
        -vector(1), vector(0), 0;
    return cross_matrix;
}


// Function to convert rotation matrix to quaternion
Eigen::Quaterniond C2q(const Eigen::Matrix3d& C) {
    double tr = C.trace();
    double p1, p2, p3, p4;
    p1 = 1 + tr;
    p2 = 1 + 2 * C(0, 0) - tr;
    p3 = 1 + 2 * C(1, 1) - tr;
    p4 = 1 + 2 * C(2, 2) - tr;
    double q1, q2, q3, q4;

    if (p1 >= p2 && p1 >= p3 && p1 >= p4)
    {
        /// p1 = max(p1, p2, p3, p4)
        q1 = 0.5 * sqrt(p1);
        q2 = (C(1, 2) - C(2, 1)) / (4 * q1);
        q3 = (C(2, 0) - C(0, 2)) / (4 * q1);
        q4 = (C(0, 1) - C(1, 0)) / (4 * q1);
    }
    else if (p2 >= p1 && p2 >= p3 && p2 >= p4)
    {
        /// p2 = max(p1, p2, p3, p4)
        q2 = 0.5 * sqrt(p2);
        q1 = (C(1, 2) - C(2, 1)) / (4 * q2);
        q3 = (C(0, 1) + C(1, 0)) / (4 * q2);
        q4 = (C(0, 2) + C(2, 0)) / (4 * q2);
    }
    else if (p3 >= p1 && p3 >= p2 && p3 >= p4)
    {
        /// p3 = max(p1, p2, p3, p4)
        q3 = 0.5 * sqrt(p3);
        q1 = (C(2, 0) - C(0, 2)) / (4 * q3);
        q2 = (C(0, 1) + C(1, 0)) / (4 * q3);
        q4 = (C(1, 2) + C(2, 1)) / (4 * q3);
    }
    else
    {
        /// p4 = max(p1, p2, p3, p4)
        q4 = 0.5 * sqrt(p4);
        q1 = (C(0, 1) - C(1, 0)) / (4 * q4);
        q2 = (C(2, 0) + C(0, 2)) / (4 * q4);
        q3 = (C(1, 2) + C(2, 1)) / (4 * q4);
    }

    if (q1 < 0)
    {
        q1 *= -1;
        q2 *= -1;
        q3 *= -1;
        q4 *= -1;
    }
    return Eigen::Quaterniond(q1, q2, q3, q4);
}

// Function to convert quaternion to rotation matrix
 Eigen::Matrix3d q2C(const Eigen::Quaterniond& q) {
     double q0, q1, q2, q3;
     q0 = q.w();
     q1 = q.x();
     q2 = q.y();
     q3 = q.z();

     Eigen::MatrixX3d M = Eigen::Matrix3d::Zero(3, 3);
     M(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
     M(0, 1) = 2 * (q1 * q2 + q0 * q3);
     M(0, 2) = 2 * (q1 * q3 - q0 * q2);
     M(1, 0) = 2 * (q1 * q2 - q0 * q3);
     M(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
     M(1, 2) = 2 * (q2 * q3 + q0 * q1);
     M(2, 0) = 2 * (q1 * q3 + q0 * q2);
     M(2, 1) = 2 * (q2 * q3 - q0 * q1);
     M(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
     return M;
}

// Function to convert quaternion to equivalent rotation vector
Eigen::Vector3d q2Phi(const Eigen::Quaterniond& q) {
     Eigen::AngleAxisd angleAxis(q);
    double angle = angleAxis.angle();
    Eigen::Vector3d axis = angleAxis.axis();
    return angle * axis;
}

// Function to convert equivalent rotation vector to rotation matrix
 Eigen::Matrix3d Phi2C(const Eigen::Vector3d& Phi) {
    double n_phi = Phi.norm();
    Eigen::Matrix3d C = Eigen::Matrix3d::Identity() + sin(n_phi) / n_phi * Skew(Phi) +
        (1 - cos(n_phi)) / (n_phi * n_phi) * Skew(Phi) * Skew(Phi);
    return C;
}

Eigen::Vector3d C2Phi(const Eigen::Matrix3d& C)
{
    return (q2Phi(C2q(C)));
}

// Function to convert equivalent rotation vector to quaternion
 Eigen::Quaterniond Phi2q(const Eigen::Vector3d& Phi) {
    double half_phi = Phi.norm();
    if (half_phi != 0) {
        return Eigen::Quaterniond(Eigen::AngleAxisd(half_phi, Phi.normalized()));
    }
    else {
        return Eigen::Quaterniond(1, 0, 0, 0);
    }
}

// Function to convert a rotation matrix to Euler angles
 Eigen::Vector3d C2Euler(const Eigen::Matrix3d& C) {
    Eigen::Vector3d Euler;
    if (abs(C(2, 0)) < 0.999) {
        Euler(1) = asin(C(1,2));  // pitch
        Euler(0) = atan2(-C(0, 2), C(2, 2));  // roll
        Euler(2) = atan2(C(1, 0), C(1, 1));  // yaw
    }
    else {
        Euler << -1, -1, -1;
    }
    if (Euler[2] < 0) {
        Euler[2] = M_PI * 2 + Euler[2];
    }
    return Euler;
}

// Function to convert Euler angles to a rotation matrix
Eigen::Matrix3d Euler2C(const Eigen::Vector3d& euler) {
    Eigen::Matrix3d C;
    C << cos(euler(2)) * cos(euler(0)) + sin(euler(2)) * sin(euler(1)) * sin(euler(0)), -sin(euler(2)) * cos(euler(0)) + cos(euler(2)) * sin(euler(1)) * sin(euler(0)), -sin(euler(0)) * cos(euler(1))
        , sin(euler(2))* cos(euler(1)), cos(euler(2))* cos(euler(1)), sin(euler(1))
        , cos(euler(2))* sin(euler(0)) - sin(euler(2)) * cos(euler(0)) * sin(euler(1)), -sin(euler(2)) * sin(euler(0)) - cos(euler(2)) * cos(euler(0)) * sin(euler(1)), cos(euler(1))* cos(euler(0));
    return C;
 }

// Function to convert Euler angles to quaternion
Eigen::Quaterniond Euler2q(const Eigen::Vector3d& euler) {
    return C2q(Euler2C(euler));
}