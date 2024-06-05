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

// Function to calculate the local gravity value
 double GRS80_g(const Eigen::Vector3d& pos) {
    double B = pos(0);
    double H = pos(2);
    double g0 = 9.7803267715 * (1 + 0.0052790414 * sin(B) * sin(B) + 0.0000232718 * sin(B) * sin(B) * sin(B) * sin(B));
    return g0 - (3.087691089e-6 - 4.397731e-9 * sin(B) * sin(B)) * H + 0.721e-12 * H * H;
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
 Eigen::Matrix3d Cross_vector(const Eigen::Vector3d& vector) {
    Eigen::Matrix3d cross_matrix;
    cross_matrix << 0, -vector(2), vector(1),
        vector(2), 0, -vector(0),
        -vector(1), vector(0), 0;
    return cross_matrix;
}


// Function to convert rotation matrix to quaternion
Eigen::Quaterniond C2q(const Eigen::Matrix3d& C) {
    return Eigen::Quaterniond(C);
}

// Function to convert quaternion to rotation matrix
 Eigen::Matrix3d q2C(const Eigen::Quaterniond& q) {
    return q.toRotationMatrix();
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
    Eigen::Matrix3d C = Eigen::Matrix3d::Identity() + sin(n_phi) / n_phi * Cross_vector(Phi) +
        (1 - cos(n_phi)) / (n_phi * n_phi) * Cross_vector(Phi) * Cross_vector(Phi);
    return C;
}

// Function to convert equivalent rotation vector to quaternion
 Eigen::Quaterniond Phi2q(const Eigen::Vector3d& Phi) {
    double half_phi = 0.5 * Phi.norm();
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
        Euler(1) = atan(-C(2, 0) / sqrt(C(2, 1) * C(2, 1) + C(2, 2) * C(2, 2)));  // pitch
        Euler(0) = atan2(C(2, 1), C(2, 2));  // roll
        Euler(2) = atan2(C(1, 0), C(0, 0));  // yaw
    }
    else {
        Euler << -1, -1, -1;
    }
    return Euler;
}

// Function to convert Euler angles to a rotation matrix
Eigen::Matrix3d Euler2C(const Eigen::Vector3d& Euler) {
    double theta = Euler(1);
    double phi = Euler(0);
    double psi = Euler(2);
    Eigen::Matrix3d C;
    C << cos(theta) * cos(psi), -cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi), sin(phi)* sin(psi) + cos(phi) * sin(theta) * cos(psi),
        cos(theta)* sin(psi), cos(phi)* cos(psi) + sin(phi) * sin(theta) * sin(psi), -sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi),
        -sin(theta), sin(phi)* cos(theta), cos(phi)* cos(theta);
    return C;
}

// Function to convert Euler angles to quaternion
Eigen::Quaterniond Euler2q(const Eigen::Vector3d& Euler) {
    Eigen::AngleAxisd rollAngle(Euler(0), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(Euler(1), Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(Euler(2), Eigen::Vector3d::UnitZ());
    Eigen::Quaterniond q = rollAngle * pitchAngle * yawAngle;
    return q;
}

// Function to update Euler angles to quaternion
Eigen::Quaterniond Update_Euler_q(const Eigen::Vector3d& E, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& pos, const Eigen::Vector3d& v, double dt) {
    double omiga_e = 7.292115e-5;
    double B = pos(0); // Convert degrees to radians

    Eigen::Vector3d omiga_ie;
    omiga_ie << omiga_e * cos(B), 0, -omiga_e * sin(B);

    Eigen::Vector3d omiga_en;
    omiga_en << v(1) / (Cal_RN(B) + pos(2)),
        -v(0) / (Cal_RM(B) + pos(2)),
        -v(1) * tan(B) / (Cal_RN(B) + pos(2));

    Eigen::Vector3d phi_k = theta1 + theta0.cross(theta1) / 12.0;
    Eigen::Vector3d zeta = (omiga_en + omiga_ie) * dt;

    Eigen::Quaterniond q_bb = Phi2q(phi_k);
    Eigen::Quaterniond q_nn = Phi2q(zeta).conjugate();
    Eigen::Quaterniond q0 = Euler2q(E);

    return (q_nn * q0) * q_bb;
}


// Function to update the position
Eigen::Vector3d Update_pos(const Eigen::Vector3d& pos0, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double dt) {
    double h = pos0(2) - (v0(2) + v1(2)) / 2.0 * dt;
    double meanh = (h + pos0(2)) / 2.0;
    double B = pos0(0) + (v0(0) + v1(0)) / (2.0 * (Cal_RM(pos0(0)) + meanh)) * dt;
    double meanB = (pos0(0) + B) / 2.0;
    double L = pos0(1) + (v0(1) + v1(1)) / ((2.0 * (Cal_RN(meanB)) + meanh) * cos(meanB)) * dt;
    return Eigen::Vector3d(B, L, h);
}

// Function to convert BLH to NE coordinates
std::vector<Eigen::Vector2d> BLH2NE(const std::vector<Eigen::Vector3d>& BLH, const Eigen::Vector3d& BLH0) {
    double B = BLH0(0);
    double L = BLH0(1);
    double h = BLH0(2);
    double R_M = Cal_RM(B);
    double R_N = Cal_RN(B);
    std::vector<Eigen::Vector2d> NE;

    for (auto blh : BLH) {
        Eigen::Vector2d ne;
        ne(0) = (blh(0) - B) * (R_M + h);
        ne(1) = (blh(1) - L) * (R_N + h) * cos(B);
        NE.push_back(ne);
    }

    return NE;
}

// Function to update attitude matrix
Eigen::Matrix3d Update_Euler_C(const Eigen::Matrix3d& E, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& pos, const Eigen::Vector3d& v, double dt) {
    double B = pos(0);
    Eigen::Vector3d omiga_ie = Eigen::Vector3d(OMEGA_E * cos(B), 0, -OMEGA_E * sin(B));
    Eigen::Vector3d omiga_en = Eigen::Vector3d(v(1) / (Cal_RN(B) + pos(2)), -v(0) / (Cal_RM(B) + pos(2)), -v(1) * tan(B) / (Cal_RN(B) + pos(2)));
    Eigen::Vector3d phi_k = theta1 + theta0.cross(theta1) / 12;
    Eigen::Vector3d zeta = (omiga_en + omiga_ie) * dt;
    Eigen::Matrix3d C_bb = Phi2C(phi_k);
    Eigen::Matrix3d C_nn = Phi2C(-zeta);
    Eigen::Matrix3d C0 = E;
    return C_nn * C0 * C_bb;
}

// Function to calculate gravity/Coriolis integral term
 Eigen::Vector3d Cal_deltv_g_cor(const Eigen::Vector3d g, const Eigen::Vector3d& omiga_ie, const Eigen::Vector3d& omiga_en, const Eigen::Vector3d& v0, double dt) {
    return (g - (2 * omiga_ie + omiga_en).cross(v0)) * dt;
}

// Function to calculate specific force integral term
 Eigen::Vector3d Cal_deltv_f(const Eigen::Vector3d& omiga_ie, const Eigen::Vector3d& omiga_en, const Eigen::Matrix3d& C, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, double dt) {
    Eigen::Vector3d zeta = (omiga_ie + omiga_en) * dt;
    Eigen::Vector3d dv_f0 = v1 + theta1.cross(v1) / 2 + (theta0.cross(v1) + v0.cross(theta1)) / 12;
    return (Eigen::Matrix3d::Identity() - 0.5 * Cross_vector(zeta)) * C * dv_f0;
}

// Function to update velocity
Eigen::Vector3d Update_velocity(const Eigen::Vector3d& v00, const Eigen::Vector3d& v0, const Eigen::Vector3d& pos, const Eigen::Vector3d& pos0, double dt, const Eigen::Vector3d& dv0, const Eigen::Vector3d& dv1, const Eigen::Matrix3d& C, const Eigen::Vector3d& theta0, const Eigen::Vector3d& theta1) {
    double B = pos(0);
    Eigen::Vector3d omiga_ie = Eigen::Vector3d(OMEGA_E * cos(B), 0, -OMEGA_E * sin(B));
    Eigen::Vector3d omiga_en = Eigen::Vector3d(v0(1) / (Cal_RN(B) + pos(2)), -v0(0) / (Cal_RM(B) + pos(2)), -v0(1) * tan(B) / (Cal_RN(B) + pos(2)));
    Eigen::Vector3d g = Eigen::Vector3d(0, 0, GRS80_g(pos));
    double B0 = pos0(0);
    Eigen::Vector3d omiga_ie0 = Eigen::Vector3d(OMEGA_E * cos(B0), 0, -OMEGA_E * sin(B0));
    Eigen::Vector3d omiga_en0 = Eigen::Vector3d(v00(1) / (Cal_RN(B0) + pos(2)), -v00(0) / (Cal_RM(B0) + pos(2)), -v00(1) * tan(B0) / (Cal_RN(B0) + pos(2)));
    Eigen::Vector3d g0 = Eigen::Vector3d(0, 0, GRS80_g(pos));
    g = Extrapol(g0, g);
    omiga_ie = Extrapol(omiga_ie0, omiga_ie);
    omiga_en = Extrapol(omiga_en0, omiga_en);
    Eigen::Vector3d dv_g_cor = Cal_deltv_g_cor(g, omiga_ie, omiga_en, Extrapol(v00, v0), dt);
    Eigen::Vector3d dv_f = Cal_deltv_f(omiga_ie, omiga_en, C, theta0, theta1, dv0, dv1, dt);
    return v0 + dv_f + dv_g_cor;
}