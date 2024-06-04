#include "comm.h"

// Function to convert degrees to radians
 double deg2rad(double degrees) {
    return degrees * DEG2RAD;
}

 Vector3d deg2rad(Vector3d& degrees)
{
    return degrees * DEG2RAD;
}

// Function to convert radians to degrees
 double rad2deg(double radians) {
    return radians * RAD2DEG;
}

 Vector3d rad2deg(Vector3d& radians)
{
    return radians * RAD2DEG;
}


// 线性外推函数
 Vector3d Extrapol(const Vector3d& v0, const Vector3d& v1) {
    // 输入 v0     历元0数据
    //      v1     历元1数据
    // 输出 result 外推历元2数据
    Vector3d result = 1.5 * v1 - 0.5 * v0;
    return result;
}

// Function to calculate the local gravity value
 double GRS80_g(const Vector3d& pos) {
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
 Matrix3d Cross_vector(const Vector3d& vector) {
    Matrix3d cross_matrix;
    cross_matrix << 0, -vector(2), vector(1),
        vector(2), 0, -vector(0),
        -vector(1), vector(0), 0;
    return cross_matrix;
}


// Function to convert rotation matrix to quaternion
 Quaterniond C2q(const Matrix3d& C) {
    return Quaterniond(C);
}

// Function to convert quaternion to rotation matrix
 Matrix3d q2C(const Quaterniond& q) {
    return q.toRotationMatrix();
}

// Function to convert quaternion to equivalent rotation vector
 Vector3d q2Phi(const Quaterniond& q) {
    AngleAxisd angleAxis(q);
    double angle = angleAxis.angle();
    Vector3d axis = angleAxis.axis();
    return angle * axis;
}

// Function to convert equivalent rotation vector to rotation matrix
 Matrix3d Phi2C(const Vector3d& Phi) {
    double n_phi = Phi.norm();
    Matrix3d C = Matrix3d::Identity() + sin(n_phi) / n_phi * Cross_vector(Phi) +
        (1 - cos(n_phi)) / (n_phi * n_phi) * Cross_vector(Phi) * Cross_vector(Phi);
    return C;
}

// Function to convert equivalent rotation vector to quaternion
 Quaterniond Phi2q(const Vector3d& Phi) {
    double half_phi = 0.5 * Phi.norm();
    if (half_phi != 0) {
        return Quaterniond(AngleAxisd(half_phi, Phi.normalized()));
    }
    else {
        return Quaterniond(1, 0, 0, 0);
    }
}

// Function to convert a rotation matrix to Euler angles
 Vector3d C2Euler(const Matrix3d& C) {
    Vector3d Euler;
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
Matrix3d Euler2C(const Vector3d& Euler) {
    double theta = Euler(1);
    double phi = Euler(0);
    double psi = Euler(2);
    Matrix3d C;
    C << cos(theta) * cos(psi), -cos(phi) * sin(psi) + sin(phi) * sin(theta) * cos(psi), sin(phi)* sin(psi) + cos(phi) * sin(theta) * cos(psi),
        cos(theta)* sin(psi), cos(phi)* cos(psi) + sin(phi) * sin(theta) * sin(psi), -sin(phi) * cos(psi) + cos(phi) * sin(theta) * sin(psi),
        -sin(theta), sin(phi)* cos(theta), cos(phi)* cos(theta);
    return C;
}

// Function to convert Euler angles to quaternion
Quaterniond Euler2q(const Vector3d& Euler) {
    AngleAxisd rollAngle(Euler(0), Vector3d::UnitX());
    AngleAxisd pitchAngle(Euler(1), Vector3d::UnitY());
    AngleAxisd yawAngle(Euler(2), Vector3d::UnitZ());
    Quaterniond q = rollAngle * pitchAngle * yawAngle;
    return q;
}

// Function to update Euler angles to quaternion
 Quaterniond Update_Euler_q(const Vector3d& E, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& pos, const Vector3d& v, double dt) {
    double omiga_e = 7.292115e-5;
    double B = pos(0); // Convert degrees to radians

    Vector3d omiga_ie;
    omiga_ie << omiga_e * cos(B), 0, -omiga_e * sin(B);

    Vector3d omiga_en;
    omiga_en << v(1) / (Cal_RN(B) + pos(2)),
        -v(0) / (Cal_RM(B) + pos(2)),
        -v(1) * tan(B) / (Cal_RN(B) + pos(2));

    Vector3d phi_k = theta1 + theta0.cross(theta1) / 12.0;
    Vector3d zeta = (omiga_en + omiga_ie) * dt;

    Quaterniond q_bb = Phi2q(phi_k);
    Quaterniond q_nn = Phi2q(zeta).conjugate();
    Quaterniond q0 = Euler2q(E);

    return (q_nn * q0) * q_bb;
}


// Function to update the position
Vector3d Update_pos(const Vector3d& pos0, const Vector3d& v0, const Vector3d& v1, double dt) {
    double h = pos0(2) - (v0(2) + v1(2)) / 2.0 * dt;
    double meanh = (h + pos0(2)) / 2.0;
    double B = pos0(0) + (v0(0) + v1(0)) / (2.0 * (Cal_RM(pos0(0)) + meanh)) * dt;
    double meanB = (pos0(0) + B) / 2.0;
    double L = pos0(1) + (v0(1) + v1(1)) / ((2.0 * (Cal_RN(meanB)) + meanh) * cos(meanB)) * dt;
    return Vector3d(B, L, h);
}

// Function to convert BLH to NE coordinates
 vector<Vector2d> BLH2NE(const vector<Vector3d>& BLH, const Vector3d& BLH0) {
    double B = BLH0(0);
    double L = BLH0(1);
    double h = BLH0(2);
    double R_M = Cal_RM(B);
    double R_N = Cal_RN(B);
    vector<Vector2d> NE;

    for (auto blh : BLH) {
        Vector2d ne;
        ne(0) = (blh(0) - B) * (R_M + h);
        ne(1) = (blh(1) - L) * (R_N + h) * cos(B);
        NE.push_back(ne);
    }

    return NE;
}

// Function to update attitude matrix
Matrix3d Update_Euler_C(const Matrix3d& E, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& pos, const Vector3d& v, double dt) {
    double B = pos(0);
    Vector3d omiga_ie = Vector3d(OMEGA_E * cos(B), 0, -OMEGA_E * sin(B));
    Vector3d omiga_en = Vector3d(v(1) / (Cal_RN(B) + pos(2)), -v(0) / (Cal_RM(B) + pos(2)), -v(1) * tan(B) / (Cal_RN(B) + pos(2)));
    Vector3d phi_k = theta1 + theta0.cross(theta1) / 12;
    Vector3d zeta = (omiga_en + omiga_ie) * dt;
    Matrix3d C_bb = Phi2C(phi_k);
    Matrix3d C_nn = Phi2C(-zeta);
    Matrix3d C0 = E;
    return C_nn * C0 * C_bb;
}

// Function to calculate gravity/Coriolis integral term
 Vector3d Cal_deltv_g_cor(const Vector3d g, const Vector3d& omiga_ie, const Vector3d& omiga_en, const Vector3d& v0, double dt) {
    return (g - (2 * omiga_ie + omiga_en).cross(v0)) * dt;
}

// Function to calculate specific force integral term
 Vector3d Cal_deltv_f(const Vector3d& omiga_ie, const Vector3d& omiga_en, const Matrix3d& C, const Vector3d& theta0, const Vector3d& theta1, const Vector3d& v0, const Vector3d& v1, double dt) {
    Vector3d zeta = (omiga_ie + omiga_en) * dt;
    Vector3d dv_f0 = v1 + theta1.cross(v1) / 2 + (theta0.cross(v1) + v0.cross(theta1)) / 12;
    return (Matrix3d::Identity() - 0.5 * Cross_vector(zeta)) * C * dv_f0;
}

// Function to update velocity
Vector3d Update_velocity(const Vector3d& v00, const Vector3d& v0, const Vector3d& pos, const Vector3d& pos0, double dt, const Vector3d& dv0, const Vector3d& dv1, const Matrix3d& C, const Vector3d& theta0, const Vector3d& theta1) {
    double B = pos(0);
    Vector3d omiga_ie = Vector3d(OMEGA_E * cos(B), 0, -OMEGA_E * sin(B));
    Vector3d omiga_en = Vector3d(v0(1) / (Cal_RN(B) + pos(2)), -v0(0) / (Cal_RM(B) + pos(2)), -v0(1) * tan(B) / (Cal_RN(B) + pos(2)));
    Vector3d g = Vector3d(0, 0, GRS80_g(pos));
    double B0 = pos0(0);
    Vector3d omiga_ie0 = Vector3d(OMEGA_E * cos(B0), 0, -OMEGA_E * sin(B0));
    Vector3d omiga_en0 = Vector3d(v00(1) / (Cal_RN(B0) + pos(2)), -v00(0) / (Cal_RM(B0) + pos(2)), -v00(1) * tan(B0) / (Cal_RN(B0) + pos(2)));
    Vector3d g0 = Vector3d(0, 0, GRS80_g(pos));
    g = Extrapol(g0, g);
    omiga_ie = Extrapol(omiga_ie0, omiga_ie);
    omiga_en = Extrapol(omiga_en0, omiga_en);
    Vector3d dv_g_cor = Cal_deltv_g_cor(g, omiga_ie, omiga_en, Extrapol(v00, v0), dt);
    Vector3d dv_f = Cal_deltv_f(omiga_ie, omiga_en, C, theta0, theta1, dv0, dv1, dt);
    return v0 + dv_f + dv_g_cor;
}