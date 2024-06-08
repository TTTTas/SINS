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

Eigen::Vector3d SQR_Mat(Eigen::Vector3d v)
{
    Eigen::Vector3d sqr;
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
    Eigen::Matrix3d C = Eigen::Matrix3d::Identity() + sin(n_phi) / n_phi * Skew(Phi) +
        (1 - cos(n_phi)) / (n_phi * n_phi) * Skew(Phi) * Skew(Phi);
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

Eigen::Quaterniond q_n2e(const Eigen::Vector3d& blh)
{
	Eigen::Quaterniond quat;

	double coslon, sinlon, coslat, sinlat;

	coslon = cos(blh[1] * 0.5);
	sinlon = sin(blh[1] * 0.5);
	coslat = cos(-M_PI * 0.25 - blh[0] * 0.5);
	sinlat = sin(-M_PI * 0.25 - blh[0] * 0.5);

	quat.w() = coslat * coslon;
	quat.x() = -sinlat * sinlon;
	quat.y() = sinlat * coslon;
	quat.z() = coslat * sinlon;

	return quat;
}

Eigen::Vector3d q_n2e_2_blh(const Eigen::Quaterniond& qne, double height)
{
	return { -2 * atan(qne.y() / qne.w()) - M_PI * 0.5, 2 * atan2(qne.z(), qne.w()), height };
}

Eigen::Matrix3d DRi(const Eigen::Vector3d& blh)
{
	Eigen::Matrix3d dri = Eigen::Matrix3d::Zero();

    Eigen::Vector2d rmn;
    rmn << Cal_RM(blh[0]), Cal_RN(blh[0]);

	dri(0, 0) = 1.0 / (rmn[0] + blh[2]);
	dri(1, 1) = 1.0 / ((rmn[1] + blh[2]) * cos(blh[0]));
	dri(2, 2) = -1;
	return dri;
}

Eigen::Matrix3d DR(const Eigen::Vector3d& blh)
{
	Eigen::Matrix3d dr = Eigen::Matrix3d::Zero();

    Eigen::Vector2d rmn;
    rmn << Cal_RM(blh[0]), Cal_RN(blh[0]);

	dr(0, 0) = rmn[0] + blh[2];
	dr(1, 1) = (rmn[1] + blh[2]) * cos(blh[0]);
	dr(2, 2) = -1;
	return dr;
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