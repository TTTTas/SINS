#include "Ins_Mech.h"

#include "cal.h"

Eigen::Vector3d Init_Yaw(std::vector<IMU*> imu_data, INS_Configure cfg)
{
    Eigen::Vector3d Atti;
    double mean_x, mean_y, mean_z;
    double meanfx = 0.0, meanfy = 0.0, meanfz = 0.0;
    mean_x = mean_y = mean_z = 0.0;
    for (auto imu: imu_data)
    {
        mean_x += imu->dtheta[0];
        mean_y += imu->dtheta[1];
        mean_z += imu->dtheta[2];
        meanfx += imu->dvel[0];
        meanfy += imu->dvel[1];
        meanfz += imu->dvel[2];
    }
    mean_x /= imu_data.size();
    mean_y /= imu_data.size();
    mean_z /= imu_data.size();
    meanfx /= imu_data.size();
    meanfy /= imu_data.size();
    meanfz /= imu_data.size();

    Eigen::Vector3d g_n;
    g_n << 0, 0, GRS80_g(cfg.gins_options.initstate.pos);
    Eigen::Vector3d v_g = g_n.normalized();
    Eigen::Vector3d omiga_n_ie(OMEGA_E * cos(Atti(0)), 0, -OMEGA_E * sin(Atti(0)));
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
    return Atti;
}

void insMech(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur)
{
    attUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    velUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    posUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
}

void attUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imupre, const IMU& imucur) {

    Eigen::Matrix3d R_theta, R_omega;
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pva_pre.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    Matrix3d R_eb_pre = pva_pre.att.cbn * get_Rot(blh[0], blh[1]);

    R_theta = Phi2C(-imucur.dtheta);
    R_omega = Phi2C(Eigen::Vector3d(0, 0, OMEGA_E * imucur.dt));

    Matrix3d R_eb_cur = R_theta * R_eb_pre * R_omega;
    
    // 姿态更新完成
    pva_cur.att.cbn = R_eb_cur * get_Rot(blh[0], blh[1]).transpose();
    pva_cur.att.qbn = C2q(pva_cur.att.cbn);
    pva_cur.att.euler = C2Euler(pva_cur.att.cbn);
}

void velUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    Eigen::Vector3d d_theta_v, d_vfe, d_vge;
    Eigen::Vector3d wie_e, wen_b;
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pva_pre.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    wie_e << 0, 0, OMEGA_E;
    Matrix3d R_be = get_Rot(blh[0], blh[1]).transpose() * pva_cur.att.cbn.transpose();
    wen_b << R_be.transpose() * wie_e;

    Eigen::Vector3d gravity_e = g_e(blh);

    // 地球自转角速度增量
    d_theta_v = imu_cur.dtheta - wen_b * imu_cur.dt;

    // 比例积分项
    d_vfe = R_be * (imu_cur.dvel - 0.5 * Skew(d_theta_v) * imu_cur.dvel);

    // 重力/科氏力补偿
    d_vge = (gravity_e - 2 * Skew(wie_e) * pva_pre.vel) * imu_cur.dt;

    // 速度更新完成
    pva_cur.vel = pva_pre.vel + d_vfe + d_vge;
}

void posUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur)
{
    pva_cur.pos = pva_pre.pos + 0.5 * (pva_pre.vel + pva_cur.vel) * imu_cur.dt;
}