#include "Ins_Mech.h"

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

void insMech(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    // 依次进行速度更新、位置更新、姿态更新, 不可调换顺序
    velUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    posUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    attUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
}

void velUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
    Eigen::Vector3d temp1, temp2, temp3;
    Eigen::Matrix3d cnn, I33 = Eigen::Matrix3d::Identity();
    Eigen::Quaterniond qne, qee, qnn;

    // 计算地理参数，子午圈半径和卯酉圈半径，地球自转角速度投影到n系, n系相对于e系转动角速度投影到n系，重力值
    Eigen::Vector2d rmrn;
    double B_pre = pva_pre.pos(0);
    rmrn << Cal_RM(pva_pre.pos(0)), Cal_RN(pva_pre.pos(0));
    Eigen::Vector3d wie_n, wen_n;
    wie_n << OMEGA_E * cos(B_pre), 0, -OMEGA_E * sin(B_pre);
    wen_n << pva_pre.vel[1] / (rmrn[1] + pva_pre.pos[2]), -pva_pre.vel[0] / (rmrn[0] + pva_pre.pos[2]),
        -pva_pre.vel[1] * tan(pva_pre.pos[0]) / (rmrn[1] + pva_pre.pos[2]);
    double gravity = GRS80_g(pva_pre.pos);

    // 旋转效应和双子样划桨效应
    temp1 = imu_cur.dtheta.cross(imu_cur.dvel) / 2;
    temp2 = imu_pre.dtheta.cross(imu_cur.dvel) / 12;
    temp3 = imu_pre.dvel.cross(imu_cur.dtheta) / 12;

    // b系比力积分项
    d_vfb = imu_cur.dvel + temp1 + temp2 + temp3;

    // 比力积分项
    temp1 = (wie_n + wen_n) * imu_cur.dt / 2;
    cnn = I33 - Skew(temp1);
    d_vfn = cnn * pva_pre.att.cbn * d_vfb;

    // 计算重力/哥式积分项
    gl << 0, 0, gravity;
    d_vgn = (gl - (2 * wie_n + wen_n).cross(pva_pre.vel)) * imu_cur.dt;

    // 得到中间时刻速度
    midvel = pva_pre.vel + (d_vfn + d_vgn) / 2;

    // 外推得到中间时刻位置
    qnn = Phi2q(temp1);
    temp2 << 0, 0, -OMEGA_E * imu_cur.dt / 2;
    qee = Phi2q(temp2);
    qne = q_n2e(pva_pre.pos);
    qne = qee * qne * qnn;
    midpos[2] = pva_pre.pos[2] - midvel[2] * imu_cur.dt / 2;
    midpos = q_n2e_2_blh(qne, midpos[2]);

    // 重新计算中间时刻的rmrn, wie_e, wen_n
    rmrn << Cal_RM(midpos(0)), Cal_RN(midpos(0));
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算n系下平均比力积分项
    temp3 = (wie_n + wen_n) * imu_cur.dt / 2;
    cnn = I33 - Skew(temp3);
    d_vfn = cnn * pva_pre.att.cbn * d_vfb;

    // 重新计算重力、哥式积分项
    gl << 0, 0, GRS80_g(midpos);
    d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * imu_cur.dt;

    // 速度更新完成
    pva_cur.vel = pva_pre.vel + d_vfn + d_vgn;
}

void posUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    Eigen::Vector3d temp1, temp2, midvel, midpos;
    Eigen::Quaterniond qne, qee, qnn;

    // 重新计算中间时刻的速度和位置
    midvel = (pva_cur.vel + pva_pre.vel) / 2;
    midpos = pva_pre.pos + DRi(pva_pre.pos) * midvel * imu_cur.dt / 2;

    // 重新计算中间时刻地理参数
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn << Cal_RM(midpos[0]), Cal_RN(midpos[0]);
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 重新计算 k时刻到k-1时刻 n系旋转矢量
    temp1 = (wie_n + wen_n) * imu_cur.dt;
    qnn = Phi2q(temp1);
    // e系转动等效旋转矢量 (k-1时刻k时刻，所以取负号)
    temp2 << 0, 0, -OMEGA_E * imu_cur.dt;
    qee = Phi2q(temp2);

    // 位置更新完成
    qne = q_n2e(pva_pre.pos);
    qne = qee * qne * qnn;
    pva_cur.pos[2] = pva_pre.pos[2] - midvel[2] * imu_cur.dt;
    pva_cur.pos = q_n2e_2_blh(qne, pva_cur.pos[2]);
}

void attUpdate(const PVA& pvapre, PVA& pvacur, const IMU& imupre, const IMU& imucur) {

    Eigen::Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
    Eigen::Vector3d temp1, midpos, midvel;

    // 重新计算中间时刻的速度和位置
    midvel = (pvapre.vel + pvacur.vel) / 2;
    qne_pre = q_n2e(pvapre.pos);
    qne_cur = q_n2e(pvacur.pos);
    temp1 = q2Phi(qne_cur.inverse() * qne_pre);
    qne_mid = qne_pre * Phi2q(temp1 / 2).inverse();
    midpos[2] = (pvacur.pos[2] + pvapre.pos[2]) / 2;
    midpos = q_n2e_2_blh(qne_mid, midpos[2]);

    // 重新计算中间时刻地理参数
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn << Cal_RM(midpos[0]), Cal_RN(midpos[0]);
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // 计算n系的旋转四元数 k-1时刻到k时刻变换
    temp1 = -(wie_n + wen_n) * imucur.dt;
    qnn = Phi2q(temp1);

    // 计算b系旋转四元数 补偿二阶圆锥误差
    temp1 = imucur.dtheta + imupre.dtheta.cross(imucur.dtheta) / 12;
    qbb = Phi2q(temp1);

    // 姿态更新完成
    pvacur.att.qbn = qnn * pvapre.att.qbn * qbb;
    pvacur.att.cbn = q2C(pvacur.att.qbn);
    pvacur.att.euler = C2Euler(pvacur.att.cbn);
}