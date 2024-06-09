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

    // ����omiga_b_ie
    Eigen::Vector3d omiga_b_ie(mean_x, mean_y, mean_z);

    // ����g_b
    Eigen::Vector3d g_b(meanfx, meanfy, meanfz);
    g_b = -g_b;

    // ����omiga_g
    Eigen::Vector3d omiga_g = g_b.normalized();

    // ����omiga_omiga
    Eigen::Vector3d omiga_omiga = (g_b.cross(omiga_b_ie)).normalized();

    // ����omiga_gomiga
    Eigen::Vector3d omiga_gomiga = (g_b.cross(omiga_b_ie).cross(g_b)).normalized();

    // ����C_n_b����
    Eigen::Matrix3d V;
    Eigen::Matrix3d O;
    V.col(0) = v_g;
    V.col(1) = v_omiga;
    V.col(2) = v_gomiga;
    O.row(0) = omiga_g.transpose();
    O.row(1) = omiga_omiga.transpose();
    O.row(2) = omiga_gomiga.transpose();

    Eigen::Matrix3d C_n_b = V * O;

    // ���㸩���ǡ�����Ǻͺ����
    Atti(0) = std::atan(-C_n_b(2, 0) / std::sqrt(C_n_b(2, 1) * C_n_b(2, 1) + C_n_b(2, 2) * C_n_b(2, 2)));
    Atti(1) = std::atan2(C_n_b(2, 1), C_n_b(2, 2));
    Atti(2) = std::atan2(C_n_b(1, 0), C_n_b(0, 0));

    std::cout << std::endl << "��ʼ��̬��:\n" << rad2deg(Atti).transpose() << std::endl;
    return Atti;
}

void insMech(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    // ���ν����ٶȸ��¡�λ�ø��¡���̬����, ���ɵ���˳��
    velUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    posUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
    attUpdate(pva_pre, pva_cur, imu_pre, imu_cur);
}

void velUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    Eigen::Vector3d d_vfb, d_vfn, d_vgn, gl, midvel, midpos;
    Eigen::Vector3d temp1, temp2, temp3;
    Eigen::Matrix3d cnn, I33 = Eigen::Matrix3d::Identity();
    Eigen::Quaterniond qne, qee, qnn;

    // ����������������Ȧ�뾶��î��Ȧ�뾶��������ת���ٶ�ͶӰ��nϵ, nϵ�����eϵת�����ٶ�ͶӰ��nϵ������ֵ
    Eigen::Vector2d rmrn;
    double B_pre = pva_pre.pos(0);
    rmrn << Cal_RM(pva_pre.pos(0)), Cal_RN(pva_pre.pos(0));
    Eigen::Vector3d wie_n, wen_n;
    wie_n << OMEGA_E * cos(B_pre), 0, -OMEGA_E * sin(B_pre);
    wen_n << pva_pre.vel[1] / (rmrn[1] + pva_pre.pos[2]), -pva_pre.vel[0] / (rmrn[0] + pva_pre.pos[2]),
        -pva_pre.vel[1] * tan(pva_pre.pos[0]) / (rmrn[1] + pva_pre.pos[2]);
    double gravity = GRS80_g(pva_pre.pos);

    // ��תЧӦ��˫��������ЧӦ
    temp1 = imu_cur.dtheta.cross(imu_cur.dvel) / 2;
    temp2 = imu_pre.dtheta.cross(imu_cur.dvel) / 12;
    temp3 = imu_pre.dvel.cross(imu_cur.dtheta) / 12;

    // bϵ����������
    d_vfb = imu_cur.dvel + temp1 + temp2 + temp3;

    // ����������
    temp1 = (wie_n + wen_n) * imu_cur.dt / 2;
    cnn = I33 - Skew(temp1);
    d_vfn = cnn * pva_pre.att.cbn * d_vfb;

    // ��������/��ʽ������
    gl << 0, 0, gravity;
    d_vgn = (gl - (2 * wie_n + wen_n).cross(pva_pre.vel)) * imu_cur.dt;

    // �õ��м�ʱ���ٶ�
    midvel = pva_pre.vel + (d_vfn + d_vgn) / 2;

    // ���Ƶõ��м�ʱ��λ��
    qnn = Phi2q(temp1);
    temp2 << 0, 0, -OMEGA_E * imu_cur.dt / 2;
    qee = Phi2q(temp2);
    qne = q_n2e(pva_pre.pos);
    qne = qee * qne * qnn;
    midpos[2] = pva_pre.pos[2] - midvel[2] * imu_cur.dt / 2;
    midpos = q_n2e_2_blh(qne, midpos[2]);

    // ���¼����м�ʱ�̵�rmrn, wie_e, wen_n
    rmrn << Cal_RM(midpos(0)), Cal_RN(midpos(0));
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // ���¼���nϵ��ƽ������������
    temp3 = (wie_n + wen_n) * imu_cur.dt / 2;
    cnn = I33 - Skew(temp3);
    d_vfn = cnn * pva_pre.att.cbn * d_vfb;

    // ���¼�����������ʽ������
    gl << 0, 0, GRS80_g(midpos);
    d_vgn = (gl - (2 * wie_n + wen_n).cross(midvel)) * imu_cur.dt;

    // �ٶȸ������
    pva_cur.vel = pva_pre.vel + d_vfn + d_vgn;
}

void posUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur) {

    Eigen::Vector3d temp1, temp2, midvel, midpos;
    Eigen::Quaterniond qne, qee, qnn;

    // ���¼����м�ʱ�̵��ٶȺ�λ��
    midvel = (pva_cur.vel + pva_pre.vel) / 2;
    midpos = pva_pre.pos + DRi(pva_pre.pos) * midvel * imu_cur.dt / 2;

    // ���¼����м�ʱ�̵������
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn << Cal_RM(midpos[0]), Cal_RN(midpos[0]);
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // ���¼��� kʱ�̵�k-1ʱ�� nϵ��תʸ��
    temp1 = (wie_n + wen_n) * imu_cur.dt;
    qnn = Phi2q(temp1);
    // eϵת����Ч��תʸ�� (k-1ʱ��kʱ�̣�����ȡ����)
    temp2 << 0, 0, -OMEGA_E * imu_cur.dt;
    qee = Phi2q(temp2);

    // λ�ø������
    qne = q_n2e(pva_pre.pos);
    qne = qee * qne * qnn;
    pva_cur.pos[2] = pva_pre.pos[2] - midvel[2] * imu_cur.dt;
    pva_cur.pos = q_n2e_2_blh(qne, pva_cur.pos[2]);
}

void attUpdate(const PVA& pvapre, PVA& pvacur, const IMU& imupre, const IMU& imucur) {

    Eigen::Quaterniond qne_pre, qne_cur, qne_mid, qnn, qbb;
    Eigen::Vector3d temp1, midpos, midvel;

    // ���¼����м�ʱ�̵��ٶȺ�λ��
    midvel = (pvapre.vel + pvacur.vel) / 2;
    qne_pre = q_n2e(pvapre.pos);
    qne_cur = q_n2e(pvacur.pos);
    temp1 = q2Phi(qne_cur.inverse() * qne_pre);
    qne_mid = qne_pre * Phi2q(temp1 / 2).inverse();
    midpos[2] = (pvacur.pos[2] + pvapre.pos[2]) / 2;
    midpos = q_n2e_2_blh(qne_mid, midpos[2]);

    // ���¼����м�ʱ�̵������
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    rmrn << Cal_RM(midpos[0]), Cal_RN(midpos[0]);
    wie_n << OMEGA_E * cos(midpos[0]), 0, -OMEGA_E * sin(midpos[0]);
    wen_n << midvel[1] / (rmrn[1] + midpos[2]), -midvel[0] / (rmrn[0] + midpos[2]),
        -midvel[1] * tan(midpos[0]) / (rmrn[1] + midpos[2]);

    // ����nϵ����ת��Ԫ�� k-1ʱ�̵�kʱ�̱任
    temp1 = -(wie_n + wen_n) * imucur.dt;
    qnn = Phi2q(temp1);

    // ����bϵ��ת��Ԫ�� ��������Բ׶���
    temp1 = imucur.dtheta + imupre.dtheta.cross(imucur.dtheta) / 12;
    qbb = Phi2q(temp1);

    // ��̬�������
    pvacur.att.qbn = qnn * pvapre.att.qbn * qbb;
    pvacur.att.cbn = q2C(pvacur.att.qbn);
    pvacur.att.euler = C2Euler(pvacur.att.cbn);
}