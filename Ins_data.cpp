#include "Ins_data.h"
#include "comm.h"
#include <fstream>

#include "Ins_Mech.h"
#include "progressbar.hpp"

#define _CRT_SECURE_NO_WARNINGS

INS_Eigen::INS_Eigen(GINSOptions& options)
{
    this->options_ = options;
    options_.print_options();
    timestamp_ = 0;

    // ����Э�������ϵͳ�������ϵͳ���״̬�����С
    MatrixXd Cov_=MatrixXd::Zero(RANK, RANK);
    Qc = MatrixXd::Zero(NOISERANK, NOISERANK);

    // ��ʼ��ϵͳ������
    auto imunoise = options_.imunoise;
    Qc.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    Qc.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    Qc.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    Qc.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();
    Qc.block(SGSTD_ID, SGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrscale_std.cwiseProduct(imunoise.gyrscale_std).asDiagonal();
    Qc.block(SASTD_ID, SASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accscale_std.cwiseProduct(imunoise.accscale_std).asDiagonal();

    // ����ϵͳ״̬(λ�á��ٶȡ���̬��IMU���)��ֵ�ͳ�ʼЭ����
	// ��ʼ��λ�á��ٶȡ���̬
    pvacur_.pos = options.initstate.pos;
    pvacur_.vel = options.initstate.vel;
    pvacur_.att.euler = options.initstate.euler;
    pvacur_.att.cbn = Euler2C(pvacur_.att.euler);
    pvacur_.att.qbn = Euler2q(pvacur_.att.euler);
    // ��ʼ��IMU���
    imuerror_ = options.initstate.imuerror;

    // ����һʱ��״̬��ͬ���ĳ�ֵ
    pvapre_ = pvacur_;

    // ��ʼ��Э����
    // initialize covariance
    ImuError imuerror_std = options.initstate_std.imuerror;
    Cov_.block(P_ID, P_ID, 3, 3) = options.initstate_std.pos.cwiseProduct(options.initstate_std.pos).asDiagonal();
    Cov_.block(V_ID, V_ID, 3, 3) = options.initstate_std.vel.cwiseProduct(options.initstate_std.vel).asDiagonal();
    Cov_.block(PHI_ID, PHI_ID, 3, 3) = options.initstate_std.euler.cwiseProduct(options.initstate_std.euler).asDiagonal();
    Cov_.block(BG_ID, BG_ID, 3, 3) = imuerror_std.gyrbias.cwiseProduct(imuerror_std.gyrbias).asDiagonal();
    Cov_.block(BA_ID, BA_ID, 3, 3) = imuerror_std.accbias.cwiseProduct(imuerror_std.accbias).asDiagonal();
    Cov_.block(SG_ID, SG_ID, 3, 3) = imuerror_std.gyrscale.cwiseProduct(imuerror_std.gyrscale).asDiagonal();
    Cov_.block(SA_ID, SA_ID, 3, 3) = imuerror_std.accscale.cwiseProduct(imuerror_std.accscale).asDiagonal();

    kf_ = new INS_KF(MatrixXd::Zero(RANK, 1), Cov_);
}

void INS_Eigen::newImuProcess()
{
    // ��ǰIMUʱ����Ϊϵͳ��ǰ״̬ʱ��,
    timestamp_ = imucur_.time;

    // ���GNSS��Ч���򽫸���ʱ������ΪGNSSʱ��
    double updatetime = gnssdata_.isvalid ? gnssdata_.time : -1;

    // �ж��Ƿ���Ҫ����GNSS����
    int res = isToUpdate(imupre_.time, imucur_.time, updatetime);

    if (res == 0) {
        // ֻ��������״̬
        insPropagation(imupre_, imucur_);
        // �ֶ�����һ��
        kf_->x_hat_ = kf_->x_hat_minus_;
        kf_->P_ = kf_->P_minus_;
    }
    else if (res == 1) {
        // GNSS���ݿ�����һ��Ԫ���ȶ���һ��Ԫ����GNSS����
        // �ֶ�Ԥ��
        kf_->x_hat_minus_ = kf_->x_hat_;
        kf_->P_minus_ = kf_->P_;
        gnssUpdate(gnssdata_);
        stateFeedback();

        pvapre_ = pvacur_;
        insPropagation(imupre_, imucur_);
    }
    else if (res == 2) {
        // GNSS���ݿ�����ǰ��Ԫ���ȶԵ�ǰIMU����״̬����
        insPropagation(imupre_, imucur_);
        gnssUpdate(gnssdata_);
        stateFeedback();
    }
    else {
        // GNSS����������IMU����֮��(�������κ�һ��), ����ǰIMU�ڲ嵽����ʱ��
        IMU midimu;
        imuInterpolate(imupre_, imucur_, updatetime, midimu);

        // ��ǰһ��IMU����״̬����
        insPropagation(imupre_, midimu);

        // ����ʱ�̽���GNSS���£�������ϵͳ״̬
        gnssUpdate(gnssdata_);
        stateFeedback();

        // �Ժ�һ��IMU����״̬����
        pvapre_ = pvacur_;
        insPropagation(midimu, imucur_);
    }

    // ���Э�������Խ���Ԫ��
    checkCov();

    // ������һʱ�̵�״̬��IMU����
    pvapre_ = pvacur_;
    imupre_ = imucur_;
}

NavState INS_Eigen::getNavState()
{
    NavState state;

    state.pos = pvacur_.pos;
    state.vel = pvacur_.vel;
    state.euler = pvacur_.att.euler;
    state.imuerror = imuerror_;

    return state;
}

void INS_Eigen::imuCompensate(IMU& imu)
{
    // ����IMU��ƫ
    imu.dtheta -= imuerror_.gyrbias * imu.dt;
    imu.dvel -= imuerror_.accbias * imu.dt;

    // ����IMU��������
    Eigen::Vector3d gyrscale, accscale;
    gyrscale = Eigen::Vector3d::Ones() + imuerror_.gyrscale;
    accscale = Eigen::Vector3d::Ones() + imuerror_.accscale;
    imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse());
    imu.dvel = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}

int INS_Eigen::isToUpdate(double imutime1, double imutime2, double updatetime) const
{
    if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR) {
        // ����ʱ�俿��imutime1
        return 1;
    }
    else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR) {
        // ����ʱ�俿��imutime2
        return 2;
    }
    else if (imutime1 < updatetime && updatetime < imutime2) {
        // ����ʱ����imutime1��imutime2֮��, ���������κ�һ��
        return 3;
    }
    else {
        // ����ʱ�䲻��imutimt1��imutime2֮�䣬�Ҳ������κ�һ��
        return 0;
    }
}

void INS_Eigen::insPropagation(IMU& imupre, IMU& imucur)
{
    // �Ե�ǰIMU����(imucur)�������, ��һIMU����(imupre)�Ѿ���������
    imuCompensate(imucur);
    // IMU״̬����(��е�����㷨)
    insMech(pvapre_, pvacur_, imupre, imucur);

    // ϵͳ������������̬������phi�����ģ��
    Eigen::MatrixXd Phi, F, Qd, G;

    // ��ʼ��Phi��(״̬ת�ƾ���)��F��Qd��(����������)��G��(����������)
    Phi.resizeLike(kf_->P_);
    F.resizeLike(kf_->P_);
    Qd.resizeLike(kf_->P_);
    G.resize(RANK, NOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // ʹ����һ��Ԫ״̬����״̬ת�ƾ���
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    double gravity;
    rmrn << Cal_RM(pvapre_.pos[0]), Cal_RN(pvapre_.pos[0]);
    gravity = GRS80_g(pvapre_.pos);
    wie_n << OMEGA_E * cos(pvapre_.pos[0]), 0, -OMEGA_E * sin(pvapre_.pos[0]);
    wen_n << pvapre_.vel[1] / (rmrn[1] + pvapre_.pos[2]), -pvapre_.vel[0] / (rmrn[0] + pvapre_.pos[2]),
        -pvapre_.vel[1] * tan(pvapre_.pos[0]) / (rmrn[1] + pvapre_.pos[2]);

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
    double rmh, rnh;

    rmh = rmrn[0] + pvapre_.pos[2];
    rnh = rmrn[1] + pvapre_.pos[2];
    accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;

    // λ�����
    temp.setZero();
    temp(0, 0) = -pvapre_.vel[2] / rmh;
    temp(0, 2) = pvapre_.vel[0] / rmh;
    temp(1, 0) = pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1) = -(pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2) = pvapre_.vel[1] / rnh;
    F.block(P_ID, P_ID, 3, 3) = temp;
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // �ٶ����
    temp.setZero();
    temp(0, 0) = -2 * pvapre_.vel[1] * OMEGA_E * cos(pvapre_.pos[0]) / rmh -
        pow(pvapre_.vel[1], 2) / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(0, 2) = pvapre_.vel[0] * pvapre_.vel[2] / rmh / rmh - pow(pvapre_.vel[1], 2) * tan(pvapre_.pos[0]) / rnh / rnh;
    temp(1, 0) = 2 * OMEGA_E * (pvapre_.vel[0] * cos(pvapre_.pos[0]) - pvapre_.vel[2] * sin(pvapre_.pos[0])) / rmh +
        pvapre_.vel[0] * pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(1, 2) = (pvapre_.vel[1] * pvapre_.vel[2] + pvapre_.vel[0] * pvapre_.vel[1] * tan(pvapre_.pos[0])) / rnh / rnh;
    temp(2, 0) = 2 * OMEGA_E * pvapre_.vel[1] * sin(pvapre_.pos[0]) / rmh;
    temp(2, 2) = -pow(pvapre_.vel[1], 2) / rnh / rnh - pow(pvapre_.vel[0], 2) / rmh / rmh +
        2 * gravity / (sqrt(rmrn[0] * rmrn[1]) + pvapre_.pos[2]);
    F.block(V_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 0) = pvapre_.vel[2] / rmh;
    temp(0, 1) = -2 * (OMEGA_E * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh);
    temp(0, 2) = pvapre_.vel[0] / rmh;
    temp(1, 0) = 2 * OMEGA_E * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1) = (pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2) = 2 * OMEGA_E * cos(pvapre_.pos[0]) + pvapre_.vel[1] / rnh;
    temp(2, 0) = -2 * pvapre_.vel[0] / rmh;
    temp(2, 1) = -2 * (OMEGA_E * cos(pvapre_.pos(0)) + pvapre_.vel[1] / rnh);
    F.block(V_ID, V_ID, 3, 3) = temp;
    F.block(V_ID, PHI_ID, 3, 3) = Skew(pvapre_.att.cbn * accel);
    F.block(V_ID, BA_ID, 3, 3) = pvapre_.att.cbn;
    F.block(V_ID, SA_ID, 3, 3) = pvapre_.att.cbn * (accel.asDiagonal());

    // ��̬���
    temp.setZero();
    temp(0, 0) = -OMEGA_E * sin(pvapre_.pos[0]) / rmh;
    temp(0, 2) = pvapre_.vel[1] / rnh / rnh;
    temp(1, 2) = -pvapre_.vel[0] / rmh / rmh;
    temp(2, 0) = -OMEGA_E * cos(pvapre_.pos[0]) / rmh - pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(2, 2) = -pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh / rnh;
    F.block(PHI_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 1) = 1 / rnh;
    temp(1, 0) = -1 / rmh;
    temp(2, 1) = -tan(pvapre_.pos[0]) / rnh;
    F.block(PHI_ID, V_ID, 3, 3) = temp;
    F.block(PHI_ID, PHI_ID, 3, 3) = -Skew(wie_n + wen_n);
    F.block(PHI_ID, BG_ID, 3, 3) = -pvapre_.att.cbn;
    F.block(PHI_ID, SG_ID, 3, 3) = -pvapre_.att.cbn * (omega.asDiagonal());

    // IMU��ƫ���ͱ�����������ģ��һ�׸�˹-����Ʒ����
    F.block(BG_ID, BG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(BA_ID, BA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SG_ID, SG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SA_ID, SA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();

    // ϵͳ������������
    G.block(V_ID, VRW_ID, 3, 3) = pvapre_.att.cbn;
    G.block(PHI_ID, ARW_ID, 3, 3) = pvapre_.att.cbn;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SG_ID, SGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SA_ID, SASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // ״̬ת�ƾ���
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // ����ϵͳ��������
    Qd = G * Qc * G.transpose() * imucur.dt;
    Qd = (Phi * Qd * Phi.transpose() + Qd) / 2;

    // EKFԤ�⴫��ϵͳЭ�����ϵͳ���״̬
    kf_->set_A(Phi);
    kf_->set_Q(Qd);
    kf_->predict();
}

void INS_Eigen::gnssUpdate(GNSS& gnssdata)
{
    // IMUλ��ת��GNSS������λ����λ��
    Eigen::Vector3d antenna_pos;
    Eigen::Matrix3d Dr, Dr_inv;
    Dr_inv = DRi(pvacur_.pos);
    Dr = DR(pvacur_.pos);
    antenna_pos = pvacur_.pos + Dr_inv * pvacur_.att.cbn * options_.antlever;

    // GNSSλ�ò�����Ϣ
    Eigen::MatrixXd dz;
    dz = Dr * (antenna_pos - gnssdata.blh);

    // ����GNSSλ�ù۲����
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos;
    H_gnsspos.resize(3, kf_->Q_.rows());
    H_gnsspos.setZero();
    H_gnsspos.block(0, P_ID, 3, 3) = Eigen::Matrix3d::Identity();
    H_gnsspos.block(0, PHI_ID, 3, 3) = Skew(pvacur_.att.cbn * options_.antlever);

    // λ�ù۲�������
    Eigen::MatrixXd R_gnsspos;
    R_gnsspos = gnssdata.std.cwiseProduct(gnssdata.std).asDiagonal();

    // EKF����Э��������״̬
    kf_->set_H(H_gnsspos);
    kf_->set_Z(dz);
    kf_->set_R(R_gnsspos);
    kf_->update();

    // GNSS����֮������Ϊ������
    gnssdata.isvalid = false;
}

void INS_Eigen::stateFeedback()
{
    Eigen::Vector3d vectemp;

    // λ������
    Eigen::Vector3d delta_r = kf_->x_hat_.block(P_ID, 0, 3, 1);
    Vector2d rmn; rmn << Cal_RM(pvapre_.pos[0]), Cal_RN(pvapre_.pos[0]);
    Eigen::Matrix3d Dr_inv = DRi(pvacur_.pos);
    pvacur_.pos -= Dr_inv * delta_r;

    // �ٶ�����
    vectemp = kf_->x_hat_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;

    // ��̬����
    vectemp = kf_->x_hat_.block(PHI_ID, 0, 3, 1);
    Eigen::Quaterniond qpn = Phi2q(vectemp);
    pvacur_.att.qbn = qpn * pvacur_.att.qbn;
    pvacur_.att.cbn = q2C(pvacur_.att.qbn);
    pvacur_.att.euler = C2Euler(pvacur_.att.cbn);

    // IMU��ƫ����
    vectemp = kf_->x_hat_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = kf_->x_hat_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // IMU������������
    vectemp = kf_->x_hat_.block(SG_ID, 0, 3, 1);
    imuerror_.gyrscale += vectemp;
    vectemp = kf_->x_hat_.block(SA_ID, 0, 3, 1);
    imuerror_.accscale += vectemp;

    // ���״̬������ϵͳ״̬��,�����״̬����
    kf_->x_hat_.setZero();
}


