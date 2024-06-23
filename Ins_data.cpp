#include "Ins_data.h"
#include "comm.h"
#include <fstream>

#include "Ins_Mech.h"
#include "progressbar.hpp"
#include "sockets.h"

#define _CRT_SECURE_NO_WARNINGS

INS_Eigen::INS_Eigen(GINSOptions& options)
{
    this->options_ = options;
    options_.print_options();
    timestamp_ = 0;

    // ����Э�������ϵͳ�������ϵͳ���״̬�����С
    MatrixXd Cov_=MatrixXd::Zero(RANK, RANK);
    Qc = MatrixXd::Zero(RANK, RANK);

    // ��ʼ��ϵͳ������
    auto imunoise = options_.imunoise;
    Qc.block(P_ID, P_ID, 3, 3) = imunoise.pos_prw.cwiseProduct(0 * imunoise.pos_prw).asDiagonal();
    Qc.block(PHI_ID, PHI_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    Qc.block(V_ID, V_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    Qc.block(BG_ID, BG_ID, 3, 3) = imunoise.gyrbias_std.cwiseProduct(0 * imunoise.gyrbias_std).asDiagonal();
    Qc.block(BA_ID, BA_ID, 3, 3) = imunoise.accbias_std.cwiseProduct(0 * imunoise.accbias_std).asDiagonal();


    // ����ϵͳ״̬(λ�á��ٶȡ���̬��IMU���)��ֵ�ͳ�ʼЭ����
	// ��ʼ��λ�á��ٶȡ���̬
    pvacur_.pos = options.initstate.pos;
    pvacur_.vel = options.initstate.vel;
    pvacur_.att.euler = options.initstate.euler;
    pvacur_.att.cbn = Euler2C(pvacur_.att.euler);
    pvacur_.att.qbn = C2q(pvacur_.att.cbn);
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

    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pvacur_.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    Eigen::Vector3d gravity_e = g_e(blh);

    // ʹ����һ��Ԫ״̬����״̬ת�ƾ���
    Eigen::Vector3d wie_e, wen_b;
    wie_e << 0, 0, OMEGA_E;
    Matrix3d R_be = get_Rot(blh[0], blh[1]).transpose() * pvacur_.att.cbn.transpose();
    wen_b << R_be.transpose() * wie_e;
    

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
	accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;
    // λ�����
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // �ٶ����
    F.block(V_ID, V_ID, 3, 3) = -2 * Skew(wie_e);
    F.block(V_ID, PHI_ID, 3, 3) = Skew(R_be * accel);
    F.block(V_ID, BA_ID, 3, 3) = R_be;

    // ��̬���
    F.block(PHI_ID, PHI_ID, 3, 3) = -Skew(wie_e);
    F.block(PHI_ID, BG_ID, 3, 3) = -R_be;

    // ϵͳ������������
    Qc.block(V_ID, V_ID, 3, 3) = R_be * options_.imunoise.acc_vrw.cwiseProduct(options_.imunoise.acc_vrw).asDiagonal() * R_be.transpose();
    Qc.block(ARW_ID,ARW_ID,3,3) = R_be * options_.imunoise.gyr_arw.cwiseProduct(options_.imunoise.gyr_arw).asDiagonal() * R_be.transpose();

    // ״̬ת�ƾ���
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // ����ϵͳ��������
    Qd = (Phi * Qc * Phi.transpose() + Qc) * imucur.dt / 2;

    // EKFԤ�⴫��ϵͳЭ�����ϵͳ���״̬
    kf_->set_A(Phi);
    kf_->set_Q(Qd);
    kf_->predict();
}

// GNSSλ�ò�����Ϣ
void INS_Eigen::gnssUpdate(GNSS& gnssdata)
{
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pvacur_.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    Matrix3d R_be = get_Rot(blh[0], blh[1]).transpose() * pvacur_.att.cbn.transpose();

    // ʹ����һ��Ԫ״̬����״̬ת�ƾ���
    // IMUλ��ת��GNSS������λ����λ��
    Eigen::Vector3d antenna_pos = pvacur_.pos + R_be * options_.antlever;

    Eigen::MatrixXd dz;
    dz = antenna_pos - gnssdata.pos;

    // ����GNSSλ�ù۲����
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos;
    H_gnsspos.resize(3, kf_->Q_.rows());
    H_gnsspos.setZero();
    H_gnsspos.block(0, P_ID, 3, 3) = Eigen::Matrix3d::Identity();
    H_gnsspos.block(0, PHI_ID, 3, 3) = Skew(R_be * options_.antlever);

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
    pvacur_.pos -= delta_r;

    // �ٶ�����
    vectemp = kf_->x_hat_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;

    // ��̬����
    vectemp = kf_->x_hat_.block(PHI_ID, 0, 3, 1);
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pvacur_.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    Matrix3d R_be = get_Rot(blh[0], blh[1]).transpose() * pvacur_.att.cbn.transpose();
    R_be = (Matrix3d::Identity(3, 3) + Skew(vectemp)) * R_be;
    pvacur_.att.cbn = R_be.transpose() * get_Rot(blh[0], blh[1]).transpose();
    pvacur_.att.qbn = C2q(pvacur_.att.cbn);
    pvacur_.att.euler = C2Euler(pvacur_.att.cbn);

    // IMU��ƫ����
    vectemp = kf_->x_hat_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = kf_->x_hat_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // ���״̬������ϵͳ״̬��,�����״̬����
    kf_->x_hat_.setZero();
}