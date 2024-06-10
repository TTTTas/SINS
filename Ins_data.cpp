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

    // 设置协方差矩阵，系统噪声阵和系统误差状态矩阵大小
    MatrixXd Cov_=MatrixXd::Zero(RANK, RANK);
    Qc = MatrixXd::Zero(NOISERANK, NOISERANK);

    // 初始化系统噪声阵
    auto imunoise = options_.imunoise;
    Qc.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    Qc.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    Qc.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    Qc.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();


    // 设置系统状态(位置、速度、姿态和IMU误差)初值和初始协方差
	// 初始化位置、速度、姿态
    pvacur_.pos = options.initstate.pos;
    pvacur_.vel = options.initstate.vel;
    pvacur_.att.euler = options.initstate.euler;
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pvacur_.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    pvacur_.att.cbe =  get_Rot(blh[0],blh[1]) * Euler2C(pvacur_.att.euler);
    pvacur_.att.qbe = C2q(pvacur_.att.cbe);
    // 初始化IMU误差
    imuerror_ = options.initstate.imuerror;

    // 给上一时刻状态赋同样的初值
    pvapre_ = pvacur_;

    // 初始化协方差
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
    // 当前IMU时间作为系统当前状态时间,
    timestamp_ = imucur_.time;

    // 如果GNSS有效，则将更新时间设置为GNSS时间
    double updatetime = gnssdata_.isvalid ? gnssdata_.time : -1;

    // 判断是否需要进行GNSS更新
    int res = isToUpdate(imupre_.time, imucur_.time, updatetime);

    if (res == 0) {
        // 只传播导航状态
        insPropagation(imupre_, imucur_);
        // 手动更新一下
        kf_->x_hat_ = kf_->x_hat_minus_;
        kf_->P_ = kf_->P_minus_;
    }
    else if (res == 1) {
        // GNSS数据靠近上一历元，先对上一历元进行GNSS更新
        // 手动预测
        kf_->x_hat_minus_ = kf_->x_hat_;
        kf_->P_minus_ = kf_->P_;
        gnssUpdate(gnssdata_);
        stateFeedback();

        pvapre_ = pvacur_;
        insPropagation(imupre_, imucur_);
    }
    else if (res == 2) {
        // GNSS数据靠近当前历元，先对当前IMU进行状态传播
        insPropagation(imupre_, imucur_);
        gnssUpdate(gnssdata_);
        stateFeedback();
    }
    else {
        // GNSS数据在两个IMU数据之间(不靠近任何一个), 将当前IMU内插到整秒时刻
        IMU midimu;
        imuInterpolate(imupre_, imucur_, updatetime, midimu);

        // 对前一半IMU进行状态传播
        insPropagation(imupre_, midimu);

        // 整秒时刻进行GNSS更新，并反馈系统状态
        gnssUpdate(gnssdata_);
        stateFeedback();

        // 对后一半IMU进行状态传播
        pvapre_ = pvacur_;
        insPropagation(midimu, imucur_);
    }

    // 检查协方差矩阵对角线元素
    checkCov();

    // 更新上一时刻的状态和IMU数据
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
    // 补偿IMU零偏
    imu.dtheta -= imuerror_.gyrbias * imu.dt;
    imu.dvel -= imuerror_.accbias * imu.dt;
}

int INS_Eigen::isToUpdate(double imutime1, double imutime2, double updatetime) const
{
    if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR) {
        // 更新时间靠近imutime1
        return 1;
    }
    else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR) {
        // 更新时间靠近imutime2
        return 2;
    }
    else if (imutime1 < updatetime && updatetime < imutime2) {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        return 3;
    }
    else {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        return 0;
    }
}

void INS_Eigen::insPropagation(IMU& imupre, IMU& imucur)
{
    // 对当前IMU数据(imucur)补偿误差, 上一IMU数据(imupre)已经补偿过了
    imuCompensate(imucur);
    // IMU状态更新(机械编排算法)
    insMech(pvapre_, pvacur_, imupre, imucur);

    // 系统噪声传播，姿态误差采用phi角误差模型
    Eigen::MatrixXd Phi, F, Qd, G;

    // 初始化Phi阵(状态转移矩阵)，F阵，Qd阵(传播噪声阵)，G阵(噪声驱动阵)
    Phi.resizeLike(kf_->P_);
    F.resizeLike(kf_->P_);
    Qd.resizeLike(kf_->P_);
    G.resize(RANK, NOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // 使用上一历元状态计算状态转移矩阵
    Eigen::Vector3d wie_e, wen_b;
    wie_e << 0, 0, OMEGA_E;
    wen_b << pvapre_.att.cbe.transpose() * wie_e;
    Eigen::Vector3d blh = get_BLH(XYZ2BLH(get_XYZ(pvapre_.pos), WGS84_e2, WGS84_a));
    blh[0] *= DEG2RAD;
    blh[1] *= DEG2RAD;
    Eigen::Vector3d gravity_e = g_e(blh);

    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
	accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;
    // 位置误差
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度误差
    F.block(V_ID, V_ID, 3, 3) = -2*Skew(wie_e) + pvapre_.att.cbe*options_.imunoise.acc_vrw.asDiagonal()* pvapre_.att.cbe.transpose();
    F.block(V_ID, PHI_ID, 3, 3) = Skew(pvapre_.att.cbe * accel);
    F.block(V_ID, BA_ID, 3, 3) = pvapre_.att.cbe;

    // 姿态误差
    F.block(PHI_ID, PHI_ID, 3, 3) = -Skew(wie_e)+ pvapre_.att.cbe * options_.imunoise.gyr_arw.asDiagonal() * pvapre_.att.cbe.transpose();
    F.block(PHI_ID, BG_ID, 3, 3) = -pvapre_.att.cbe;

    // 系统噪声驱动矩阵
    G.block(V_ID, VRW_ID, 3, 3) = pvapre_.att.cbe;
    G.block(PHI_ID, ARW_ID, 3, 3) = pvapre_.att.cbe;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 状态转移矩阵
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // 计算系统传播噪声
    Qd = G * Qc * G.transpose() * imucur.dt;
    Qd = (Phi * Qd * Phi.transpose() + Qd) * imucur.dt / 2;

    // EKF预测传播系统协方差和系统误差状态
    kf_->set_A(Phi);
    kf_->set_Q(Qd);
    kf_->predict();
}

void INS_Eigen::gnssUpdate(GNSS& gnssdata)
    // GNSS位置测量新息
{
    // IMU位置转到GNSS天线相位中心位置
    Eigen::Vector3d antenna_pos = pvacur_.pos + pvacur_.att.cbe * options_.antlever;

    Eigen::MatrixXd dz;
    dz = antenna_pos - gnssdata.blh;

    // 构造GNSS位置观测矩阵
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos;
    H_gnsspos.resize(3, kf_->Q_.rows());
    H_gnsspos.setZero();
    H_gnsspos.block(0, P_ID, 3, 3) = Eigen::Matrix3d::Identity();
    H_gnsspos.block(0, PHI_ID, 3, 3) = Skew(pvacur_.att.cbe * options_.antlever);

    // 位置观测噪声阵
    Eigen::MatrixXd R_gnsspos;
    R_gnsspos = gnssdata.std.cwiseProduct(gnssdata.std).asDiagonal();

    // EKF更新协方差和误差状态
    kf_->set_H(H_gnsspos);
    kf_->set_Z(dz);
    kf_->set_R(R_gnsspos);
    kf_->update();

    // GNSS更新之后设置为不可用
    gnssdata.isvalid = false;
}

void INS_Eigen::stateFeedback()
{
    Eigen::Vector3d vectemp;

    // 位置误差反馈
    Eigen::Vector3d delta_r = kf_->x_hat_.block(P_ID, 0, 3, 1);
    pvacur_.pos -= delta_r;

    // 速度误差反馈
    vectemp = kf_->x_hat_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;

    // 姿态误差反馈
    vectemp = kf_->x_hat_.block(PHI_ID, 0, 3, 1);
    Eigen::Quaterniond qpn = Phi2q(vectemp);
    pvacur_.att.qbe = qpn * pvacur_.att.qbe;
    pvacur_.att.cbe = q2C(pvacur_.att.qbe);
    pvacur_.att.euler = C2Euler(pvacur_.att.cbe);

    // IMU零偏误差反馈
    vectemp = kf_->x_hat_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = kf_->x_hat_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // 误差状态反馈到系统状态后,将误差状态清零
    kf_->x_hat_.setZero();
}