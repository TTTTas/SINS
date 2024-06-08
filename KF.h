#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <map>

#include "Configure.h"
#include <vector>
using namespace std;
using namespace Eigen;

#define MODE_SPP 1
#define MODE_RTK 2

class KalmanFilter
{
public:
    KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance)
        : x_hat_(initial_state), P_(initial_covariance) {}

    virtual void predict()
    {
        x_hat_minus_ = A_ * x_hat_;
        P_minus_ = A_ * P_ * A_.transpose() + Q_;
    }

    virtual void update()
    {
        MatrixXd S = H_ * P_minus_ * H_.transpose() + R_;
        K_ = P_minus_ * H_.transpose() * S.inverse();
        x_hat_ = x_hat_minus_ + K_ * (z_ - H_ * x_hat_minus_);
        P_ = (MatrixXd::Identity(P_.rows(), P_.cols()) - K_ * H_) * P_minus_;
    }

    virtual void set_A(MatrixXd A) { A_ = A; }
    virtual void set_H(MatrixXd H) { H_ = H; }
    virtual void set_Q(MatrixXd Q) { Q_ = Q; }
    virtual void set_R(MatrixXd R) { R_ = R; }
    virtual void set_Z(MatrixXd Z) { z_ = Z; }
    virtual void setState(MatrixXd state) { x_hat_ = state; }
    virtual void reset()
    {
        x_hat_ = MatrixXd::Zero(x_hat_.rows(), x_hat_.cols());
        P_ = MatrixXd::Zero(P_.rows(), P_.cols());
    }

    MatrixXd getState() const { return x_hat_; }
    MatrixXd getState_minus() const { return x_hat_minus_; }

    MatrixXd A_; // 状态转移矩阵
    MatrixXd Q_; // 过程噪声协方差
    MatrixXd H_; // 观测矩阵
    MatrixXd R_; // 测量噪声协方差
    MatrixXd x_hat_; // 状态估计
    MatrixXd P_;     // 估计误差协方差
    MatrixXd x_hat_minus_; // 预测状态估计
    MatrixXd P_minus_;     // 预测估计误差协方差
    MatrixXd K_; // 卡尔曼增益
    MatrixXd z_; //观测值
};

class GNSS_KF : public KalmanFilter
{
public:
    GNSS_KF(MatrixXd initial_state, MatrixXd initial_covariance);
    void predict();
    void update();
    void set_A(MatrixXd A);
    void set_H(MatrixXd H, int mode);
    void set_Q(MatrixXd Q);
    void set_R(MatrixXd R, int mode);
    void set_Z(MatrixXd Z);
    void setState(MatrixXd state);
    void reset();
    MatrixXd getState() const;
    MatrixXd getState_minus();

    MatrixXd A_; // 状态转移矩阵
    MatrixXd Q_; // 过程噪声协方差
    MatrixXd H_; // 观测矩阵
    MatrixXd R_; // 测量噪声协方差
    MatrixXd x_hat_; // 状态估计
    MatrixXd P_;     // 估计误差协方差
    MatrixXd x_hat_minus_; // 预测状态估计
    MatrixXd P_minus_;     // 预测估计误差协方差
    MatrixXd K_; // 卡尔曼增益
    MatrixXd z_ ; //观测值
};

MatrixXd getA(double T, GNSS_Configure cfg);

MatrixXd getH(MatrixXd B);

MatrixXd getQ(double T, GNSS_Configure cfg);

MatrixXd getR(double ROW);

MatrixXd getz(MatrixXd l_P,MatrixXd l_V);

int RTK_getA(vector<int> old_prn, map<int, int>& new_prn, MatrixXd& A, int old_ref);

class INS_KF:public KalmanFilter
{
public:
    INS_KF(MatrixXd initial_state, MatrixXd initial_covariance);
    void predict();
    void update();
    void set_A(MatrixXd A);
    void set_H(MatrixXd H);
    void set_Q(MatrixXd Q);
    void set_R(MatrixXd R);
    void set_Z(MatrixXd Z);
    void setState(MatrixXd state);
    void reset();
    MatrixXd getState() const;
    MatrixXd getState_minus();

    MatrixXd A_; // 状态转移矩阵
    MatrixXd Q_; // 过程噪声协方差
    MatrixXd H_; // 观测矩阵
    MatrixXd R_; // 测量噪声协方差
    MatrixXd x_hat_; // 状态估计
    MatrixXd P_;     // 估计误差协方差
    MatrixXd x_hat_minus_; // 预测状态估计
    MatrixXd P_minus_;     // 预测估计误差协方差
    MatrixXd K_; // 卡尔曼增益
    MatrixXd z_; //观测值
};