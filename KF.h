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
    KalmanFilter(MatrixXd initial_state, MatrixXd initial_covariance);
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

    MatrixXd A_; // ״̬ת�ƾ���
    MatrixXd Q_; // ��������Э����
    MatrixXd H_; // �۲����
    MatrixXd R_; // ��������Э����
    MatrixXd x_hat_; // ״̬����
    MatrixXd P_;     // �������Э����
    MatrixXd x_hat_minus_; // Ԥ��״̬����
    MatrixXd P_minus_;     // Ԥ��������Э����
    MatrixXd K_; // ����������
    MatrixXd z_ ; //�۲�ֵ
};

MatrixXd getA(double T, Configure cfg);

MatrixXd getH(MatrixXd B);

MatrixXd getQ(double T, Configure cfg);

MatrixXd getR(double ROW);

MatrixXd getz(MatrixXd l_P,MatrixXd l_V);

int RTK_getA(vector<int> old_prn, map<int, int>& new_prn, MatrixXd& A, int old_ref);
