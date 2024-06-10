#pragma once
#include <vector>
#include <string>
#include <cmath>
#include "INS_types.h"

#include "cal.h"
#include "comm.h"

class INS_Eigen
{
public:
    explicit INS_Eigen(GINSOptions& options);

    ~INS_Eigen() = default;

    /**
     * @brief ����µ�IMU���ݣ�(��)����IMU���
     * @param [in] imu        �µ�IMUԭʼ����
     * @param [in] compensate �Ƿ񲹳�IMU���
     * */
    void addImuData(const IMU& imu, bool compensate = false) {

        imupre_ = imucur_;
        imucur_ = imu;

        if (compensate) {
            imuCompensate(imucur_);
        }
    }

    /**
     * @brief ����µ�GNSS����
     * @param [in] gnss �µ�GNSS����
     * */
    void addGnssData(const GNSS& gnss) {
        gnssdata_ = gnss;
    }

    /**
     * @brief �����µ�IMU����
     * */
    void newImuProcess();

    /**
     * @brief �ڲ�������ʽ��IMU���ݵ�ָ��ʱ��
     * @param [in]     imu1      ǰһʱ��IMU����
     * @param [in,out] imu2      ��ǰʱ��IMU����
     * @param [in]     timestamp �����ڲ嵽��ʱ��
     * @param [in,out] midimu    ����ڲ�ʱ�̵�IMU����
     * */
    static void imuInterpolate(const IMU& imu1, IMU& imu2, const double timestamp, IMU& midimu) {

        if (imu1.time > timestamp || imu2.time < timestamp) {
            return;
        }

        double lamda = (timestamp - imu1.time) / (imu2.time - imu1.time);

        midimu.time = timestamp;
        midimu.dtheta = imu2.dtheta * lamda;
        midimu.dvel = imu2.dvel * lamda;
        midimu.dt = timestamp - imu1.time;

        imu2.dtheta = imu2.dtheta - midimu.dtheta;
        imu2.dvel = imu2.dvel - midimu.dvel;
        imu2.dt = imu2.dt - midimu.dt;
    }

    /**
     * @brief ��ȡ��ǰʱ��
     * */
    double timestamp() const {
        return timestamp_;
    }

    /**
     * @brief ��ȡ��ǰIMU״̬
     * */
    NavState getNavState();

    /**
     * @brief ��ȡ��ǰ״̬Э����
     * */
    Eigen::MatrixXd getCovariance() {
        return kf_->P_;
    }

private:
    /**
     * @brief ��ǰIMU������IMU������
     * @param [in,out] imu ��Ҫ������IMU����
     * */
    void imuCompensate(IMU& imu);

    /**
     * @brief �ж��Ƿ���Ҫ����,�Լ�������һʱ��ϵͳ״̬
     * @param [in] imutime1   ��һIMU״̬ʱ��
     * @param [in] imutime2   ��ǰIMU״̬ʱ��
     * @param [in] updatetime ״̬���µ�ʱ��
     * @return 0: ����Ҫ����
     *         1: ��Ҫ������һIMU״̬
     *         2: ��Ҫ���µ�ǰIMU״̬
     *         3: ��Ҫ��IMU�����ڲ嵽״̬����ʱ��
     * */
    int isToUpdate(double imutime1, double imutime2, double updatetime) const;

    /**
     * @brief ����INS״̬����(IMU��е�����㷨), ������IMU״̬ת�ƾ����������
     * @param [in,out] imupre ǰһʱ��IMU����
     * @param [in,out] imucur ��ǰʱ��IMU����
     * */
    void insPropagation(IMU& imupre, IMU& imucur);

    /**
     * @brief ʹ��GNSSλ�ù۲����ϵͳ״̬
     * @param [in,out] gnssdata
     * */
    void gnssUpdate(GNSS& gnssdata);

    /**
     * @brief ʹ��ODO���ٸ��� 
     */
    void ODOUpdate();

    /**
     * @brief ʹ��ZUPT����
     */
    void ZUPTUpdate();

    /**
     * @brief �������״̬����ǰ״̬
     * */
    void stateFeedback();

    /**
     * @brief ���Э����Խ���Ԫ���Ƿ�Ϊ��
     * */
    void checkCov() {

        for (int i = 0; i < RANK; i++) {
            if (kf_->P_(i, i) < 0) {
                std::cout << "Covariance is negative at " << std::setprecision(10) << timestamp_ << " !" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

private:
    GINSOptions options_;

    double timestamp_;

    // ����ʱ�������IMU״̬�͹۲���Ϣ���С��������Ϊ���߶���
    const double TIME_ALIGN_ERR = 0.001;

    // IMU��GNSSԭʼ����
    IMU imupre_;
    IMU imucur_;
    GNSS gnssdata_;

    // IMU״̬��λ�á��ٶȡ���̬��IMU��
    PVA pvacur_;
    PVA pvapre_;
    ImuError imuerror_;

    // Kalman�˲����
    INS_KF* kf_;
    MatrixXd Qc;

    const int RANK = 15;
    const int NOISERANK = 12;

    // ״̬ID������ID
    enum State_ID { P_ID = 0, V_ID = 3, PHI_ID = 6, BG_ID = 9, BA_ID = 12};
    enum Noise_ID { VRW_ID = 0, ARW_ID = 3, BGSTD_ID = 6, BASTD_ID = 9};
};
