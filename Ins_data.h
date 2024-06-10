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
     * @brief 添加新的IMU数据，(不)补偿IMU误差
     * @param [in] imu        新的IMU原始数据
     * @param [in] compensate 是否补偿IMU误差
     * */
    void addImuData(const IMU& imu, bool compensate = false) {

        imupre_ = imucur_;
        imucur_ = imu;

        if (compensate) {
            imuCompensate(imucur_);
        }
    }

    /**
     * @brief 添加新的GNSS数据
     * @param [in] gnss 新的GNSS数据
     * */
    void addGnssData(const GNSS& gnss) {
        gnssdata_ = gnss;
    }

    /**
     * @brief 处理新的IMU数据
     * */
    void newImuProcess();

    /**
     * @brief 内插增量形式的IMU数据到指定时刻
     * @param [in]     imu1      前一时刻IMU数据
     * @param [in,out] imu2      当前时刻IMU数据
     * @param [in]     timestamp 给定内插到的时刻
     * @param [in,out] midimu    输出内插时刻的IMU数据
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
     * @brief 获取当前时间
     * */
    double timestamp() const {
        return timestamp_;
    }

    /**
     * @brief 获取当前IMU状态
     * */
    NavState getNavState();

    /**
     * @brief 获取当前状态协方差
     * */
    Eigen::MatrixXd getCovariance() {
        return kf_->P_;
    }

private:
    /**
     * @brief 当前IMU误差补偿到IMU数据中
     * @param [in,out] imu 需要补偿的IMU数据
     * */
    void imuCompensate(IMU& imu);

    /**
     * @brief 判断是否需要更新,以及更新哪一时刻系统状态
     * @param [in] imutime1   上一IMU状态时间
     * @param [in] imutime2   当前IMU状态时间
     * @param [in] updatetime 状态更新的时间
     * @return 0: 不需要更新
     *         1: 需要更新上一IMU状态
     *         2: 需要更新当前IMU状态
     *         3: 需要将IMU进行内插到状态更新时间
     * */
    int isToUpdate(double imutime1, double imutime2, double updatetime) const;

    /**
     * @brief 进行INS状态更新(IMU机械编排算法), 并计算IMU状态转移矩阵和噪声阵
     * @param [in,out] imupre 前一时刻IMU数据
     * @param [in,out] imucur 当前时刻IMU数据
     * */
    void insPropagation(IMU& imupre, IMU& imucur);

    /**
     * @brief 使用GNSS位置观测更新系统状态
     * @param [in,out] gnssdata
     * */
    void gnssUpdate(GNSS& gnssdata);

    /**
     * @brief 使用ODO轮速更新 
     */
    void ODOUpdate();

    /**
     * @brief 使用ZUPT修正
     */
    void ZUPTUpdate();

    /**
     * @brief 反馈误差状态到当前状态
     * */
    void stateFeedback();

    /**
     * @brief 检查协方差对角线元素是否都为正
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

    // 更新时间对齐误差，IMU状态和观测信息误差小于它则认为两者对齐
    const double TIME_ALIGN_ERR = 0.001;

    // IMU和GNSS原始数据
    IMU imupre_;
    IMU imucur_;
    GNSS gnssdata_;

    // IMU状态（位置、速度、姿态和IMU误差）
    PVA pvacur_;
    PVA pvapre_;
    ImuError imuerror_;

    // Kalman滤波相关
    INS_KF* kf_;
    MatrixXd Qc;

    const int RANK = 15;
    const int NOISERANK = 12;

    // 状态ID和噪声ID
    enum State_ID { P_ID = 0, V_ID = 3, PHI_ID = 6, BG_ID = 9, BA_ID = 12};
    enum Noise_ID { VRW_ID = 0, ARW_ID = 3, BGSTD_ID = 6, BASTD_ID = 9};
};
