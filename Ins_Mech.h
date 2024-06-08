#pragma once
#include "comm.h"
#include "INS_types.h"
#include "Configure.h"
#include <Eigen/Dense>

/**
 * @brief  航向角初始化 
 * @param  [in]  imu_data 用于初始化的IMU数据
 * @param  [in]  配置表
 * @return [out] 初始航向欧拉角
 */
Eigen::Vector3d Init_Yaw(std::vector<IMU*> imu_data, INS_Configure cfg);

/**
 * @brief INS机械编排算法, 利用IMU数据进行速度、位置和姿态更新
 * @param [in]     pva_pre 上一时刻状态
 * @param [in,out] pva_cur 输出当前时刻状态
 * @param [in]     imu_pre, imu_cur imudata
 * */
void insMech(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
     * @breif 位置更新
     *        position update
     * */
void posUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
 * @breif 速度更新
 *        velocity update
 * */
void velUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
 * @breif 姿态更新
 *        attitude update
 * */
void attUpdate(const PVA& pvapre, PVA& pvacur, const IMU& imupre, const IMU& imucur);