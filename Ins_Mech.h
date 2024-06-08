#pragma once
#include "comm.h"
#include "INS_types.h"
#include "Configure.h"
#include <Eigen/Dense>

/**
 * @brief  ����ǳ�ʼ�� 
 * @param  [in]  imu_data ���ڳ�ʼ����IMU����
 * @param  [in]  ���ñ�
 * @return [out] ��ʼ����ŷ����
 */
Eigen::Vector3d Init_Yaw(std::vector<IMU*> imu_data, INS_Configure cfg);

/**
 * @brief INS��е�����㷨, ����IMU���ݽ����ٶȡ�λ�ú���̬����
 * @param [in]     pva_pre ��һʱ��״̬
 * @param [in,out] pva_cur �����ǰʱ��״̬
 * @param [in]     imu_pre, imu_cur imudata
 * */
void insMech(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
     * @breif λ�ø���
     *        position update
     * */
void posUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
 * @breif �ٶȸ���
 *        velocity update
 * */
void velUpdate(const PVA& pva_pre, PVA& pva_cur, const IMU& imu_pre, const IMU& imu_cur);

/**
 * @breif ��̬����
 *        attitude update
 * */
void attUpdate(const PVA& pvapre, PVA& pvacur, const IMU& imupre, const IMU& imucur);