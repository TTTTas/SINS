#pragma once
#include "decode_binary.h"
#include "transform.h"
#include "Configure.h"
#include "GNSS_data.h"

using namespace std;

/*地球椭球相关参数*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926

/*code IDs*/
#define UNKOWN 0
/*GPS*/
#define CODE_L1C 1
#define CODE_L2P 2
#define CODE_L2W 3
#define CODE_L5Q 4
#define CODE_L1L 5
#define CODE_L2S 6
/*BDS*/
#define CODE_L2I 7
#define CODE_L7I 8
#define CODE_L6I 9
#define CODE_L1P 10
#define CODE_L5P 11

/*FREQUNCY*/ // MHz
/*GPS*/
#define Freq_L1 1575.42
#define Freq_L2 1227.60
#define Freq_L5 1176.45
/*BDS*/
#define Freq_B1 1561.098
#define Freq_B1_C 1575.42
#define Freq_B2 1207.14
#define Freq_B2_a 1176.45
#define Freq_B3 1268.52

/*Hopefiled*/
#define H0 0
#define T0 288.16
#define P0 1013.25
#define RH0 0.5

/*粗差探测阈值*/
#define GF_THRESH 0.0131
#define MW_THRESH 5

#define Method_LS 1
#define Method_KF 2

/*平方数*/
double SQR(double x);
/*模长*/
double Len(XYZ *pos);

double Len(XYZ* pos1, XYZ* pos2);

/*计算DOP值*/
double Cal_PDOP(MatrixXd Qxx);
/*最小二乘计算*/
unsigned int Cal_LEAST_SQR(MatrixXd B, MatrixXd l, MatrixXd P, MatrixXd *Qxx, MatrixXd &x, double *thegma, double *DOP);
/*解码频点*/
unsigned int decode_SYN(int sys, int signal);
/*解码频率*/
double CODE2FREQ(int code);

// 钟差改正
double CORRECT_CLK(double t, EPHEMERIS *eph);

// TGD计算
double TGD(EPHEMERIS *e, double f, int sys);

// 星历位置
unsigned int SAT_POS_CAL(double t, EPHEMERIS *eph, XYZ *xyz, double &clk, double dt, int SYS);

// 卫星高度角计算
double Ele_Angle(XYZ SatPos, XYZ RcvPos, int sys);

// Hopefiled对流层改正(m)
double Hopefield(double E, double H);

// Hopefiled对流层改正(m)
double Hopefield(XYZ SatPos, XYZ RcvPos, int sys);

/*Klobuchar模型改正电离层*/
double Klobuchar(XYZ RcvPos, double E, double A, double alpha[4], double beta[4], double UT, double code, int sys);

/*网口下位置解算*/
unsigned int LS_SPV(DATA_SET* data, GNSS_Configure cfg);

unsigned int setup_LS(DATA_SET* data, GNSS_Configure cfg, int sys);

unsigned int setup_KF(DATA_SET* data, GNSS_Configure cfg, int sys);

double get_measure(Satellate* sate, GNSS_Configure cfg, EPHEMERIS* eph, double &p);

double P_height1(double ELEV);

double P_height2(double ELEV);

// SPP单点定位KF
unsigned int KF_SPV(DATA_SET* data, double dt_e, GNSS_Configure cfg);

unsigned int Select_Common_Sates(DATA_SET* rove, DATA_SET* base, RTK_DATA* rtk, GNSS_Configure cfg);

int get_mat(MatrixXd& B, MatrixXd& L, MatrixXd& Q, GPSTIME* time, XYZ* RcvPos, vector<Satellate*>sates, int sys, GNSS_Configure cfg, RTK_DATA* rtk, vector<int>& PRN_used);

void RemoveRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);

void RemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

unsigned int RTK_Set(DATA_SET* base, RTK_DATA* rtk, GNSS_Configure cfg, int sys, int method);

unsigned int RTK_Solve(DATA_SET* base, RTK_DATA* rtk, GNSS_Configure cfg);

MatrixXd Cal_V(RTK_DATA* rtk, GNSS_Configure cfg, XYZ* rove_pos);