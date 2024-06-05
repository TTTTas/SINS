#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <iomanip>
#include <io.h>
#include <direct.h>
#include <thread>
#include <mutex>

#include "transform.h"
#include "read.h"
#include "cal.h"
#include "sockets.h"
#include "KF.h"
#include "data.h"
#include "Configure.h"

using namespace std;
using namespace Eigen;

#define _CRT_SECURE_NO_WARNINGS
#define TIME_SPAN_THRESH 0.001

int main()
{
	/* 配置 */
	Configure CfgInfo;
	
	/* 结果存储变量 */
	DATA_SET* rove_data = new DATA_SET(CfgInfo);
	DATA_SET* base_data = new DATA_SET(CfgInfo);
	vector<OBS_DATA*> rove_datas;
	vector<OBS_DATA*> base_datas;
	RTK_DATA* RTK;
	/* 初始判断变量 */
	double dt_epoch = 1; // 文件流历元间时间差
	double temp_t = 0;
	FILE* DATA_Fobs;   // log文件指针
	FILE* DATA_base_Fobs;   // log文件指针
	FILE* Pos_Fobs;    // pos文件指针
	FILE* KF_Fobs;
	
	/* 网口输入数据相关变量 */
	int lenR, lenD;
	int lenR_base, lenD_base;
	unsigned char curbuff[MAXRAWLEN];
	unsigned char curbuff_base[MAXRAWLEN];
	lenD = 0;
	lenD_base = 0;
	unsigned char decBuff[2 * MAXRAWLEN];
	unsigned char decBuff_base[2 * MAXRAWLEN];

	SOCKET NetGps_1;
	SOCKET NetGps_2;
	int ReadFlag;
	int ReadFlag_Base;

	int choice = 0;

	printf("请选择输入方式\n1. 文件\t2. 网口\t3. RTK文件\t4. RTK实时网口\n");
	std::cin >> choice;
	
	/* 获取文件生成时间 */
	time_t nowtime;
	time(&nowtime); // 获取1970年1月1日0点0分0秒到现在经过的秒数
	tm p;
	localtime_s(&p, &nowtime); // 将秒数转换为本地时间,年从1900算起,需要+1900,月为0-11,所以要+1
	string filetime = to_string(p.tm_year + 1900) + "_" + to_string(p.tm_mon + 1) + "_" + to_string(p.tm_mday) + "_" + to_string(p.tm_hour) + "_" + to_string(p.tm_min) + "_" + to_string(p.tm_sec);
	string logpath = "D:\\GitHub\\Integrated-Navigation\\data\\logs\\";
	createDirectory(logpath);
	logpath += filetime + string(".log");
	string logpath_base = "D:\\GitHub\\Integrated-Navigation\\data\\logs_base\\";
	createDirectory(logpath_base);
	logpath_base += filetime + string(".log");
	string pospath = "D:\\GitHub\\Integrated-Navigation\\data\\Pos\\";
	string KFpath = "D:\\GitHub\\Integrated-Navigation\\data\\KF\\";
	createDirectory(pospath);
	createDirectory(KFpath);
	pospath += filetime + string(".pos");
	KFpath += filetime + string(".kf");
	CfgInfo.ObsDatFile_1 = logpath.c_str();
	CfgInfo.ObsDatFile_2 = logpath_base.c_str();
	CfgInfo.ResDatFile = pospath.c_str();
	CfgInfo.KFDatFile = KFpath.c_str();

	FILE* file;
	char load_filepath[200];
	char base_filepath[200];
	FILE* base_file;
	FILE* rove_file;

	bool base_flag = true;
	bool rove_flag = true;

	double judge = 1000000;
	switch (choice)
	{
	case 1:
		std::cout << "请输入文件路径" << endl;
		std::cin >> load_filepath;
		decodefile(rove_data, CfgInfo, load_filepath);
		break;

	case 2:
		if (OpenSocket(NetGps_1, CfgInfo.NetIP_1, CfgInfo.NetPort_1) == false)
		{
			printf("The ip %s was not opened\n", CfgInfo.NetIP_1);
			return 0;
		}

		if ((DATA_Fobs = fopen(CfgInfo.ObsDatFile_1, "wb")) == NULL)
		{
			printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile_1);
			exit(0);
		}

		if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
		{
			printf("The pos file %s was not opened\n", CfgInfo.ResDatFile);
			exit(0);
		}

		if ((KF_Fobs = fopen(CfgInfo.KFDatFile, "w")) == NULL)
		{
			printf("The kf file %s was not opened\n", "KF.pos");
			exit(0);
		}

		while (1)
		{
			Sleep(980);

			if ((lenR = recv(NetGps_1, (char*)curbuff, MAXRAWLEN, 0)) > 0) // 读取数据
			{
				printf("%5d\n", lenR);
				fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // 记录二进制数据流到文件中

				if ((lenD + lenR) > 2 * MAXRAWLEN)
					lenD = 0;

				memcpy(decBuff + lenD, curbuff, lenR); // 缓存拼接
				lenD += lenR;

				ReadFlag = decodestream(rove_data, decBuff, lenD); // 解码

				if (ReadFlag != 1)
				{
					printf("Data acquisition and decode failed \n");
				}
				else
				{
					if (rove_data->LS_first)
						dt_epoch = 1;
					else
						dt_epoch = rove_data->OBSTIME->SecOfWeek - temp_t;

					if (dt_epoch == 0)
						break;

					temp_t = rove_data->OBSTIME->SecOfWeek;
					rove_data->DetectOut(CfgInfo, dt_epoch);
					rove_data->Sate_pos_pre(rove_data->OBSTIME->SecOfWeek, CfgInfo);
					if (CfgInfo.SPP_LS_used)
					{
						if (LS_SPV(rove_data, CfgInfo))
							rove_data->LS_first = false;

						std::cout << "LS" << endl;
						rove_data->LS_print(CfgInfo);              // 输出至控制台
						rove_data->LS_Filewrite(Pos_Fobs, CfgInfo); // 输出至文件
					}

					if (CfgInfo.SPP_KF_used)
					{
						if (KF_SPV(rove_data, dt_epoch, CfgInfo))
							rove_data->KF_first = false;

						std::cout << "KF" << endl;
						rove_data->KF_Print(KF_Fobs, CfgInfo);
					}

					rove_data->reset();
				}
			}
			else
			{
				printf("NO MESSAGES IN!\n");
			}
		}

		break;

	case 3:
		base_data = new DATA_SET(CfgInfo);
		RTK = new RTK_DATA(CfgInfo);
		std::cout << "基站文件" << endl;
		std::cin >> base_filepath;
		std::cout << "流动站文件" << endl;
		std::cin >> load_filepath;
		if ((base_file = fopen(base_filepath, "rb")) == NULL)
		{
			printf("未能打开文件\n");
			return -1;
		}
		if ((rove_file = fopen(load_filepath, "rb")) == NULL)
		{
			printf("未能打开文件\n");
			return -1;
		}
		if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
		{
			printf("The pos file %s was not opened\n", CfgInfo.ResDatFile);
			exit(0);
		}
		if ((KF_Fobs = fopen(CfgInfo.KFDatFile, "w")) == NULL)
		{
			printf("The kf file %s was not opened\n", "KF.pos");
			exit(0);
		}
		while (!feof(rove_file)&&!feof(base_file))
		{
			if (rove_flag && (lenR = fread(curbuff, 1, 10000, rove_file)) > 0)
			{
				printf("%5d\n", lenR);

				if ((lenD + lenR) > 2 * MAXRAWLEN)
					lenD = 0;

				memcpy(decBuff + lenD, curbuff, lenR); // 缓存拼接
				lenD += lenR;

				ReadFlag = decodestream(rove_data, decBuff, lenD); // 解码

				if (ReadFlag != 1)
				{
					printf("Rove Data acquisition and decode failed \n");
				}
			}
			
			if (base_flag && (lenR_base = fread(curbuff_base, 1, 10000, base_file)) > 0)
			{
				printf("%5d\n", lenR_base);

				if ((lenD_base + lenR_base) > 2 * MAXRAWLEN)
					lenD_base = 0;

				memcpy(decBuff_base + lenD_base, curbuff_base, lenR_base); // 缓存拼接
				lenD_base += lenR_base;

				ReadFlag_Base = decodestream(base_data, decBuff_base, lenD_base); // 解码

				if (ReadFlag_Base != 1)
				{
					printf("Base Data acquisition and decode failed \n");
				}
			}
			
			judge = rove_data->range->OBS_TIME->SecOfWeek - base_data->range->OBS_TIME->SecOfWeek;
			if (judge > TIME_SPAN_THRESH)
			{
				base_flag = true;
				rove_flag = false;
				base_data->reset();
				continue;
			}
			else if (judge < -TIME_SPAN_THRESH)
			{
				base_flag = false;
				rove_flag = true;
				rove_data->reset();
				continue;
			}
			else if (ReadFlag_Base != 1 && ReadFlag != 1)
			{
				base_data->reset();
				rove_data->reset();
				continue;
			}
			else
			{
				base_flag = true;
				rove_flag = true;
			}
			std::cout << rove_data->range->OBS_TIME->SecOfWeek << "\t" << base_data->range->OBS_TIME->SecOfWeek << endl;
			if (rove_data->LS_first)
				dt_epoch = 1;
			else
				dt_epoch = rove_data->OBSTIME->SecOfWeek - temp_t;

			if (dt_epoch == 0)
				break;
			if(dt_epoch > 5)
			{
				MatrixXd temp_x = RTK->KF->x_hat_.topRows(3);
				RTK->KF->x_hat_ = temp_x;
				RTK->KF->P_ = MatrixXd::Identity(3, 3) * SQR(30);
				RTK->Fix_N = MatrixXd::Zero(0, 1);
			}
			temp_t = rove_data->OBSTIME->SecOfWeek;
			rove_data->Sate_pos_pre(rove_data->OBSTIME->SecOfWeek, CfgInfo);
			rove_data->DetectOut(CfgInfo, dt_epoch);
			base_data->Sate_pos_pre(base_data->OBSTIME->SecOfWeek, CfgInfo);
			base_data->DetectOut(CfgInfo, dt_epoch);
			if (LS_SPV(rove_data, CfgInfo))
				rove_data->LS_first = false;
			if (LS_SPV(base_data, CfgInfo))
				base_data->LS_first = false;
			if (rove_data->LS_result == Success_Solve && base_data->LS_result == Success_Solve)
			{
				Select_Common_Sates(rove_data, base_data, RTK, CfgInfo);
				RTK_Solve(base_data, RTK, CfgInfo);
				if (CfgInfo.RTK_LS_used)
					RTK->LS_Output(Pos_Fobs, CfgInfo);
				if (CfgInfo.RTK_KF_used)
				{
					RTK->KF_Output(KF_Fobs, CfgInfo);
					RTK->KF->reset();
				}
				std::cout << endl;
			}
			rove_data->reset();
			base_data->reset();
			RTK->reset();
		}
		break;

	case 4:
		base_data = new DATA_SET(CfgInfo);
		RTK = new RTK_DATA(CfgInfo);
		if (OpenSocket(NetGps_1, CfgInfo.NetIP_1, CfgInfo.NetPort_1) == false)
		{
			printf("The ip %s was not opened\n", CfgInfo.NetIP_1);
			return 0;
		}
		if (OpenSocket(NetGps_2, CfgInfo.NetIP_2, CfgInfo.NetPort_2) == false)
		{
			printf("The ip %s was not opened\n", CfgInfo.NetIP_2);
			return 0;
		}
		if ((DATA_Fobs = fopen(CfgInfo.ObsDatFile_1, "wb")) == NULL)
		{
			printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile_1);
			exit(0);
		}
		if ((DATA_base_Fobs = fopen(CfgInfo.ObsDatFile_2, "wb")) == NULL)
		{
			printf("The obs file %s was not opened\n", CfgInfo.ObsDatFile_2);
			exit(0);
		}
		if ((Pos_Fobs = fopen(CfgInfo.ResDatFile, "w")) == NULL)
		{
			printf("The pos file %s was not opened\n", CfgInfo.ResDatFile);
			exit(0);
		}
		if ((KF_Fobs = fopen(CfgInfo.KFDatFile, "w")) == NULL)
		{
			printf("The kf file %s was not opened\n", CfgInfo.KFDatFile);
			exit(0);
		}
		while (true)
		{
			Sleep(800);
			if ((lenR = recv(NetGps_1, (char*)curbuff, MAXRAWLEN, 0)) > 0) // 读取数据
			{
				printf("%5d\n", lenR);
				fwrite(curbuff, sizeof(unsigned char), lenR, DATA_Fobs); // 记录二进制数据流到文件中

				if ((lenD + lenR) > 2 * MAXRAWLEN)
					lenD = 0;

				memcpy(decBuff + lenD, curbuff, lenR); // 缓存拼接
				lenD += lenR;

				ReadFlag = decodestream(rove_data, decBuff, lenD); // 解码

				if (ReadFlag != 1)
				{
					printf("Data acquisition and decode failed \n");
				}
				else
				{
					rove_datas.push_back(rove_data->range);
					rove_data->range = new OBS_DATA();
				}
			}
			if ((lenR_base = recv(NetGps_2, (char*)curbuff_base, MAXRAWLEN, 0)) > 0)
			{
				printf("%5d\n", lenR_base);
				fwrite(curbuff_base, sizeof(unsigned char), lenR_base, DATA_base_Fobs); // 记录二进制数据流到文件中

				if ((lenD_base + lenR_base) > 2 * MAXRAWLEN)
					lenD_base = 0;

				memcpy(decBuff_base + lenD_base, curbuff_base, lenR_base); // 缓存拼接
				lenD_base += lenR_base;

				ReadFlag_Base = decodestream(base_data, decBuff_base, lenD_base); // 解码

				if (ReadFlag_Base != 1)
				{
					printf("Base Data acquisition and decode failed \n");
				}
				else
				{
					base_datas.push_back(base_data->range);
					base_data->range = new OBS_DATA();
				}
			}
			do
			{
				if (rove_datas.empty() || base_datas.empty())
					break;
				judge = rove_datas[0]->OBS_TIME->SecOfWeek - base_datas[0]->OBS_TIME->SecOfWeek;
				if (judge > TIME_SPAN_THRESH)
				{
					delete base_datas[0];
					base_datas.erase(base_datas.begin());
				}
				else if (judge < -TIME_SPAN_THRESH)
				{
					delete rove_datas[0];
					rove_datas.erase(rove_datas.begin());
				}
			}
			while (abs(judge)>TIME_SPAN_THRESH);
			if (rove_datas.empty() || base_datas.empty())
				continue;
			std::cout << "Rove: " << rove_datas.size() << "\t" << rove_datas[0]->OBS_TIME->SecOfWeek << "\nBase: " << base_datas.size() << "\t" << base_datas[0]->OBS_TIME->SecOfWeek << endl;
			*rove_data->range = *rove_datas[0];
			*rove_data->OBSTIME = *rove_datas[0]->OBS_TIME;
			*base_data->range = *base_datas[0];
			*base_data->OBSTIME = *base_datas[0]->OBS_TIME;
			if (rove_data->LS_first)
				dt_epoch = 1;
			else
				dt_epoch = rove_data->OBSTIME->SecOfWeek - temp_t;

			if (dt_epoch == 0)
				break;
			if (dt_epoch > 5)
			{
				MatrixXd temp_x = RTK->KF->x_hat_.topRows(3);
				RTK->KF->x_hat_ = temp_x;
				RTK->KF->P_ = MatrixXd::Identity(3, 3) * SQR(30);
				RTK->Fix_N = MatrixXd::Zero(0, 1);
			}
			temp_t = rove_data->OBSTIME->SecOfWeek;
			rove_data->DetectOut(CfgInfo, dt_epoch);
			rove_data->Sate_pos_pre(rove_data->OBSTIME->SecOfWeek, CfgInfo);
			base_data->DetectOut(CfgInfo, dt_epoch);
			base_data->Sate_pos_pre(base_data->OBSTIME->SecOfWeek, CfgInfo);
			if (LS_SPV(rove_data, CfgInfo))
				rove_data->LS_first = false;
			if (LS_SPV(base_data, CfgInfo))
				base_data->LS_first = false;
			if (rove_data->LS_result == Success_Solve && base_data->LS_result == Success_Solve)
			{
				Select_Common_Sates(rove_data, base_data, RTK, CfgInfo);
				RTK_Solve(base_data, RTK, CfgInfo);
				if (CfgInfo.RTK_LS_used)
					RTK->LS_Output(Pos_Fobs, CfgInfo);
				if (CfgInfo.RTK_KF_used)
				{
					RTK->KF_Output(KF_Fobs, CfgInfo);
					RTK->KF->reset();
				}
				std::cout << endl;
			}
			rove_data->reset();
			base_data->reset();
			delete rove_datas[0];
			rove_datas.erase(rove_datas.begin());
			delete base_datas[0];
			base_datas.erase(base_datas.begin());
			RTK->reset();
		}
		break;
	default:
		break;
	}
	std::system("pause");
	return 0;
}