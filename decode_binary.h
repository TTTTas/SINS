#pragma once
#pragma region header
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <iterator>
#include "transform.h"
#include "Configure.h"
#include "GNSS_data.h"
using namespace std;

#define MAXRAWLEN 40960		 // ����ȡ���ݳ���
#define MAXNUM 8			 // ������
#define POLYCRC32 0xEDB88320 // CRC32У�������
#define OEM4SYNC1 0xAA		 /* oem7/6/4 message start sync code 1 */
#define OEM4SYNC2 0x44		 /* oem7/6/4 message start sync code 2 */
#define OEM4SYNC3 0x12		 /* oem7/6/4 message start sync code 3 */
#define OEM4HLEN 28			 /* oem7/6/4 message header length (bytes) */

/* message IDs */
#define ID_RANGE 43			 /* oem7/6/4 range measurement */
#define ID_GPSEPHEMERIS 7	 /* oem7/6 decoded gps ephemeris */
#define ID_BDSEPHEMERIS 1696 /* oem7/6 decoded bds ephemeris */

/*system IDs*/
#define SYS_GPS 0
#define SYS_BDS 4

#pragma endregion

#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t *)(p)))
static uint16_t U2(uint8_t *p)
{
	uint16_t u;
	memcpy(&u, p, 2);
	return u;
} // 2�ֽ�����
static uint32_t U4(uint8_t *p)
{
	uint32_t u;
	memcpy(&u, p, 4);
	return u;
} // 4�ֽ�����
static int32_t I4(uint8_t *p)
{
	int32_t i;
	memcpy(&i, p, 4);
	return i;
} // 4�ֽ�����
static float R4(uint8_t *p)
{
	float r;
	memcpy(&r, p, 4);
	return r;
} // 4�ֽ�float
static double R8(uint8_t *p)
{
	double r;
	memcpy(&r, p, 8);
	return r;
} // 8�ֽ�double
/*���ұ�ʶ��AA4412*/
int check_syn(unsigned char *buff, unsigned char data);
/*����У����*/
unsigned int CRC32(const unsigned char *buff, int len);
/*����У����*/
unsigned int UCRC32(const unsigned char *buff, int len);
/*����ϵͳ*/
unsigned short decode_SYS(unsigned char *buff);

/*
 * α��
 * α�ྫ��
 * �ز�
 * �ز�����
 * ������
 * �����
 * ����
 * ����
 * ͨ��״̬
 */
/*�������ǵ�Ƶ��۲�����*/
unsigned int decode_RANGE_STAT(unsigned char *buff, Satellate *sate);
/*����������������Ƶ��۲�����*/
unsigned int decode_RANGE(unsigned char *buff, int num, OBS_DATA *obs);
/*���뵥��GPS��������*/
unsigned int decode_GPSEPH_STAT(unsigned char *buff, EPHEMERIS *epoch);
/*���뵥��BDS��������*/
unsigned int decode_BDSEPH_STAT(unsigned char *buff, EPHEMERIS *bdse);

/*����*/
int decodestream(DATA_SET* result, unsigned char Buff[], int& d);

int decodefile(DATA_SET* result, GNSS_Configure cfg, const char* path);
