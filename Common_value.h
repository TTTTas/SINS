#pragma once
#define GPS_SAT_QUAN 32
#define BDS_SAT_QUAN 63

/*解算结果*/
#define UN_Solve 0
#define Success_Solve 1
#define Epoch_Loss -1
#define OBS_DATA_Loss -2
#define Set_UP_B_fail -3

#define MAXRAWLEN 40960		 // 最大读取数据长度
#define MAXNUM 8			 // 波段数
#define POLYCRC32 0xEDB88320 // CRC32校验码参数
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

/*地球椭球相关参数*/
#define WGS84_e2 0.0066943799013
#define WGS84_a 6378137.0
#define CGCS2000_e2 0.00669438002290
#define CGCS2000_a 6378137.0
#define Pi 3.1415926

/*大地常数*/
#define WGS84_GM 3.986005E+14
#define CGCS2000_GM 3.986004418E+14
#define omiga_earth 7.2921151467E-05
/*光速*/
#define velocity_c 2.99792458E8

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
#define L1 1575.42
#define L2 1227.60
#define L5 1176.45
/*BDS*/
#define B1 1561.098
#define B1_C 1575.42
#define B2 1207.14
#define B2_a 1176.45
#define B3 1268.52

/*Hopefiled*/
#define H0 0
#define T0 288.16
#define P0 1013.25
#define RH0 0.5

/*粗差探测阈值*/
#define GF_THRESH 0.05
#define MW_THRESH 10