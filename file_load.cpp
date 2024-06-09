#include "file_load.h"

#include <unordered_map>

#include "cal.h"
#include "progressbar.hpp"
#include <absl/strings/str_split.h>


IMU read_line_imu_txt(const std::string& line)
{
    vector<double>data_;
    IMU imu;
    std::stringstream lineStream(line);
    lineStream >> imu.time >> imu.dtheta[0] >> imu.dtheta[1] >> imu.dtheta[2];
    lineStream >> imu.dvel[0] >> imu.dvel[1] >> imu.dvel[2];
    imu.dt = 0;
	return imu;
}

IMU read_line_imu(const std::string& line)
{
    double g_x, g_y, g_z, a_x, a_y, a_z;
    g_x = g_y = g_z = a_x = a_y = a_z = 0.0;
    std::stringstream lineStream(line);
    std::string temp;

    std::string prefix;
    char comma, semicolon;
    int unknown;
    double timestamp1, timestamp2;


    // 忽略 "%RAWIMUSA,"
    lineStream.ignore(10);

    // 读取 prefix 和 timestamp1
    lineStream >> unknown >> comma >> timestamp1 >> semicolon;
    // 读取 prefix 和 timestamp2
    lineStream >> unknown >> comma >> timestamp2 >> comma;

    // 读取 a_x, a_y, a_z, g_x, g_y, g_z
    lineStream >> unknown >> comma >> a_z >> comma >> a_y >> comma >> a_x >> comma;
    lineStream >> g_z >> comma >> g_y >> comma >> g_x;

    IMU imu;
    imu.time = timestamp1;
    imu.dtheta << -g_y * GYR_SCALE, g_x* GYR_SCALE, -g_z * GYR_SCALE;
    imu.dvel << -a_y * ACC_SCALE, a_x* ACC_SCALE, -a_z * ACC_SCALE;
    imu.dt = 0;
    imu.odovel = 0;
    return imu;
}

int read_imu_asc(vector<IMU*>& Imu, INS_Configure cfg)
{
    std::ifstream inputFile(cfg.Imu_path);
    std::string line;
    int val = 0;
    if (!inputFile) {
        std::cerr << "无法打开文件!" << std::endl;
        return 1;
    }

    // 获取文件总大小（字节）
    inputFile.seekg(0, std::ios::end);
    std::streamsize totalSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // 初始化进度条
    progressbar bar(totalSize);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");
    std::streamsize bytesRead = 0;

    std::printf("读取文件:\n");
    std::streamsize bytes_last = 0;
    double rate = cfg.Samp_rate;
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue;
        }
        else
        {
            val = 1;
            IMU* imu = new IMU();
            *imu = read_line_imu(line);
            Imu.push_back(imu);
            Imu.back()->dt = 1 / rate;
            if(Imu.size()>1)
            {
                double dt_ = Imu.back()->time - Imu[Imu.size() - 2]->time;
                if (dt_ < 2 / rate)Imu.back()->dt = dt_;
            }

            // 更新已读取字节数
            bytesRead = inputFile.tellg();
            // 计算进度百分比并更新进度条
            if(bytesRead-bytes_last>0.005*totalSize)
            {
                bar.update(bytesRead);
                bytes_last = bytesRead;
            }
        }

    }
    return val;
}

vector<double> read_line_gnss(const std::string& line)
{
    std::vector<double> datas;
    string d;
    std::istringstream iss(line);
    while(iss>>d)
    {
        datas.push_back(stod(d));
    }

    return datas;
}

int read_GNSS(vector<GNSS*>& Gnss, INS_Configure cfg)
{
    std::ifstream inputFile(cfg.GNSS_path);
    std::string line;
    int val = 0;
    if (!inputFile) {
        std::cerr << "无法打开文件!" << std::endl;
        return 1;
    }

    // 获取文件总大小（字节）
    inputFile.seekg(0, std::ios::end);
    std::streamsize totalSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // 初始化进度条
    progressbar bar(totalSize);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");
    std::streamsize bytesRead = 0;

    std::printf("读取文件:\n");
    std::streamsize bytes_last = 0;
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue;
        }
        else
        {
            val = 1;
            GNSS* gnss = new GNSS();
            vector<double> data = read_line_gnss(line);
            gnss->time = data[0];
            if(data.size()==7)
            {
                gnss->blh << deg2rad(data[1]), deg2rad(data[2]), data[3];
                gnss->std << data[4], data[4], data[5];
                gnss->vel << 0, 0, 0;
                gnss->vel_std << 0, 0, 0;
                cfg.use_GNSS_vel = false;
            }
            else
            {
                gnss->blh << data[1], data[2], data[3];
                gnss->vel << data[4], data[4], data[5];
                gnss->std << data[7], data[8], data[9];
                gnss->vel_std << data[10], data[11], data[12];
            }
            Gnss.push_back(gnss);

            // 更新已读取字节数
            bytesRead = inputFile.tellg();
            // 计算进度百分比并更新进度条
            if (bytesRead - bytes_last > 0.005 * totalSize)
            {
                bar.update(bytesRead);
                bytes_last = bytesRead;
            }
        }

    }
    return val;
}

bool isAllSpaces(const std::string& str) {
    return std::all_of(str.begin(), str.end(), [](char c) {
        return std::isspace(static_cast<unsigned char>(c));
        });
}

Satellate* read_OBS_line(string line, vector<pair<string, int>> pairs, vector<string> signals)
{
    Satellate* sate = new Satellate();
    int sys = 0;
    if (line[0] == 'G')sys = SYS_GPS;
    else if (line[0] == 'C')sys = SYS_BDS;
    sate->SYS = sys;
    sate->PRN = stoi(line.substr(1, 2));
    sate->Phase_NUM = signals.size();
    for(int i=0;i<sate->Phase_NUM;i++)
    {
        int type = 0;
        if(sys==SYS_GPS)
        {
            switch (signals[i][0])
            {
            case '1':type = CODE_L1C; break;
            case '2':type = CODE_L2P; break;
            case '5':type = CODE_L5Q; break;
            default:type = UNKOWN; break;
            }
        }
        if(sys==SYS_BDS)
        {
            switch (signals[i][0])
            {
            case '1':type = CODE_L1P; break;
            case '2':type = CODE_L2I; break;
            case '5':type = CODE_L5P; break;
            case '6':type = CODE_L6I; break;
            case '7':type = CODE_L7I; break;
            default:type = UNKOWN; break;
            }
            break;
        }
        
        sate->SYG_TYPE[i] = type;
    }
    int index = 3;
    for(auto p:pairs)
    {
        double value = 0;
        string value_string = line.substr(index, 14);
        if (!isAllSpaces(value_string))value = stod(value_string);
        if (!value&&sate->SYG_TYPE[p.second]!=-1)
        {
            sate->Phase_NUM--;
            sate->SYG_TYPE[p.second] = -1;
        }
        switch (p.first[0])
        {
        case 'C':
            sate->PSERA[p.second] = value;
            break;
        case 'L':
            sate->PHASE[p.second] = value;
            break;
        case 'D':
            sate->DOPPLER[p.second] = value;
            break;
        case 'S':
            sate->SNR[p.second] = value;
            break;
        default:break;
        }
        index += 16;
    }


    return sate;
}


int read_OBS_Rnx(vector<OBS_DATA*>* obs, const char* path, XYZ* appro_pos)
{
    std::ifstream inputFile(path);
    std::string line;
    int val = 0;
    if (!inputFile) {
        std::cerr << "无法打开文件!" << std::endl;
        return val;
    }

    // 获取文件总大小（字节）
    inputFile.seekg(0, std::ios::end);
    std::streamsize totalSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // 初始化进度条
    progressbar bar(totalSize);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");
    std::streamsize bytesRead = 0;

    std::printf("读取Rinex观测值文件:\n");
    std::streamsize bytes_last = 0;

    do
    {
        std::getline(inputFile, line);
    } while (line.find("APPROX POSITION XYZ")==string::npos);

    
    appro_pos->X = stod(line.substr(0, 14));
    appro_pos->Y = stod(line.substr(14, 28));
    appro_pos->Z = stod(line.substr(28, 52));

    do
    {
        getline(inputFile, line);
    } while (line.find("SYS / # / OBS TYPES")==string::npos);


    vector<pair<string, int>> gps_pairs;
    vector<pair<string, int>> bds_pairs;
    vector<string> gps_signal;
    vector<string> bds_signal;
    while (line.find("SYS / # / OBS TYPES")!=string::npos)
    {
	    if(line[0]=='G')
	    {
            int index = 6;
            string signal_type;
	    	do
            {
                signal_type = line.substr(index + 1, 3);
                vector<string>::iterator it = find(gps_signal.begin(), gps_signal.end(), signal_type.substr(1, 2));
                if(it!=gps_signal.end())
                {
                    gps_pairs.push_back(pair<string, int>(signal_type, it - gps_signal.begin()));
                }
                else
                {
                    gps_signal.push_back(signal_type.substr(1, 2));
                    gps_pairs.push_back(pair<string, int>(signal_type, gps_signal.size() - 1));
                }
                index += 4;
            } while (line[index + 2] != ' ');
	    }
        else if (line[0] == 'C')
        {
            int index = 6;
            string signal_type;
            do
            {
                signal_type = line.substr(index + 1, 3);
                vector<string>::iterator it = find(bds_signal.begin(), bds_signal.end(), signal_type.substr(1, 2));
                if (it != bds_signal.end())
                {
                    bds_pairs.push_back(pair<string, int>(signal_type, it - bds_signal.begin()));
                }
                else
                {
                    bds_signal.push_back(signal_type.substr(1, 2));
                    bds_pairs.push_back(pair<string, int>(signal_type, bds_signal.size() - 1));
                }
                index += 4;
            } while (line[index + 2] != ' ');
        }

        getline(inputFile, line);
    }

    while (line.find("END OF HEADER") == string::npos)getline(inputFile, line);

    // 更新已读取字节数
    bytesRead = inputFile.tellg();
	bar.update(bytesRead);
	bytes_last = bytesRead;

    // 此行为文件头结束行
    while(std::getline(inputFile,line))
    {
        UTC time;
        int Sate_Num;
        time.Year = stoi(line.substr(2, 4));
        time.Month = stoi(line.substr(7, 2));
        time.Day = stoi(line.substr(10, 2));
        time.Hour = stoi(line.substr(13, 2));
        time.Min = stoi(line.substr(16, 2));
        time.Sec = int(stod(line.substr(19, 10)));
        Sate_Num = stoi(line.substr(32, 3));

        OBS_DATA* obs_data = new OBS_DATA();
        *obs_data->OBS_TIME = MJD2GPSTIME(UTC2MJD(time));
        for(int i=0;i<Sate_Num;i++)
        {
            getline(inputFile, line);
            if (line[0] == 'G')obs_data->GPS_SATE.push_back(read_OBS_line(line, gps_pairs, gps_signal));
            if (line[0] == 'C')obs_data->BDS_SATE.push_back(read_OBS_line(line, bds_pairs, bds_signal));
        }
        obs_data->Sate_Num = obs_data->GPS_SATE.size() + obs_data->BDS_SATE.size();
        obs->push_back(obs_data);
        // 更新已读取字节数
        bytesRead = inputFile.tellg();
        // 计算进度百分比并更新进度条
        if (bytesRead - bytes_last > 0.005 * totalSize)
        {
            bar.update(bytesRead);
            bytes_last = bytesRead;
        }
    }

    inputFile.close();
    val = 1;
    return val;
}

int read_EPHEMERIS_One(EPHEMERIS* eph, ifstream& fpr, int prn)
{
    string line;
    //第一行
    getline(fpr, line);
    UTC time;
    time.Year = stoi(line.substr(1, 4));
    time.Month = stoi(line.substr(6, 2));
    time.Day = stoi(line.substr(9, 2));
    time.Hour = stoi(line.substr(12, 2));
    time.Min = stoi(line.substr(15, 2));
    time.Sec = stoi(line.substr(18, 2));
    GPSTIME gpst = MJD2GPSTIME(UTC2MJD(time));
    if (eph->PRN != 0 && ((gpst.Week - eph->toe_wn) * 604800 + gpst.SecOfWeek - eph->toe_tow) < 0)
        return 1;
    eph->toe_wn = gpst.Week;
    eph->toe_tow = gpst.SecOfWeek;
    eph->PRN = prn;
    eph->a_f0 = stod(line.substr(23, 19));
    eph->a_f1 = stod(line.substr(42, 19));
    eph->a_f2 = stod(line.substr(61, 19));
    //第二行
    getline(fpr, line);
    eph->IODE1 = stod(line.substr(4, 19));
    eph->Crs = stod(line.substr(23, 19));
    eph->delt_n = stod(line.substr(42, 19));
    eph->M0 = stod(line.substr(61, 19));

    //第三行
    getline(fpr, line);
    eph->Cuc = stod(line.substr(4, 19));
    eph->e = stod(line.substr(23, 19));
    eph->Cus = stod(line.substr(42, 19));
    eph->sqrt_A = stod(line.substr(61, 19));

    //第四行
    getline(fpr, line);
    eph->toe_tow = stod(line.substr(4, 19));
    eph->Cic = stod(line.substr(23, 19));
    eph->Omiga0 = stod(line.substr(42, 19));
    eph->Cis = stod(line.substr(61, 19));

    //第五行
    getline(fpr, line);
    eph->i0 = stod(line.substr(4, 19));
    eph->Crc = stod(line.substr(23, 19));
    eph->omiga = stod(line.substr(42, 19));
    eph->dot_Omiga = stod(line.substr(61, 19));

    //第六行
    getline(fpr, line);
    eph->dot_i = stod(line.substr(4, 19));
    eph->toe_wn = stod(line.substr(42, 19));

    //第七行
    getline(fpr, line);
    eph->T_GD1 = stod(line.substr(42, 19));
    eph->IODC = stod(line.substr(61, 19));
    //第八行
    getline(fpr, line);

    return 1;
}

int read_EPHEMERIS_Rnx(EPHEMERIS** gps_eph, EPHEMERIS** bds_eph, const char* path)
{
    std::ifstream inputFile(path, ios::in || ios::binary);
    std::string line;
    int val = 0;
    if (!inputFile) {
        std::cerr << "无法打开文件!" << std::endl;
        return val;
    }

    // 获取文件总大小（字节）
    inputFile.seekg(0, std::ios::end);
    std::streamsize totalSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    // 初始化进度条
    progressbar bar(totalSize);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");
    std::streamsize bytesRead = 0;

    std::printf("读取Rinex星历文件:\n");
    std::streamsize bytes_last = 0;

    do
    {
        std::getline(inputFile, line);
    } while (line.find("END OF HEADER") == string::npos);

    char sys=0;
    int prn=0;
    char ch1, ch2;
    while(inputFile.peek() != EOF)
    {
        inputFile.get(sys);
        inputFile.get(ch1); inputFile.get(ch2);
        prn = (ch1 - '0') * 10 + ch2 - '0';
        switch (sys)
        {
        case 'G':read_EPHEMERIS_One(gps_eph[prn], inputFile, prn); break;
        case 'C':read_EPHEMERIS_One(bds_eph[prn], inputFile, prn); break;
        case 'S':for (int i = 0; i < 4; i++)getline(inputFile, line); break;
        case 'R':for (int i = 0; i < 4; i++)getline(inputFile, line); break;
        case 'E':for (int i = 0; i < 8; i++)getline(inputFile, line); break;
        case 'J':for (int i = 0; i < 8; i++)getline(inputFile, line); break;
        case 'I':for (int i = 0; i < 8; i++)getline(inputFile, line); break;
        default:break;
        }
        // 更新已读取字节数
        bytesRead = inputFile.tellg();
        // 计算进度百分比并更新进度条
        if (bytesRead - bytes_last > 0.005 * totalSize)
        {
            bar.update(bytesRead);
            bytes_last = bytesRead;
        }
    }
}