#include <iostream>
#include <Eigen/Dense>

#include "comm.h"
#include "Filerloader.h"
#include "Filesaver.h"
#include "GNSS_data.h"
#include "Ins_data.h"
#include <absl/time/clock.h>
#include "progressbar.hpp"

void writeNavResult(double time, NavState& navstate, FileSaver& navfile, FileSaver& imuerrfile);
void writeSTD(double time, Eigen::MatrixXd& cov, FileSaver& stdfile);

int main()
{
    INS_Configure cfg;
    auto ts = absl::Now();
    // imu�������ã����ݴ�������
    double starttime, endtime;
    starttime = cfg.start_time;
    endtime = cfg.end_time;

    // ����GNSS�ļ���IMU�ļ�
    GnssFileLoader gnssfile(cfg.GNSS_path);
    ImuFileLoader imufile(cfg);

    // ����GIEngine
    INS_Eigen giengine(cfg.gins_options);

    // ��������ļ�
    // navfile: gnssweek(1) + time(1) + pos(3) + vel(3) + euler angle(3) = 11
    // imuerrfile: time(1) + gyrbias(3) + accbias(3) + gyrscale(3) + accscale(3) = 13
    // stdfile: time(1) + pva_std(9) + imubias_std(6) + imuscale_std(6) = 22
    int nav_columns = 11, imuerr_columns = 13, std_columns = 22;
    FileSaver navfile(cfg.Out_Folder + "/Nav_result.nav", nav_columns);
    FileSaver imuerrfile(cfg.Out_Folder + "/IMU_ERR.txt", imuerr_columns);
    FileSaver stdfile(cfg.Out_Folder + "/STD.txt", std_columns);

    // ����ļ��Ƿ���ȷ��
    // check if these files are all opened
    if (!gnssfile.isOpen() || !imufile.isOpen() || !navfile.isOpen() || !imuerrfile.isOpen() || !stdfile.isOpen()) {
        std::cout << "Failed to open data file!" << std::endl;
        return -1;
    }

    // ��鴦��ʱ��
    // check process time
    if (endtime < 0) {
        endtime = imufile.endtime();
    }
    if (endtime > 604800 || starttime < imufile.starttime() || starttime > endtime) {
        std::cout << "Process time ERROR!" << std::endl;
        return -1;
    }

    // ���ݶ���
    IMU imu_cur;
    do {
        imu_cur = imufile.next();
    } while (imu_cur.time < starttime);

    GNSS gnss;
    do {
        gnss = gnssfile.next();
    } while (gnss.time <= starttime);

    // ���IMU���ݵ�GIEngine�У�����IMU���
    giengine.addImuData(imu_cur, true);

    // ���GNSS���ݵ�GIEngine
    giengine.addGnssData(gnss);

    // ���ڱ��洦����
    // used to save processing results
    double timestamp;
    NavState navstate;
    Eigen::MatrixXd cov;

    // ������ʾ�������
    // used to display processing progress
    int percent = 0, lastpercent = 0;
    double interval = endtime - starttime;
    // ��ʼ��������
    progressbar bar(interval);
    bar.set_todo_char(" ");
    bar.set_done_char("��");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");

    cout << "��ʼ������ϵ�������:\n";
    bar.update(0);
    while (true) {
        // ��ǰIMU״̬ʱ������GNSSʱ��ʱ����ȡ������µ�GNSS���ݵ�GIEngine
        if (gnss.time < imu_cur.time && !gnssfile.isEof()) {
            gnss = gnssfile.next();
            giengine.addGnssData(gnss);
        }

        // ��ȡ������µ�IMU���ݵ�GIEngine
        imu_cur = imufile.next();
        if (imu_cur.time > endtime || imufile.isEof()) {
            break;
        }
        giengine.addImuData(imu_cur);

        // �����µ�IMU����
        giengine.newImuProcess();

        // ��ȡ��ǰʱ�䣬IMU״̬��Э����
        timestamp = giengine.timestamp();
        navstate = giengine.getNavState();
        cov = giengine.getCovariance();

        // ���洦����
        writeNavResult(timestamp, navstate, navfile, imuerrfile);
        writeSTD(timestamp, cov, stdfile);

        // ��ʾ�����չ
        percent = (imu_cur.time - starttime) / interval * 100;
        if (percent - lastpercent >= 0.5) {
            bar.update(imu_cur.time - starttime);
            lastpercent = percent;
        }

    }

    // �رմ򿪵��ļ�
    imufile.close();
    gnssfile.close();
    navfile.close();
    imuerrfile.close();
    stdfile.close();

    // �������
    auto te = absl::Now();
    std::cout << std::endl << std::endl << "KF-GINS Process Finish! ";
    std::cout << "From " << starttime << " s to " << endtime << " s, total " << interval << " s!" << std::endl;
    std::cout << "Cost " << absl::ToDoubleSeconds(te - ts) << " s in total" << std::endl;

	std::system("pause");

	return 0;
}

/**
 * @brief ���浼�������IMU����ת��Ϊ���õ�λ
 * */
void writeNavResult(double time, NavState& navstate, FileSaver& navfile, FileSaver& imuerrfile) {

    std::vector<double> result;

    // ���浼�����
    result.clear();
    result.push_back(0);
    result.push_back(time);
    result.push_back(navstate.pos[0]);
    result.push_back(navstate.pos[1]);
    result.push_back(navstate.pos[2]);
    result.push_back(navstate.vel[0]);
    result.push_back(navstate.vel[1]);
    result.push_back(navstate.vel[2]);
    result.push_back(navstate.euler[0] * RAD2DEG);
    result.push_back(navstate.euler[1] * RAD2DEG);
    result.push_back(navstate.euler[2] * RAD2DEG);
    navfile.dump(result);

    // ����IMU���
    auto imuerr = navstate.imuerror;
    result.clear();
    result.push_back(time);
    result.push_back(imuerr.gyrbias[0] * RAD2DEG * 3600);
    result.push_back(imuerr.gyrbias[1] * RAD2DEG * 3600);
    result.push_back(imuerr.gyrbias[2] * RAD2DEG * 3600);
    result.push_back(imuerr.accbias[0] * 1e5);
    result.push_back(imuerr.accbias[1] * 1e5);
    result.push_back(imuerr.accbias[2] * 1e5);
    imuerrfile.dump(result);
}

/**
 * @brief �����׼���ת��Ϊ���õ�λ
 * */
void writeSTD(double time, Eigen::MatrixXd& cov, FileSaver& stdfile) {

    std::vector<double> result;

    result.clear();
    result.push_back(time);
    // ����λ�á��ٶȡ���̬��׼��
    for (int i = 0; i < 6; i++) {
        result.push_back(sqrt(cov(i, i)));
    }
    for (int i = 6; i < 9; i++) {
        result.push_back(sqrt(cov(i, i)) * RAD2DEG);
    }

    // ����IMU����׼��
    for (int i = 9; i < 12; i++) {
        result.push_back(sqrt(cov(i, i)) * RAD2DEG * 3600);
    }
    for (int i = 12; i < 15; i++) {
        result.push_back(sqrt(cov(i, i)) * 1e5);
    }
    stdfile.dump(result);
}