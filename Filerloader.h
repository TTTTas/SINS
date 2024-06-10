#pragma once
#include <corecrt_io.h>
#include <iosfwd>
#include <string>

#include "Configure.h"
#include "file_load.h"
#include "INS_types.h"

class ImuFileLoader
{
public:
    ImuFileLoader() = delete;
    ImuFileLoader(INS_Configure cfg) {
        filefp_.open(cfg.Imu_path);

        dt_ = 1.0 / (double)cfg.Samp_rate;

        imu_.time = 0;
    }

    const IMU& next() {
        imu_pre_ = imu_;

        string line;
        getline(filefp_, line);
        imu_ = read_line_imu(line);

        double dt = imu_.time - imu_pre_.time;
        if (dt < 0.1) {
            imu_.dt = dt;
        }
        else {
            imu_.dt = dt_;
        }
        return imu_;
    }

    double starttime() {

        double starttime;
        std::streampos sp = filefp_.tellg();

        filefp_.seekg(0, std::ios_base::beg);
        string line;
        getline(filefp_, line);
        starttime = read_line_imu(line).time;
        filefp_.seekg(sp, std::ios_base::beg);
        return starttime;
    }

    double endtime() {

        double endtime = -1;
        std::streampos sp = filefp_.tellg();

        filefp_.seekg(-2, std::ios_base::end);
        char byte = 0;
        auto pos = filefp_.tellg();
        do {
            pos -= 1;
            filefp_.seekg(pos);
            filefp_.read(&byte, 1);
        } while (byte != '\n');
        string line;
        getline(filefp_, line);
        endtime = read_line_imu(line).time;
        filefp_.seekg(sp, std::ios_base::beg);
        return endtime;
    }

    void close() {
        filefp_.close();
    }

    bool isOpen() {
        return filefp_.is_open();
    }

    bool isEof() {
        return filefp_.eof();
    }

private:
    std::fstream filefp_;
    double dt_;
    IMU imu_, imu_pre_;
};

class GnssFileLoader
{
public:
    GnssFileLoader() = delete;
    explicit GnssFileLoader(const string& filename) {
        filefp_.open(filename);
    }

    const GNSS& next() {
        string line;
        getline(filefp_, line);
        if (line.empty())
        {
            gnss_.isvalid = false;
            return gnss_;
        }
        data = read_line_gnss(line);
        gnss_.time = data[0];
        gnss_.isvalid = true;
        if (data.size() == 7)
        {
            gnss_.blh << data[1], data[2], data[3];
            gnss_.std << data[4], data[4], data[5];
            gnss_.vel << 0, 0, 0;
            gnss_.vel_std << 0, 0, 0;
        }
        else
        {
            gnss_.blh << data[1], data[2], data[3];
            gnss_.vel << data[4], data[4], data[5];
            gnss_.std << data[7], data[8], data[9];
            gnss_.vel_std << data[10], data[11], data[12];
        }

        return gnss_;
    }

    void close() {
        filefp_.close();
    }

    bool isOpen() {
        return filefp_.is_open();
    }

    bool isEof() {
        return filefp_.eof();
    }

private:
    std::fstream filefp_;
    GNSS gnss_;
    vector<double> data;
};