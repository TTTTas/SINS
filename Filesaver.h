#pragma once
#include <fstream>
#include <string>
#include <vector>

#include <absl/strings/str_format.h>

class FileSaver
{
public:
    FileSaver() = default;
    FileSaver(const std::string& filename, int columns) {
        open(filename, columns);
    }

    bool open(const std::string& filename, int columns) {
        filefp_.open(filename, std::ios_base::out);

        columns_ = columns;

        return filefp_.is_open();
    }

    void dump(const std::vector<double>& data) {
        dump_(data);
    }

    void dumpn(const std::vector<std::vector<double>>& data) {
        for (const auto& k : data) {
            dump_(k);
        }
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
    void dump_(const std::vector<double>& data)
    {
        std::string line;

        constexpr absl::string_view format = "%-15.9lf ";

        line = absl::StrFormat(format, data[0]);
        for (size_t k = 1; k < data.size(); k++) {
            absl::StrAppendFormat(&line, format, data[k]);
        }

        filefp_ << line << "\n";
    }

private:
    std::fstream filefp_;
    int columns_;
};



