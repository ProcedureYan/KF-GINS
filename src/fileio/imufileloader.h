/*
 * OB_GINS: An Optimization-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Hailiang Tang
 *    Contact : thl@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef IMUFILELOADER_H
#define IMUFILELOADER_H

#include "common/types.h"
#include "fileloader.h"

class ImuFileLoader : public FileLoader {

public:
    ImuFileLoader() = delete;
    // 构造函数，接受文件名、列数和采样率（默认200 Hz）作为参数
    ImuFileLoader(const string &filename, int columns, int rate = 200) {
        open(filename, columns, FileLoader::TEXT);// 调用基类的 open 函数打开文件

        dt_ = 1.0 / (double) rate;// 计算采样周期 dt（单位为秒）

        imu_.time = 0;// 初始化 IMU 数据的时间戳
    }

    // 返回下一条 IMU 数据记录
    const IMU &next() {
        imu_pre_ = imu_; // 记录上一条 IMU 数据

        data_ = load(); // 加载新的一条数据

        imu_.time = data_[0]; // 将时间戳设置为数据的第一列值
        memcpy(imu_.dtheta.data(), &data_[1], 3 * sizeof(double)); // 复制角增量数据到 imu_.dtheta
        memcpy(imu_.dvel.data(), &data_[4], 3 * sizeof(double)); // 复制速度增量数据到 imu_.dvel

        double dt = imu_.time - imu_pre_.time; // 计算与上一次数据之间的时间差
        if (dt < 0.1) {
            imu_.dt = dt; // 如果时间差小于0.1秒，使用该时间差作为dt
        } else {
            imu_.dt = dt_; // 否则使用默认采样周期 dt_
        }

        // 增量形式，如果列数超过7，处理额外的增量速度数据
        if (columns_ > 7) {
            imu_.odovel = data_[7] * imu_.dt; // 计算车轮速度增量
        }

        return imu_;
    }

    // 返回文件中数据的开始时间
    double starttime() {

        double starttime;
        std::streampos sp = filefp_.tellg(); // 保存当前文件指针位置

        filefp_.seekg(0, std::ios_base::beg); // 移动文件指针到文件开头
        starttime = load().front(); // 读取第一条数据的时间
        filefp_.seekg(sp, std::ios_base::beg); // 恢复文件指针位置
        return starttime;
    }

    // 返回文件中数据的结束时间
    double endtime() {

        double endtime    = -1; // 初始化结束时间
        std::streampos sp = filefp_.tellg(); // 保存当前文件指针位置

        if (filetype_ == TEXT) {// 如果文件类型是文本格式
            filefp_.seekg(-2, std::ios_base::end); // 从文件末尾向前移动指针
            char byte = 0;
            auto pos  = filefp_.tellg();
            do {
                pos -= 1;
                filefp_.seekg(pos);
                filefp_.read(&byte, 1);
            } while (byte != '\n'); // 逐字节向前读取，直到找到换行符
        } else { // 如果文件类型是二进制格式
            filefp_.seekg(-columns_ * sizeof(double), std::ios_base::end); // 定位到最后一条记录
        }
        endtime = load().front(); // 读取最后一条数据的时间
        filefp_.seekg(sp, std::ios_base::beg); // 恢复文件指针位置
        return endtime;
    }

private:
    double dt_; // 采样周期

    IMU imu_, imu_pre_; // 当前和上一条 IMU 数据
    vector<double> data_; // 用于存储从文件加载的数据
};

#endif // IMUFILELOADER_H
