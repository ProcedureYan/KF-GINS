/*
 * KF-GINS: An EKF-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Liqiang Wang
 *    Contact : wlq@whu.edu.cn
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

#include "common/earth.h"
#include "common/rotation.h"
#include <deque>

#include "gi_engine.h"
#include "insmech.h"

GIEngine::GIEngine(GINSOptions &options) {

    this->options_ = options;
    options_.print_options();
    timestamp_ = 0;

    // 设置协方差矩阵，系统噪声阵和系统误差状态矩阵大小
    // resize covariance matrix, system noise matrix, and system error state matrix
    Cov_.resize(RANK, RANK);
    Qc_.resize(NOISERANK, NOISERANK);
    dx_.resize(RANK, 1);
    Cov_.setZero();
    Qc_.setZero();
    dx_.setZero();

    // 初始化系统噪声阵
    // initialize noise matrix
    auto imunoise                   = options_.imunoise;
    Qc_.block(ARW_ID, ARW_ID, 3, 3) = imunoise.gyr_arw.cwiseProduct(imunoise.gyr_arw).asDiagonal();
    Qc_.block(VRW_ID, VRW_ID, 3, 3) = imunoise.acc_vrw.cwiseProduct(imunoise.acc_vrw).asDiagonal();
    Qc_.block(BGSTD_ID, BGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrbias_std.cwiseProduct(imunoise.gyrbias_std).asDiagonal();
    Qc_.block(BASTD_ID, BASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accbias_std.cwiseProduct(imunoise.accbias_std).asDiagonal();
    Qc_.block(SGSTD_ID, SGSTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.gyrscale_std.cwiseProduct(imunoise.gyrscale_std).asDiagonal();
    Qc_.block(SASTD_ID, SASTD_ID, 3, 3) =
        2 / imunoise.corr_time * imunoise.accscale_std.cwiseProduct(imunoise.accscale_std).asDiagonal();

    // 设置系统状态(位置、速度、姿态和IMU误差)初值和初始协方差
    // set initial state (position, velocity, attitude and IMU error) and covariance
    initialize(options_.initstate, options_.initstate_std);

    pvacur_.att.cbv = Rotation::euler2matrix(options_.odo_euler);
}

void GIEngine::initialize(const NavState &initstate, const NavState &initstate_std) {

    // 初始化位置、速度、姿态
    // initialize position, velocity and attitude
    pvacur_.pos       = initstate.pos;
    pvacur_.vel       = initstate.vel;
    pvacur_.att.euler = initstate.euler;
    pvacur_.att.cbn   = Rotation::euler2matrix(pvacur_.att.euler);
    pvacur_.att.qbn   = Rotation::euler2quaternion(pvacur_.att.euler);
    // 初始化IMU误差
    // initialize imu error
    imuerror_ = initstate.imuerror;

    // 给上一时刻状态赋同样的初值
    // set the same value to the previous state
    pvapre_ = pvacur_;

    // 初始化协方差
    // initialize covariance
    ImuError imuerror_std            = initstate_std.imuerror;
    Cov_.block(P_ID, P_ID, 3, 3)     = initstate_std.pos.cwiseProduct(initstate_std.pos).asDiagonal();
    Cov_.block(V_ID, V_ID, 3, 3)     = initstate_std.vel.cwiseProduct(initstate_std.vel).asDiagonal();
    Cov_.block(PHI_ID, PHI_ID, 3, 3) = initstate_std.euler.cwiseProduct(initstate_std.euler).asDiagonal();
    Cov_.block(BG_ID, BG_ID, 3, 3)   = imuerror_std.gyrbias.cwiseProduct(imuerror_std.gyrbias).asDiagonal();
    Cov_.block(BA_ID, BA_ID, 3, 3)   = imuerror_std.accbias.cwiseProduct(imuerror_std.accbias).asDiagonal();
    Cov_.block(SG_ID, SG_ID, 3, 3)   = imuerror_std.gyrscale.cwiseProduct(imuerror_std.gyrscale).asDiagonal();
    Cov_.block(SA_ID, SA_ID, 3, 3)   = imuerror_std.accscale.cwiseProduct(imuerror_std.accscale).asDiagonal();
}

void GIEngine::newImuProcess() {

    // 当前IMU时间作为系统当前状态时间,
    // set current IMU time as the current state time
    timestamp_ = imucur_.time;

    // 如果GNSS有效，则将更新时间设置为GNSS时间
    // set update time as the gnss time if gnssdata is valid
    double updatetime = gnssdata_.isvalid ? gnssdata_.time : -1;

    // 判断具体的更新过程
    // determine if we should do GNSS update
    int res = isToUpdate(imupre_.time, imucur_.time, updatetime);

    if (res == 0) {
        // 只传播导航状态（只有IMU更新）
        // only propagate navigation state
        insPropagation(imupre_, imucur_);
    } else if (res == 1) {
        // GNSS数据靠近上一历元，先对上一历元进行GNSS更新
        // gnssdata is near to the previous imudata, we should firstly do gnss update
        gnssUpdate(gnssdata_);
        stateFeedback();

        pvapre_ = pvacur_;
        insPropagation(imupre_, imucur_);
    } else if (res == 2) {
        // GNSS数据靠近当前历元，先对当前IMU进行状态传播
        // gnssdata is near current imudata, we should firstly propagate navigation state
        insPropagation(imupre_, imucur_);
        // if(if_add_odo == true){
        //     odo_ = GetOdoVel(odo_for_calculate, timestamp_);
        //     if(odo_.isvalid){
        //         ODONHCUpdate(odo_, imucur_);
        //     }
        //     if_add_odo = false;
        // }
        gnssUpdate(gnssdata_);
        stateFeedback();
    } else {
        // GNSS数据在两个IMU数据之间(不靠近任何一个), 将当前IMU内插到整秒时刻
        // gnssdata is between the two imudata, we interpolate current imudata to gnss time
        IMU midimu;
        imuInterpolate(imupre_, imucur_, updatetime, midimu);

        // 对前一半IMU进行状态传播
        // propagate navigation state for the first half imudata
        insPropagation(imupre_, midimu);
        // if(if_add_odo == true){
        //     odo_ = GetOdoVel(odo_for_calculate, midimu.time);
        //     if(odo_.isvalid){
        //         ODONHCUpdate(odo_, imucur_);
        //     }
        //     if_add_odo = false;
        // }
        // 整秒时刻进行GNSS更新，并反馈系统状态
        // do GNSS position update at the whole second and feedback system states
        gnssUpdate(gnssdata_);
        stateFeedback();

        // 对后一半IMU进行状态传播
        // propagate navigation state for the second half imudata
        pvapre_ = pvacur_;
        insPropagation(midimu, imucur_);
    }

    if(if_add_odo == true){
        odo_ = GetOdoVel(odo_for_calculate, timestamp_);
        if(odo_.isvalid){
            ODONHCUpdate(odo_, imucur_);
            stateFeedback();
        }
        if_add_odo = false;
    }

    // 检查协方差矩阵对角线元素
    // check diagonal elements of current covariance matrix
    checkCov();

    // 更新上一时刻的状态和IMU数据
    // update system state and imudata at the previous epoch
    pvapre_ = pvacur_;
    imupre_ = imucur_;
}

void GIEngine::imuCompensate(IMU &imu) {

    // 补偿IMU零偏
    // compensate the imu bias
    imu.dtheta -= imuerror_.gyrbias * imu.dt;
    imu.dvel -= imuerror_.accbias * imu.dt;

    // 补偿IMU比例因子
    // compensate the imu scale
    Eigen::Vector3d gyrscale, accscale;
    gyrscale   = Eigen::Vector3d::Ones() + imuerror_.gyrscale;
    accscale   = Eigen::Vector3d::Ones() + imuerror_.accscale;
    imu.dtheta = imu.dtheta.cwiseProduct(gyrscale.cwiseInverse());
    imu.dvel   = imu.dvel.cwiseProduct(accscale.cwiseInverse());
}

void GIEngine::insPropagation(IMU &imupre, IMU &imucur) {

    // 对当前IMU数据(imucur)补偿误差, 上一IMU数据(imupre)已经补偿过了
    // compensate imu error to 'imucur', 'imupre' has been compensated
    imuCompensate(imucur);
    // IMU状态更新(机械编排算法)
    // update imustate(mechanization)
    INSMech::insMech(pvapre_, pvacur_, imupre, imucur);

    // 系统噪声传播，姿态误差采用phi角误差模型
    // system noise propagate, phi-angle error model for attitude error
    Eigen::MatrixXd Phi, F, Qd, G;

    // 初始化Phi阵(状态转移矩阵)，F阵(Phi的微分形式)，Qd阵(传播噪声阵)，G阵(噪声驱动阵)
    // initialize Phi (state transition), F matrix, Qd(propagation noise) and G(noise driven) matrix
    Phi.resizeLike(Cov_);
    F.resizeLike(Cov_);
    Qd.resizeLike(Cov_);
    G.resize(RANK, NOISERANK);
    Phi.setIdentity();
    F.setZero();
    Qd.setZero();
    G.setZero();

    // 使用上一历元状态计算状态转移矩阵
    // compute state transition matrix using the previous state
    Eigen::Vector2d rmrn;
    Eigen::Vector3d wie_n, wen_n;
    double gravity;

    // 计算经度和纬度的子午圈曲率半径和卯酉圈曲率半径
    rmrn  = Earth::meridianPrimeVerticalRadius(pvapre_.pos[0]);
    // 计算地球重力加速度
    gravity = Earth::gravity(pvapre_.pos);
    // 地球自转角速度在导航系下的投影
    wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, - WGS84_WIE * sin(pvapre_.pos[0]);
    // 地球旋转角速度加上由于地球曲率和速度产生的角速度
    wen_n << pvapre_.vel[1] / (rmrn[1] + pvapre_.pos[2]), -pvapre_.vel[0] / (rmrn[0] + pvapre_.pos[2]),
        -pvapre_.vel[1] * tan(pvapre_.pos[0]) / (rmrn[1] + pvapre_.pos[2]);

    // 临时矩阵和向量，用于存储计算中的加速度和角速度
    Eigen::Matrix3d temp;
    Eigen::Vector3d accel, omega;
    double rmh, rnh;
    // 计算在当前高度的经度和纬度方向的曲率半径
    rmh   = rmrn[0] + pvapre_.pos[2];
    rnh   = rmrn[1] + pvapre_.pos[2];
    // 当前IMU的加速度和角速度
    accel = imucur.dvel / imucur.dt;
    omega = imucur.dtheta / imucur.dt;

    // 位置误差
    // position error
    temp.setZero();
    temp(0, 0)                = -pvapre_.vel[2] / rmh;
    temp(0, 2)                = pvapre_.vel[0] / rmh;
    temp(1, 0)                = pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                = -(pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                = pvapre_.vel[1] / rnh;
    // 填充状态矩阵F中的位置误差项，并初始化速度误差项
    F.block(P_ID, P_ID, 3, 3) = temp;
    F.block(P_ID, V_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 速度误差
    // velocity error
    temp.setZero();
    temp(0, 0) = -2 * pvapre_.vel[1] * WGS84_WIE * cos(pvapre_.pos[0]) / rmh -
                 pow(pvapre_.vel[1], 2) / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(0, 2) = pvapre_.vel[0] * pvapre_.vel[2] / rmh / rmh - pow(pvapre_.vel[1], 2) * tan(pvapre_.pos[0]) / rnh / rnh;
    temp(1, 0) = 2 * WGS84_WIE * (pvapre_.vel[0] * cos(pvapre_.pos[0]) - pvapre_.vel[2] * sin(pvapre_.pos[0])) / rmh +
                 pvapre_.vel[0] * pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(1, 2) = (pvapre_.vel[1] * pvapre_.vel[2] + pvapre_.vel[0] * pvapre_.vel[1] * tan(pvapre_.pos[0])) / rnh / rnh;
    temp(2, 0) = 2 * WGS84_WIE * pvapre_.vel[1] * sin(pvapre_.pos[0]) / rmh;
    temp(2, 2) = -pow(pvapre_.vel[1], 2) / rnh / rnh - pow(pvapre_.vel[0], 2) / rmh / rmh +
                 2 * gravity / (sqrt(rmrn[0] * rmrn[1]) + pvapre_.pos[2]);
    F.block(V_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 0)                  = pvapre_.vel[2] / rmh;
    temp(0, 1)                  = -2 * (WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh);
    temp(0, 2)                  = pvapre_.vel[0] / rmh;
    temp(1, 0)                  = 2 * WGS84_WIE * sin(pvapre_.pos[0]) + pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh;
    temp(1, 1)                  = (pvapre_.vel[2] + pvapre_.vel[0] * tan(pvapre_.pos[0])) / rnh;
    temp(1, 2)                  = 2 * WGS84_WIE * cos(pvapre_.pos[0]) + pvapre_.vel[1] / rnh;
    temp(2, 0)                  = -2 * pvapre_.vel[0] / rmh;
    temp(2, 1)                  = -2 * (WGS84_WIE * cos(pvapre_.pos(0)) + pvapre_.vel[1] / rnh);
    F.block(V_ID, V_ID, 3, 3)   = temp;
    F.block(V_ID, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvapre_.att.cbn * accel);
    F.block(V_ID, BA_ID, 3, 3)  = pvapre_.att.cbn;
    F.block(V_ID, SA_ID, 3, 3)  = pvapre_.att.cbn * (accel.asDiagonal());

    // 姿态误差
    // attitude error
    temp.setZero();
    temp(0, 0) = -WGS84_WIE * sin(pvapre_.pos[0]) / rmh;
    temp(0, 2) = pvapre_.vel[1] / rnh / rnh;
    temp(1, 2) = -pvapre_.vel[0] / rmh / rmh;
    temp(2, 0) = -WGS84_WIE * cos(pvapre_.pos[0]) / rmh - pvapre_.vel[1] / rmh / rnh / pow(cos(pvapre_.pos[0]), 2);
    temp(2, 2) = -pvapre_.vel[1] * tan(pvapre_.pos[0]) / rnh / rnh;
    F.block(PHI_ID, P_ID, 3, 3) = temp;
    temp.setZero();
    temp(0, 1)                    = 1 / rnh;
    temp(1, 0)                    = -1 / rmh;
    temp(2, 1)                    = -tan(pvapre_.pos[0]) / rnh;
    F.block(PHI_ID, V_ID, 3, 3)   = temp;
    F.block(PHI_ID, PHI_ID, 3, 3) = -Rotation::skewSymmetric(wie_n + wen_n);
    F.block(PHI_ID, BG_ID, 3, 3)  = -pvapre_.att.cbn;
    F.block(PHI_ID, SG_ID, 3, 3)  = -pvapre_.att.cbn * (omega.asDiagonal());

    // IMU零偏误差和比例因子误差，建模成一阶高斯-马尔科夫过程
    // imu bias error and scale error, modeled as the first-order Gauss-Markov process
    F.block(BG_ID, BG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(BA_ID, BA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SG_ID, SG_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();
    F.block(SA_ID, SA_ID, 3, 3) = -1 / options_.imunoise.corr_time * Eigen::Matrix3d::Identity();

    // 系统噪声驱动矩阵
    // system noise driven matrix
    G.block(V_ID, VRW_ID, 3, 3)    = pvapre_.att.cbn;
    G.block(PHI_ID, ARW_ID, 3, 3)  = pvapre_.att.cbn;
    G.block(BG_ID, BGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(BA_ID, BASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SG_ID, SGSTD_ID, 3, 3) = Eigen::Matrix3d::Identity();
    G.block(SA_ID, SASTD_ID, 3, 3) = Eigen::Matrix3d::Identity();

    // 状态转移矩阵
    // compute the state transition matrix
    Phi.setIdentity();
    Phi = Phi + F * imucur.dt;

    // 计算系统传播噪声
    // compute system propagation noise
    Qd = G * Qc_ * G.transpose() * imucur.dt;
    Qd = (Phi * Qd * Phi.transpose() + Qd) / 2;

    // EKF预测传播系统协方差和系统误差状态
    // do EKF predict to propagate covariance and error state
    EKFPredict(Phi, Qd);
}

void GIEngine::gnssUpdate(GNSS &gnssdata) {

    // IMU位置转到GNSS天线相位中心位置
    // convert IMU position to GNSS antenna phase center position
    Eigen::Vector3d antenna_pos; // 定义天线相位中心位置的向量
    Eigen::Matrix3d Dr, Dr_inv; // 定义旋转矩阵和其逆矩阵
    Dr_inv      = Earth::DRi(pvacur_.pos); // 计算当前位置的逆旋转矩阵
    Dr          = Earth::DR(pvacur_.pos);  // 计算当前位置的旋转矩阵
    // 计算天线相位中心位置，结合IMU的位置、旋转矩阵和天线杆臂
    antenna_pos = pvacur_.pos + Dr_inv * pvacur_.att.cbn * options_.antlever;

    std::cout << antenna_pos << std::endl;
    std::cout << gnssdata.blh<<std::endl;
    std::cout << std::endl;

    // GNSS位置测量新息
    // compute GNSS position innovation
    Eigen::MatrixXd dz; // 定义测量新息矩阵
    dz = Dr * (antenna_pos - gnssdata.blh); // 计算天线相位中心与GNSS测量位置之间的差异，并进行旋转

    // 构造GNSS位置观测矩阵
    // construct GNSS position measurement matrix
    Eigen::MatrixXd H_gnsspos; // 定义观测矩阵
    H_gnsspos.resize(3, Cov_.rows()); // 设置矩阵大小
    H_gnsspos.setZero(); // 初始化矩阵为零
    H_gnsspos.block(0, P_ID, 3, 3)   = Eigen::Matrix3d::Identity(); // 设置位置部分为单位矩阵
    // 设置旋转部分为天线杆臂的反对称矩阵
    H_gnsspos.block(0, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvacur_.att.cbn * options_.antlever);

    // 位置观测噪声阵
    // construct measurement noise matrix
    Eigen::MatrixXd R_gnsspos; // 定义噪声矩阵
    R_gnsspos = gnssdata.std.cwiseProduct(gnssdata.std).asDiagonal(); // 根据GNSS数据的标准差构造对角噪声矩阵

    // EKF更新协方差和误差状态
    // do EKF update to update covariance and error state
    EKFUpdate(dz, H_gnsspos, R_gnsspos);

    // GNSS更新之后设置为不可用
    // Set GNSS invalid after update
    gnssdata.isvalid = false;
}

int GIEngine::isToUpdate(double imutime1, double imutime2, double updatetime) const {

    if (abs(imutime1 - updatetime) < TIME_ALIGN_ERR) {
        // 更新时间靠近imutime1
        // updatetime is near to imutime1
        return 1;
    } else if (abs(imutime2 - updatetime) <= TIME_ALIGN_ERR) {
        // 更新时间靠近imutime2
        // updatetime is near to imutime2
        return 2;
    } else if (imutime1 < updatetime && updatetime < imutime2) {
        // 更新时间在imutime1和imutime2之间, 但不靠近任何一个
        // updatetime is between imutime1 and imutime2, but not near to either
        return 3;
    } else {
        // 更新时间不在imutimt1和imutime2之间，且不靠近任何一个
        // updatetime is not bewteen imutime1 and imutime2, and not near to either.
        return 0;
    }
}

void GIEngine::EKFPredict(Eigen::MatrixXd &Phi, Eigen::MatrixXd &Qd) {

    assert(Phi.rows() == Cov_.rows());
    assert(Qd.rows() == Cov_.rows());

    // 传播系统协方差和误差状态
    // propagate system covariance and error state
    Cov_ = Phi * Cov_ * Phi.transpose() + Qd;
    dx_  = Phi * dx_;
}

void GIEngine::EKFUpdate(Eigen::MatrixXd &dz, Eigen::MatrixXd &H, Eigen::MatrixXd &R) {

    // 断言检查，确保观测矩阵和协方差矩阵的维度匹配
    assert(H.cols() == Cov_.rows()); // H的列数应等于Cov_的行数(观测矩阵 H 的列数等于状态协方差矩阵 Cov_ 的行数)
    assert(dz.rows() == H.rows()); // dz的行数应等于H的行数(观测差 dz 的行数与观测矩阵 H 的行数相等)
    assert(dz.rows() == R.rows()); // dz的行数应等于R的行数（观测差 dz 的行数与观测噪声协方差矩阵 R 的行数相等）
    assert(dz.cols() == 1); // dz应为列向量

    // 计算Kalman增益
    // Compute Kalman Gain
    auto temp         = H * Cov_ * H.transpose() + R;
    Eigen::MatrixXd K = Cov_ * H.transpose() * temp.inverse();

    // 更新系统误差状态和协方差
    // update system error state and covariance
    Eigen::MatrixXd I;
    I.resizeLike(Cov_);
    I.setIdentity();
    I = I - K * H; // 更新I以反映当前的卡尔曼增益

    // 如果每次更新后都进行状态反馈，则更新前dx_一直为0，下式可以简化为：dx_ = K * dz;
    // if state feedback is performed after every update, dx_ is always zero before the update
    // the following formula can be simplified as : dx_ = K * dz;
    
    dx_  = dx_ + K * (dz - H * dx_); // 更新误差状态dx_

    Cov_ = I * Cov_ * I.transpose() + K * R * K.transpose(); // 更新协方差矩阵
}

void GIEngine::stateFeedback() {

    Eigen::Vector3d vectemp;

    // 位置误差反馈
    // posisiton error feedback
    Eigen::Vector3d delta_r = dx_.block(P_ID, 0, 3, 1);
    Eigen::Matrix3d Dr_inv  = Earth::DRi(pvacur_.pos);
    pvacur_.pos -= Dr_inv * delta_r;

    // 速度误差反馈
    // velocity error feedback
    vectemp = dx_.block(V_ID, 0, 3, 1);
    pvacur_.vel -= vectemp;

    // std::cout << pvacur_.vel<<std::endl;
    // std::cout << std::endl;

    // 姿态误差反馈
    // attitude error feedback
    vectemp                = dx_.block(PHI_ID, 0, 3, 1);
    Eigen::Quaterniond qpn = Rotation::rotvec2quaternion(vectemp);
    pvacur_.att.qbn        = qpn * pvacur_.att.qbn;
    pvacur_.att.cbn        = Rotation::quaternion2matrix(pvacur_.att.qbn);
    pvacur_.att.euler      = Rotation::matrix2euler(pvacur_.att.cbn);

    // IMU零偏误差反馈
    // IMU bias error feedback
    vectemp = dx_.block(BG_ID, 0, 3, 1);
    imuerror_.gyrbias += vectemp;
    vectemp = dx_.block(BA_ID, 0, 3, 1);
    imuerror_.accbias += vectemp;

    // IMU比例因子误差反馈
    // IMU sacle error feedback
    vectemp = dx_.block(SG_ID, 0, 3, 1);
    imuerror_.gyrscale += vectemp;
    vectemp = dx_.block(SA_ID, 0, 3, 1);
    imuerror_.accscale += vectemp;

    // 误差状态反馈到系统状态后,将误差状态清零
    // set 'dx' to zero after feedback error state to system state
    dx_.setZero();
}

NavState GIEngine::getNavState() {

    NavState state;

    state.pos      = pvacur_.pos;
    state.vel      = pvacur_.vel;
    state.euler    = pvacur_.att.euler;
    state.imuerror = imuerror_;

    return state;
}

ODO GIEngine::GetOdoVel(std::deque<ODO> &odoraw, double time){
    ODO tmp_odo;
    tmp_odo.odo_vel = 0;
    tmp_odo.isvalid = false;
    if (odoraw.size() < 5){
        std::cout << "WARN: too little data to get odovel!" << std::endl;
        return tmp_odo;
    }

    double firstvel = 0;
    double secondvel = 0;
    double firsttime = 0;
    double secondtime = 0;
    int firstindex = 0;
    int secondindex = 0;

    // Loop through the rows of odoraw
    for (int i = 0; i < odoraw.size(); ++i) {
        if (odoraw[i].time <= time) {
            // Update firstvel and firsttime
            firstvel = (firstvel * firstindex + odoraw[i].odo_vel) / (firstindex + 1);
            firsttime = (firsttime * firstindex + odoraw[i].time) / (firstindex + 1);
            firstindex++;
        } else {
            // Update secondvel and secondtime
            secondvel = (secondvel * secondindex + odoraw[i].odo_vel) / (secondindex + 1);
            secondtime = (secondtime * secondindex + odoraw[i].time) / (secondindex + 1);
            secondindex++;
        }
    }

    // Compute odovel
    tmp_odo.odo_vel = firstvel + (time - firsttime) / (secondtime - firsttime) * (secondvel - firstvel);
    tmp_odo.isvalid = true;

    return tmp_odo;
}

// Matlab 版本的navstate其中一部分就是gi_engine类的成员变量pvacur_
// kf.P 和 kf.x分别是gi_engine类成员变量Cov_和dx_
void GIEngine::ODONHCUpdate(ODO odocur, IMU imucur){
    Eigen::Vector3d odonhc_vel(odocur.odo_vel, 0, 0);
    Eigen::Vector3d wib_b, wie_n, wen_n, win_n, wnb_b, vel_pre;

    wib_b = imucur.dtheta / imucur.dt;
    // 地球自转角速度在导航系下的投影
    wie_n << WGS84_WIE * cos(pvapre_.pos[0]), 0, - WGS84_WIE * sin(pvapre_.pos[0]);
    // 地球旋转角速度加上由于地球曲率和速度产生的角速度
    Eigen::Vector2d rmrn = Earth::meridianPrimeVerticalRadius(pvapre_.pos[0]);
    wen_n << pvapre_.vel[1] / (rmrn[1] + pvapre_.pos[2]), -pvapre_.vel[0] / (rmrn[0] + pvapre_.pos[2]),
        -pvapre_.vel[1] * tan(pvapre_.pos[0]) / (rmrn[1] + pvapre_.pos[2]);
    win_n = wie_n + wen_n;
    wnb_b = wib_b - pvacur_.att.cbn.transpose() * win_n;
    vel_pre = pvacur_.att.cbv * (pvacur_.att.cbn.transpose() * pvapre_.vel + Rotation::skewSymmetric(wnb_b) * options_.odolever); // 注意：“cbv” 和 “odolever”
    
    // std::cout << pvapre_.vel <<std::endl;
    // std::cout << odonhc_vel <<std::endl;
    // std::cout << std::endl;

    Eigen::MatrixXd dz_v = vel_pre + odonhc_vel;
    Eigen::MatrixXd dz =pvacur_.att.cbn * pvacur_.att.cbv.transpose() * dz_v;

    // std::cout << dz <<std::endl;
    // std::cout << std::endl;
    // dz.setZero();
    /*measurement equation and noise*/

    /*TODO: add measurement equation and noise matrix here!!*/

    // 构造车轮里程计速度观测矩阵
    Eigen::MatrixXd H_odovel; // 定义观测矩阵
    H_odovel.resize(3, Cov_.rows()); // 设置矩阵大小
    H_odovel.setZero(); // 初始化矩阵为零
    H_odovel.block(0, V_ID, 3 , 3) = Eigen::Matrix3d::Identity();
    // 设置旋转部分为轮速计杆臂的反对称矩阵
    H_odovel.block(0, PHI_ID, 3, 3) = Rotation::skewSymmetric(pvacur_.att.cbn * options_.odolever); 
    
    // 速度观测噪声阵
    // construct measurement noise matrix
    Eigen::Matrix3d R_odovel; // 定义噪声矩阵
    R_odovel.resize(3, 3);
    R_odovel.setZero();
    R_odovel(0, 0) = 0.1;
    R_odovel(1, 1) = 0.1;
    R_odovel(2, 2) = 0.1;

    auto temp         = H_odovel * Cov_ * H_odovel.transpose() + R_odovel;
    Eigen::MatrixXd K = Cov_ * H_odovel.transpose() * temp.inverse();

    Eigen::MatrixXd I;
    I.resizeLike(Cov_);
    I.setIdentity();
    I = I - K * H_odovel; // 更新I以反映当前的卡尔曼增益
    dx_  = dx_ + K * (dz - H_odovel * dx_); // 更新误差状态dx_
    Cov_ = I * Cov_ * I.transpose() + K * R_odovel * K.transpose(); // 更新协方差矩阵
}