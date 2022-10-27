/**
 * @file time_evolution.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 一次时间步长的更新，包括密度、速度、质量总能的更新和时间步长的选取
 * @version v1.01
 * @date 2022-10-07
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

/**
 * @brief Get the density
 * 
 * @param xi 节点的参考坐标横坐标
 * @param eta 节点的参考坐标纵坐标
 * @param r 节点所在单元的编号
 * @return double 
 */
double get_density(double xi, double eta, int r);

/**
 * @brief 得到单元速度的右端矩阵
 * 
 * @param i 单元的行编号
 * @param j 单元的列编号
 * @return double** 
 */
double ** velocity_matrix(int i, int j);

/**
 * @brief 得到单元质量总能的右端矩阵
 * 
 * @param i 单元的行编号
 * @param j 单元的列编号
 * @return double* 
 */
double * energy_matrix(int i, int j);

/**
 * @brief RK第一步计算
 * 
 * @param un n时刻原始向量
 * @param M 系数矩阵
 * @param Rn 右端矩阵
 * @param timeStep 时间步长
 * @param d 向量维数
 * @return double* 
 * 
 * @note 
 * - 输入的方程为M du/dt = R
 * - 输入的M在计算逆矩阵时会变成单位矩阵，若要避免此问题需另设一个指针
 */
double * RK_step1(double * un, double ** M, double * Rn, double timeStep, int d);

/**
 * @brief RK第二步计算
 * 
 * @param un n时刻原始向量
 * @param us1 第一次RK结果
 * @param M 系数矩阵
 * @param Rs1 第一次RK后重新计算的右端矩阵
 * @param timeStep 时间步长
 * @param d 向量维数
 * @return double* 
 * 
 *  * @note 
 * - 输入的方程为M du/dt = R
 * - 输入的M在计算逆矩阵时会变成单位矩阵，若要避免此问题需另设一个指针
 */
double * RK_step2(double * un, double * us1, double ** M, double * Rs1, double timeStep, int d);

/**
 * @brief 选择时间步长dt，
 * 使得每个单元内声速运动距离都不超过单元的某种长度度量
 * 
 * @param dtn 之前的时间步长
 * @return double
 * 
 */
double choose_dt(double dtn);


#endif