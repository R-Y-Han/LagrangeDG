/**
 * @file nodal_solver.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 定义节点解法器求解节点速度和corner force
 * @version v1.01
 * @date 2022-10-05
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef NODAL_SOLVER_H
#define NODAL_SOLVER_H

/**
 * @brief 计算第(i,j)个节点处每个相邻单元在该点的速度重构值的平均
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double* 
 */
double * u_average(int i, int j);

/**
 * @brief 计算边界节点的速度
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double* 
 */
double * nodal_velocity_boundary(int i, int j);

/**
 * @brief 计算节点速度
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double *
 */
double * nodal_velocity(int i, int j);

/**
 * @brief 计算第(i,j)个节点在第k个相邻单元中的两个corner force
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @param r 所考虑的单元编号
 * @return double** 
 * 
 * @note 
 * - 注意在第k个单元内有两个corner force
 * - 注意法向量的方向是否指向单元外部
 */
double ** corner_force(int i, int j, int r);

#endif