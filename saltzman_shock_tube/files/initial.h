/**
 * @file initial.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 生成网格并记录所有初始信息
 * @version v1.01
 * @date 2022-10-04
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef INITIAL_H
#define INITIAL_H

/**
 * @brief 生成初始网格划分：
 * - 分配内存
 * - 记录网格的编号
 * - 记录节点的编号
 * - 记录节点物理坐标
 * - 记录网格顶点编号、顶点物理坐标
 */
void generatemesh();

/**
 * @brief 对第(i,j)个单元寻找相邻单元，
 * 从(i,j-1)开始逆时针数第一个存在的单元作为第一个项。
 * 存储一维序号。
 * 
 * @param i 网格的二维行编号
 * @param j 网格的二维列编号
 */
void element_findneighbor(int i, int j);

/**
 * @brief 
 * - 对第(i,j)个节点寻找相邻的网格，
 * 从(i,j)个网格逆时针数第一个存在的网格开始记录，逆时针顺序。
 * - 对第(i,j)个节点寻找相邻的节点，
 * 从(i-1,j)开始逆时针算起第一个存在的节点开始记录，逆时针顺序。
 * - 计算各边的（部分）边长和外法向量
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 */
void node_findneighbor(int i, int j);

/**
 * @brief 计算第(i,j)个网格的质心
 * 
 * @param i 网格行编号
 * @param j 网格列编号
 * @return double* 质心的参考单元坐标
 * @note 返回的是指针，存储时是变量，记得删除指针清空内存
 */
double * mass_center(int i, int j);

/**
 * @brief 计算第(i,j)个网格的质量矩阵
 * 
 * @param i 网格行编号
 * @param j 网格列编号
 * @return double** 质量矩阵
 */
double ** mass_matrix(int i, int j);

/**
 * @brief 初始化各个变量，生成网格并记录初始信息
 * 
 */
void initial();



#endif