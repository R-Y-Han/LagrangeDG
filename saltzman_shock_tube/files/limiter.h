/**
 * @file limiter.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 
 * @version 0.1
 * @date 2022-10-14
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Reights Reserved.
 * 
 */

#ifndef LIMITER_H
#define LIMITER_H

/**
 * @brief 对标量进行限制
 * 
 * @param i 单元行编号
 * @param j 单元列编号
 * @param up 限制前的值
 * @param alpha 限制的系数
 * @return double 
 */
double tau_limiter(int i, int j, double up, double alpha);

double ux_limiter(int i, int j, double up, double alpha);

double uy_limiter(int i, int j, double up, double alpha);
#endif