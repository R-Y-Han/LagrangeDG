/**
 * @file plot.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 绘制plt文件，生成网格和密度、速度、内能的分布图
 * @version v1.01
 * @date 2022-10-09
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef PLOT_H
#define PLOT_H

extern const char* fn_mesh;   /**< 网格的文件名路径*/
extern const char* fn_pre;  /**< 网格压力的文件名路径*/
extern const char* fn_interene; /**< 网格内能的文件名路径*/
extern const char* fn_ux;

/**
 * @brief 绘制网格图像
 * 
 */
void plotmesh();

/**
 * @brief 绘制网格压力
 * 
 */
void plotpressure();

/**
 * @brief 绘制网格内能
 * 
 */
void plotinternalenergy();

/**
 * @brief 绘制网格速度
 * 
 */
void plotux();

void plotrho1d();
#endif