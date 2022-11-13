/**
 * @file plot.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief ����plt�ļ�������������ܶȡ��ٶȡ����ܵķֲ�ͼ
 * @version v1.01
 * @date 2022-10-09
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef PLOT_H
#define PLOT_H

extern const char* fn_mesh;   /**< ������ļ���·��*/
extern const char* fn_pre;  /**< ����ѹ�����ļ���·��*/
extern const char* fn_interene; /**< �������ܵ��ļ���·��*/
extern const char* fn_ux;

/**
 * @brief ��������ͼ��
 * 
 */
void plotmesh();

/**
 * @brief ��������ѹ��
 * 
 */
void plotpressure();

/**
 * @brief ������������
 * 
 */
void plotinternalenergy();

/**
 * @brief ���������ٶ�
 * 
 */
void plotux();

void plotrho1d();
#endif