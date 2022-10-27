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
 * @brief �Ա�����������
 * 
 * @param i ��Ԫ�б��
 * @param j ��Ԫ�б��
 * @param up ����ǰ��ֵ
 * @param alpha ���Ƶ�ϵ��
 * @return double 
 */
double tau_limiter(int i, int j, double up, double alpha);

double ux_limiter(int i, int j, double up, double alpha);

double uy_limiter(int i, int j, double up, double alpha);
#endif