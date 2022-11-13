/**
 * @file time_evolution.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief һ��ʱ�䲽���ĸ��£������ܶȡ��ٶȡ��������ܵĸ��º�ʱ�䲽����ѡȡ
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
 * @param xi �ڵ�Ĳο����������
 * @param eta �ڵ�Ĳο�����������
 * @param r �ڵ����ڵ�Ԫ�ı��
 * @return double 
 */
double get_density(double xi, double eta, int r);

/**
 * @brief �õ���Ԫ�ٶȵ��Ҷ˾���
 * 
 * @param i ��Ԫ���б��
 * @param j ��Ԫ���б��
 * @return double** 
 */
double ** velocity_matrix(int i, int j);

/**
 * @brief �õ���Ԫ�������ܵ��Ҷ˾���
 * 
 * @param i ��Ԫ���б��
 * @param j ��Ԫ���б��
 * @return double* 
 */
double * energy_matrix(int i, int j);

/**
 * @brief RK��һ������
 * 
 * @param un nʱ��ԭʼ����
 * @param M ϵ������
 * @param Rn �Ҷ˾���
 * @param timeStep ʱ�䲽��
 * @param d ����ά��
 * @return double* 
 * 
 * @note 
 * - ����ķ���ΪM du/dt = R
 * - �����M�ڼ��������ʱ���ɵ�λ������Ҫ���������������һ��ָ��
 */
double * RK_step1(double * un, double ** M, double * Rn, double timeStep, int d);

/**
 * @brief RK�ڶ�������
 * 
 * @param un nʱ��ԭʼ����
 * @param us1 ��һ��RK���
 * @param M ϵ������
 * @param Rs1 ��һ��RK�����¼�����Ҷ˾���
 * @param timeStep ʱ�䲽��
 * @param d ����ά��
 * @return double* 
 * 
 *  * @note 
 * - ����ķ���ΪM du/dt = R
 * - �����M�ڼ��������ʱ���ɵ�λ������Ҫ���������������һ��ָ��
 */
double * RK_step2(double * un, double * us1, double ** M, double * Rs1, double timeStep, int d);

/**
 * @brief ѡ��ʱ�䲽��dt��
 * ʹ��ÿ����Ԫ�������˶����붼��������Ԫ��ĳ�ֳ��ȶ���
 * 
 * @param dtn ֮ǰ��ʱ�䲽��
 * @return double
 * 
 */
double choose_dt(double dtn);


#endif