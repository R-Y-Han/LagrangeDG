/**
 * @file config.h
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief Configuration for 2D Taylor Green vortex.
 * Define some parameters and functions shared by all parts of the program.
 * @version 1.01
 * @date 2022-10-12
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
using namespace std;

/** @name �����������
 * @{
*/

const int dim = 3;  /**< ����������*/

const int n = 10;    /**< ��������ʷִ�С��0-n*/
const int m = 10;    /**< ���������ʷִ�С��0-m*/

const double X = 1; /**< �������������*/
const double Y = 1; /**< ��������������*/

const double hx = X / (double) n;   /**< spatial step in x*/
const double hy = Y / (double) m;   /**< spatial step in y*/

const double T = 0.1;   /**< end time*/

extern double dt;  /**< time step*/

const double c_CFL = 0.1;  /**< CFL����ϵ��*/

const double ref_xi[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};   /**< �ο���Ԫ���ĸ����㣬�����½�����ʱ��*/

/** @}�����������*/

/** @name ����Ĳ����ͱ���
 * @{
 */

const double gamma = 5 /(double) 3; /**< specific heat*/

double EOS(double rho, double e);   /**< EOS*/

extern double p;   /**< pressure*/
extern double rho; /**< density*/

/** @struct Omega
 * @brief ����Ԫ��Ķ���
 * - ������������м������Ͳ���������
 * - ����Ԫ�Ͳο���Ԫ��Ϊͬһ����Ԫ��������ͬһ���ṹ����
 * - ����һά��źͶ�ά��ŵĹ�ϵ��q = m*i + j
 */
struct Omega{
    int i, j;   /**< ��������꣬��i�У���j��*/
    int q;  /**< ����ı�ţ��������ϣ�������������q = m*i + j */
    int vertex[4];  /**< ������ĸ����㣬�����½�Ϊ��һ������ʱ��˳��*/
    double vx[4];   /**< �����ĸ�������������������*/
    double vy[4];   /**< �����ĸ��������������������*/
    double vx0[4];   /**< �����ĸ�����ĳ�ʼ������������꣬���ڼ����ܶ�*/
    double vy0[4];   /**< �����ĸ�����ĳ�ʼ�������������꣬���ڼ����ܶ�*/
    vector<int> neighbor_element;   /**< ��������ڵ�Ԫ*/
    double xi_c;    /**< �������ĵĺ����꣬�ο��ռ�*/
    double eta_c;   /**< �������ĵ������꣬�ο��ռ�*/
    double mass;    /**< ���������*/

    double phi_x(double xi, double eta);    /**< �Ӳο��ռ䵽����ռ��ӳ�䣬x����*/
    double phi_y(double xi, double eta);    /**< �Ӳο��ռ䵽����ռ��ӳ�䣬y����*/
    double Psi(int i, double xi, double eta);  /**< ����Ļ�����*/
    double Psi_xi(int i, double xi, double eta);    /**< ��������xi��*/
    double Psi_eta(int i, double xi, double eta);   /**< ��������eta��*/
    double ** M;    /**< �������������*/
    double Jacobi(double xi, double eta);   /**< ����(xi,eta)����Jacobi���������ʽ*/
    double Jacobi_0(double xi, double eta);   /**< �����ʼʱ��(xi,eta)����Jacobi���������ʽ*/

    double * uxlast;    /**< element velocity vector x component for t_n*/
    double * uxnext;    /**< element velocity vector x component for t_{n+1}*/
    double * uylast;    /**< element velocity vector y component for t_n*/
    double * uynext;    /**< element velocity vector y component for t_{n+1}*/

    double * taulast;  /**< specific energy for t_n */
    double * taunext; /**< specific energy for t_{n+1}*/

};

/** @struct node
 * @brief ����ڵ���Ķ���
 * - ��������ڵ�����м������Ͳ���������
 * - һά��źͶ�ά��ŵĹ�ϵ��q = (m+1)* i + j
 */
struct node{
    int i, j;   /**< �ڵ�Ĳο����꣬��i�е�j��*/
    int q;  /**< �ڵ��һά��ţ� q = (m+1)* i + j */
    int boundary;   /**< ��¼�Ƿ�Ϊ�߽�ڵ㣬boundary = 0,���ǣ�=1���ǣ�else,����ڵ�*/
    double x, y;    /**< �ڵ����������*/
    double xtemp, ytemp;    /**< ��ʱ�洢��������*/
    double x0, y0;  /**< �ڵ�ĳ�ʼ�������꣬���ڼ����ܶ�*/
    double upstarx;  /**< �ڵ��ٶ�x����*/
    double upstary;  /**< �ڵ��ٶ�x����*/
    vector<int > neighbor_element;  /**< �ڵ���������񣬴����Ͻǿ�ʼ��ʱ��*/
    vector<int > neighbor_node; /**< �ڵ�����ڽڵ�*/
    vector<double > a;  /**< �ڵ㵽���ڽڵ�Ĳ��ֱ߳�*/
    vector<double > nx; /**< �ڵ㵽���ڽڵ���ϵ��ⷨ����x���������ע�ⶨ��������Ϊ��������ϵ��*/
    vector<double > ny; /**< �ڵ㵽���ڽڵ���ϵ��ⷨ����y���������ע�ⶨ��������Ϊ��������ϵ��*/
};

extern Omega o[n][m];  /**< ����n * m������*/
extern node point[n+1][m+1];   /**< ����(n+1)*(m+1)���ڵ�*/
/** @} ��������ͱ���*/


/** @name ��ѧ��⹤��
 * @{
 */
const double pi = 3.1415927;    /**< Բ����*/

const int gd = 9;
extern double Gausspoint_xi[gd];  /**< Gauss����ڵ�����꣬[-1,1]x[1,1]*/
extern double Gausspoint_eta[gd];  /**< Gauss����ڵ�������*/
extern double Gaussweight[gd];   /**< Gauss���Ȩ��*/

/**
 * @brief �κ�����������ο��ռ䵽����ռ��ͬ��ӳ��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp(int i, double xi, double eta);

/**
 * @brief �κ�����xi��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp_xi(int i, double xi, double eta);

/**
 * @brief �κ�����eta��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp_eta(int i, double xi, double eta);

/**
 * @brief ����a��b��������
 * 
 * @param ax a�����������
 * @param ay a������������
 * @param bx b�����������
 * @param by b������������
 * @return double ab�䳤��
 */
double length(double ax, double ay, double bx, double by);

/**
 * @brief ����ab�ߵĵ�λ������
 * 
 * @param ax a�����������
 * @param ay a������������
 * @param bx b�����������
 * @param by b������������
 * @return double* ����ռ�ĵ�λ������
 * @note 
 * -ab��������aָ��b���ͷ��صķ�����n�����������ϵ
 * -ע���ڽڵ�ṹ���д洢����double�������ﷵ�ص���ָ�룬ÿ�μ������������ڴ�
 */
double * normal(double ax, double ay, double bx, double by);
/** @} */


/**
 * @name ��ʼ����
 * @brief �����ֵ��������
 * @{
 */

/**
 * @brief ��ʼ�ܶȳ�
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_rho(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���x�������
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_ux(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���x�󵼵�x����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_ux_x(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���x�󵼵�y����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_ux_y(double x, double y);

/**
 * @brief ��ʼ�ٶȳ�y�������
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_uy(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���y�󵼵�x����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_uy_x(double x, double y);

/**
 * @brief ��ʼ�ٶȳ���y�󵼵�y����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_uy_y(double x, double y);

/**
 * @brief ��ʼѹ����
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_e(double x, double y);

/**
 * @brief ��ʼѹ������x��ƫ��
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_e_x(double x, double y);

/**
 * @brief ��ʼѹ������y��ƫ��
 * 
 * @param x ���������
 * @param y ����������
 * @return double 
 */
double ini_e_y(double x, double y);

double ini_p(double x, double y);


/** @} */



/** @} */

#endif
