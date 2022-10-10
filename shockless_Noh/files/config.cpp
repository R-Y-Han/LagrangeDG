/**
 * @file config.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief config.h�ļ���һЩ�����Ķ��壬����������κ�����
 * @version v1.01
 * @date 2022-10-3
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <cmath>
#include "config.h"


double dt = 1000;  /**< time step*/

/**
 * @brief ��������״̬����
 * 
 * @param rho �ܶȣ���ֵ��
 * @param e ���ܣ���ֵ��
 * @return double p ѹ��
 */
double EOS(double rho, double e)
{
    double p;
    p = (gamma - 1) * rho * e;
    return p;
}

double p;   /**< pressure*/
double rho; /**< density*/


/**
 * @brief �Ӳο��ռ䵽�����ռ��ӳ��x����
 * 
 * @param xi �ο��ռ������
 * @param eta �ο��ռ�������
 * @return double �����ռ������
 */
double Omega::phi_x(double xi, double eta)
{
    double xtemp;
    xtemp = 0;
    for (int i=0; i< 4; i++)
    {
        xtemp = xtemp + this->vx[i] * bp(i,xi,eta);
    }
    return xtemp;
}

/**
 * @brief �Ӳο��ռ䵽�����ռ��ӳ��y����
 * 
 * @param xi �ο��ռ������
 * @param eta �ο��ռ�������
 * @return double �����ռ�������
 */
double Omega::phi_y(double xi, double eta)
{
    double ytemp;
    ytemp = 0;
    for (int i=0; i<4; i++)
    {
        ytemp = ytemp + this->vy[i] * bp(i,xi,eta);
    }
    return ytemp;
}

/**
 * @brief ÿ����Ԫ�Ļ��������������ĵ�Taylorչ���õ�
 * 
 * @param i ������ѡ��
 * @param xi �����ο��ռ������
 * @param eta �����ο��ռ�������
 * @return double ���ؽ����
 *  - i = 0, Psi = 1
 *  - i = 1, Psi = xi - xi_c
 *  - i = 2, Psi = eta - eta_c
 *  - otherwise, no return (only 2nd-order).
 */
double Omega::Psi(int i, double xi, double eta)
{
    double ans;
    if ( i == 0)
    {
        ans = 1;
    }
    else if ( i == 1)
    {
        ans = (xi - this->xi_c);
    }
    else if ( i == 2)
    {
        ans = (eta - this->eta_c);
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief ��������xi��
 * 
* @param i ������ѡ��
 * @param xi �����ο��ռ������
 * @param eta �����ο��ռ�������
 * @return double ���ؽ����
 *  - i = 0, Psi = 0
 *  - i = 1, Psi = 1
 *  - i = 2, Psi = 0
 *  - otherwise, no return (only 2nd-order).
 */
double Omega::Psi_xi(int i, double xi, double eta)
{
    double ans;
    if ( i == 0)
    {
        ans = 0;
    }
    else if ( i == 1)
    {
        ans = 1;
    }
    else if ( i == 2)
    {
        ans = 0;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief ��������eta��
 * 
* @param i ������ѡ��
 * @param xi �����ο��ռ������
 * @param eta �����ο��ռ�������
 * @return double ���ؽ����
 *  - i = 0, Psi = 0
 *  - i = 1, Psi = 0
 *  - i = 2, Psi = 1
 *  - otherwise, no return (only 2nd-order).
 */
double Omega::Psi_eta(int i, double xi, double eta)
{
    double ans;
    if ( i == 0)
    {
        ans = 0;
    }
    else if ( i == 1)
    {
        ans = 0;
    }
    else if ( i == 2)
    {
        ans = 1;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief ��(xi,eta)����Jacobi����ʽ
 * 
 * @param xi �����ο��ռ������
 * @param eta �����ο��ռ�������
 * @return double |J|
 */
double Omega::Jacobi(double xi, double eta)
{
    double ans;
    double a, b, c, d;
    
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    for ( int i = 0; i < 4; i++)
    {
        a = a + this->vx[i] * bp_xi(i,xi,eta);
        b = b + this->vx[i] * bp_eta(i,xi,eta);
        c = c + this->vy[i] * bp_xi(i,xi,eta);
        d = d + this->vy[i] * bp_eta(i,xi,eta);
    }
    ans = a * d - b * c;
    return ans;
}

/**
 * @brief ���ʼʱ��(xi,eta)����Jacobi����ʽ
 * 
 * @param xi �����ο��ռ������
 * @param eta �����ο��ռ�������
 * @return double |J|
 */
double Omega::Jacobi_0(double xi, double eta)
{
    double ans;
    double a, b, c, d;
    
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    for ( int i = 0; i < 4; i++)
    {
        a = a + this->vx0[i] * bp_xi(i,xi,eta);
        b = b + this->vx0[i] * bp_eta(i,xi,eta);
        c = c + this->vy0[i] * bp_xi(i,xi,eta);
        d = d + this->vy0[i] * bp_eta(i,xi,eta);
    }
    ans = a * d - b * c;
    return ans;
}

Omega o[n][m];  /**< ����n * m������*/
node point[n+1][m+1];   /**< ����(n+1)*(m+1)���ڵ�*/


double Gausspoint_xi[4] = { - 1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3), - 1 / sqrt(3)};  /**< Gauss����ڵ�����꣬[-1,1]x[1,1]*/
double Gausspoint_eta[4] = { - 1 / sqrt(3), - 1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3)};  /**< Gauss����ڵ�������*/
double Gaussweight[4] = {1, 1, 1, 1};   /**< Gauss���Ȩ��*/

/**
 * @brief �κ�����������ο��ռ䵽�����ռ��ͬ��ӳ��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = (1 + xit * xi) * ( 1 + etat * eta) / 4;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief �κ�����xi��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp_xi(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = xit * ( 1 + etat * eta) / 4;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief �κ�����eta��
 * 
 * @param i ѡ������ i = 0��1��2��3
 * @param xi �ο��ռ����������
 * @param eta �ο��ռ����������
 * @return double 
 * - i=0�������½ǣ���ʱ��ת
 */
double bp_eta(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = (1 + xit * xi) * etat / 4;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief ����a��b��������
 * 
 * @param ax a������������
 * @param ay a������������
 * @param bx b������������
 * @param by b������������
 * @return double ab�䳤��
 */
double length(double ax, double ay, double bx, double by)
{
    double ans;
    ans = (ax - bx) * (ax - bx);
    ans = ans + (ay - by ) * ( ay - by);
    ans = sqrt(ans);
    return ans;
}

/**
 * @brief ����ab�ߵĵ�λ������
 * 
 * @param ax a������������
 * @param ay a������������
 * @param bx b������������
 * @param by b������������
 * @return double* �����ռ�ĵ�λ������
 * @note 
 * -ab��������aָ��b���ͷ��صķ�����n�����������ϵ
 * -ע���ڽڵ�ṹ���д洢����double�������ﷵ�ص���ָ�룬ÿ�μ������������ڴ�
 */
double * normal(double ax, double ay, double bx, double by)
{
    double * n = new double [2];
    double norm;
    norm = length(ax,ay,bx,by);
    n[0] = -(by - ay) / norm;
    n[1] = (bx - ax) / norm;
    return n;
}


/**
 * @brief ��ʼ�ܶȳ�
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_rho(double x, double y)
{
    return 1;
}

/**
 * @brief ��ʼ�ٶȳ���x�������
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_ux(double x, double y)
{
    return -x;
}

/**
 * @brief ��ʼ�ٶȳ���x�󵼵�x����
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_ux_x(double x, double y)
{
    return  -1;
}

/**
 * @brief ��ʼ�ٶȳ���x�󵼵�y����
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_ux_y(double x, double y)
{
    return 0;
}

/**
 * @brief ��ʼ�ٶȳ�y�������
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_uy(double x, double y)
{
    return -y;
}

/**
 * @brief ��ʼ�ٶȳ���y�󵼵�x����
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_uy_x(double x, double y)
{
    return 0;
}

/**
 * @brief ��ʼ�ٶȳ���y�󵼵�y����
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_uy_y(double x, double y)
{
    return -1;
}

/**
 * @brief ��ʼ���ܳ�
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_e(double x, double y)
{
    return 1;
}

/**
 * @brief ��ʼ���ܳ���x��ƫ��
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_e_x(double x, double y)
{
    return 0;
}

/**
 * @brief ��ʼ���ܳ���y��ƫ��
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_e_y(double x, double y)
{
    return 0;
}

/**
 * @brief ��ʼѹ����
 * 
 * @param x ����������
 * @param y ����������
 * @return double 
 */
double ini_p(double x, double y)
{
    double ans;
    ans = ( gamma - 1 ) * ini_rho(x,y) * ini_e(x,y);
    return ans;
}

/**
 * @brief ����(i,j)tʱ���ܶȽ�����
 * 
 * @param i �����б��
 * @param j �����б��
 * @param xit �ο��ڵ������
 * @param etat �ο��ڵ�������
 * @param t ʱ��
 * @return double 
 */
double ana_rho(int i, int j, double xit, double etat, double t)
{
    double ans, xt, yt;
    xt = 0;
    yt = 0;
    for (int k=0; k<4; k++)
    {
        xt = xt + o[i][j].vx0[k] * bp(k,xit,etat);
        yt = yt + o[i][j].vy0[k] * bp(k,xit,etat);
    }//xt,ytΪ��ʼʱ��(xit,etat)��������
    ans = ini_rho(xt,yt) / pow((1-t), 2);
    return ans;
}

/**
 * @brief �ٶ�x�������������
 * 
 * @param i �����б��
 * @param j �����б��
 * @param xit �ο��ڵ������
 * @param etat �ο��ڵ�������
 * @param t ʱ��
 * @return double 
 */
double ana_ux(int i, int j, double xit, double etat, double t)
{
    double ans, xt, yt;
    xt = 0;
    yt = 0;
    for (int k=0; k<4; k++)
    {
        xt = xt + o[i][j].vx0[k] * bp(k,xit,etat);
        yt = yt + o[i][j].vy0[k] * bp(k,xit,etat);
    }//xt,ytΪ��ʼʱ��(xit,etat)��������
    ans = ini_ux(xt,yt);
    return ans;
}

/**
 * @brief �ٶ�y�������������
 * 
 * @param i �����б��
 * @param j �����б��
 * @param xit �ο��ڵ������
 * @param etat �ο��ڵ�������
 * @param t ʱ��
 * @return double 
 */
double ana_uy(int i, int j, double xit, double etat, double t)
{
    double ans, xt, yt;
    xt = 0;
    yt = 0;
    for (int k=0; k<4; k++)
    {
        xt = xt + o[i][j].vx0[k] * bp(k,xit,etat);
        yt = yt + o[i][j].vy0[k] * bp(k,xit,etat);
    }//xt,ytΪ��ʼʱ��(xit,etat)��������
    ans = ini_uy(xt,yt);
    return ans;
}

/**
 * @brief ���ܽ�����
 * 
 * @param i �����б��
 * @param j �����б��
 * @param xit �ο��ڵ������
 * @param etat �ο��ڵ�������
 * @param t ʱ��
 * @return double 
 */
double ana_e(int i, int j, double xit, double etat, double t)
{
    double ans, xt, yt;
    xt = 0;
    yt = 0;
    for (int k=0; k<4; k++)
    {
        xt = xt + o[i][j].vx0[k] * bp(k,xit,etat);
        yt = yt + o[i][j].vy0[k] * bp(k,xit,etat);
    }//xt,ytΪ��ʼʱ��(xit,etat)��������
    ans = pow( (1-t), 2 * (gamma - 1) );
    ans = ini_e(xt,yt) / ans;
    return ans;
}