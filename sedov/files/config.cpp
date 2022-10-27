/**
 * @file config.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief config.h文件中一些函数的定义，如基函数，形函数等
 * @version v1.01
 * @date 2022-10-3
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include <cmath>
#include "config.h"


double dt = 1000;  /**< time step*/

/**
 * @brief 理想气体状态方程
 * 
 * @param rho 密度（点值）
 * @param e 内能（点值）
 * @return double p 压力
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
 * @brief 从参考空间到物理空间的映射x分量
 * 
 * @param xi 参考空间横坐标
 * @param eta 参考空间纵坐标
 * @return double 物理空间横坐标
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
 * @brief 从参考空间到物理空间的映射y分量
 * 
 * @param xi 参考空间横坐标
 * @param eta 参考空间纵坐标
 * @return double 物理空间纵坐标
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
 * @brief 每个单元的基函数，由在质心的Taylor展开得到
 * 
 * @param i 基函数选择
 * @param xi 变量参考空间横坐标
 * @param eta 变量参考空间纵坐标
 * @return double 返回结果：
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
 * @brief 基函数对xi求导
 * 
 * @param i 基函数选择
 * @param xi 变量参考空间横坐标
 * @param eta 变量参考空间纵坐标
 * @return double 返回结果：
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
 * @brief 基函数对eta求导
 * 
* @param i 基函数选择
 * @param xi 变量参考空间横坐标
 * @param eta 变量参考空间纵坐标
 * @return double 返回结果：
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
 * @brief 求(xi,eta)处的Jacobi行列式
 * 
 * @param xi 变量参考空间横坐标
 * @param eta 变量参考空间纵坐标
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
 * @brief 求初始时刻(xi,eta)处的Jacobi行列式
 * 
 * @param xi 变量参考空间横坐标
 * @param eta 变量参考空间纵坐标
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

Omega o[n][m];  /**< 生成n * m个网格*/
node point[n+1][m+1];   /**< 生成(n+1)*(m+1)个节点*/

double gp = sqrt(3.0/5.0);
double Gausspoint_xi[gd] = {-gp, gp, gp, -gp, 0, gp, 0, -gp,0};  /**< Gauss求积节点横坐标，[-1,1]x[1,1]*/
double Gausspoint_eta[gd] = {-gp, -gp, gp, gp, -gp, 0, gp, 0,0};  /**< Gauss求积节点纵坐标*/
double Gaussweight[gd] = {25/81.0, 25/81.0, 25/81.0, 25/81.0, 40/81.0, 40/81.0, 40/81.0, 40/81.0, 64/81.0};   /**< Gauss求积权重*/

/**
 * @brief 形函数，用于求参考空间到物理空间的同胚映射
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = (1 + xit * xi) * ( 1 + etat * eta) / 4.0;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief 形函数对xi求导
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp_xi(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = xit * ( 1 + etat * eta) / 4.0;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief 形函数对eta求导
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp_eta(int i, double xi, double eta)
{
    double ans;
    if ( i < 4)
    {
        double xit, etat;
        xit = ref_xi[i][0];
        etat = ref_xi[i][1];
        ans = (1 + xit * xi) * etat / 4.0;
    }
    else{
        ans = 0;
    }
    return ans;
}

/**
 * @brief 计算a，b两点间距离
 * 
 * @param ax a的物理横坐标
 * @param ay a的物理纵坐标
 * @param bx b的物理横坐标
 * @param by b的物理纵坐标
 * @return double ab间长度
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
 * @brief 计算ab边的单位法向量
 * 
 * @param ax a的物理横坐标
 * @param ay a的物理纵坐标
 * @param bx b的物理横坐标
 * @param by b的物理纵坐标
 * @return double* 物理空间的单位法向量
 * @note 
 * -ab向量（从a指向b）和返回的法向量n组成右手坐标系
 * -注意在节点结构体中存储的是double，在这里返回的是指针，每次计算完后需清空内存
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
 * @brief 初始密度场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_rho(double x, double y)
{
    double ans;
    ans = 1;
    return ans;
}

/**
 * @brief 初始速度场的x方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux(double x, double y)
{
    double ans;
    ans = 0;
    ans = -x;
    return ans;
}

/**
 * @brief 初始速度场对x求导的x分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux_x(double x, double y)
{
    double ans;
    ans = 0;
    ans = -1;
    return  ans;
}

/**
 * @brief 初始速度场对x求导的y分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux_y(double x, double y)
{
    double ans;
    ans = 0;
    return ans;
}

/**
 * @brief 初始速度场y方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy(double x, double y)
{
    double ans;
    ans = 0;
    ans = -y;
    return ans;
}

/**
 * @brief 初始速度场对y求导的x分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy_x(double x, double y)
{
    double ans;
    ans = 0;
    return ans;
}

/**
 * @brief 初始速度场对y求导的y分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy_y(double x, double y)
{
    double ans;
    ans = 0;
    ans = -1;
    return ans;
}

/**
 * @brief 初始压力场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_p(double x, double y)
{
    double ans;
    ans = 1e-6;
    return ans;
}

/**
 * @brief 初始压力场对x求偏导
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_p_x(double x, double y)
{
    double ans;
    ans = 0;
    return ans;
}

/**
 * @brief 初始压力场对y求偏导
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_p_y(double x, double y)
{
    double ans;
    ans = 0;
    return ans;
}
