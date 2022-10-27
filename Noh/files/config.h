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

/** @name 计算区域参数
 * @{
*/

const int dim = 3;  /**< 基函数个数*/

const int n = 10;    /**< 网格横向剖分大小，0-n*/
const int m = 10;    /**< 网格纵向剖分大小，0-m*/

const double X = 1; /**< 计算区域横坐标*/
const double Y = 1; /**< 计算区域纵坐标*/

const double hx = X / (double) n;   /**< spatial step in x*/
const double hy = Y / (double) m;   /**< spatial step in y*/

const double T = 0.1;   /**< end time*/

extern double dt;  /**< time step*/

const double c_CFL = 0.1;  /**< CFL条件系数*/

const double ref_xi[4][2] = {{-1,-1}, {1,-1}, {1,1}, {-1,1}};   /**< 参考单元的四个顶点，从左下角起逆时针*/

/** @}计算区域参数*/

/** @name 流体的参数和变量
 * @{
 */

const double gamma = 5 /(double) 3; /**< specific heat*/

double EOS(double rho, double e);   /**< EOS*/

extern double p;   /**< pressure*/
extern double rho; /**< density*/

/** @struct Omega
 * @brief 网格单元类的定义
 * - 包含网格的所有几何量和部分物理量
 * - 物理单元和参考单元视为同一个单元，保存在同一个结构体中
 * - 网格一维编号和二维编号的关系：q = m*i + j
 */
struct Omega{
    int i, j;   /**< 网格的坐标，第i列，第j行*/
    int q;  /**< 网格的编号，从下往上，从左往右数，q = m*i + j */
    int vertex[4];  /**< 网格的四个顶点，以左下角为第一个，逆时针顺序*/
    double vx[4];   /**< 网格四个顶点的物理坐标横坐标*/
    double vy[4];   /**< 网格四个顶点的物理坐标纵坐标*/
    double vx0[4];   /**< 网格四个顶点的初始物理坐标横坐标，用于计算密度*/
    double vy0[4];   /**< 网格四个顶点的初始物理坐标纵坐标，用于计算密度*/
    vector<int> neighbor_element;   /**< 网格的相邻单元*/
    double xi_c;    /**< 网格质心的横坐标，参考空间*/
    double eta_c;   /**< 网格质心的纵坐标，参考空间*/
    double mass;    /**< 网格的质量*/

    double phi_x(double xi, double eta);    /**< 从参考空间到物理空间的映射，x分量*/
    double phi_y(double xi, double eta);    /**< 从参考空间到物理空间的映射，y分量*/
    double Psi(int i, double xi, double eta);  /**< 网格的基函数*/
    double Psi_xi(int i, double xi, double eta);    /**< 基函数对xi求导*/
    double Psi_eta(int i, double xi, double eta);   /**< 基函数对eta求导*/
    double ** M;    /**< 网格的质量矩阵*/
    double Jacobi(double xi, double eta);   /**< 计算(xi,eta)出的Jacobi矩阵的行列式*/
    double Jacobi_0(double xi, double eta);   /**< 计算初始时刻(xi,eta)出的Jacobi矩阵的行列式*/

    double * uxlast;    /**< element velocity vector x component for t_n*/
    double * uxnext;    /**< element velocity vector x component for t_{n+1}*/
    double * uylast;    /**< element velocity vector y component for t_n*/
    double * uynext;    /**< element velocity vector y component for t_{n+1}*/

    double * taulast;  /**< specific energy for t_n */
    double * taunext; /**< specific energy for t_{n+1}*/

};

/** @struct node
 * @brief 网格节点类的定义
 * - 包含网格节点的所有几何量和部分物理量
 * - 一维序号和二维序号的关系：q = (m+1)* i + j
 */
struct node{
    int i, j;   /**< 节点的参考坐标，第i行第j列*/
    int q;  /**< 节点的一维编号， q = (m+1)* i + j */
    int boundary;   /**< 记录是否为边界节点，boundary = 0,不是，=1，是，else,虚拟节点*/
    double x, y;    /**< 节点的物理坐标*/
    double xtemp, ytemp;    /**< 暂时存储部分数据*/
    double x0, y0;  /**< 节点的初始物理坐标，用于计算密度*/
    double upstarx;  /**< 节点速度x分量*/
    double upstary;  /**< 节点速度x分量*/
    vector<int > neighbor_element;  /**< 节点的相邻网格，从右上角开始逆时针*/
    vector<int > neighbor_node; /**< 节点的相邻节点*/
    vector<double > a;  /**< 节点到相邻节点的部分边长*/
    vector<double > nx; /**< 节点到相邻节点边上的外法向量x方向分量（注意定向：切向法向为右手坐标系）*/
    vector<double > ny; /**< 节点到相邻节点边上的外法向量y方向分量（注意定向：切向法向为右手坐标系）*/
};

extern Omega o[n][m];  /**< 生成n * m个网格*/
extern node point[n+1][m+1];   /**< 生成(n+1)*(m+1)个节点*/
/** @} 流体参数和变量*/


/** @name 数学求解工具
 * @{
 */
const double pi = 3.1415927;    /**< 圆周率*/

const int gd = 9;
extern double Gausspoint_xi[gd];  /**< Gauss求积节点横坐标，[-1,1]x[1,1]*/
extern double Gausspoint_eta[gd];  /**< Gauss求积节点纵坐标*/
extern double Gaussweight[gd];   /**< Gauss求积权重*/

/**
 * @brief 形函数，用于求参考空间到物理空间的同胚映射
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp(int i, double xi, double eta);

/**
 * @brief 形函数对xi求导
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp_xi(int i, double xi, double eta);

/**
 * @brief 形函数对eta求导
 * 
 * @param i 选择函数， i = 0，1，2，3
 * @param xi 参考空间变量横坐标
 * @param eta 参考空间变量纵坐标
 * @return double 
 * - i=0代表左下角，逆时针转
 */
double bp_eta(int i, double xi, double eta);

/**
 * @brief 计算a，b两点间距离
 * 
 * @param ax a的物理横坐标
 * @param ay a的物理纵坐标
 * @param bx b的物理横坐标
 * @param by b的物理纵坐标
 * @return double ab间长度
 */
double length(double ax, double ay, double bx, double by);

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
double * normal(double ax, double ay, double bx, double by);
/** @} */


/**
 * @name 初始条件
 * @brief 定义初值条件函数
 * @{
 */

/**
 * @brief 初始密度场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_rho(double x, double y);

/**
 * @brief 初始速度场的x方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux(double x, double y);

/**
 * @brief 初始速度场对x求导的x分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux_x(double x, double y);

/**
 * @brief 初始速度场对x求导的y分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_ux_y(double x, double y);

/**
 * @brief 初始速度场y方向分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy(double x, double y);

/**
 * @brief 初始速度场对y求导的x分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy_x(double x, double y);

/**
 * @brief 初始速度场对y求导的y分量
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_uy_y(double x, double y);

/**
 * @brief 初始压力场
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_e(double x, double y);

/**
 * @brief 初始压力场对x求偏导
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_e_x(double x, double y);

/**
 * @brief 初始压力场对y求偏导
 * 
 * @param x 物理横坐标
 * @param y 物理纵坐标
 * @return double 
 */
double ini_e_y(double x, double y);

double ini_p(double x, double y);


/** @} */



/** @} */

#endif
