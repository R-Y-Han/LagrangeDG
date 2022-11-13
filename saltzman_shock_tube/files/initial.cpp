/**
 * @file initial.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 生成网格并记录所有初始信息
 * @version v1.01
 * @date 2022-10-04
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include "initial.h"
#include "config.h"
#include <cmath>

#include <iostream>
using namespace std;


/**
 * @brief 生成初始网格划分
 * - 分配内存
 * - 记录网格的编号
 * - 记录节点的编号
 * - 记录节点物理坐标
 * - 记录网格顶点编号、顶点物理坐标
 */
void generatemesh()
{
    int i, j;

    //下面记录节点编号，判断是否为边界节点，记录物理坐标
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            //记录编号
            point[i][j].i = i;
            point[i][j].j = j;
            point[i][j].q = (m+1) * i + j;
            
            //判断是否为边界节点
            point[i][j].boundary = 0;
            if ( i == 0 || i == n)
            {
                point[i][j].boundary = 1;
            }
            if (j == 0 || j == m)
            {
                point[i][j].boundary = 1;
            }

            //记录物理坐标
            point[i][j].x = hx * i;
            point[i][j].y = hy * j;
            
            point[i][j].x = point[i][j].x + (0.1 - point[i][j].y) * sin(pi * point[i][j].x);
            
            point[i][j].x0 = point[i][j].x;
            point[i][j].y0 = point[i][j].y;
        }
    }

    //下面记录网格编号
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            //记录网格编号
            o[i][j].i = i;
            o[i][j].j = j;
            o[i][j].q = m * i + j;
        }
    }

    //下面记录网格顶点
   
    int *ti = new int [4];  //记录顶点逆时针的行编号
    int *tj = new int [4];  //记录顶点逆时针的列编号
    //第一个顶点为左下角的点，逆时针旋转

    for (i = 0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            int tempi, tempj;  //分别表示顶点的行编号和列编号

            ti[0] = i;
            ti[1] = i+1;
            ti[2] = i+1;
            ti[3] = i;
            tj[0] = j;
            tj[1] = j;
            tj[2] = j+1;
            tj[3] = j+1;
            
            for (int k=0; k<4; k++)
            {
                tempi = ti[k];
                tempj = tj[k]; 
                o[i][j].vertex[k] = point[tempi][tempj].q;
                o[i][j].vx[k] = point[tempi][tempj].x;
                o[i][j].vy[k] = point[tempi][tempj].y;

                o[i][j].vx0[k] = point[tempi][tempj].x0;
                o[i][j].vy0[k] = point[tempi][tempj].y0;
            }
        }
    }

    delete[] ti;
    delete[] tj;

    return ;
}

/**
 * @brief 对第(i,j)个单元寻找相邻单元，
 * 从(i,j-1)开始逆时针数第一个存在的单元作为第一个项
 * 存储一维序号
 * 
 * @param i 网格的二维行编号
 * @param j 网格的二维列编号
 */
void element_findneighbor(int i, int j)
{
    o[i][j].neighbor_element.clear();

    int ti;  //代表相邻网格的行编号
    int tj;  //代表相邻网格的列编号
    int num;    //代表相邻网格的一维编号
    
    //判断(i,j-1)单元是否存在
    if ( j == 0) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i;
        tj = j-1;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //判断(i+1,j)单元是否存在
    if ( i == n - 1 ) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i+1;
        tj = j;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //判断(i,j+1)单元是否存在
    if ( j == m-1) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i;
        tj = j+1;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //判断(i-1,j)单元是否存在
    if ( i == 0) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i-1;
        tj = j;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    return ;
}

/**
 * @brief 
 * - 对第(i,j)个节点寻找相邻的网格
 * 从(i,j)个网格逆时针数第一个存在的网格开始记录，逆时针顺序
 * - 对第(i,j)个节点寻找相邻的节点
 * 从(i-1,j)开始逆时针算起第一个存在的节点开始记录，逆时针顺序
 * - 计算各边的（部分）边长和外法向量
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 */
void node_findneighbor(int i, int j)
{
    int ti, tj; //代表相邻网格或节点的行编号和列编号
    int num;    //代表相邻网格或节点的一维编号
    
    //**************先记录相邻网格****************//
    point[i][j].neighbor_element.clear();
    
    //判断网格(i,j)是否存在
    if ( i == n)
    {
        //不存在
        ;
    }
    else{
        ti = i;
        if (j == m)
        {
            //节点为右上角
            ;
        }
        else{
            tj = j;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //判断网格(i-1,j)是否存在
    if ( i == 0)
    {
        //不存在
        ;
    }
    else{
        ti = i-1;
        if (j == m)
        {
            //不存在
            ;
        }
        else{
            tj = j;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //判断网格(i-1,j-1)是否存在
    if ( i == 0)
    {
        //不存在
        ;
    }
    else{
        ti = i-1;
        if (j == 0)
        {
            //不存在
            ;
        }
        else{
            tj = j-1;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //判断网格(i,j-1)是否存在
    if ( i == n)
    {
        //不存在
        ;
    }
    else{
        ti = i;
        if (j == 0)
        {
            //不存在
            ;
        }
        else{
            tj = j-1;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }
    //*****************相邻网格记录完毕******************//

    //*****************下面记录相邻节点******************//
    point[i][j].neighbor_node.clear();

    //判断(i,j-1)节点是否存在
    if ( j == 0) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i;
        tj = j-1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //判断(i+1,j)节点是否存在
    if ( i == n) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i+1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //判断(i,j+1)节点是否存在
    if ( j == m) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i;
        tj = j+1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //判断(i-1,j)节点是否存在
    if ( i == 0) 
    {
        //不存在
        ;
    }
    else{
        //存在
        ti = i-1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }
    //**********相邻节点记录完毕******************//

    //**********下面计算到相邻节点的长度和外法向**********//
    point[i][j].a.clear();
    point[i][j].nx.clear();
    point[i][j].ny.clear();
    for (int k=0; k < point[i][j].neighbor_node.size(); k++)
    {
        int itemp, jtemp;   //相邻节点的行列编号

        itemp = point[i][j].neighbor_node[k] / (m+1);
        jtemp = point[i][j].neighbor_node[k] % (m+1);

        double ax, ay, bx, by;  //a代表本节点，b代表相邻节点

        ax = point[i][j].x;
        ay = point[i][j].y;
        bx = point[itemp][jtemp].x;
        by = point[itemp][jtemp].y;

        //计算到邻点的部分长度
        double l;  //长度
        l = length(ax,ay,bx,by);

        l = l / 2.0; //除以2是因为计算中每个节点被分配到一半的边长
        point[i][j].a.push_back(l);

        //计算法向
        double * ntemp;    //法向
        ntemp = normal(ax,ay,bx,by);
        point[i][j].nx.push_back(ntemp[0]);
        point[i][j].ny.push_back(ntemp[1]);

        delete[] ntemp;
    }

    return ;
}

/**
 * @brief 计算第(i,j)个网格的质心，同时计算每个网格的质量
 * 
 * @param i 网格行编号
 * @param j 网格列编号
 * @return double* 质心的参考单元坐标
 * @note 返回的是指针，存储时是变量，记得删除指针清空内存
 */
double * mass_center(int i, int j)
{
    double * ce = new double [2];

    int k;

    double masst;    //单元质量
    double mf_xi, mf_eta;   //质量矩

    //********下面利用Gauss积分计算质量和质量矩************//
    double rhot, jt;    //记录密度和Jacobi行列式的点值
    double xt, yt;  //Gauss点的物理坐标
    masst = 0;
    mf_xi = 0;
    mf_eta = 0;
    for (k = 0; k < gd; k++)
    {
        xt = o[i][j].phi_x(Gausspoint_xi[k],Gausspoint_eta[k]);
        yt = o[i][j].phi_y(Gausspoint_xi[k],Gausspoint_eta[k]);

        rhot = ini_rho(xt,yt);
        jt = o[i][j].Jacobi(Gausspoint_xi[k],Gausspoint_eta[k]);

        masst = masst + rhot * jt * Gaussweight[k];
        mf_xi = mf_xi + rhot * Gausspoint_xi[k] * jt * Gaussweight[k];
        mf_eta = mf_eta + rhot * Gausspoint_eta[k] * jt * Gaussweight[k];
    }
    o[i][j].mass = masst;

    //**********下面计算质心坐标***********//

    ce[0] = mf_xi / masst;  //质心参考空间横坐标
    ce[1] = mf_eta / masst; //质心参考空间纵坐标

    return ce;
}

/**
 * @brief 计算第(i,j)个网格的质量矩阵
 * 
 * @param i 网格行编号
 * @param j 网格列编号
 * @return double** 质量矩阵
 */
double ** mass_matrix(int i, int j)
{
    int k, l;
    double rhot, jt, xt, yt, temp;
    double ** mma = new double * [dim];
    for ( k = 0; k < dim; k++)
    {
        mma[k] = new double [dim];
    }

    for (k = 0; k<dim; k++)
    {
        for (l=0; l<dim; l++)
        {
            mma[k][l] = 0;
            for (int r = 0; r <gd; r++)
            {
                xt = o[i][j].phi_x(Gausspoint_xi[r], Gausspoint_eta[r]);
                yt = o[i][j].phi_y(Gausspoint_xi[r], Gausspoint_eta[r]);

                rhot = ini_rho(xt,yt);

                jt = o[i][j].Jacobi(Gausspoint_xi[r], Gausspoint_eta[r]);
                
                temp = o[i][j].Psi(k,Gausspoint_xi[r],Gausspoint_eta[r]);
                temp = temp * o[i][j].Psi(l,Gausspoint_xi[r],Gausspoint_eta[r]);
                
                mma[k][l] = mma[k][l] + rhot * temp * jt * Gaussweight[r];
            }
        }
    }

    return mma;
}

/**
 * @brief 初始化各个变量，生成网格并记录初始信息
 * 
 */
void initial()
{
    int i, j, k;
    //***********下面生成网格划分**********//
    generatemesh();

    //***********下面对每个网格寻找相邻网格**********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            element_findneighbor(i,j);
        }
    }

    //***********下面对每个节点寻找相邻节点***********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            node_findneighbor(i,j);
        }
    }


    //************下面计算每个网格的质量和质心，并计算质量矩阵*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            //下面计算质心
            double * ce;
            ce = mass_center(i,j);
            o[i][j].xi_c = ce[0];
            o[i][j].eta_c = ce[1];

            delete[] ce;

            //下面计算质量矩阵
            o[i][j].M = mass_matrix(i,j);
        }
    }

    //***********下面给每个单元赋初值***************//
    double xt, yt; //质心的物理坐标
    double x_xi, x_eta, y_xi, y_eta;    //Jacobi矩阵的四个分量
    double taux, tauy;
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            o[i][j].uxlast = new double [dim];
            o[i][j].uxnext = new double [dim];
            o[i][j].uylast = new double [dim];
            o[i][j].uynext = new double [dim];
            o[i][j].taulast = new double [dim];
            o[i][j].taunext = new double [dim];

            xt = o[i][j].phi_x(o[i][j].xi_c,o[i][j].eta_c);
            yt = o[i][j].phi_y(o[i][j].xi_c,o[i][j].eta_c);
           
            x_xi = 0;
            x_eta = 0;
            y_xi = 0;
            y_eta = 0;
            for (k = 0; k< 4; k++)
            {
                x_xi = x_xi + o[i][j].vx[k] * bp_xi(k,o[i][j].xi_c,o[i][j].eta_c);
                x_eta = x_eta + o[i][j].vx[k] * bp_eta(k,o[i][j].xi_c,o[i][j].eta_c);
                y_xi = y_xi + o[i][j].vy[k] * bp_xi(k,o[i][j].xi_c,o[i][j].eta_c);
                y_eta = y_eta + o[i][j].vy[k] * bp_eta(k,o[i][j].xi_c,o[i][j].eta_c);
            }

            //********下面初始化速度场************//
            o[i][j].uxlast[0] = ini_ux(xt,yt);
            o[i][j].uylast[0] = ini_uy(xt,yt);  //在质心的值

            o[i][j].uxlast[1] = ini_ux_x(xt,yt) * x_xi + ini_ux_y(xt,yt) * y_xi;    //质心处ux对xi求偏导
            o[i][j].uylast[1] = ini_uy_x(xt,yt) * x_xi + ini_uy_y(xt,yt) * y_xi;    //质心处uy对xi求偏导

            o[i][j].uxlast[2] = ini_ux_x(xt,yt) * x_eta + ini_ux_y(xt,yt) * y_eta;  //质心处ux对eta求偏导
            o[i][j].uylast[2] = ini_uy_x(xt,yt) * x_eta + ini_uy_y(xt,yt) * y_eta;  //质心处uy对eta求偏导
            
            //*********下面初始化总能量场*********//
            double rho;
            rho = ini_rho(xt,yt);
            o[i][j].rholast = rho;
            o[i][j].rho0 = rho;
            o[i][j].taulast[0] = ini_p(xt,yt) / ( rho * (gamma - 1)) + 0.5 * (ini_ux(xt,yt) * ini_ux(xt,yt) + ini_uy(xt,yt) * ini_uy(xt,yt));
            
            taux = ini_p_x(xt,yt) / ( rho * (gamma - 1)) + ini_ux(xt,yt) * ini_ux_x(xt,yt) + ini_uy(xt,yt) * ini_uy_x(xt,yt);
            tauy = ini_p_y(xt,yt) / ( rho * (gamma - 1)) + ini_ux(xt,yt) * ini_ux_y(xt,yt) + ini_uy(xt,yt) * ini_uy_y(xt,yt);
            
            o[i][j].taulast[1] = taux * x_xi + tauy * y_xi;
            o[i][j].taulast[2] = taux * x_eta + tauy * y_eta;
        }
    }

    return ;
}