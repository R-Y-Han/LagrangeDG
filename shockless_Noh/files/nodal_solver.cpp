/**
 * @file nodal_solver.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 定义节点解法器求解节点速度和corner force
 * @version v1.01
 * @date 2022-10-05
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include "config.h"
#include "nodal_solver.h"
#include <cmath>

#include <iostream>
using namespace std;

/**
 * @brief 计算第(i,j)个节点处每个相邻单元在该点的速度重构值的平均
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double* 
 */
double * u_average(int i, int j)
{
    double * uavg = new double [2];
    int k, l, r, temp;
    double xit, etat, norm;

    uavg[0] = 0;
    uavg[1] = 0;
    for (r = 0; r < point[i][j].neighbor_element.size(); r++)
    {
        //************对第k个相邻单元执行*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;    //节点第r个相邻单元的行列编号
        
        //********下面判断节点(i,j)在单元(k,l)中是第几个顶点*******//
        for (temp=0; temp<4; temp++)
        {
            if (o[k][l].vertex[temp] == point[i][j].q)
            {
                //此时循环指标temp刚好就是节点(i,j)在(k,l)单元中的局部位置
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];

        //**********下面计算该单元内节点处速度重构值*******//
        for (temp = 0; temp < dim; temp++)
        {
            uavg[0] = uavg[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
            uavg[1] = uavg[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
        }
    }

    //********此时uavg记录了所有单元在该点处速度重构值的总和********//
    uavg[0] = uavg[0] / (double) point[i][j].neighbor_element.size();
    uavg[1] = uavg[1] / (double) point[i][j].neighbor_element.size();
    //********uavg为节点处速度重构值的平均****************//

    return uavg;
}

/**
 * @brief 计算节点速度
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double *
 */
double * nodal_velocity(int i, int j)
{
    double* ustar = new double [2];

    //**************先判断是否为边界点*********//
    if (point[i][j].boundary == 1)
    {
        //******是边界点*******//
        ustar[0] = - point[i][j].x0;
        ustar[1] = - point[i][j].y0;
        //边界速度不变
        return ustar;
    }

    //*************不是边界节点***************//

    int r, temp, k, l, node1, node2, p1, p2;
    double xit, etat, nx1, ny1, nx2, ny2, a1, a2;
    double mu, topx, topy, bottom;
    topx = 0;
    topy = 0;
    bottom = 0;
    //***************对每个corner进行计算***********//
    for (r=0; r<point[i][j].neighbor_element.size(); r++)
    {
        //***********确定corner所在单元编号*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;
        
        //************确定节点在单元内的位置和参考空间坐标***********//
        for (temp=0; temp < 4; temp++)
        {
            if (o[k][l].vertex[temp] == point[i][j].q)
            {
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];

        //***********确定相邻节点的编号***********//
        if (temp == 3)
        {
            node1 = o[k][l].vertex[0];
        }
        else{
            node1 = o[k][l].vertex[temp + 1];
        }
        
        if (temp == 0)
        {
            node2 = o[k][l].vertex[3];
        }
        else{
            node2 = o[k][l].vertex[temp - 1];
        }
        //node1, node2为相邻顶点的编号

        //*********确定相邻节点和该节点的关系*********//
        for (p1 = 0; p1 < point[i][j].neighbor_node.size(); p1++)
        {
            if (point[i][j].neighbor_node[p1] == node1)
            {
                break;
            }
        }

        for (p2 = 0; p2 < point[i][j].neighbor_node.size(); p2++)
        {
            if (point[i][j].neighbor_node[p2] == node2)
            {
                break;
            }
        }
        //p1, p2为(i,j)指向相邻顶点的指针

        //***********记录单位外法向和相应长度********//
        a1 = point[i][j].a[p1];
        a2 = point[i][j].a[p2];

        nx1 = - point[i][j].nx[p1];
        ny1 = - point[i][j].ny[p1]; //注意此处加了负号，为使指向外向

        nx2 = point[i][j].nx[p2];
        ny2 = point[i][j].ny[p2];

        //*******下面计算第r个corner中的相应分量********//
        double xt, yt;
        xt = o[k][l].phi_x(xit,etat);
        yt = o[k][l].phi_y(xit,etat);

        //********下面计算密度**********//
        rho = ini_rho(xt,yt) * o[k][l].Jacobi_0(xit,etat);
        rho = rho / o[k][l].Jacobi(xit,etat);

        //********下面计算节点处速度重构值********//
        double * uc = new double [2];
        uc[0] = 0;
        uc[1] = 0;
        for (temp = 0; temp < dim; temp++)
        {
            uc[0] = uc[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
            uc[1] = uc[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
        }

        //********下面计算内能**********//
        
        double ve;
        ve = 0;
        for (temp = 0; temp < dim; temp++)
        {
            ve = ve + o[k][l].taulast[temp] * o[k][l].Psi(temp,xit,etat);
        }
        ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        
        //**************下面计算阻抗*************//
        mu = (gamma - 1) * ve;
        mu = sqrt(mu);  //此时mu为声速
        mu = rho * mu;
        
        //*************下面计算激波单位方向*********//
        double * e;
        e = u_average(i,j);
        e[0] = e[0] - uc[0];
        e[1] = e[1] - uc[1];
        
        double norm;
        norm = e[0] * e[0] + e[1] * e[1];
        norm = sqrt(norm);
        if (norm == 0)
        {
            e[0] = 0;
            e[1] = 0;
        }
        else{
            e[0] = e[0] / norm;
            e[1] = e[1] / norm;
        }
        
        //*********下面计算分子分母各分量********//
        p = EOS (rho, ve);
        
        double mid;
        //*********下面计算第一条边的量***********//
        mid = nx1 * e[0] + ny1 * e[1];
        mid = mu * abs(mid) * a1;

        topx = topx + mid * uc[0] + a1 * nx1 * p;
        topy = topy + mid * uc[1] + a1 * ny1 * p;
        bottom = bottom + mid;

        //************下面计算第二条边的量*******//
        mid = nx2 * e[0] + ny2 * e[1];
        mid = mu * abs(mid) * a2;

        topx = topx + mid * uc[0] + a2 * nx2 * p;
        topy = topy + mid * uc[1] + a2 * ny2 * p;
        bottom = bottom + mid;

        delete[] uc;
        delete[] e;
    }
    if (bottom == 0)
    {
        double * utemp = u_average(i,j);
        ustar[0] = utemp[0];
        ustar[1] = utemp[1];
        delete[] utemp;
    }
    else{
        ustar[0] = topx / bottom;
        ustar[1] = topy / bottom;
    }

    return ustar;
}

/**
 * @brief 计算第(i,j)个节点在第k个相邻单元中的两个corner force
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @param r 所考虑的单元编号
 * @return double** 
 * 
 * @note 
 * - 注意在第k个单元内有两个corner force
 * - 注意法向量的方向是否指向单元外部
 */
double ** corner_force(int i, int j, int r)
{
    int k, l, temp, node1, node2, p1, p2;
    double xit, etat, a1, a2, nx1, ny1, nx2, ny2, mu;

    double ** f = new double * [2]; //在每个角落内有两个corner force
    for (temp=0; temp < 2; temp++)
    {
        f[temp] = new double [2];
    }

    //***********确定corner所在单元编号*********//
    k = r / m;
    l = r % m;
        
    //************确定节点在单元内的位置和参考空间坐标***********//
    for (temp=0; temp < 4; temp++)
    {
        if (o[k][l].vertex[temp] == point[i][j].q)
        {
            break;
        }
    }
    xit = ref_xi[temp][0];
    etat = ref_xi[temp][1];

    //***********确定相邻节点的编号***********//
    //此时temp记录节点在单元中的位置
    if (temp == 3)
    {
        node1 = o[k][l].vertex[0];
    }
    else{
        node1 = o[k][l].vertex[temp + 1];
    }
        
    if (temp == 0)
    {
        node2 = o[k][l].vertex[3];
    }
    else{
        node2 = o[k][l].vertex[temp - 1];
    }
    //node1, node2为相邻顶点的编号

    //*********确定相邻节点和该节点的关系*********//
    for (p1 = 0; p1 < point[i][j].neighbor_node.size(); p1++)
    {
        if (point[i][j].neighbor_node[p1] == node1)
        {
            break;
        }
    }

    for (p2 = 0; p2 < point[i][j].neighbor_node.size(); p2++)
    {
        if (point[i][j].neighbor_node[p2] == node2)
        {
            break;
        }
    }
    //p1, p2为(i,j)指向相邻顶点的指针

    //***************下面记录每条边上的外法向量和长度************//
    a1 = point[i][j].a[p1];
    a2 = point[i][j].a[p2];

    nx1 = - point[i][j].nx[p1];
    ny1 = - point[i][j].ny[p1];
    nx2 = point[i][j].nx[p2];
    ny2 = point[i][j].ny[p2];

    //*******下面计算第r个corner中的相应分量********//
    double xt, yt;
    xt = o[k][l].phi_x(xit,etat);
    yt = o[k][l].phi_y(xit,etat);

    //********下面计算密度**********//
    rho = ini_rho(xt,yt) * o[k][l].Jacobi_0(xit,etat);
    rho = rho / o[k][l].Jacobi(xit,etat);

    //********下面计算节点处速度重构值********//
    double * uc = new double [2];
    uc[0] = 0;
    uc[1] = 0;
    for (temp = 0; temp < dim; temp++)
    {
        uc[0] = uc[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
        uc[1] = uc[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
    }

    //********下面计算内能**********//
        
    double ve;
    ve = 0;
    for (temp = 0; temp < dim; temp++)
    {
        ve = ve + o[k][l].taulast[temp] * o[k][l].Psi(temp,xit,etat);
    }
    ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        
    //**************下面计算阻抗*************//
    mu = (gamma - 1) * ve;
    mu = sqrt(mu);  //此时mu为声速
    mu = rho * mu;
        
    //*************下面计算激波单位方向*********//
    double * e;
    e = u_average(i,j);
    e[0] = e[0] - uc[0];
    e[1] = e[1] - uc[1];
        
    double norm;
    norm = e[0] * e[0] + e[1] * e[1];
    norm = sqrt(norm);
    if (norm == 0)
    {
        e[0] = 0;
        e[1] = 0;
    }
    else{
        e[0] = e[0] / norm;
        e[1] = e[1] / norm;
    }

    //*********下面计算压强********//
    p = EOS (rho, ve);

    //***********下面计算两条边上的corner force*********//
    double mid;
    
    mid = nx1 * e[0] + ny1 * e[1];
    mid = mu * abs(mid) * a1;
    f[0][0] = - a1 * nx1 * p + mid * ( point[i][j].upstarx - uc[0] );
    f[0][1] = - a1 * ny1 * p + mid * ( point[i][j].upstary - uc[1] );

    mid = nx2 * e[0] + ny2 * e[1];
    mid = mu * abs(mid) * a2;

    f[1][0] = - a2 * nx2 * p + mid * ( point[i][j].upstarx - uc[0] );
    f[1][1] = - a2 * ny2 * p + mid * ( point[i][j].upstary - uc[1] );

    delete[] uc;
    delete[] e;
    return f;
}