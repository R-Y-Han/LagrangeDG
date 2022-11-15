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
#include "limiter.h"
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
        
            int qt;
        //********下面判断节点(i,j)在单元(k,l)中是第几个顶点*******//
        if (j==0)
            {
                if (l == m-1)
                {
                    qt = point[i][m].q;
                }
                else{
                    qt = point[i][j].q;
                }
            }
            else if (j == m)
            {
                if (l == 0)
                {
                    qt = point[i][0].q;
                }
                else{
                    qt = point[i][j].q;
                }
            }
            else{
                qt = point[i][j].q;
            }
        for (temp=0; temp<4; temp++)
        {
            if (o[k][l].vertex[temp] == qt)
            {
                //此时循环指标temp刚好就是节点(i,j)在(k,l)单元中的局部位置
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];
        //**********下面计算该单元内节点处速度重构值*******//
        double uxt, uyt;
        uxt = 0;
        uyt = 0;
        for (temp = 0; temp < dim; temp++)
        {
            uxt = uxt + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
            uyt = uyt + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
        }
        //uxt = ux_limiter(k,l,uxt,0.6);
        //uyt = uy_limiter(k,l,uyt,0.6);
        uavg[0] = uavg[0] + uxt;
        uavg[1] = uavg[1] + uyt;
    }
    //********此时uavg记录了所有单元在该点处速度重构值的总和********//
    uavg[0] = uavg[0] / (double) point[i][j].neighbor_element.size();
    uavg[1] = uavg[1] / (double) point[i][j].neighbor_element.size();
    //********uavg为节点处速度重构值的平均****************//

    return uavg;
}

/**
 * @brief 计算边界节点的速度
 * 
 * @param i 节点的行编号
 * @param j 节点的列编号
 * @return double* 
 */
double * nodal_velocity_boundary(int i, int j)
{
    double *ustar = new double [2];
    
    //********先求解四个顶点*******//
    if (i == 0 || i == n)
    {
        ustar[0] = 0;
        ustar[1] = 0;
        return ustar;
    }

    //********边上的点*********//
    int k, l, r, temp;
    int node1, node2, bnode1, bnode2, p1, p2;
    double nstarx, nstary;
    double xit, etat, nx1, ny1, nx2, ny2, a1, a2;
    double mu;
    //********下面先找出边界的外法向*******//
    for (r=0; r<point[i][j].neighbor_node.size(); r++)
    {
        k = point[i][j].neighbor_node[r] / (m+1);
        l = point[i][j].neighbor_node[r] % (m+1);
        if (point[k][l].boundary == 1)
        {
            bnode1 = point[k][l].q;
            break;
        }
    }
    for (r=0; r<point[i][j].neighbor_node.size(); r++)
    {
        k = point[i][j].neighbor_node[r] / (m+1);
        l = point[i][j].neighbor_node[r] % (m+1);
        if (point[k][l].boundary == 1)
        {
            if (point[k][l].q != bnode1)
            {
                bnode2 = point[k][l].q;
                break;
            }
        }
    }
    
    nstarx = 0;
    nstary = 0;
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

        nx1 = - point[i][j].nx[p1];
        ny1 = - point[i][j].ny[p1]; //注意此处加了负号，为使指向外向

        nx2 = point[i][j].nx[p2];
        ny2 = point[i][j].ny[p2];

        int kt, lt;
        kt = node1 / (m+1);
        lt = node1 % (m+1);
        if (point[kt][lt].boundary == 1)
        {
            nstarx = nstarx + nx1;
            nstary = nstary + ny1;
        }
        kt = node2 / (m+1);
        lt = node2 % (m+1);
        if (point[kt][lt].boundary == 1)
        {
            nstarx = nstarx + nx2;
            nstary = nstary + ny2;
        }
    }
    nstarx = nstarx / 2.0;
    nstary = nstary / 2.0;


    double ** M = new double * [3];
    double ** M_inv = new double * [3];
    double * R = new double [3];
    for (r=0; r<3; r++)
    {
        M[r] = new double [3];
        M_inv[r] = new double [3];
        R[r] = 0;
    }
    for (k=0; k<3; k++)
    {
        for (l=0; l<3; l++)
        {
            M[k][l] = 0;
            if (k == l)
            {
                M_inv[k][l] = 1;
            }
            else{
                M_inv[k][l] = 0;
            }
        }
    }
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
        double xt_v, yt_v;
        xt_v = o[k][l].phi_x0(o[k][l].xi_c, o[k][l].eta_c);
        yt_v = o[k][l].phi_y0(o[k][l].xi_c, o[k][l].eta_c);

        //********下面计算密度**********//
        double rho_vertex;
        rho_vertex = ini_rho(xt_v,yt_v) * o[k][l].Jacobi_0(o[k][l].xi_c, o[k][l].eta_c)
                                        / o[k][l].Jacobi(o[k][l].xi_c, o[k][l].eta_c);
        
        //********下面计算节点处速度和总能重构值********//
        double * uc = new double [2];
        double ve;
        ve=0;
        uc[0] = 0;
        uc[1] = 0;
        for (temp = 0; temp < dim; temp++)
        {
            uc[0] = uc[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
            uc[1] = uc[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
            ve = ve + o[k][l].taulast[temp] * o[k][l].Psi(temp,xit,etat);
        }
        if (dim == 3)
        {
            uc[0] = ux_limiter(k,l,uc[0],ux_alpha);
            uc[1] = uy_limiter(k,l,uc[1],uy_alpha);
            ve = tau_limiter(k,l,ve,tau_alpha);
        }
        //********下面计算内能**********//
        ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        ve = max(1e-12,ve);
        double p_vertex;
        p_vertex = EOS(rho_vertex,ve);

        //**************下面计算阻抗*************//
        double rho_center;
        double xtc, ytc;
        xtc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
        ytc = o[k][l].phi_y0(o[k][l].xi_c,o[k][l].eta_c);
        rho_center = ini_rho(xtc,ytc) * o[k][l].Jacobi_0(o[k][l].xi_c,o[k][l].eta_c)
                                      / o[k][l].Jacobi(o[k][l].xi_c,o[k][l].eta_c);

        double p_center, ve_center;
        ve_center = o[k][l].uxlast[0]*o[k][l].uxlast[0] + o[k][l].uylast[0] * o[k][l].uylast[0];
        ve_center = o[k][l].taulast[0] - ve_center * 0.5;
        p_center = EOS(rho_center,ve_center);
        p_center = max(1e-12, p_center);
        mu = gamma * p_center / rho_center;
        mu = sqrt(mu);  //此时mu为声速
        mu = rho_center * mu;
        
        //*********下面计算分子分母各分量********//
        
        double mid;
        //*********下面计算第一条边的量***********//

        M[0][0] = M[0][0] + mu * a1 * nx1 * nx1;
        M[0][1] = M[0][1] + mu * a1 * nx1 * ny1;
        M[1][0] = M[1][0] + mu * a1 * ny1 * nx1;
        M[1][1] + M[1][1] + mu * a1 * ny1 * ny1;

        mid = nx1 * uc[0] + ny1 * uc[1];
        
        R[0] = R[0] + a1 * nx1 * p_vertex + mu * a1 * nx1 * mid;
        R[1] = R[1] + a1 * ny1 * p_vertex + mu * a1 * ny1 * mid;

        //************下面计算第二条边的量*******//
        
        M[0][0] = M[0][0] + mu * a2 * nx2 * nx2;
        M[0][1] = M[0][1] + mu * a2 * nx2 * ny2;
        M[1][0] = M[1][0] + mu * a2 * ny2 * nx2;
        M[1][1] + M[1][1] + mu * a2 * ny2 * ny2;

        mid = nx2 * uc[0] + ny2 * uc[2];
        
        R[0] = R[0] + a2 * nx2 * p_vertex + mu * a2 * nx2 * mid;
        R[1] = R[1] + a2 * ny2 * p_vertex + mu * a2 * ny2 * mid;

        delete[] uc;
    }

    M[0][2] = nstarx;
    M[1][2] = nstary;
    M[2][0] = nstarx;
    M[2][1] = nstary;
    M[2][2] = 0;
    R[2] = 0;

    for (k=1; k<3; k++)
    {
        M[k][0] = M[k][0] / M[0][0];
    }
    for (k=1; k<3; k++)
    {
        for (l=k; l<3; l++)
        {
            for (r=0; r<k; r++)
            {
                M[k][l] = M[k][l] - M[k][r] * M[r][l];
            }
        }
        for (l= k+1; l<3; l++)
        {
            for (r=0; r<k; r++)
            {
                M[l][k] = M[l][k] - M[l][r] * M[r][k];
            }
            M[l][k] = M[l][k] / M[k][k];
        }

        for (l=0; l<3; l++)
        {
            R[k] = R[k] - M[k][l] * R[l];
        }
    }

    R[2] = R[2] / M[2][2];
    for (k=1; k>=0; k--)
    {
        for (l=2; l>k; l--)
        {
            R[k] = R[k] - M[k][l] * R[l];
        }
        R[k] = R[k] / M[k][k];
    }

    ustar[0] = R[0];
    ustar[1] = R[1];
    
    for (r=0; r<3; r++)
    {
        delete[] M[r];
        delete[] M_inv[r];
    }
    delete[] M;
    delete[] M_inv;
    delete[] R;
    return ustar;
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
    double* ustar;
    double ** Mpn = new double * [2];
    Mpn[0] = new double [2];
    Mpn[1] = new double [2];
    for (int k=0; k<2; k++)
    {
        for (int l=0; l<2; l++)
        {
            Mpn[k][l] = 0;
        }
    }

    //**************先判断是否为边界点*********//
    if (point[i][j].boundary == 1)
    {
        //******是边界点*******//
        ustar = nodal_velocity_boundary(i,j);
        return ustar;
    }

    //*************不是边界节点***************//
    ustar = new double [2];
    int r, temp, k, l, node1, node2, p1, p2;
    double xit, etat, nx1, ny1, nx2, ny2, a1, a2;
    double mu, mid;
    double rhs[2];
    rhs[0] = 0;
    rhs[1] = 0;
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
        double xt_v, yt_v;
        xt_v = o[k][l].phi_x(xit,etat);
        yt_v = o[k][l].phi_y(xit,etat);

        //********下面计算密度**********//
        double rho_vertex;
        double xt0, yt0;
        xt0 = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
        yt0 = o[k][l].phi_y0(o[k][l].xi_c,o[k][l].eta_c);
        rho_vertex = ini_rho(xt0,yt0) * o[k][l].Jacobi_0(o[k][l].xi_c,o[k][l].eta_c)
                                      / o[k][l].Jacobi(o[k][l].xi_c,o[k][l].eta_c);

        //********下面计算节点处速度和总能重构值********//
        double * uc = new double [2];
        double ve;
        ve=0;
        uc[0] = 0;
        uc[1] = 0;
        for (temp = 0; temp < dim; temp++)
        {
            uc[0] = uc[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
            uc[1] = uc[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
            ve = ve + o[k][l].taulast[temp] * o[k][l].Psi(temp,xit,etat);
        }
        if (dim == 3)
        {
            uc[0] = ux_limiter(k,l,uc[0],0.6);
            uc[1] = uy_limiter(k,l,uc[1],0.6);
            ve = tau_limiter(k,l,ve,0.6);
        }
        

        //********下面计算内能**********//
        ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        ve = max(1e-9, ve);
        double p_vertex;
        p_vertex = EOS(rho_vertex,ve);

        //**************下面计算阻抗*************//
        //mu = (gamma - 1) * ve;
        //mu = sqrt(mu);  //此时mu为声速
        //mu = rho * mu;
        double rho_center;
        double xtc, ytc;
        xtc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
        ytc = o[k][l].phi_y0(o[k][l].xi_c,o[k][l].eta_c);
        rho_center = ini_rho(xtc,ytc) * o[k][l].Jacobi_0(o[k][l].xi_c,o[k][l].eta_c)
                                      / o[k][l].Jacobi(o[k][l].xi_c,o[k][l].eta_c);

        double p_center, ve_center;
        ve_center = o[k][l].uxlast[0]*o[k][l].uxlast[0] + o[k][l].uylast[0] * o[k][l].uylast[0];
        ve_center = o[k][l].taulast[0] - ve_center * 0.5;
        p_center = EOS(rho_center,ve_center);
        p_center = max(1e-12, p_center);
        mu = gamma * p_center / rho_center;
        mu = sqrt(mu);  //此时mu为声速
        mu = rho_center * mu;
        
        //**********下面计算第一条边分量*******//
        mid = nx1 * uc[0] + ny1 * uc[1];
        rhs[0] = rhs[0] + a1 * p_vertex * nx1 + mu * a1 * mid * nx1;
        rhs[1] = rhs[1] + a1 * p_vertex * ny1 + mu * a1 * mid * ny1;

        Mpn[0][0] = Mpn[0][0] + mu * a1 * nx1 * nx1;
        Mpn[0][1] = Mpn[0][1] + mu * a1 * ny1 * nx1;
        Mpn[1][0] = Mpn[1][0] + mu * a1 * nx1 * ny1;
        Mpn[1][1] = Mpn[1][1] + mu * a1 * ny1 * ny1;


        //*********下面计算第二条边分量********//
        mid = nx2 * uc[0] + ny2 * uc[1];
        rhs[0] = rhs[0] + a2 * p_vertex * nx2 + mu * a2 * mid * nx2;
        rhs[1] = rhs[1] + a2 * p_vertex * ny2 + mu * a2 * mid * ny2;

        Mpn[0][0] = Mpn[0][0] + mu * a2 * nx2 * nx2;
        Mpn[0][1] = Mpn[0][1] + mu * a2 * ny2 * nx2;
        Mpn[1][0] = Mpn[1][0] + mu * a2 * nx2 * ny2;
        Mpn[1][1] = Mpn[1][1] + mu * a2 * ny2 * ny2;
        
        delete[] uc;
    }

    ustar[0] = 0;
    ustar[1] = 0;
    double detMpn;
    detMpn = Mpn[0][0] * Mpn[1][1] - Mpn[0][1] * Mpn[1][0];
    
    ustar[0] = Mpn[1][1] * rhs[0] - Mpn[0][1] * rhs[1];
    ustar[0] = ustar[0] / detMpn;
    
    ustar[1] = -Mpn[1][0] * rhs[0] + Mpn[0][0] * rhs[1];
    ustar[1] = ustar[1] / detMpn;

    delete[] Mpn[0];
    delete[] Mpn[1];
    delete[] Mpn;
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

    //********下面计算密度**********//
    double rho_vertex;
    double xt0, yt0;
    xt0 = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
    yt0 = o[k][l].phi_y0(o[k][l].xi_c,o[k][l].eta_c);
    rho_vertex = ini_rho(xt0,yt0) * o[k][l].Jacobi_0(o[k][l].xi_c,o[k][l].eta_c)
                                  / o[k][l].Jacobi(o[k][l].xi_c,o[k][l].eta_c);
    

    //********下面计算节点处速度重构值********//
    double * uc = new double [2];
    double ve;
    ve = 0;
    uc[0] = 0;
    uc[1] = 0;
    for (temp = 0; temp < dim; temp++)
    {
        uc[0] = uc[0] + o[k][l].uxlast[temp] * o[k][l].Psi(temp,xit,etat);
        uc[1] = uc[1] + o[k][l].uylast[temp] * o[k][l].Psi(temp,xit,etat);
        ve = ve + o[k][l].taulast[temp] * o[k][l].Psi(temp,xit,etat);
    }
    if (dim == 3)
    {
        uc[0] = ux_limiter(k,l,uc[0],0.6);
        uc[1] = uy_limiter(k,l,uc[1],0.6);
        ve = tau_limiter(k,l,ve,0.6);
    }
        
    //********下面计算内能**********//
    ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
    ve = max(1e-9, ve);
    double p_vertex;
    p_vertex = EOS(rho_vertex,ve);

    //**************下面计算阻抗*************//
        double rho_center;
        double xtc, ytc;
        xtc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
        ytc = o[k][l].phi_y0(o[k][l].xi_c,o[k][l].eta_c);
        rho_center = ini_rho(xtc,ytc) * o[k][l].Jacobi_0(o[k][l].xi_c,o[k][l].eta_c)
                                      / o[k][l].Jacobi(o[k][l].xi_c,o[k][l].eta_c);

        double p_center, ve_center;
        ve_center = o[k][l].uxlast[0]*o[k][l].uxlast[0] + o[k][l].uylast[0] * o[k][l].uylast[0];
        ve_center = o[k][l].taulast[0] - ve_center * 0.5;
        p_center = EOS(rho_center,ve_center);
        p_center = max(1e-12, p_center);
        mu = gamma * p_center / rho_center;
        mu = sqrt(mu);  //此时mu为声速
        mu = rho_center * mu;
    
    //***********下面计算两条边上的corner force*********//
    double mid;
    mid = (point[i][j].upstarx - uc[0]) * nx1 + (point[i][j].upstary - uc[1]) * ny1;

    f[0][0] = - a1 * nx1 * p_vertex + a1 * mu * mid * nx1;
    f[0][1] = - a1 * ny1 * p_vertex + a1 * mu * mid * ny1;

    mid = (point[i][j].upstarx - uc[0]) * nx2 + (point[i][j].upstary - uc[1]) * ny2;

    f[1][0] = - a2 * nx2 * p_vertex + a2 * mu * mid * nx2;
    f[1][1] = - a2 * ny2 * p_vertex + a2 * mu * mid * ny2;

    delete[] uc;
    return f;
}