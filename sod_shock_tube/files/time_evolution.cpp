/**
 * @file time_evolution.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 一次时间步长的更新，包括密度、速度、质量总能的更新和时间步长的选取
 * @version v1.01
 * @date 2022-10-07
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
 * 
 */

#include "time_evolution.h"
#include "config.h"
#include "nodal_solver.h"
#include "limiter.h"
#include <cmath>

#include <iostream>
using namespace std;


/**
 * @brief 得到单元速度的右端矩阵
 * 
 * @param i 单元的行编号
 * @param j 单元的列编号
 * @return double** 
 */
double ** velocity_matrix(int i, int j)
{
    double ** R;
    int k, l, r, temp;
    double xit, etat, rhot, jt, et, mid;
    R = new double * [dim];
    for (temp = 0; temp < dim; temp++)
    {
        R[temp] = new double [2];
    }
    for ( k=0; k<dim; k++)
    {
        for (l=0; l<2; l++)
        {
            R[k][l] = 0;
        }
    }

    //*********下面计算边界通量*********//
  
    //*********对每个顶点计算***********//
    for (r = 0; r<4; r++)
    {
        k = o[i][j].vertex[r] / (m+1);
        l = o[i][j].vertex[r] % (m+1);  //节点的行列编号

        xit = ref_xi[r][0];
        etat = ref_xi[r][1];    //节点的参考空间坐标

        double ** f = corner_force(k,l,o[i][j].q);

        for ( temp = 0; temp < dim; temp++)
        {
            R[temp][0] = R[temp][0] + o[i][j].Psi(temp,xit,etat) * (f[0][0] + f[1][0]);
            R[temp][1] = R[temp][1] + o[i][j].Psi(temp,xit,etat) * (f[0][1] + f[1][1]);
        }
        
        for (temp=0; temp<2; temp++)
        {
            delete[] f[temp];
        }
        delete[] f;
    }

    //**********下面计算体积分**********//
    for (temp = 0; temp < gd; temp++)
    {
        xit = Gausspoint_xi[temp];
        etat = Gausspoint_eta[temp];
        
        //rhot = o[i][j].rholast;
        double xt0, yt0;
        xt0 = o[i][j].phi_x0(xit,etat);
        yt0 = o[i][j].phi_y0(xit,etat);
        if (xt0 == 0.5)
        {
            double xc;
            xc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
            xt0 = xt0 - 1e-9 * (xt0 - xc);
        }
        rhot = ini_rho(xt0,yt0) * o[i][j].Jacobi_0(xit,etat)
                    / o[i][j].Jacobi(xit,etat);
        jt = o[i][j].Jacobi(xit,etat);
        
        //*********下面计算压强********//
        double* ut = new double [2];
        et = 0;
        ut[0] = 0;
        ut[1] = 0;
        for (r=0; r<dim; r++)
        {
            et = et + o[i][j].taulast[r] * o[i][j].Psi(r,xit,etat);
            ut[0] = ut[0] + o[i][j].uxlast[r] * o[i][j].Psi(r,xit,etat);
            ut[1] = ut[1] + o[i][j].uylast[r] * o[i][j].Psi(r,xit,etat);
        }
        et = tau_limiter(i,j,et,0.6);
        ut[0] = ux_limiter(i,j,ut[0],0.6);
        ut[1] = uy_limiter(i,j,ut[1],0.6);
        et = et - (ut[0] * ut[0] + ut[1] * ut[1]) * 0.5;
        et = max(1e-12, et);
        double p_vertex;
        p_vertex = EOS(rhot,et);

        delete[] ut;

        //***********下面计算J的逆矩阵**********//
        double** J_inv = new double * [2];
        for (r=0; r<2; r++)
        {
            J_inv[r] = new double [2];
        }
        J_inv[0][0] = 0;
        J_inv[0][1] = 0;
        J_inv[1][0] = 0;
        J_inv[1][1] = 0;
        for (r=0; r<4; r++)
        {
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,xit,etat);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,xit,etat);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,xit,etat);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,xit,etat);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;

        for (r=0; r<dim; r++)
        {
            mid = o[i][j].Psi_xi(r,xit,etat) * J_inv[0][0]
                + o[i][j].Psi_eta(r,xit,etat) * J_inv[1][0];
            mid = - mid * jt * p_vertex;
            R[r][0] = R[r][0] - mid * Gaussweight[temp];

            mid = o[i][j].Psi_xi(r,xit,etat) * J_inv[0][1]
                + o[i][j].Psi_eta(r,xit,etat) * J_inv[1][1];
            mid = - mid * jt * p_vertex;
            R[r][1] = R[r][1] - mid * Gaussweight[temp];
        }

        for (r=0; r<2; r++)
        {
            delete[] J_inv[r];
        }
        delete[] J_inv;
    }

    double ** Rt = new double * [2];
    Rt[0] = new double [dim];
    Rt[1] = new double [dim];
    for (r=0; r<dim; r++)
    {
        Rt[0][r] = R[r][0];
        Rt[1][r] = R[r][1];
    }
    
    for (r = 0; r<dim; r++)
    {
        delete[] R[r];
    }
    delete[] R;

    return Rt;
}

/**
 * @brief 得到单元质量总能的右端矩阵
 * 
 * @param i 单元的行编号
 * @param j 单元的列编号
 * @return double* 
 */
double * energy_matrix(int i, int j)
{
    double * R = new double [dim];
    int k, l, r, temp;
    double xit, etat, rhot, jt, et, mid;
    for (r=0; r<dim; r++)
    {
        R[r] = 0;
    }

    //*********下面计算边界通量*********//
    //*********对每个顶点计算***********//
    for (r = 0; r<4; r++)
    {
        k = o[i][j].vertex[r] / (m+1);
        l = o[i][j].vertex[r] % (m+1);  //节点的行列编号

        xit = ref_xi[r][0];
        etat = ref_xi[r][1];    //节点的参考空间坐标

        double ** f = corner_force(k,l,o[i][j].q);

        for ( temp = 0; temp < dim; temp++)
        {
            mid = f[0][0] * point[k][l].upstarx + f[0][1] * point[k][l].upstary
                + f[1][0] * point[k][l].upstarx + f[1][1] * point[k][l].upstary;
            R[temp] = R[temp] + o[i][j].Psi(temp,xit,etat) * mid;
        }
        
        for (temp=0; temp<2; temp++)
        {
            delete[] f[temp];
        }
        delete[] f;
    }
    
    //**********下面计算体积分**********//
    for (temp = 0; temp < gd; temp++)
    {
        xit = Gausspoint_xi[temp];
        etat = Gausspoint_eta[temp];
        
        //rhot = o[i][j].rholast;
        double xt0, yt0;
        xt0 = o[i][j].phi_x0(xit,etat);
        yt0 = o[i][j].phi_y0(xit,etat);
        if (xt0 == 0.5)
        {
            double xc;
            xc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
            xt0 = xt0 - 1e-9 * (xt0 - xc);
        }
        rhot = ini_rho(xt0,yt0) * o[i][j].Jacobi_0(xit,etat)
                    / o[i][j].Jacobi(xit,etat);
        jt = o[i][j].Jacobi(xit,etat);

        //*********下面计算压强********//
        double* ut = new double [2];
        et = 0;
        ut[0] = 0;
        ut[1] = 0;
        for (r=0; r<dim; r++)
        {
            et = et + o[i][j].taulast[r] * o[i][j].Psi(r,xit,etat);
            ut[0] = ut[0] + o[i][j].uxlast[r] * o[i][j].Psi(r,xit,etat);
            ut[1] = ut[1] + o[i][j].uylast[r] * o[i][j].Psi(r,xit,etat);
        }
        et = tau_limiter(i,j,et,0.6);
        ut[0] = ux_limiter(i,j,ut[0],0.6);
        ut[1] = uy_limiter(i,j,ut[1],0.6);
        et = et - (ut[0] * ut[0] + ut[1] * ut[1]) * 0.5;
        et = max(1e-12, et);
        double p_vertex;
        p_vertex = EOS(rhot,et);

        //***********下面计算J的逆矩阵**********//
        double** J_inv = new double * [2];
        for (r=0; r<2; r++)
        {
            J_inv[r] = new double [2];
        }

        J_inv[0][0] = 0;
        J_inv[0][1] = 0;
        J_inv[1][0] = 0;
        J_inv[1][1] = 0;
        for (r=0; r<4; r++)
        {
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,xit,etat);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,xit,etat);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,xit,etat);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,xit,etat);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;

        double t1, t2;
        for (r=0; r<dim; r++)
        {
            t1 = o[i][j].Psi_xi(r,xit,etat) * J_inv[0][0]
               + o[i][j].Psi_eta(r,xit,etat) * J_inv[0][1];
            t2 = o[i][j].Psi_xi(r,xit,etat) * J_inv[1][0]
               + o[i][j].Psi_eta(r,xit,etat) * J_inv[1][1];
            mid = t1 * ut[0] + t2 * ut[1];
            mid = - mid * jt * p_vertex;
            R[r] = R[r] - mid * Gaussweight[temp];
        }

        for (r=0; r<2; r++)
        {
            delete[] J_inv[r];
        }
        delete[] J_inv;
        delete[] ut;
    }


    return R;
}

/**
 * @brief RK第一步计算
 * 
 * @param un n时刻原始向量
 * @param M 系数矩阵
 * @param Rn 右端矩阵
 * @param timeStep 时间步长
 * @param d 向量维数
 * @return double* 
 * 
 * @note 
 * - 输入的方程为M du/dt = R
 * - 输入的M在计算逆矩阵时会变成单位矩阵，若要避免此问题需另设一个指针
 */
double * RK_step1(double * un, double ** M, double * Rn, double timeStep, int d)
{
    double * us1 = new double [d];
    int i, j, k, l;
    double temp;

    //*************下面计算M的逆矩阵**********//
    double ** M_inv = new double * [d];
    for (i=0; i<d; i++)
    {
        M_inv[i] = new double [d];
    }
    for (i=0; i<d; i++)
    {
        for (j=0; j<d; j++)
        {
            M_inv[i][j] = 0;
        }
    }

    if (d == 1)
    {
        M_inv[0][0] = 1.0 / M[0][0];
    }
    if (d == 3)
    {
        M_inv[0][0] = 1.0 / M[0][0];
        double det;
        det = M[1][1] * M[2][2] - M[1][2] * M[2][1];
        M_inv[1][1] = M[2][2] / det;
        M_inv[1][2] = - M[1][2] / det;
        M_inv[2][1] = - M[2][1] / det;
        M_inv[2][2] = M[1][1] / det;
    }
    

    for (i=0; i<d; i++)
    {
        us1[i] = 0;
    }
    for (i=0; i<d; i++)
    {
        for (j=0; j<d; j++)
        {
            us1[i] = us1[i] + M_inv[i][j] * Rn[j];
        }
    }
    for (i=0; i<d; i++)
    {
        us1[i] = us1[i] * timeStep + un[i];
    }


    for (i=0; i<d; i++)
    {
        delete[] M_inv[i];
    }
    delete[] M_inv;

    return us1;
}

/**
 * @brief RK第二步计算
 * 
 * @param un n时刻原始向量
 * @param us1 第一次RK结果
 * @param M 系数矩阵
 * @param Rs1 第一次RK后重新计算的右端矩阵
 * @param timeStep 时间步长
 * @param d 向量维数
 * @return double* 
 * 
 *  * @note 
 * - 输入的方程为M du/dt = R
 * - 输入的M在计算逆矩阵时会变成单位矩阵，若要避免此问题需另设一个指针
 */
double * RK_step2(double * un, double * us1, double ** M, double * Rs1, double timeStep, int d)
{
    double * us2 = new double [d];
    int i, j, k, l;
    double temp;

    //*************下面计算M的逆矩阵**********//
    double ** M_inv = new double * [d];
    for (i=0; i<d; i++)
    {
        M_inv[i] = new double [d];
    }
    for (i=0; i<d; i++)
    {
        for (j=0; j<d; j++)
        {
            M_inv[i][j] = 0;
        }
    }

    if (d == 1)
    {
        M_inv[0][0] = 1.0 / M[0][0];
    }
    if (d == 3)
    {
        M_inv[0][0] = 1.0 / M[0][0];
        double det;
        det = M[1][1] * M[2][2] - M[1][2] * M[2][1];
        M_inv[1][1] = M[2][2] / det;
        M_inv[1][2] = - M[1][2] / det;
        M_inv[2][1] = - M[2][1] / det;
        M_inv[2][2] = M[1][1] / det;
    }

    for (i=0; i<d; i++)
    {
        us2[i] = 0;
    }
    for (i=0; i<d; i++)
    {
        for (j=0; j<d; j++)
        {
            us2[i] = us2[i] + M_inv[i][j] * Rs1[j];
        }   
    }

    for (i=0; i<d; i++)
    {
        us2[i] = 0.5 * timeStep * us2[i] + 0.5 * un[i] + 0.5 * us1[i];
    }

    for (i=0; i<d; i++)
    {
        delete[] M_inv[i];
    }
    delete[] M_inv;

    return us2;
}

/**
 * @brief 选择时间步长dt，
 * 使得每个单元内声速运动距离都不超过单元的某种长度度量
 * 
 * @param dtn 之前的时间步长
 * @return double
 * 
 */
double choose_dt(double dtn)
{
    double timeStep;

    int i, j, k, l, r;
    double et, len, c, temp;

    //************对每个单元计算CFL时间步长**********//
    timeStep = 10000;
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            //*****声速取为质心的声速*********//
            et = o[i][j].uxlast[0] * o[i][j].uxlast[0];
            et = et + o[i][j].uylast[0] * o[i][j].uylast[0];
            et = o[i][j].taulast[0] - 0.5 * et;
            
            c = (gamma - 1) * et;
            c = sqrt(c);
            
            //*******下面计算最短的边长********//
            len = 100000;
            for (r=0; r<4; r++)
            {
                k = o[i][j].vertex[r] / (m+1);
                l = o[i][j].vertex[r] % (m+1);

                int node1, node2, p1, p2;
                if (r == 3)
                {
                    node1 = o[i][j].vertex[0];
                }
                else{
                    node1 = o[i][j].vertex[r+1];
                }

                if (r == 0)
                {
                    node2 = o[i][j].vertex[3];
                }
                else{
                    node2 = o[i][j].vertex[r-1];
                }

                for (p1 = 0; p1 < point[k][l].neighbor_node.size(); p1++)
                {
                    if (node1 == point[k][l].neighbor_node[p1])
                    {
                        break;
                    }
                }
                for (p2 = 0; p2 < point[k][l].neighbor_node.size(); p2++)
                {
                    if (node2 == point[k][l].neighbor_node[p2])
                    {
                        break;
                    }
                }
                
                temp = point[k][l].a[p1] * 2;
                if (temp < len)
                {
                    len = temp;
                }
                temp = point[k][l].a[p2] * 2;
                if (temp < len)
                {
                    len = temp;
                }
            }

            temp = len / c;

            if (temp < timeStep)
            {
                timeStep = temp;
            }
        }
    }

    timeStep = timeStep * c_CFL;

    if (timeStep > dtn)
    {
        timeStep = 1.01 * dtn;
    }

    return timeStep;
}