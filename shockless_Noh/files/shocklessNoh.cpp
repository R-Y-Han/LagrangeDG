/**
* @mainpage  A Lagrangian Discontinuous Galerkin Hydrodynamic Method
* <table>
* <tr><th>Project  <td> Code implementation for A Lagrangian Discontinuous Galerkin Hydrodynamic Method
* <tr><th>Author   <td> R.Y. Han
* <tr><th>Institution   <td> IAPCM
* </table>
* @section   Project details
* Code the method from paper A Lagrangian Discontinuous Galerkin Hydrodynamic Method on several testproblems.
* 
**********************************************************************************
*/
/** 
* @file [shocklessNoh.cpp]
* @brief Implementation for 2D shockless Noh test problem.
* @author R.Y. Han,  IAPCM
* @email: hanruoyu21@gscaep.ac.cn
* @date 9.29.2022
* @version v1.01
*
* @details To compute the states in the star region, and identify two nonlinear waves.   \n 
* Plot the result at a given finish time T, in a region I = [a,b].             \n 
* @htmlonly 
* <span style="font-weight: bold">History</span> 
* @endhtmlonly 
* Project|Version|Auther|Date|Describe
* -------------|------|----|------|-------- 
* A Lagrangian \n Discontinuous Galerkin \n Hydrodynamic Method |V1.01|R.Y. Han|9.29.2022| 2D shockless Noh
* @date 2022-05-14
* 
* @copyright Copyright (c) 2022 R.Y. Han. All Rights Reserved.
*/ 

#include <iostream>
#include <cmath>
#include "config.h"
#include "initial.h"
#include "nodal_solver.h"
#include "time_evolution.h"
#include "plot.h"
#include <fstream>

using namespace std;

/**
 * @brief 进行一步时间推进
 * 
 */
void onetimestep();

void onetimestep2();

/**
 * @brief 计算t时刻内能的L2误差
 * 
 * @param t 终止时间
 * @return double 
 */
double inter_ene_difference(double t);

int main()
{
    int i, j;
    double time = 0;
    initial();

    while(time < T)
    {
        dt = choose_dt(dt);
        if (time + dt > T)
        {
            dt = T - time;
        }
        time = time + dt;
        onetimestep2();
        cout<<time<<endl;
    }

    double norm;
    norm = inter_ene_difference(T);
    cout<<"norm="<<norm<<endl;

    const char* fn = "F:\\C++Codes\\LagrangeDG\\shockless_Noh\\output\\norm.txt";
    fstream f;
    f.open(fn, ios::out | ios::app);
    f<<"n="<<n<<"\t"<<"m="<<m<<"\t"<<"CFL number="<<c_CFL<<endl;
    f<<"norm="<<norm<<endl<<endl;
    f.close();


    plotmesh();
    plotpressure();
    plotinternalenergy();
    system("pause");
}

void onetimestep()
{
    int i, j, k, l, r;
    
    //********下面进行第一次RK计算**********//
    //*********下面对每个节点求解节点速度*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            utemp = nodal_velocity(i,j);
            point[i][j].upstarx = utemp[0];
            point[i][j].upstary = utemp[1];
            delete[] utemp;
        }
    }

    //********下面进行第一次单元速度和质量总能的计算*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            double ** Rvs1;
            double * Res1;
            Rvs1 = velocity_matrix(i,j);
            Res1 = energy_matrix(i,j);

            double ** Mtemp = new double * [dim];
            for (k=0; k<dim; k++)
            {
                Mtemp[k] = new double [dim];
            }

            double * uxs1;
            double * uys1;
            double * taus1;
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uxs1 = RK_step1(o[i][j].uxlast,Mtemp,Rvs1[0],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uys1 = RK_step1(o[i][j].uylast,Mtemp,Rvs1[1],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            taus1 = RK_step1(o[i][j].taulast,Mtemp,Res1,dt,dim);

            for (k=0; k<dim; k++)
            {
                o[i][j].uxnext[k] = o[i][j].uxlast[k];
                o[i][j].uynext[k] = o[i][j].uylast[k];
                o[i][j].taunext[k] = o[i][j].taulast[k];

                o[i][j].uxlast[k] = uxs1[k];
                o[i][j].uylast[k] = uys1[k];
                o[i][j].taulast[k] = taus1[k];
            }


            for (k=0; k<dim; k++)
            {
                delete[] Mtemp[k];
            }
            delete[] Rvs1[0];
            delete[] Rvs1[1];
            delete[] Rvs1;
            delete[] Res1;
            delete[] Mtemp;
            delete[] uxs1;
            delete[] uys1;
            delete[] taus1;

        }
    }
    //注意此时tn时刻的值保存在next中，第一步RK的值保存在last中
    
    //**********下面进行第一次节点坐标的计算*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double* xlast = new double [2];
            xlast[0] = point[i][j].x;
            xlast[1] = point[i][j].y;
            double* uplast = new double [2];
            uplast[0] = point[i][j].upstarx;
            uplast[1] = point[i][j].upstary;

            double ** Mtemp = new double * [2];
            Mtemp[0] = new double [2];
            Mtemp[1] = new double [2];
            Mtemp[0][0] = 1;
            Mtemp[0][1] = 0;
            Mtemp[1][0] = 0;
            Mtemp[1][1] = 1;
            
            double* xnext;
            xnext = RK_step1(xlast,Mtemp,uplast,dt,2);
            point[i][j].xtemp = xnext[0];
            point[i][j].ytemp = xnext[1];

            delete[] Mtemp[0];
            delete[] Mtemp[1];
            delete[] Mtemp;
            delete[] xlast;
            delete[] xnext;
            delete[] uplast;
        }
    }
    //此时xy仍保存n时刻位置，temp保存第一次RK的位置


    //**********下面进行第二次RK计算***********//
    //*********下面对每个节点求解节点速度*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            //注意此时last指针记录的是第一步RK的结果
            utemp = nodal_velocity(i,j);
            point[i][j].upstarx = utemp[0];
            point[i][j].upstary = utemp[1];
            delete[] utemp;
        }
    }

    //********下面进行第二次单元速度和质量总能的计算*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            double ** Rvs2;
            double * Res2;
            Rvs2 = velocity_matrix(i,j);
            Res2 = energy_matrix(i,j);

            double ** Mtemp = new double * [dim];
            for (k=0; k<dim; k++)
            {
                Mtemp[k] = new double [dim];
            }

            double * uxs2;
            double * uys2;
            double * taus2;
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uxs2 = RK_step2(o[i][j].uxnext,o[i][j].uxlast,Mtemp,Rvs2[0],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uys2 = RK_step2(o[i][j].uynext,o[i][j].uylast,Mtemp,Rvs2[1],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            taus2 = RK_step2(o[i][j].taunext,o[i][j].taulast,Mtemp,Res2,dt,dim);

            for (k=0; k<dim; k++)
            {
                o[i][j].uxlast[k] = uxs2[k];
                o[i][j].uylast[k] = uys2[k];
                o[i][j].taulast[k] = taus2[k];
            }


            for (k=0; k<dim; k++)
            {
                delete[] Mtemp[k];
            }
            delete[] Rvs2[0];
            delete[] Rvs2[1];
            delete[] Rvs2;
            delete[] Res2;
            delete[] Mtemp;
            delete[] uxs2;
            delete[] uys2;
            delete[] taus2;

        }
    }
    //注意此时n+1时刻的值被保存在last中，next中无关紧要

    //**********下面进行第二次节点坐标的计算*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double* xlast = new double [2];
            xlast[0] = point[i][j].x;
            xlast[1] = point[i][j].y;
            double* xmid = new double [2];
            xmid[0] = point[i][j].xtemp;
            xmid[1] = point[i][j].ytemp;
            double * uplast = new double [2];
            uplast[0] = point[i][j].upstarx;
            uplast[1] = point[i][j].upstary;

            double ** Mtemp = new double * [2];
            Mtemp[0] = new double [2];
            Mtemp[1] = new double [2];
            Mtemp[0][0] = 1;
            Mtemp[0][1] = 0;
            Mtemp[1][0] = 0;
            Mtemp[1][1] = 1;
            
            double* xnext;
            xnext = RK_step2(xlast,xmid,Mtemp,uplast,dt,2);
            point[i][j].x = xnext[0];
            point[i][j].y = xnext[1];

            delete[] Mtemp[0];
            delete[] Mtemp[1];
            delete[] Mtemp;
            delete[] xlast;
            delete[] xmid;
            delete[] xnext;
            delete[] uplast;
        }
    }

    //**********下面更新节点到相邻点的长度和外法向量****//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double ax, ay, bx, by;
            ax = point[i][j].x;
            ay = point[i][j].y;
            for (r=0; r<point[i][j].neighbor_node.size(); r++)
            {
                k = point[i][j].neighbor_node[r] / (m + 1);
                l = point[i][j].neighbor_node[r] % (m + 1);
                
                bx = point[k][l].x;
                by = point[k][l].y;

                point[i][j].a[r] = length(ax,ay,bx,by) * 0.5;

                double * ntemp;
                ntemp = normal(ax,ay,bx,by);
                point[i][j].nx[r] = ntemp[0];
                point[i][j].ny[r] = ntemp[1];
                delete[] ntemp;                
            }
        }
    }

    //***********更新单元顶点位置*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (r=0; r<4; r++)
            {
                k = o[i][j].vertex[r] / (m + 1);
                l = o[i][j].vertex[r] % (m + 1);

                o[i][j].vx[r] = point[k][l].x;
                o[i][j].vy[r] = point[k][l].y;
            }
        }
    }

    return ;
}

double inter_ene_difference(double t)
{
    double ans, e_ana, e_r, temp;
    int i, j, k, l;
    double xit, etat, uxt, uyt, taut, jt;

    ans = 0;
    //*********下面计算每个单元内的L2误差**********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            temp = 0;
            //计算每个Gauss节点处的误差
            for (k=0; k<4; k++)
            {
                xit = Gausspoint_xi[k];
                etat = Gausspoint_eta[k];
                //********计算理论值*********//
                e_ana = ana_e(i,j,xit,etat,t);

                //**********计算重构值********//
                e_r = 0;
                uxt = 0;
                uyt = 0;
                taut = 0;
                for (l=0; l<dim; l++)
                {
                    uxt = uxt + o[i][j].uxlast[l] * o[i][j].Psi(l,xit,etat);
                    uyt = uyt + o[i][j].uylast[l] * o[i][j].Psi(l,xit,etat);
                    taut = taut + o[i][j].taulast[l] * o[i][j].Psi(l,xit,etat);
                }
                e_r = uxt * uxt + uyt * uyt;
                e_r = taut - 0.5 * e_r;
                
                jt = o[i][j].Jacobi(xit,etat);
                temp = temp + Gaussweight[k] * (e_ana - e_r) * (e_ana - e_r) * jt;
            }
            ans = ans + temp;
        }
    }
    ans = sqrt(ans);

    return ans;
}




void onetimestep2()
{
    int i, j, k, l, r;
    
    //********下面进行第一次RK计算**********//
    //*********下面对每个节点求解节点速度*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            utemp = nodal_velocity(i,j);
            point[i][j].upstarx = utemp[0];
            point[i][j].upstary = utemp[1];
            delete[] utemp;
        }
    }

    //********下面进行第一次单元速度和质量总能的计算*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            double ** Rvs1;
            double * Res1;
            Rvs1 = velocity_matrix(i,j);
            Res1 = energy_matrix(i,j);

            double ** Mtemp = new double * [dim];
            for (k=0; k<dim; k++)
            {
                Mtemp[k] = new double [dim];
            }

            double * uxs1;
            double * uys1;
            double * taus1;
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uxs1 = RK_step1(o[i][j].uxlast,Mtemp,Rvs1[0],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uys1 = RK_step1(o[i][j].uylast,Mtemp,Rvs1[1],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            taus1 = RK_step1(o[i][j].taulast,Mtemp,Res1,dt,dim);

            for (k=0; k<dim; k++)
            {
                o[i][j].uxnext[k] = o[i][j].uxlast[k];
                o[i][j].uynext[k] = o[i][j].uylast[k];
                o[i][j].taunext[k] = o[i][j].taulast[k];

                o[i][j].uxlast[k] = uxs1[k];
                o[i][j].uylast[k] = uys1[k];
                o[i][j].taulast[k] = taus1[k];
            }


            for (k=0; k<dim; k++)
            {
                delete[] Mtemp[k];
            }
            delete[] Rvs1[0];
            delete[] Rvs1[1];
            delete[] Rvs1;
            delete[] Res1;
            delete[] Mtemp;
            delete[] uxs1;
            delete[] uys1;
            delete[] taus1;

        }
    }
    //注意此时tn时刻的值保存在next中，第一步RK的值保存在last中
    
    //**********下面进行第一次节点坐标的计算*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double* xlast = new double [2];
            xlast[0] = point[i][j].x;
            xlast[1] = point[i][j].y;
            double* uplast = new double [2];
            uplast[0] = point[i][j].upstarx;
            uplast[1] = point[i][j].upstary;

            double ** Mtemp = new double * [2];
            Mtemp[0] = new double [2];
            Mtemp[1] = new double [2];
            Mtemp[0][0] = 1;
            Mtemp[0][1] = 0;
            Mtemp[1][0] = 0;
            Mtemp[1][1] = 1;
            
            double* xnext;
            xnext = RK_step1(xlast,Mtemp,uplast,dt,2);

            point[i][j].xtemp = point[i][j].x;
            point[i][j].ytemp = point[i][j].y;
            point[i][j].x = xnext[0];
            point[i][j].y = xnext[1];

            delete[] Mtemp[0];
            delete[] Mtemp[1];
            delete[] Mtemp;
            delete[] xlast;
            delete[] xnext;
            delete[] uplast;
        }
    }
    //此时temp保存n时刻位置，x,y保存第一次RK的位置

    //**********下面更新节点到相邻点的长度和外法向量****//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double ax, ay, bx, by;
            ax = point[i][j].x;
            ay = point[i][j].y;
            for (r=0; r<point[i][j].neighbor_node.size(); r++)
            {
                k = point[i][j].neighbor_node[r] / (m + 1);
                l = point[i][j].neighbor_node[r] % (m + 1);
                
                bx = point[k][l].x;
                by = point[k][l].y;

                point[i][j].a[r] = length(ax,ay,bx,by) * 0.5;

                double * ntemp;
                ntemp = normal(ax,ay,bx,by);
                point[i][j].nx[r] = ntemp[0];
                point[i][j].ny[r] = ntemp[1];
                delete[] ntemp;                
            }
        }
    }

    //***********更新单元顶点位置*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (r=0; r<4; r++)
            {
                k = o[i][j].vertex[r] / (m + 1);
                l = o[i][j].vertex[r] % (m + 1);

                o[i][j].vx[r] = point[k][l].x;
                o[i][j].vy[r] = point[k][l].y;
            }
        }
    }


    //**********下面进行第二次RK计算***********//
    //*********下面对每个节点求解节点速度*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            //注意此时last指针记录的是第一步RK的结果
            utemp = nodal_velocity(i,j);
            point[i][j].upstarx = utemp[0];
            point[i][j].upstary = utemp[1];
            delete[] utemp;
        }
    }

    //********下面进行第二次单元速度和质量总能的计算*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            double ** Rvs2;
            double * Res2;
            Rvs2 = velocity_matrix(i,j);
            Res2 = energy_matrix(i,j);

            double ** Mtemp = new double * [dim];
            for (k=0; k<dim; k++)
            {
                Mtemp[k] = new double [dim];
            }

            double * uxs2;
            double * uys2;
            double * taus2;
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uxs2 = RK_step2(o[i][j].uxnext,o[i][j].uxlast,Mtemp,Rvs2[0],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            uys2 = RK_step2(o[i][j].uynext,o[i][j].uylast,Mtemp,Rvs2[1],dt,dim);
            for (k=0; k<dim; k++)
            {
                for (l=0; l<dim; l++)
                {
                    Mtemp[k][l] = o[i][j].M[k][l];
                }
            }
            taus2 = RK_step2(o[i][j].taunext,o[i][j].taulast,Mtemp,Res2,dt,dim);

            for (k=0; k<dim; k++)
            {
                o[i][j].uxlast[k] = uxs2[k];
                o[i][j].uylast[k] = uys2[k];
                o[i][j].taulast[k] = taus2[k];
            }


            for (k=0; k<dim; k++)
            {
                delete[] Mtemp[k];
            }
            delete[] Rvs2[0];
            delete[] Rvs2[1];
            delete[] Rvs2;
            delete[] Res2;
            delete[] Mtemp;
            delete[] uxs2;
            delete[] uys2;
            delete[] taus2;

        }
    }
    //注意此时n+1时刻的值被保存在last中，next中无关紧要

    //**********下面进行第二次节点坐标的计算*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double* xlast = new double [2];
            xlast[0] = point[i][j].xtemp;
            xlast[1] = point[i][j].ytemp;
            double* xmid = new double [2];
            xmid[0] = point[i][j].x;
            xmid[1] = point[i][j].y;
            double * uplast = new double [2];
            uplast[0] = point[i][j].upstarx;
            uplast[1] = point[i][j].upstary;

            double ** Mtemp = new double * [2];
            Mtemp[0] = new double [2];
            Mtemp[1] = new double [2];
            Mtemp[0][0] = 1;
            Mtemp[0][1] = 0;
            Mtemp[1][0] = 0;
            Mtemp[1][1] = 1;
            
            double* xnext;
            xnext = RK_step2(xlast,xmid,Mtemp,uplast,dt,2);
            point[i][j].x = xnext[0];
            point[i][j].y = xnext[1];

            delete[] Mtemp[0];
            delete[] Mtemp[1];
            delete[] Mtemp;
            delete[] xlast;
            delete[] xmid;
            delete[] xnext;
            delete[] uplast;
        }
    }

    //**********下面更新节点到相邻点的长度和外法向量****//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double ax, ay, bx, by;
            ax = point[i][j].x;
            ay = point[i][j].y;
            for (r=0; r<point[i][j].neighbor_node.size(); r++)
            {
                k = point[i][j].neighbor_node[r] / (m + 1);
                l = point[i][j].neighbor_node[r] % (m + 1);
                
                bx = point[k][l].x;
                by = point[k][l].y;

                point[i][j].a[r] = length(ax,ay,bx,by) * 0.5;

                double * ntemp;
                ntemp = normal(ax,ay,bx,by);
                point[i][j].nx[r] = ntemp[0];
                point[i][j].ny[r] = ntemp[1];
                delete[] ntemp;                
            }
        }
    }

    //***********更新单元顶点位置*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (r=0; r<4; r++)
            {
                k = o[i][j].vertex[r] / (m + 1);
                l = o[i][j].vertex[r] % (m + 1);

                o[i][j].vx[r] = point[k][l].x;
                o[i][j].vy[r] = point[k][l].y;
            }
        }
    }

    return ;
}
