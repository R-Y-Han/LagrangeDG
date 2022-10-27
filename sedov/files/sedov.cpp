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
* @file [sodShock.cpp]
* @brief Implementation for 2D Taylor-Green vortex test problem.
* @author R.Y. Han,  IAPCM
* @email: hanruoyu21@gscaep.ac.cn
* @date 10.12.2022
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
 * @brief ����һ��ʱ���ƽ�
 * 
 */
void onetimestep();

void onetimestep2();

int main()
{
    int i, j, k;
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

        onetimestep();
        cout<<time<<endl;
    }



    plotmesh();
    plotpressure();
    plotinternalenergy();
    system("pause");
}

void onetimestep()
{
    int i, j, k, l, r;
    
    //********������е�һ��RK����**********//
    //*********�����ÿ���ڵ����ڵ��ٶ�*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            utemp = nodal_velocity(i,j);
            delete[] utemp;
        }
    }

    //**********������е�һ�νڵ�����ļ���*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            point[i][j].xtemp = point[i][j].x;
            point[i][j].ytemp = point[i][j].y;
            
            point[i][j].x = point[i][j].xtemp + dt * point[i][j].upstarx;
            point[i][j].y = point[i][j].ytemp + dt * point[i][j].upstary;
        }
    }
    //��ʱtemp����nʱ��λ�ã�x,y�����һ��RK��λ��

    //**********������½ڵ㵽���ڵ�ĳ��Ⱥ��ⷨ����****//
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

    //***********���µ�Ԫ����λ��*********//
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

    //********������е�һ�ε�Ԫ�ٶȺ��������ܵļ���*********//
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
    //ע���ʱtnʱ�̵�ֵ������next�У���һ��RK��ֵ������last��


    //**********������еڶ���RK����***********//
    //*********�����ÿ���ڵ����ڵ��ٶ�*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            double * utemp;
            //ע���ʱlastָ���¼���ǵ�һ��RK�Ľ��
            utemp = nodal_velocity(i,j);
            point[i][j].upstarx = utemp[0];
            point[i][j].upstary = utemp[1];
            delete[] utemp;
        }
    }

     //**********������еڶ��νڵ�����ļ���*********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            point[i][j].x = 0.5 * point[i][j].xtemp + 0.5 * point[i][j].x + 0.5 * dt * point[i][j].upstarx;
            point[i][j].y = 0.5 * point[i][j].ytemp + 0.5 * point[i][j].y + 0.5 * dt * point[i][j].upstary;
        }
    }

    //**********������½ڵ㵽���ڵ�ĳ��Ⱥ��ⷨ����****//
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

    //***********���µ�Ԫ����λ��*********//
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



    //********������еڶ��ε�Ԫ�ٶȺ��������ܵļ���*********//
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
    //ע���ʱn+1ʱ�̵�ֵ��������last�У�next���޹ؽ�Ҫ

    return ;
}
