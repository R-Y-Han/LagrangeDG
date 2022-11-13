/**
 * @file plot.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 绘制plt文件，生成网格和密度、速度、内能的分布图
 * @version v1.01
 * @date 2022-10-09
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Rights Reseved.
 * 
 */

#include "plot.h"
#include "config.h"
#include <fstream>
using namespace std;


const char* fn_mesh = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\mesh.plt";
const char* fn_pre = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\pressure.plt";
const char* fn_interene = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\e.plt";
const char* fn_ux = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\ux.plt";

/**
 * @brief 绘制网格图像
 * 
 */
void plotmesh()
{
    int i, j, k;
    remove(fn_mesh);
    ofstream f(fn_mesh);
    f<<"VARIABLES = X, Y"<<endl;
    f<<"ZONE N = "<<(n + 1) * (m + 1)<<", E = "<<n*m<<", F = FEPOINT, ET = QUADRILATERAL"<<endl;
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            f<<"\t"<<point[i][j].x<<"\t"<<point[i][j].y<<endl;
        }
    }

    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (k=0; k<4; k++)
            {
                f<<"\t"<<o[i][j].vertex[k]+1;
            }
            f<<endl;
        }
    }

    f.close();
    return ;
}

/**
 * @brief 绘制网格压力
 * 
 */
void plotpressure()
{
    int i, j, k;
    double xt, yt, temp;
    remove(fn_pre);
    ofstream f(fn_pre);
    f<<"VARIABLES = X, P"<<endl;

    for (i=0; i<n; i++)
    {
        f<<point[i][0].x<<"\t";
            double rho, p;
            rho = o[i][0].rholast;
            temp = o[i][0].uxlast[0] * o[i][0].uxlast[0]
                 + o[i][0].uylast[0] * o[i][0].uylast[0];
            temp = o[i][0].taulast[0] - 0.5 * temp;
            p = EOS(rho,temp);
            f<<p<<endl;
    }

    f.close();
    return ;
}

/**
 * @brief 绘制网格内能
 * 
 */
void plotinternalenergy()
{
    int i, j, k;
    double temp;
    remove(fn_interene);
    ofstream f(fn_interene);
    f<<"VARIABLES = X, Y, e"<<endl;
    f<<"ZONE N = "<<(n+1)*(m+1)<<", E = "<<n*m<<endl;
    f<<"ZONETYPE = FEQUADRILATERAL"<<endl;
    f<<"DATAPACKING = BLOCK"<<endl;
    f<<"varlocation=([3]=cellcentered)"<<endl;
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            f<<point[i][j].x<<endl;
        }
    }
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            f<<point[i][j].y<<endl;
        }
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            temp = o[i][j].uxlast[0] * o[i][j].uxlast[0];
            temp = temp + o[i][j].uylast[0] * o[i][j].uylast[0];
            temp = o[i][j].taulast[0] - 0.5 * temp;
            f<<temp<<endl;
        }
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (k=0; k<4; k++)
            {
                f<<"\t"<<o[i][j].vertex[k]+1;
            }
            f<<endl;
        }
    }

    f.close();
    return ;
}

/**
 * @brief 绘制网格压力
 * 
 */
void plotux()
{
    int i, j, k;
    double xt, yt, temp;
    remove(fn_ux);
    ofstream f(fn_ux);
    f<<"VARIABLES = X, Y, ux"<<endl;
    f<<"ZONE N = "<<(n+1)*(m+1)<<", E = "<<n*m<<endl;
    f<<"ZONETYPE = FEQUADRILATERAL"<<endl;
    f<<"DATAPACKING = BLOCK"<<endl;
    f<<"varlocation=([3]=cellcentered)"<<endl;
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            f<<point[i][j].x<<endl;
        }
    }
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            f<<point[i][j].y<<endl;
        }
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            temp = o[i][j].uxlast[0];
            f<<temp<<endl;
        }
    }
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            for (k=0; k<4; k++)
            {
                f<<"\t"<<o[i][j].vertex[k]+1;
            }
            f<<endl;
        }
    }

    f.close();
    return ;
}

void plote1d()
{
    int i, j, k;
    double xt, yt, temp;
    const char* fn_e1d = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\e1d.plt";
    remove(fn_e1d);
    ofstream f(fn_e1d);
    f<<"VARIABLES = X, e"<<endl;

    for (i=0; i<n; i++)
    {
        f<<point[i][0].x<<"\t";
            temp = o[i][0].uxlast[0] * o[i][0].uxlast[0]
                 + o[i][0].uylast[0] * o[i][0].uylast[0];
            temp = o[i][0].taulast[0] - 0.5 * temp;
            f<<temp<<endl;
    }

    f.close();
    return ;
}


void plotux1d()
{
    int i, j, k;
    double xt, yt, temp;
    const char* fn_ux1d = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\ux1d.plt";
    remove(fn_ux1d);
    ofstream f(fn_ux1d);
    f<<"VARIABLES = X, numerical_velocity"<<endl;

    for (i=0; i<n; i++)
    {
        f<<point[i][0].x<<"\t";
            temp = o[i][m/2].uxlast[0];
            f<<temp<<endl;
    }

    f.close();
    return ;
}

void plotrho1d()
{
    int i, j, k;
    double xt, yt, temp;
    const char* fn_rho1d = "F:\\C++Code\\LagrangeDG\\sod_shock_tube\\output\\rho1d.plt";
    remove(fn_rho1d);
    ofstream f(fn_rho1d);
    f<<"VARIABLES = X, density"<<endl;

    for (i=0; i<n; i++)
    {
        f<<point[i][0].x<<"\t";
            double rho;
            rho = o[i][m/2].rholast;
            f<<rho<<endl;
    }

    f.close();
    return ;
}