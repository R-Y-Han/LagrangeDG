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


const char* fn_mesh = "F:\\C++Code\\LagrangeDG\\sedov\\output\\mesh.plt";
const char* fn_pre = "F:\\C++Code\\LagrangeDG\\sedov\\output\\pressure.plt";
const char* fn_interene = "F:\\C++Code\\LagrangeDG\\sedov\\output\\e.plt";

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
    double xt0, yt0, temp;
    remove(fn_pre);
    ofstream f(fn_pre);
    f<<"VARIABLES = X, Y, P"<<endl;
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
            
            xt0 = o[i][j].phi_x0(o[i][j].xi_c,o[i][j].eta_c);
            yt0 = o[i][j].phi_y0(o[i][j].xi_c,o[i][j].eta_c);
            double rho, p;
            rho = ini_rho(xt0,yt0) * o[i][j].Jacobi_0(o[i][j].xi_c,o[i][j].eta_c)
                                   / o[i][j].Jacobi(o[i][j].xi_c,o[i][j].eta_c);
            temp = o[i][j].uxlast[0] * o[i][j].uxlast[0]
                 + o[i][j].uylast[0] * o[i][j].uylast[0];
            temp = o[i][j].taulast[0] - 0.5 * temp;
            p = EOS(rho,temp);
            f<<p<<endl;
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

void plotrho1d()
{
    int i, j, k;
    double xt, yt, temp;
    const char* fn_e1d = "F:\\C++Code\\LagrangeDG\\sedov\\output\\rho1d.plt";
    remove(fn_e1d);
    ofstream f(fn_e1d);
    f<<"VARIABLES = X, density"<<endl;

    for (i=n/2; i<n; i++)
    {
        f<<point[i][0].x<<"\t";
        double rho;
        double xc, yc;
        xc = o[i][m/2].phi_x0(o[i][m/2].xi_c,o[i][m/2].eta_c);
        yc = o[i][m/2].phi_y0(o[i][m/2].xi_c,o[i][m/2].eta_c);
        rho = ini_rho(xc,yc) * o[i][m/2].Jacobi_0(o[i][m/2].xi_c,o[i][m/2].eta_c)
                             / o[i][m/2].Jacobi(o[i][m/2].xi_c,o[i][m/2].eta_c);
        f<<rho<<endl;
    }

    f.close();
    return ;
}