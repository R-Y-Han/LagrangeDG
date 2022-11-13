/**
 * @file limiter.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief 
 * @version 0.1
 * @date 2022-10-14
 * 
 * @copyright Copyright (c) 2022 R.Y. Han. All Reights Reserved.
 * 
 */

#include "config.h"
#include "limiter.h"
#include <cmath>

/**
 * @brief 对标量进行限制
 * 
 * @param i 单元行编号
 * @param j 单元列编号
 * @param up 限制前的值
 * @param alpha 限制的系数
 * @return double 
 */
double tau_limiter(int i, int j, double up, double alpha)
{
    double ul;
    double phi;
    double beta, ma, om, phiinva, temp;
    double um;
    int k, l, r, s;
        
    if (up - o[i][j].taulast[0] > 1e-6)
    {
        um = o[i][j].taulast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].taulast[0] > um)
            {
                um = o[k][l].taulast[0];
            }
        }

        phi = alpha * (um - o[i][j].taulast[0]) / (up - o[i][j].taulast[0]);
        phi = min(1.0,phi);
    }
    else if (up - o[i][j].taulast[0] < -1e-6)
    {
        um = o[i][j].taulast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].taulast[0] < um)
            {
                um = o[k][l].taulast[0];
            }
        }

        phi = alpha * (um - o[i][j].taulast[0]) / (up - o[i][j].taulast[0]);
        phi = min(1.0,phi);
    }
    else{
        phi = 1;
    }

    om = 0.2;
    ma = sqrt(o[i][j].mass / o[i][j].rholast);
    beta = min(1.0,om * ma);
    phiinva = 1;
    double divu, jt;
    if (1)
    {
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        divu = o[i][j].uxlast[1] * J_inv[0][0] + o[i][j].uxlast[2] * J_inv[1][0]
             + o[i][j].uylast[1] * J_inv[0][1] + o[i][j].uylast[2] * J_inv[1][1];
        
        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    for (s=0; s<o[i][j].neighbor_element.size(); s++)
    {
        k = o[i][j].neighbor_element[s] / m;
        l = o[i][j].neighbor_element[s] % m;

        double temp, jt;
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        temp = o[k][l].uxlast[1] * J_inv[0][0] + o[k][l].uxlast[2] * J_inv[1][0]
             + o[k][l].uylast[1] * J_inv[0][1] + o[k][l].uylast[2] * J_inv[1][1];

        temp = temp / (divu + 1e-6);

        if (temp < phiinva)
        {
            phiinva = temp;
        }

        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    
    phiinva = beta * max(0.0, phiinva) + (1-beta);
    phi = phiinva * phi;

    ul = o[i][j].taulast[0] + phi * (up - o[i][j].taulast[0]);
    return ul;
}

double ux_limiter(int i, int j, double up, double alpha)
{
    double ul;
    double phi;
    double beta, ma, om, phiinva, temp;
    double um;
    int k, l, r, s;
        
    if (up - o[i][j].uxlast[0]> 1e-6)
    {
        um = o[i][j].uxlast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uxlast[0] > um)
            {
                um = o[k][l].uxlast[0];
            }
        }

        phi = alpha * (um - o[i][j].uxlast[0]) / (up - o[i][j].uxlast[0]);
        phi = min(1.0,phi);
    }
    else if (up - o[i][j].uxlast[0]<-1e-6)
    {
        um = o[i][j].uxlast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uxlast[0] < um)
            {
                um = o[k][l].uxlast[0];
            }
        }

        phi = alpha * (um - o[i][j].uxlast[0]) / (up - o[i][j].uxlast[0]);
        phi = min(1.0,phi);
    }
    else{
        phi = 1;
    }

    om = 0.2;
    ma = sqrt(o[i][j].mass / o[i][j].rholast);
    beta = min(1.0,om * ma);
    phiinva = 1;
    double divu, jt;
    if (1)
    {
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        divu = o[i][j].uxlast[1] * J_inv[0][0] + o[i][j].uxlast[2] * J_inv[1][0]
             + o[i][j].uylast[1] * J_inv[0][1] + o[i][j].uylast[2] * J_inv[1][1];
        
        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    for (s=0; s<o[i][j].neighbor_element.size(); s++)
    {
        k = o[i][j].neighbor_element[s] / m;
        l = o[i][j].neighbor_element[s] % m;

        double temp, jt;
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        temp = o[k][l].uxlast[1] * J_inv[0][0] + o[k][l].uxlast[2] * J_inv[1][0]
             + o[k][l].uylast[1] * J_inv[0][1] + o[k][l].uylast[2] * J_inv[1][1];

        temp = temp / (divu + 1e-6);

        if (temp < phiinva)
        {
            phiinva = temp;
        }

        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    
    phiinva = beta * max(0.0, phiinva) + (1-beta);
    phi = phiinva * phi;

    ul = o[i][j].uxlast[0] + phi * (up - o[i][j].uxlast[0]);
    return ul;
}

double uy_limiter(int i, int j, double up, double alpha)
{
    double ul;
    double phi;
    double beta, ma, om, phiinva, temp;
    double um;
    int k, l, r, s;
        
    if (up - o[i][j].uylast[0] > 1e-6)
    {
        um = o[i][j].uylast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uylast[0] > um)
            {
                um = o[k][l].uylast[0];
            }
        }

        phi = alpha * (um - o[i][j].uylast[0]) / (up - o[i][j].uylast[0]);
        phi = min(1.0,phi);
    }
    else if (up - o[i][j].uylast[0] < -1e-6)
    {
        um = o[i][j].uylast[0];
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uylast[0] < um)
            {
                um = o[k][l].uylast[0];
            }
        }

        phi = alpha * (um - o[i][j].uylast[0]) / (up - o[i][j].uylast[0]);
        phi = min(1.0,phi);
    }
    else{
        phi = 1;
    }

    om = 0.2;
    ma = sqrt(o[i][j].mass / o[i][j].rholast);
    beta = min(1.0,om * ma);
    phiinva = 1;
    double divu, jt;
    if (1)
    {
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[i][j].xi_c,o[i][j].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[i][j].xi_c,o[i][j].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        divu = o[i][j].uxlast[1] * J_inv[0][0] + o[i][j].uxlast[2] * J_inv[1][0]
             + o[i][j].uylast[1] * J_inv[0][1] + o[i][j].uylast[2] * J_inv[1][1];
        
        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    for (s=0; s<o[i][j].neighbor_element.size(); s++)
    {
        k = o[i][j].neighbor_element[s] / m;
        l = o[i][j].neighbor_element[s] % m;

        double temp, jt;
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
            J_inv[0][0] = J_inv[0][0] + o[i][j].vx[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[0][1] = J_inv[0][1] + o[i][j].vx[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][0] = J_inv[1][0] + o[i][j].vy[r] * bp_xi(r,o[k][l].xi_c,o[k][l].eta_c);
            J_inv[1][1] = J_inv[1][1] + o[i][j].vy[r] * bp_eta(r,o[k][l].xi_c,o[k][l].eta_c);
        }
        J_inv[0][0] = J_inv[1][1] / jt;
        J_inv[0][1] = - J_inv[0][1] / jt;
        J_inv[1][0] = - J_inv[1][0] / jt;
        J_inv[1][1] = J_inv[0][0] / jt;
        temp = o[k][l].uxlast[1] * J_inv[0][0] + o[k][l].uxlast[2] * J_inv[1][0]
             + o[k][l].uylast[1] * J_inv[0][1] + o[k][l].uylast[2] * J_inv[1][1];

        temp = temp / (divu + 1e-6);

        if (temp < phiinva)
        {
            phiinva = temp;
        }

        delete[] J_inv[0];
        delete[] J_inv[1];
        delete[] J_inv;
    }
    
    phiinva = beta * max(0.0, phiinva) + (1-beta);
    phi = phiinva * phi;

    ul = o[i][j].uylast[0] + phi * (up - o[i][j].uylast[0]);
    return ul;
}