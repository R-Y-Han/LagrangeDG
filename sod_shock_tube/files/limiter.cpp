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
    double um;
    int k, l, r;
        
    if (up > o[i][j].taulast[0])
    {
        um = 0;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].taulast[0] > um)
            {
                um = o[k][l].taulast[0];
            }
        }

        ul = alpha * (um - o[i][j].taulast[0]) / (up - o[i][j].taulast[0]);
        ul = min(1.0,ul);
    }
    else if (up < o[i][j].taulast[0])
    {
        um = 10000;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].taulast[0] < um)
            {
                um = o[k][l].taulast[0];
            }
        }
        ul = alpha * (um - o[i][j].taulast[0]) / (up - o[i][j].taulast[0]);
        ul = min(1.0,ul);
    }
    else{
        ul = 1;
    }

    ul = o[i][j].taulast[0] + ul * (up - o[i][j].taulast[0]);
    return ul;
}

double ux_limiter(int i, int j, double up, double alpha)
{
    double ul;
    double um;
    int k, l, r;
        
    if (up > o[i][j].uxlast[0])
    {
        um = 0;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uxlast[0] > um)
            {
                um = o[k][l].uxlast[0];
            }
        }

        ul = alpha * (um - o[i][j].uxlast[0]) / (up - o[i][j].uxlast[0]);
        ul = min(1.0,ul);
    }
    else if (up < o[i][j].uxlast[0])
    {
        um = 10000;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uxlast[0] < um)
            {
                um = o[k][l].uxlast[0];
            }
        }
        ul = alpha * (um - o[i][j].uxlast[0]) / (up - o[i][j].uxlast[0]);
        ul = min(1.0,ul);
    }
    else{
        ul = 1;
    }

    ul = o[i][j].uxlast[0] + ul * (up - o[i][j].uxlast[0]);
    return ul;
}

double uy_limiter(int i, int j, double up, double alpha)
{
    double ul;
    double um;
    int k, l, r;
        
    if (up > o[i][j].uylast[0])
    {
        um = 0;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uylast[0] > um)
            {
                um = o[k][l].uylast[0];
            }
        }

        ul = alpha * (um - o[i][j].uylast[0]) / (up - o[i][j].uylast[0]);
        ul = min(1.0,ul);
    }
    else if (up < o[i][j].uylast[0])
    {
        um = 10000;
        for (r=0; r<o[i][j].neighbor_element.size(); r++)
        {
            k = o[i][j].neighbor_element[r] / m;
            l = o[i][j].neighbor_element[r] % m;

            if (o[k][l].uylast[0] < um)
            {
                um = o[k][l].uylast[0];
            }
        }
        ul = alpha * (um - o[i][j].uylast[0]) / (up - o[i][j].uylast[0]);
        ul = min(1.0,ul);
    }
    else{
        ul = 1;
    }

    ul = o[i][j].uylast[0] + ul * (up - o[i][j].uylast[0]);
    return ul;
}