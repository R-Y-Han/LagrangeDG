/**
 * @file initial.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief �������񲢼�¼���г�ʼ��Ϣ
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
 * @brief ���ɳ�ʼ���񻮷�
 * - �����ڴ�
 * - ��¼����ı��
 * - ��¼�ڵ�ı��
 * - ��¼�ڵ���������
 * - ��¼���񶥵��š�������������
 */
void generatemesh()
{
    int i, j;

    //�����¼�ڵ��ţ��ж��Ƿ�Ϊ�߽�ڵ㣬��¼��������
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            //��¼���
            point[i][j].i = i;
            point[i][j].j = j;
            point[i][j].q = (m+1) * i + j;
            
            //�ж��Ƿ�Ϊ�߽�ڵ�
            point[i][j].boundary = 0;
            if ( i == 0 || i == n)
            {
                point[i][j].boundary = 1;
            }
            if (j == 0 || j == m)
            {
                point[i][j].boundary = 1;
            }

            //��¼��������
            point[i][j].x = hx * i;
            point[i][j].y = hy * j;
            
            point[i][j].x = point[i][j].x + (0.1 - point[i][j].y) * sin(pi * point[i][j].x);
            
            point[i][j].x0 = point[i][j].x;
            point[i][j].y0 = point[i][j].y;
        }
    }

    //�����¼������
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            //��¼������
            o[i][j].i = i;
            o[i][j].j = j;
            o[i][j].q = m * i + j;
        }
    }

    //�����¼���񶥵�
   
    int *ti = new int [4];  //��¼������ʱ����б��
    int *tj = new int [4];  //��¼������ʱ����б��
    //��һ������Ϊ���½ǵĵ㣬��ʱ����ת

    for (i = 0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            int tempi, tempj;  //�ֱ��ʾ������б�ź��б��

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
 * @brief �Ե�(i,j)����ԪѰ�����ڵ�Ԫ��
 * ��(i,j-1)��ʼ��ʱ������һ�����ڵĵ�Ԫ��Ϊ��һ����
 * �洢һά���
 * 
 * @param i ����Ķ�ά�б��
 * @param j ����Ķ�ά�б��
 */
void element_findneighbor(int i, int j)
{
    o[i][j].neighbor_element.clear();

    int ti;  //��������������б��
    int tj;  //��������������б��
    int num;    //�������������һά���
    
    //�ж�(i,j-1)��Ԫ�Ƿ����
    if ( j == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j-1;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i+1,j)��Ԫ�Ƿ����
    if ( i == n - 1 ) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i+1;
        tj = j;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i,j+1)��Ԫ�Ƿ����
    if ( j == m-1) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j+1;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    //�ж�(i-1,j)��Ԫ�Ƿ����
    if ( i == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i-1;
        tj = j;
        num = o[ti][tj].q;
        o[i][j].neighbor_element.push_back(num);
    }

    return ;
}

/**
 * @brief 
 * - �Ե�(i,j)���ڵ�Ѱ�����ڵ�����
 * ��(i,j)��������ʱ������һ�����ڵ�����ʼ��¼����ʱ��˳��
 * - �Ե�(i,j)���ڵ�Ѱ�����ڵĽڵ�
 * ��(i-1,j)��ʼ��ʱ�������һ�����ڵĽڵ㿪ʼ��¼����ʱ��˳��
 * - ������ߵģ����֣��߳����ⷨ����
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 */
void node_findneighbor(int i, int j)
{
    int ti, tj; //�������������ڵ���б�ź��б��
    int num;    //�������������ڵ��һά���
    
    //**************�ȼ�¼��������****************//
    point[i][j].neighbor_element.clear();
    
    //�ж�����(i,j)�Ƿ����
    if ( i == n)
    {
        //������
        ;
    }
    else{
        ti = i;
        if (j == m)
        {
            //�ڵ�Ϊ���Ͻ�
            ;
        }
        else{
            tj = j;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //�ж�����(i-1,j)�Ƿ����
    if ( i == 0)
    {
        //������
        ;
    }
    else{
        ti = i-1;
        if (j == m)
        {
            //������
            ;
        }
        else{
            tj = j;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //�ж�����(i-1,j-1)�Ƿ����
    if ( i == 0)
    {
        //������
        ;
    }
    else{
        ti = i-1;
        if (j == 0)
        {
            //������
            ;
        }
        else{
            tj = j-1;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }

    //�ж�����(i,j-1)�Ƿ����
    if ( i == n)
    {
        //������
        ;
    }
    else{
        ti = i;
        if (j == 0)
        {
            //������
            ;
        }
        else{
            tj = j-1;
            num = o[ti][tj].q;
            point[i][j].neighbor_element.push_back(num);
        }
    }
    //*****************���������¼���******************//

    //*****************�����¼���ڽڵ�******************//
    point[i][j].neighbor_node.clear();

    //�ж�(i,j-1)�ڵ��Ƿ����
    if ( j == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j-1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i+1,j)�ڵ��Ƿ����
    if ( i == n) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i+1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i,j+1)�ڵ��Ƿ����
    if ( j == m) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i;
        tj = j+1;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }

    //�ж�(i-1,j)�ڵ��Ƿ����
    if ( i == 0) 
    {
        //������
        ;
    }
    else{
        //����
        ti = i-1;
        tj = j;
        num = point[ti][tj].q;
        point[i][j].neighbor_node.push_back(num);
    }
    //**********���ڽڵ��¼���******************//

    //**********������㵽���ڽڵ�ĳ��Ⱥ��ⷨ��**********//
    point[i][j].a.clear();
    point[i][j].nx.clear();
    point[i][j].ny.clear();
    for (int k=0; k < point[i][j].neighbor_node.size(); k++)
    {
        int itemp, jtemp;   //���ڽڵ�����б��

        itemp = point[i][j].neighbor_node[k] / (m+1);
        jtemp = point[i][j].neighbor_node[k] % (m+1);

        double ax, ay, bx, by;  //a�����ڵ㣬b�������ڽڵ�

        ax = point[i][j].x;
        ay = point[i][j].y;
        bx = point[itemp][jtemp].x;
        by = point[itemp][jtemp].y;

        //���㵽�ڵ�Ĳ��ֳ���
        double l;  //����
        l = length(ax,ay,bx,by);

        l = l / 2.0; //����2����Ϊ������ÿ���ڵ㱻���䵽һ��ı߳�
        point[i][j].a.push_back(l);

        //���㷨��
        double * ntemp;    //����
        ntemp = normal(ax,ay,bx,by);
        point[i][j].nx.push_back(ntemp[0]);
        point[i][j].ny.push_back(ntemp[1]);

        delete[] ntemp;
    }

    return ;
}

/**
 * @brief �����(i,j)����������ģ�ͬʱ����ÿ�����������
 * 
 * @param i �����б��
 * @param j �����б��
 * @return double* ���ĵĲο���Ԫ����
 * @note ���ص���ָ�룬�洢ʱ�Ǳ������ǵ�ɾ��ָ������ڴ�
 */
double * mass_center(int i, int j)
{
    double * ce = new double [2];

    int k;

    double masst;    //��Ԫ����
    double mf_xi, mf_eta;   //������

    //********��������Gauss���ּ���������������************//
    double rhot, jt;    //��¼�ܶȺ�Jacobi����ʽ�ĵ�ֵ
    double xt, yt;  //Gauss�����������
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

    //**********���������������***********//

    ce[0] = mf_xi / masst;  //���Ĳο��ռ������
    ce[1] = mf_eta / masst; //���Ĳο��ռ�������

    return ce;
}

/**
 * @brief �����(i,j)���������������
 * 
 * @param i �����б��
 * @param j �����б��
 * @return double** ��������
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
 * @brief ��ʼ�������������������񲢼�¼��ʼ��Ϣ
 * 
 */
void initial()
{
    int i, j, k;
    //***********�����������񻮷�**********//
    generatemesh();

    //***********�����ÿ������Ѱ����������**********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            element_findneighbor(i,j);
        }
    }

    //***********�����ÿ���ڵ�Ѱ�����ڽڵ�***********//
    for (i=0; i<=n; i++)
    {
        for (j=0; j<=m; j++)
        {
            node_findneighbor(i,j);
        }
    }


    //************�������ÿ����������������ģ���������������*********//
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            //�����������
            double * ce;
            ce = mass_center(i,j);
            o[i][j].xi_c = ce[0];
            o[i][j].eta_c = ce[1];

            delete[] ce;

            //���������������
            o[i][j].M = mass_matrix(i,j);
        }
    }

    //***********�����ÿ����Ԫ����ֵ***************//
    double xt, yt; //���ĵ���������
    double x_xi, x_eta, y_xi, y_eta;    //Jacobi������ĸ�����
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

            //********�����ʼ���ٶȳ�************//
            o[i][j].uxlast[0] = ini_ux(xt,yt);
            o[i][j].uylast[0] = ini_uy(xt,yt);  //�����ĵ�ֵ

            o[i][j].uxlast[1] = ini_ux_x(xt,yt) * x_xi + ini_ux_y(xt,yt) * y_xi;    //���Ĵ�ux��xi��ƫ��
            o[i][j].uylast[1] = ini_uy_x(xt,yt) * x_xi + ini_uy_y(xt,yt) * y_xi;    //���Ĵ�uy��xi��ƫ��

            o[i][j].uxlast[2] = ini_ux_x(xt,yt) * x_eta + ini_ux_y(xt,yt) * y_eta;  //���Ĵ�ux��eta��ƫ��
            o[i][j].uylast[2] = ini_uy_x(xt,yt) * x_eta + ini_uy_y(xt,yt) * y_eta;  //���Ĵ�uy��eta��ƫ��
            
            //*********�����ʼ����������*********//
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