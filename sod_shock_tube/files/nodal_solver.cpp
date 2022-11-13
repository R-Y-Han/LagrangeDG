/**
 * @file nodal_solver.cpp
 * @author R.Y. Han (hanruoyu21@gscaep.ac.cn)
 * @brief ����ڵ�ⷨ�����ڵ��ٶȺ�corner force
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
 * @brief �����(i,j)���ڵ㴦ÿ�����ڵ�Ԫ�ڸõ���ٶ��ع�ֵ��ƽ��
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
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
        //************�Ե�k�����ڵ�Ԫִ��*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;    //�ڵ��r�����ڵ�Ԫ�����б��
        
            int qt;
        //********�����жϽڵ�(i,j)�ڵ�Ԫ(k,l)���ǵڼ�������*******//
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
                //��ʱѭ��ָ��temp�պþ��ǽڵ�(i,j)��(k,l)��Ԫ�еľֲ�λ��
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];
        //**********�������õ�Ԫ�ڽڵ㴦�ٶ��ع�ֵ*******//
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
    //********��ʱuavg��¼�����е�Ԫ�ڸõ㴦�ٶ��ع�ֵ���ܺ�********//
    uavg[0] = uavg[0] / (double) point[i][j].neighbor_element.size();
    uavg[1] = uavg[1] / (double) point[i][j].neighbor_element.size();
    //********uavgΪ�ڵ㴦�ٶ��ع�ֵ��ƽ��****************//

    return uavg;
}

/**
 * @brief ����߽�ڵ���ٶ�
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 * @return double* 
 */
double * nodal_velocity_boundary(int i, int j)
{
    double *ustar = new double [2];
    
    //********������ĸ�����*******//
    if (i == 0 || i == n)
    {
        ustar[0] = 0;
        ustar[1] = 0;
        return ustar;
    }

    //********���ϵĵ�*********//
    int k, l, r, temp;
    int node1, node2, bnode1, bnode2, p1, p2;
    double nstarx, nstary;
    double xit, etat, nx1, ny1, nx2, ny2, a1, a2;
    double mu, topx, topy, bottom;
    //********�������ҳ��߽���ⷨ��*******//
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
        //***********ȷ��corner���ڵ�Ԫ���*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;

        //************ȷ���ڵ��ڵ�Ԫ�ڵ�λ�úͲο��ռ�����***********//
        for (temp=0; temp < 4; temp++)
        {
            if (o[k][l].vertex[temp] == point[i][j].q)
            {
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];

        //***********ȷ�����ڽڵ�ı��***********//
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
        //node1, node2Ϊ���ڶ���ı��

        //*********ȷ�����ڽڵ�͸ýڵ�Ĺ�ϵ*********//
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
        //p1, p2Ϊ(i,j)ָ�����ڶ����ָ��

        nx1 = - point[i][j].nx[p1];
        ny1 = - point[i][j].ny[p1]; //ע��˴����˸��ţ�Ϊʹָ������

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


    double b_pre =0;
    topx = 0;
    topy = 0;
    bottom = 0;
    //***************��ÿ��corner���м���***********//
    for (r=0; r<point[i][j].neighbor_element.size(); r++)
    {
        //***********ȷ��corner���ڵ�Ԫ���*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;

        //************ȷ���ڵ��ڵ�Ԫ�ڵ�λ�úͲο��ռ�����***********//
        for (temp=0; temp < 4; temp++)
        {
            if (o[k][l].vertex[temp] == point[i][j].q)
            {
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];

        //***********ȷ�����ڽڵ�ı��***********//
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
        //node1, node2Ϊ���ڶ���ı��

        //*********ȷ�����ڽڵ�͸ýڵ�Ĺ�ϵ*********//
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
        //p1, p2Ϊ(i,j)ָ�����ڶ����ָ��

        //***********��¼��λ�ⷨ�����Ӧ����********//
        a1 = point[i][j].a[p1];
        a2 = point[i][j].a[p2];

        nx1 = - point[i][j].nx[p1];
        ny1 = - point[i][j].ny[p1]; //ע��˴����˸��ţ�Ϊʹָ������

        nx2 = point[i][j].nx[p2];
        ny2 = point[i][j].ny[p2];

        //*******��������r��corner�е���Ӧ����********//
        double xt_v, yt_v;
        xt_v = o[k][l].phi_x(xit,etat);
        yt_v = o[k][l].phi_y(xit,etat);

        //********��������ܶ�**********//
        double rho;
        rho = o[k][l].rholast;

        //********�������ڵ㴦�ٶȺ������ع�ֵ********//
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
        //uc[0] = ux_limiter(k,l,uc[0],0.6);
        //uc[1] = uy_limiter(k,l,uc[1],0.6);
        //ve = tau_limiter(k,l,ve,0.6);
        //********�����������**********//
        ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        ve = max(1e-9,ve);
        //**************��������迹*************//
        mu = (gamma - 1) * ve;
        mu = sqrt(mu);  //��ʱmuΪ����
        mu = rho * mu;
        
        //*************������㼤����λ����*********//
        double * e;
        e = u_average(i,j);

        e[0] = e[0] - uc[0];
        e[1] = e[1] - uc[1];
        
        double norm;
        norm = e[0] * e[0] + e[1] * e[1];
        norm = sqrt(norm);
        if (norm < 1e-9)
        {
            e[0] = 1e-9;
            e[1] = 1e-9;
        }
        else{
            e[0] = e[0] / norm;
            e[1] = e[1] / norm;
        }
        //*********���������ӷ�ĸ������********//
        double p_vertex;
        p_vertex = EOS (rho, ve);
        
        double mid;
        //*********��������һ���ߵ���***********//
        mid = nx1 * e[0] + ny1 * e[1];
        mid = mu * abs(mid) * a1;

        
        topx = topx + mid * uc[0] + a1 * nx1 * p_vertex;
        topy = topy + mid * uc[1] + a1 * ny1 * p_vertex;
        bottom = bottom + mid;

        //************�������ڶ����ߵ���*******//
        mid = nx2 * e[0] + ny2 * e[1];
        mid = mu * abs(mid) * a2;

        
        topx = topx + mid * uc[0] + a2 * nx2 * p_vertex;
        topy = topy + mid * uc[1] + a2 * ny2 * p_vertex;
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
        b_pre = topx * nstarx + topy * nstary;
        ustar[0] = (topx - b_pre * nstarx) / bottom;
        ustar[1] = (topy - b_pre * nstary) / bottom;
    }
    return ustar;
}


/**
 * @brief ����ڵ��ٶ�
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
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

    //**************���ж��Ƿ�Ϊ�߽��*********//
    if (i == 0 || i == n)
    {
        //******�Ǳ߽��*******//
        //ustar = nodal_velocity_boundary(i,j);
        ustar = new double [2];
        ustar[0] = 0;
        ustar[1] = 0;
        return ustar;
    }

    //*************���Ǳ߽�ڵ�***************//
    ustar = new double [2];
    int r, temp, k, l, node1, node2, p1, p2;
    double xit, etat, nx1, ny1, nx2, ny2, a1, a2;
    double mu, mid;
    double rhs[2];
    rhs[0] = 0;
    rhs[1] = 0;
    //***************��ÿ��corner���м���***********//
    for (r=0; r<point[i][j].neighbor_element.size(); r++)
    {
        //***********ȷ��corner���ڵ�Ԫ���*********//
        k = point[i][j].neighbor_element[r] / m;
        l = point[i][j].neighbor_element[r] % m;

        //************ȷ���ڵ��ڵ�Ԫ�ڵ�λ�úͲο��ռ�����***********//
        int qt;
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
        for (temp=0; temp < 4; temp++)
        {
            if (o[k][l].vertex[temp] == qt)
            {
                break;
            }
        }
        xit = ref_xi[temp][0];
        etat = ref_xi[temp][1];

        //***********ȷ�����ڽڵ�ı��***********//
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
        //node1, node2Ϊ���ڶ���ı��

        //*********ȷ�����ڽڵ�͸ýڵ�Ĺ�ϵ*********//
        for (p1 = 0; p1 < point[i][j].neighbor_node.size(); p1++)
        {
            if (j==0)
            {
                if (l == m-1)
                {
                    qt = point[i][m].neighbor_node[p1];
                }
                else{
                    qt = point[i][j].neighbor_node[p1];
                }
            }
            else if (j == m)
            {
                if (l == 0)
                {
                    qt = point[i][0].neighbor_node[p1];
                }
                else{
                    qt = point[i][j].neighbor_node[p1];
                }
            }
            else{
                qt = point[i][j].neighbor_node[p1];
            }
            if (qt == node1)
            {
                break;
            }
        }
        

        for (p2 = 0; p2 < point[i][j].neighbor_node.size(); p2++)
        {
            if (j==0)
            {
                if (l == m-1)
                {
                    qt = point[i][m].neighbor_node[p2];
                }
                else{
                    qt = point[i][j].neighbor_node[p2];
                }
            }
            else if (j == m)
            {
                if (l == 0)
                {
                    qt = point[i][0].neighbor_node[p2];
                }
                else{
                    qt = point[i][j].neighbor_node[p2];
                }
            }
            else{
                qt = point[i][j].neighbor_node[p2];
            }
            if (qt == node2)
            {
                break;
            }
        }
        //p1, p2Ϊ(i,j)ָ�����ڶ����ָ��

        //***********��¼��λ�ⷨ�����Ӧ����********//
        a1 = point[i][j].a[p1];
        a2 = point[i][j].a[p2];

        nx1 = - point[i][j].nx[p1];
        ny1 = - point[i][j].ny[p1]; //ע��˴����˸��ţ�Ϊʹָ������

        nx2 = point[i][j].nx[p2];
        ny2 = point[i][j].ny[p2];

        //*******��������r��corner�е���Ӧ����********//
        double xt_v, yt_v;
        xt_v = o[k][l].phi_x(xit,etat);
        yt_v = o[k][l].phi_y(xit,etat);

        //********��������ܶ�**********//
        double rho_vertex;
        //rho_vertex = o[k][l].rholast;
        double xt0, yt0;
        xt0 = o[k][l].phi_x0(xit,etat);
        yt0 = o[k][l].phi_y0(xit,etat);
        if (xt0 == 0.5)
        {
            double xc;
            xc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
            xt0 = xt0 - 1e-9 * (xt0 - xc);
        }
        rho_vertex = ini_rho(xt0,yt0) * o[k][l].Jacobi_0(xit,etat)
                    / o[k][l].Jacobi(xit,etat);

        //********�������ڵ㴦�ٶȺ������ع�ֵ********//
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
        

        //********�����������**********//
        ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
        ve = max(1e-9, ve);
        double p_vertex;
        p_vertex = EOS(rho_vertex,ve);

        //**************��������迹*************//
        //mu = (gamma - 1) * ve;
        //mu = sqrt(mu);  //��ʱmuΪ����
        //mu = rho * mu;
        double rho_center;
        rho_center = o[k][l].rholast;

        double p_center, ve_center;
        ve_center = o[k][l].uxlast[0]*o[k][l].uxlast[0] + o[k][l].uylast[0] * o[k][l].uylast[0];
        ve_center = o[k][l].taulast[0] - ve_center * 0.5;
        p_center = EOS(rho_center,ve_center);
        p_center = max(1e-12, p_center);
        mu = gamma * p_center / rho_center;
        mu = sqrt(mu);  //��ʱmuΪ����
        mu = rho_center * mu;
        
        //**********��������һ���߷���*******//
        mid = nx1 * uc[0] + ny1 * uc[1];
        rhs[0] = rhs[0] + a1 * p_vertex * nx1 + mu * a1 * mid * nx1;
        rhs[1] = rhs[1] + a1 * p_vertex * ny1 + mu * a1 * mid * ny1;

        Mpn[0][0] = Mpn[0][0] + mu * a1 * nx1 * nx1;
        Mpn[0][1] = Mpn[0][1] + mu * a1 * ny1 * nx1;
        Mpn[1][0] = Mpn[1][0] + mu * a1 * nx1 * ny1;
        Mpn[1][1] = Mpn[1][1] + mu * a1 * ny1 * ny1;


        //*********�������ڶ����߷���********//
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
 * @brief �����(i,j)���ڵ��ڵ�k�����ڵ�Ԫ�е�����corner force
 * 
 * @param i �ڵ���б��
 * @param j �ڵ���б��
 * @param r �����ǵĵ�Ԫ���
 * @return double** 
 * 
 * @note 
 * - ע���ڵ�k����Ԫ��������corner force
 * - ע�ⷨ�����ķ����Ƿ�ָ��Ԫ�ⲿ
 */
double ** corner_force(int i, int j, int r)
{
    int k, l, temp, node1, node2, p1, p2;
    double xit, etat, a1, a2, nx1, ny1, nx2, ny2, mu;

    double ** f = new double * [2]; //��ÿ��������������corner force
    for (temp=0; temp < 2; temp++)
    {
        f[temp] = new double [2];
    }

    //***********ȷ��corner���ڵ�Ԫ���*********//
    k = r / m;
    l = r % m;

    int qt;
    //************ȷ���ڵ��ڵ�Ԫ�ڵ�λ�úͲο��ռ�����***********//
    for (temp=0; temp < 4; temp++)
    {
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
        if (o[k][l].vertex[temp] == qt)
        {
            break;
        }
    }
    xit = ref_xi[temp][0];
    etat = ref_xi[temp][1];

    //***********ȷ�����ڽڵ�ı��***********//
    //��ʱtemp��¼�ڵ��ڵ�Ԫ�е�λ��
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
    //node1, node2Ϊ���ڶ���ı��

    //*********ȷ�����ڽڵ�͸ýڵ�Ĺ�ϵ*********//
    for (p1 = 0; p1 < point[i][j].neighbor_node.size(); p1++)
    {
            if (j==0)
            {
                if (l == m-1)
                {
                    qt = point[i][m].neighbor_node[p1];
                }
                else{
                    qt = point[i][j].neighbor_node[p1];
                }
            }
            else if (j == m)
            {
                if (l == 0)
                {
                    qt = point[i][0].neighbor_node[p1];
                }
                else{
                    qt = point[i][j].neighbor_node[p1];
                }
            }
            else{
                qt = point[i][j].neighbor_node[p1];
            }
        if (qt == node1)
        {
            break;
        }
    }

    for (p2 = 0; p2 < point[i][j].neighbor_node.size(); p2++)
    {
            if (j==0)
            {
                if (l == m-1)
                {
                    qt = point[i][m].neighbor_node[p2];
                }
                else{
                    qt = point[i][j].neighbor_node[p2];
                }
            }
            else if (j == m)
            {
                if (l == 0)
                {
                    qt = point[i][0].neighbor_node[p2];
                }
                else{
                    qt = point[i][j].neighbor_node[p2];
                }
            }
            else{
                qt = point[i][j].neighbor_node[p2];
            }
        if (qt == node2)
        {
            break;
        }
    }
    //p1, p2Ϊ(i,j)ָ�����ڶ����ָ��

    //***************�����¼ÿ�����ϵ��ⷨ�����ͳ���************//
    a1 = point[i][j].a[p1];
    a2 = point[i][j].a[p2];

    nx1 = - point[i][j].nx[p1];
    ny1 = - point[i][j].ny[p1];
    nx2 = point[i][j].nx[p2];
    ny2 = point[i][j].ny[p2];

    //*******��������r��corner�е���Ӧ����********//
    double xt_v, yt_v;
    xt_v = o[k][l].phi_x(xit,etat);
    yt_v = o[k][l].phi_y(xit,etat);

    //********��������ܶ�**********//
    double rho_vertex;
    //rho_vertex = o[k][l].rholast;
    double xt0, yt0;
    xt0 = o[k][l].phi_x0(xit,etat);
    yt0 = o[k][l].phi_y0(xit,etat);
    if (xt0 == 0.5)
    {
        double xc;
        xc = o[k][l].phi_x0(o[k][l].xi_c,o[k][l].eta_c);
        xt0 = xt0 - 1e-9 * (xt0 - xc);
    }
    rho_vertex = ini_rho(xt0,yt0) * o[k][l].Jacobi_0(xit,etat)
                / o[k][l].Jacobi(xit,etat);
    

    //********�������ڵ㴦�ٶ��ع�ֵ********//
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
        
    //********�����������**********//
    ve = ve - 0.5 * (uc[0] * uc[0] + uc[1] * uc[1]);
    ve = max(1e-9, ve);
    double p_vertex;
    p_vertex = EOS(rho_vertex,ve);

    //**************��������迹*************//
    double rho_center;
        rho_center = o[k][l].rholast;

        double p_center, ve_center;
        ve_center = o[k][l].uxlast[0]*o[k][l].uxlast[0] + o[k][l].uylast[0] * o[k][l].uylast[0];
        ve_center = o[k][l].taulast[0] - ve_center * 0.5;
        p_center = EOS(rho_center,ve_center);
        p_center = max(1e-12, p_center);
        mu = gamma * p_center / rho_center;
        mu = sqrt(mu);  //��ʱmuΪ����
        mu = rho_center * mu;
    
    //***********��������������ϵ�corner force*********//
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