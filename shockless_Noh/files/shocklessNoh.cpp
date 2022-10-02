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
