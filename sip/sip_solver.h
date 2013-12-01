/* This program analyses the sample average approximation method for solving the
   stochastic knapsack problem (version inspired from the paper by Kleywegt,
   Shapiro, and Homem-de-Mello on the same topic). The analysis is conducted by
   solving successive sample average approximation problems with increasing
   sample sizes, and with predefined stopping rule.
   To accelerate the process, we implement a dynamic programming procedure to
   optimize the sample average knapsack problem.

Author: Péguy Pierre-Louis
N.B. a portion of this code (normal cdf, random number generator) is due to
the works of others. Credits are explicitly given where credit is due.
*/

/*============================================================================*/
/*============================================================================*/

#ifndef __SIPSOLVER_H__
#define __SIPSOLVER_H__

#include "defs.h"

/*============================================================================*/
/*============================================================================*/
/*Below is the declaration of the NORMP funbction. This function implements the
  Normal CDF.
  */


void NORMP( double Z, double *P, double *Q, double *PDF );


/*============================================================================*/
/*============================================================================*/
/*============================ Function f_bar_n===============================*/
/* This function returns the value of the objective of the current SAA problem
   at a given (candidate) solution.
   */
double f_bar_n( const int x[], double W[][SET], int n, int start );


/*============================================================================*/
/*============================================================================*/
/*========================== Function f_scenario =============================*/
/* This function returns the value of the objective function for 1 scenario only.
   That is, a SAA problem with one sample point. Used in the variance calculus.
   */
double f_scenario( const int x[], double W[][SET], int scenarioNum );


/*============================================================================*/
/*============================================================================*/
/*================Function to manage Active combinations======================*/
/* This function helps make the dynamic programming algorithm efficient.
   It keeps track of what combinations are no longer promising. See the
   description of the dynamic programming for more clarification on promising
   combinations.
   */
void Set_Aktiv( int Aktiv[][POWSET], double f[], int s );


/*============================================================================*/
/*============================================================================*/
/*======================== Function to solve SIP_n ===========================*/
/* This function implements the dynamic programming that solves the SAA knapsack
   problem. We arbitrarily order all the elements. Each stage corresponds to a \
   particular element (ii.e. stage ii corresponds to element ii => decision to be
   made == whether to include element ii or not in the knapsack. At stage ii, we
   can choose up to element ii to put in the knapsack.). The optimal solution is
   found when the last stage is solved. Non-promising states, together with
   their descendants, are eliminated in advance (at each stage). At each stage,
   we keep "active" only the best state among all comparable states to be
   included in the next stage.
   */
void Solve_SIP( int n_k, int start, double W[][SET], int x_sol[], double& z );


// (Almost) Pierre's original Solve_SIP code.  Used here for debugging only
void Solve_SIP_Orig( int n_k, int start, double W[][SET], int x_sol[], double& z );

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=========== Function to Evaluate True Objective Value At any X =============*/
/* This function evaluates the "true" (stochastic) objective value at any
   arbitrary vector X.
   */
double obj_func( const int x[] );


/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=======This Function finds the next Combination of Xi=======================*/
/* This function is used in finding the true optimal solution (to the stochastic
   problem -- not the SAA) through complete explicit enumeration.
   */
void Find_Next( int x[], int last );


/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*====================Function to Find True Optimal===========================*/
/*
   This function finds the true optimal solution of the stochastic problem
   through explicit enumeration, i.e. brute force.
   */
void Optimal( int x_star[], double& z_star );


/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*============= Calculates Sample Variance for x_hat and x_gap ===============*/
/*
   Cf. Güzin's paper
   */
double S_variance( const int x_hat[], int x_gap[], const int m_k, int start, double W2[][SET] );


/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*===================== Function to estimate Gap =============================*/
/*
   This function calculates the gap estimate (Cf. Güzin's paper), checks for the
   stopping criterion, and increases the sample size accordingly.
   */
double Gap_estimate( const int x_hat[], const int n_k, double& s_std, double W2[][SET] );

#endif
