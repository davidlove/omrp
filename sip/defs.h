#ifndef __DEFS_H__
#define __DEFS_H__

/*File Inclusion*/
//#include "rng.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <unistd.h>

/*Defining min and max functions for scalars*/
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*namespace for input-output*/
using namespace std;

/*Constant Definition*/
#define  SET     10          // Number of elements to choose from
#define  C       4           // Cost for each unit of the knapsack size exceeded
#define  KNAPSACK_SIZE       100         // Size of knapsack
#define  COUNT   100000      // Large Constant
#define  PI      3.14159265  // Pi
#define  Nk      11          // # of different n_0 values
#define  Mk      11          // # of different m_0 values
#define  Dk      5           // # of different d values
#define  KFN     25          // Resampling Frequency For Quality Assesment problem
#define  KFM     25          // # of different d values

#define CAP      100          // Maximum increase of sample size when appropriate
#define ALPHA    0.1         // Alpha
#define Z_ALPHA  1.282       // Z_ALPHA
#define REPLICAT 50          // Number of Replications
#define RP       2

#define POWSET   1024        // Cardinality of the power set



#define  FALSE 0
#define  TRUE  1
#define  PRECISION 10000

/*============================================================================*/
/*============================================================================*/

const double rewards[SET]  = {10.1, 11.1, 12.1, 13.1, 14.1, 15.1, 16.1, 17.1, 18.1, 19.1};
const double mu[SET]       = {20.1, 21.1, 22.1, 23.1, 24.1, 25.1, 26.1, 27.1, 28.1, 29.1};
const double sigma[SET] = {5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5};

const int n[Nk] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
const int m[Mk] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
const double d[Dk] = {0.05, 0.04, 0.03, 0.02, 0.01};

#endif
