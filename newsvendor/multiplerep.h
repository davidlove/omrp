#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <search.h>        //to use qsort --> quicksort function
#include <unistd.h>
#include <string>
#include "circ_sort.h"

using namespace std;
using std::string;

#define DEBUG_DEMAND 0  // True to debug demand generation
#define DEBUG_GAP    0  // True to debug gap-estimation 
#define DEBUG_COVER  0  // True to debug coverage estimation
#define DEBUG_OBJ    0  // True to debug optimality gap determination
#define DEBUG        DEBUG_DEMAND || DEBUG_GAP || DEBUG_COVER || DEBUG_OBJ

#define QUICKSORT        0
#define BUBBLESORT       0
#define PEBBLEBUBBLESORT 0
#define INSERTIONSORT    1

// calculate exact solution value
double FindObjValue( double xhat, double price, double cost, double a, double b, ofstream& debugFile );

// Calculate the expected profit for some candidate solution
double FindExpectedProfit( double xhat, double price, double cost, double a, double b, ofstream& debugfile );

// Calculate the expected profit for some candidate solution
double FindSampleProfit( double xhat, double xi, double price, double cost, ofstream& debugfile );

// generate the demaind
void GenerateDemand( circ_sort& demand, const int sampleSize, const double lowerLimit, const double upperLimit, ofstream& debugFile );

// Find the gap estimate for a given batch
double FindGapEstimate( circ_sort& demand, const int sampleSize, const int discreteQuantile, const double xnstar, const double xhat, const int r, const int c, ofstream& debugFile );

// Find the coverage probability of the gap estimates
double FindCoverageProbability( const double* gapest, const double* gbar, const int rep, const double ng_nonoverlap, const int ng, const double za, const double OptGap, const int degreesFreedom, ofstream& debugFile, double* var = NULL, double* gapbar = NULL );

//this is the required Compare function for quicksort
int Compare( const void *arg1, const void *arg2 );

//generates truncated normal variates
double TruncNormal ( double mean1, double stdev1 ); 

//generates Uniform(a,b) variates
double Uniform ( double a, double b ); 

//calculates combined sample variance 
double UnusedCalcVar ( double *y1, double *y2, double sum1, double sum2, int iss, int r ); 

//calculates plain old sample variance 
double CalculateVariance ( const double* y1, double ybar1, int size, int repCounter, int degreesFreedom ); 
