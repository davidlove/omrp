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

/*File Inclusion*/
#include <cmath>
#include <stdio.h>
// #include "rng.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <time.h>

#include "randomc.h"                   // define classes for random number generators

#define RANDOM_GENERATOR TRandomMersenne
#include "stocc.h"                     // define random library classes



/*Defining min and max functions for scalars*/
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/*namespace for input-output*/
using namespace std;

/*Constant Definition*/
#define  SET     10          // Number of elements to choose from
#define  c       4           // Cost for each unit of the knapsack size exceeded
#define  q       100         // Size of knapsack
#define  Count   100000      // Large Constant
#define  PI      3.14159265  // Pi
#define  n_1     500          // Initial Sample Size for SAA problem
#define  m_1     500          // Initial Sample size for gap estimate problem
#define  eps     0.426733    // Epsillon

#define CAP      100          // Maximum increase of sample size when appropriate
#define alpha    0.1         // Alpha
#define Z_alpha  1.282       // Z_alpha
#define Replicat 50          // Number of Replications
#define rp       2

#define PowSet   1024        // Cardinality of the power set



#define  FALSE 0
#define  TRUE  1
#define  Precision 10000

/*============================================================================*/
/*============================================================================*/

double rewards[SET]  = {10.1, 11.1, 12.1, 13.1, 14.1, 15.1, 16.1, 17.1, 18.1, 19.1};
double mu[SET]       = {20.1, 21.1, 22.1, 23.1, 24.1, 25.1, 26.1, 27.1, 28.1, 29.1};
double sigma[SET] = {5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5};

int x_gap[SET];                  // Stores solution of gap problem
double z_gap;                    // objective value of gap problem
int x_star[SET];                 // True Solution
double z_star;                   // True optimal objective
double F_x_starN[Replicat];    // obj. value at true solution of each SAA
double Fbar_x_hat[Replicat];     // obj. value at SAA solution
int Jumps;


//vector<double> v2(Count);        // Stores large number of i.i.d. random numbers
//vector<double> v(Count);         // Stores large number of i.i.d. random numbers

double W[Count][SET];            // Stores random weights
double W2[Count][SET];           // Stores random weights

double f[PowSet];                // For use in dynamic programming algorithm
double delta[Count][PowSet];     // Cf. Dynamic Programming Function (Solve_SIP)
int Active[SET][PowSet];         // Idem
/*============================================================================*/
/*============================================================================*/
/*Below is the declaration of the NORMP funbction. This function implements the
  Normal CDF.
*/
double zero = 0.0, one = 1.0, half = 0.5;
double ltone = 7.0, utzero = 18.66, con = 1.28;

double P1, P2, Q, X;

void NORMP(double Z, double *P, double *Q, double *PDF) {

/*	Normal distribution probabilities accurate to 1.e-15.
 	Z = no. of standard deviations from the mean.
	P, Q = probabilities to the left & right of Z.   P + Q = 1.
        PDF = the probability density.
        Based upon algorithm 5666 for the error function, from:
        Hart, J.F. et al, 'Computer Approximations', Wiley 1968

        Programmer: Alan Miller

        Latest revision - 30 March 1986       */

    double P0,P1,P2,P3,P4,P5,P6;
    double Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    double CUTOFF,EXPNTL,ROOT2PI,ZABS;

	P0=220.2068679123761; P1=221.2135961699311; P2=112.0792914978709;
    P3= 33.91286607838300; P4= 6.373962203531650;
    P5=  0.7003830644436881; P6= 0.3526249659989109E-01;
    Q0=440.4137358247522; Q1=793.8265125199484; Q2=637.3336333788311;
   	Q3=296.5642487796737; Q4= 86.78073220294608; Q5=16.06417757920695;
    Q6=  1.755667163182642; Q7=0.8838834764831844E-01;
    CUTOFF= 7.071; ROOT2PI=2.506628274631001;

	ZABS = fabs(Z);

//	|Z| > 37.0

	if (ZABS > 37.0) {
	  *PDF = zero;
	  if (Z >  zero) {
	    *P = one;
	    *Q = zero;
      }
	  else {
	    *P = zero;
	    *Q = one;
	  }
	  return;
	}

//  |Z| <= 37.0

	EXPNTL = exp(-half*ZABS*ZABS);
	*PDF = EXPNTL/ROOT2PI;

//  |Z| < CUTOFF = 10/sqrt(2)

	if (ZABS < CUTOFF)
	  *P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS +
     	       P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS  +
     	       Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + Q0);

//  |Z| >= CUTOFF

	else
	  *P = *PDF/(ZABS + one/(ZABS + 2.0/(ZABS + 3.0/(ZABS + 4.0/
     		(ZABS + 0.65)))));

	if (Z < zero)
	  *Q = one - *P;
	else {
	  *Q = *P;
	  *P = one - *Q;
	}
}
/*============================================================================*/
/*============================================================================*/
/*============================ Function f_bar_n===============================*/
/* This function returns the value of the objective of the current SAA problem
   at a given (candidate) solution.
*/
double f_bar_n(const int x[], double W[][SET], int n, int start){
    double Sum1 = 0;
    double Sum2 = 0;

    for (int i = 0; i < SET; i++){
        Sum1 = rewards[i]*x[i] + Sum1;
    }

    double Temp = 0;
    for (int j = start; j < start + n; j++){
        for (int i = 0; i < SET; i ++)
             Temp =  Temp + W[j][i]*x[i];
        Temp = Temp - q;
        Sum2 = Sum2 + max(0, Temp);
        Temp = 0;
    }

    return (Sum1 - double(c*Sum2/n));
}
/*============================================================================*/
/*============================================================================*/
/*========================== Function f_scenario =============================*/
/* This function returns the value of the objective function for 1 scenario only.
   That is, a SAA problem with one sample point. Used in the variance calculus.
*/
double f_scenario(const int x[], double W[][SET], int j){
    double Sum1 = 0;
    double Sum2 = 0;

    for (int i = 0; i < SET; i++){
        Sum1 = rewards[i]*x[i] + Sum1;
    }

    double Temp = 0;
    for (int i = 0; i < SET; i ++)
        Temp =  Temp + W[j][i]*x[i];
    Temp = Temp - q;
    Sum2 = Sum2 + max(0, Temp);

    return (Sum1 - double(c*Sum2));
}
/*============================================================================*/
/*============================================================================*/
/*================Function to manage Active combinations======================*/
/* This function helps make the dynamic programming algorithm efficient.
   It keeps track of what combinations are no longer promising. See the
   description of the dynamic programming for more clarification on promising
   combinations.
*/
void Set_Aktiv(int Aktiv[][PowSet], double f[], int s){

     for(int i = 1; i < pow(2, s + 1); i = i + 1){
         if (!Aktiv[s][i]){
             Aktiv[s + 1][2*i] = 0;
             Aktiv[s + 1][2*i + 1] = 0;
         }
         else{
           if(i%2 == 1){
                if (f[i] < f[i - 1]){
                    Aktiv[s + 1][2*i] = 0;
                    Aktiv[s + 1][2*i + 1] = 0;
                }
           }
         }
     }
}
/*============================================================================*/
/*============================================================================*/
/*======================== Function to solve SIP_n ===========================*/
/* This function implements the dynamic programming that solves the SAA knapsack
   problem. We arbitrarily order all the elements. Each stage corresponds to a \
   particular element (i.e. stage i corresponds to element i => decision to be
   made == whether to include element i or not in the knapsack. At stage i, we
   can choose up to element i to put in the knapsack.). The optimal solution is
   found when the last stage is solved. Non-promising states, together with
   their descendants, are eliminated in advance (at each stage). At each stage,
   we keep "active" only the best state among all comparable states to be
   included in the next stage.
*/
void Solve_SIP(int n_k, int start, double W[][SET], int x_sol[], double& z){

     /*Initializations*/
     for(int i = 0; i < pow(2, SET); i++){
		 f[i] = 0;
         for (int k = 0; k < SET; k++)
              Active[k][i] = 1;
         for(int j = start; j < start + n_k; j++)
              delta[j][i] = - q;
     }

     /*Stage Zero Computations*/
     f[0] = 0;
     f[1] = rewards[0];
     for(int j = start; j < start + n_k; j++){
         delta[j][1] = delta[j][1] + W[j][0];
         if (W[j][0] > q)
             f[1] = f[1] - (double(c)/double(n_k))*(W[j][0] - q);
     }

     /*Inactivate "inferior" combinations*/
     Set_Aktiv(Active, f, 0);

     /*Continue with procedure up to stage SET - 2*/
     for(int s = 1; s < SET - 1; s++){
         for(int i = int(pow(2, s + 1)) - 2; i >= 0; i = i - 2){
             if(Active[s][i]){
               f[i] = f[i/2];
               f[i + 1] = f[i] + rewards[s];
               for(int j = start; j < start + n_k; j++){
                   delta[j][i] = delta[j][i/2];
                   delta[j][i + 1] = delta[j][i] + W[j][s];
                   if((delta[j][i + 1] > 0) && (delta[j][i] <= 0))
                      f[i + 1] = f[i + 1] - (double(c)/double(n_k))*delta[j][i + 1];
                   if((delta[j][i + 1] > 0) && (delta[j][i] > 0))
                      f[i + 1] = f[i + 1] - (double(c)/double(n_k))*W[j][s];
               }
             }
         }
         Set_Aktiv(Active, f, s);
     }

     /*Last Stage -- And finding optimal solution in decimal*/
     int sol = 0;
     z = 0;
     for(int i = int(pow(2, SET)) - 2; i >= 0; i = i - 2){
             if(Active[SET - 1][i]){
               f[i] = f[i/2];
               if(f[i] > z){
                  sol = i;
                  z = f[i];
               }
             }
             if(Active[SET - 1][i + 1]){
                f[i + 1] = f[i] + rewards[s];
                for(int j = start; j < start + n_k; j++){
                   delta[j][i] = delta[j][i/2];
                   delta[j][i + 1] = delta[j][i] + W[j][s];
                   if((delta[j][i + 1] > 0) && (delta[j][i] <= 0))
                      f[i + 1] = f[i + 1] - (double(c)/double(n_k))*delta[j][i + 1];
                   if((delta[j][i + 1] > 0) && (delta[j][i] > 0))
                      f[i + 1] = f[i + 1] - (double(c)/double(n_k))*W[j][s];
                }
                if(f[i + 1] > z){
                  sol = i + 1;
                  z = f[i + 1];
               }
             }
     }

     /*Converting optimal solution to 0-1 format*/
     for(int i = 0; i < SET; i++)
        x_sol[i] = 0;

     int k = SET - 1;
     while(sol > 1){
        x_sol[k] = sol % 2;
        sol = int(floor(double(sol/2)));
		k--;
     }
     x_sol[k] = 1;

}
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=========== Function to Evaluate True Objective Value At any X =============*/
/* This function evaluates the "true" (stochastic) objective value at any
   arbitrary vector X.
*/
double obj_func(const int x[]){
   double mu_x = 0;
   double sigma_x = 0;
   double prob;
   double z_temp;

   for (int i = 0; i < SET; i++)
         mu_x = mu_x + mu[i]*x[i];
   mu_x = mu_x - q;

   double Sum = 0;
   for (int i = 0; i < SET; i++)
        Sum = Sum + pow(sigma[i], 2)*pow(x[i], 2);
   sigma_x = sqrt(Sum);

   if (sigma_x != 0)
       prob = double(mu_x/sigma_x);
   else
       prob = double(mu_x/0.00001);


   NORMP(prob, &P1, &P2, &Q);
   double Sum1 = 0;
   for (int i = 0; i < SET; i++)
       Sum1 = Sum1 + rewards[i]*x[i];
   z_temp = Sum1  - c*(mu_x*P1 + (sigma_x/(sqrt(2*PI)*exp(0.5*pow(prob, 2)))));

   return z_temp;

}
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=======This Function finds the next Combination of Xi=======================*/
/* This function is used in finding the true optimal solution (to the stochastic
   problem -- not the SAA) through complete explicit enumeration.
*/
void Find_Next(int x[], int last){
    if (x[last] == 0)
        x[last] = 1;
    else{
        x[last] = 0;
        last = last - 1;
        if(!(last < 0))
           Find_Next(x, last);
    }
}

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*====================Function to Find True Optimal===========================*/
/*
   This function finds the true optimal solution of the stochastic problem
   through explicit enumeration, i.e. brute force.
*/
void Optimal(int x_star[], double& z_star){
 int x[SET];
 double temp;
 z_star = -100000;


 for (int i = 0; i < SET; i++)
     x[i] = 0;

 for (int j = 0; j < pow(2, SET); j++){
      temp = obj_func(x);
      if ( temp > z_star){
          z_star = temp;
          for (int i = 0; i < SET; i++)
               x_star[i] = x[i];
      }
      Find_Next(x, SET - 1);

 }

}
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*============= Calculates Sample Variance for x_hat and x_gap ===============*/
/*
   Cf. Güzin's paper
*/
double S_variance(const int x_hat[], int x_gap[], const int m_k, int start, double W2[][SET]){
      double sum = 0;
      for (int j = start; j < start + m_k; j++){
           sum = sum + pow(f_scenario(x_hat, W2, j) - f_scenario(x_gap, W2, j) - f_bar_n(x_hat, W2, m_k, start) + f_bar_n(x_gap, W2, m_k, start) ,2);
      }

      return double(sum/(m_k - 1));

}
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*===================== Function to estimate Gap =============================*/
/*
  This function calculates the gap estimate (Cf. Güzin's paper), checks for the
  stopping criterion, and increases the sample size accordingly.
*/
void Gap_estimate(const int x_hat[], int& m_k, double& mu_hat, double W2[][SET], int& stop){

     int start = 0;
     int m = Count/rp;
     mu_hat = 0.0;
     double s_var = 0.0;

     for(int i = 0; i < rp; i++){
        start = i*m;
        Solve_SIP(m_k, start, W2, x_gap, z_gap);
        mu_hat = mu_hat + z_gap - f_bar_n(x_hat, W2, m_k, start);
        s_var = s_var + S_variance(x_hat, x_gap, m_k, start, W2);
     }

     mu_hat = double(mu_hat/rp);
     s_var = double(s_var/rp);


     if ( mu_hat + (1 + Z_alpha*sqrt(s_var))/sqrt(double(m_k)) <= eps ){
         stop = 1;
     }
     else{

           m_k = m_k + 1;
        /* 
		   int old_mk = m_k;
		 
           double b;
           double delta;

	       b = sqrt(s_var)*Z_alpha + 1;

           delta = pow(b, 2) + 4*eps*m_k*mu_hat;

           m_k = min(CAP + m_k, ceil((pow(b, 2) + delta + 2*b*sqrt(delta))/(4*pow(eps, 2))));
		   Jumps = Jumps + (m_k >= old_mk + CAP);

        */
           cout<<"Test new m_k = "<<m_k<<endl;
            
     }

     cout<<"Test Output: Gap Estimate = "<<mu_hat<<endl;
}
/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=========================== Main Function ==================================*/
int main(){
    /* Defining/Declaring Necessary Parameters*/

    time_t seconds1, seconds2, total_time;     // Keep Track of Time



    int k[Replicat];                          // Iterations
    int n_k[Replicat];                        // SAA sample size
    int m_k[Replicat];                        // Gap Problem sample size
    double mu_hat[Replicat];                  // Gap for each replication
    double gap[Replicat];
    double CI[Replicat];
	int ref[Replicat];                        // Number of times we refresh

    int coverage = 0;                         // Coverage Probability


    int x_hat[Replicat][SET];                 // Cand. Solutions
    double z_hat[Replicat];                   // Obj. Values

    int x_temp[SET];                          // Temp Variables
    double z_temp;
    
	/*Save Output to File*/
   ofstream OutFile;
   OutFile.open("newran_f2rp_1_gsn_n0500_refreshed25.txt");

   OutFile<<setw(3)<<"REP"<<setw(5)<<"n_k"<<setw(5)<<"K"<<setw(11)<<"GAP"<<setw(15)<<"Estimate"<<setw(15)<<"Solution"
          <<setw(11)<<"Z_nk(x^)"<<setw(14)<<"Z(x_n*)"<<setw(5)<<"Jumps"<<endl<<endl;
    	
    //Initializations
    for (int i = 0; i < Replicat; i++){
        k[i] = 0;
		ref[i] = 0;
        n_k[i] = n_1;
        m_k[i] = m_1;
	}

    int stop = 0;                             // Boolean Variable
    
    Optimal(x_star, z_star);                  // True Optimals
    cout<<"Test OutPut ; z* = "<<z_star<<endl;
    cout<<"X Star = ";
    for (int i = 0 ; i < SET; i++)
        cout<<x_star[i];
    cout<<endl;

    int r_seed = 1;

    seconds1 = time (NULL);                   // Start Time

    for(int u = 0; u < Replicat; u++){
      Jumps = 0;
		                                      
     /* RNG y = 1313 + 171*u; //  12345 + 121*u;
      for (int i = 0; i < SET; i++){          // Generating Random Variates
         y.normal(v, mu[i], sigma[i]);
         for (int j = 0; j < Count; j++) {
              W[j][i] = v[j];
         }
      }


      for (int i = 0; i < SET; i++){          // Generating Random Variates
         y.normal(v2, mu[i], sigma[i]);
         for (int j = 0; j < Count; j++){
              W2[j][i] = v2[j];
         }
      }*/
      
	  
	  int seed = 1234567890 + 12345678*u;  //int32 seed = time(0);                // random seed
      StochasticLib1 sto(seed);            // make instance of random library

	  for (int i = 0; i < SET; i++){       // Generating Random Variates
		  for (int j = 0; j < Count; j++){
              W[j][i] = sto.Normal(mu[i], sigma[i]);
			  W2[j][i] = sto.Normal(mu[i], sigma[i]);
		  }         
      }
      

      while (!stop){
          k[u]++;
          cout<<"Iteration Test Candidate: "<<u<<"."<<k[u]<<endl;
	      Solve_SIP(n_k[u], 0, W, x_temp, z_temp);
		  /*Refreshing*/
		 if ((k[u]%25)==0){
			  ref[u]++;
			  seed = 1234567890 + 1313*k[u] + 12345678*u;
			  StochasticLib1 sto(seed);
			  for (int i = 0; i < SET; i++){          
                   for (int j = 0; j < Count; j++){
                        W2[j][i] = sto.Normal(mu[i], sigma[i]);
				   }
		      } 
              
			  /*RNG y = 1313 + 17*k[u];
			  for (int i = 0; i < SET; i++){          // Generating Random Variates
                   y.normal(v2, mu[i], sigma[i]);
                   for (int j = 0; j < Count; j++){
                        W2[j][i] = v2[j];
                   }
              } */
		  } 

          Gap_estimate(x_temp, m_k[u], mu_hat[u], W2, stop);
          n_k[u] = m_k[u];

      }

      for (int i = 0; i < SET; i++){
         x_hat[u][i] = x_temp[i];
      }
      z_hat[u] = z_temp;
	  F_x_starN[u] = obj_func(x_temp);

      gap[u] = z_star - F_x_starN[u];
      if(gap[u] < eps)
         coverage++;

      stop = 0;

	  OutFile<<setw(3)<< u + 1<<setw(5)<<n_k[u]<<setw(5)<<k[u]<<setw(11)<<gap[u]<<setw(15)<<mu_hat[u];
      OutFile<<setw(6);
      for (int t = 0; t < SET; t++)
            OutFile<<x_hat[u][t];
      OutFile<<setw(11)<<z_hat[u]<<setw(13)<<F_x_starN[u]<<setw(5)<<Jumps<<endl;
    }

    seconds2 = time (NULL);                       // Stopping Time
    total_time = seconds2 - seconds1;             // Total Time


   /* Printing Out Results */

   /* Averages */
   double avg_It = 0.0;
   double avg_gap = 0.0;
   double avg_gap_est = 0.0;
   double avg_Nt = 0.0;

   
   for (int u = 0; u < Replicat; u++){
       avg_It = avg_It + k[u];
       avg_gap = avg_gap + gap[u];
       avg_gap_est = avg_gap_est + mu_hat[u];
       avg_Nt = avg_Nt + n_k[u];
   }
   avg_It = avg_It/Replicat;
   avg_gap = avg_gap/Replicat;
   avg_gap_est = avg_gap_est/Replicat;
   avg_Nt = avg_Nt/Replicat;

   OutFile<<endl<<endl;
   OutFile<<"True Solution: X = ["<<x_star[0];
   for (int i = 1; i < SET; i++)
       OutFile<<", "<<x_star[i];
   OutFile<<"]"<<endl;
   OutFile<<"True Optimal Value = "<<z_star<<endl;
   double percent = double(coverage)/Replicat;
   OutFile<<"Coverage P = "<<percent<<endl;
   OutFile<<"Average Gap Estimate = "<<avg_gap_est<<endl;
   OutFile<<"Average True Gap = "<<avg_gap<<endl;
   OutFile<<"Average N_T = "<<avg_Nt<<endl;
   OutFile<<"Average Stopping Time = "<<avg_It<<endl;
   OutFile<<"Total Time in secs = "<<total_time<<endl;

   OutFile.close();

   cout<<"===================The End==========================="<<endl;

   getchar();

   return 0;
}
/*============================= The End ======================================*/


