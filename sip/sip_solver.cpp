#include "sip_solver.h"

void NORMP( double Z, double *P, double *Q, double *PDF )
{

   /*	Normal distribution probabilities accurate to 1.e-15.
        Z = no. of standard deviations from the mean.
        P, Q = probabilities to the left & right of Z.   P + Q = 1.
        PDF = the probability density.
        Based upon algorithm 5666 for the error function, from:
        Hart, J.F. et al, 'Computer Approximations', Wiley 1968

   Programmer: Alan Miller

   Latest revision - 30 March 1986       */
   double zero = 0.0, one = 1.0, half = 0.5;
   double ltone = 7.0, utzero = 18.66, con = 1.28;

   double P0, P1, P2, P3, P4, P5, P6;
   double Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
   double CUTOFF, EXPNTL, ROOT2PI, ZABS;

   P0 = 220.2068679123761;
   P1 = 221.2135961699311;
   P2 = 112.0792914978709;
   P3 = 33.91286607838300;
   P4 = 6.373962203531650;
   P5 =  0.7003830644436881;
   P6 = 0.3526249659989109E-01;
   Q0 = 440.4137358247522;
   Q1 = 793.8265125199484;
   Q2 = 637.3336333788311;
   Q3 = 296.5642487796737;
   Q4 = 86.78073220294608;
   Q5 = 16.06417757920695;
   Q6 =  1.755667163182642;
   Q7 = 0.8838834764831844E-01;
   CUTOFF = 7.071;
   ROOT2PI = 2.506628274631001;

   ZABS = fabs( Z );

   //	|Z| > 37.0

   if ( ZABS > 37.0 )
   {
      *PDF = zero;
      if ( Z >  zero )
      {
         *P = one;
         *Q = zero;
      }
      else
      {
         *P = zero;
         *Q = one;
      }
      return;
   }

   //  |Z| <= 37.0

   EXPNTL = exp( -half * ZABS * ZABS );
   *PDF = EXPNTL / ROOT2PI;

   //  |Z| < CUTOFF = 10/sqrt(2)

   if ( ZABS < CUTOFF )
      *P = EXPNTL * ( ( ( ( ( ( P6 * ZABS + P5 ) * ZABS + P4 ) * ZABS + P3 ) * ZABS +
                          P2 ) * ZABS + P1 ) * ZABS + P0 ) / ( ( ( ( ( ( ( Q7 * ZABS + Q6 ) * ZABS  +
                                Q5 ) * ZABS + Q4 ) * ZABS + Q3 ) * ZABS + Q2 ) * ZABS + Q1 ) * ZABS + Q0 );

   //  |Z| >= CUTOFF

   else
      *P = *PDF / ( ZABS + one / ( ZABS + 2.0 / ( ZABS + 3.0 / ( ZABS + 4.0 /
                                   ( ZABS + 0.65 ) ) ) ) );

   if ( Z < zero )
      *Q = one - *P;
   else
   {
      *Q = *P;
      *P = one - *Q;
   }
}



double f_bar_n( const int x[], double W[][SET], int n, int start )
{
   double totalRewards = 0;
   double totalOverKnapsackSize = 0;

   for ( int ii = 0; ii < SET; ii++ )
   {
      totalRewards = rewards[ii] * x[ii] + totalRewards;
   }

   double itemsWeight;
   for ( int jj = start; jj < start + n; jj++ )
   {
      itemsWeight = 0;
      for ( int ii = 0; ii < SET; ii ++ )
         itemsWeight =  itemsWeight + W[jj][ii] * x[ii];
      itemsWeight = itemsWeight - KNAPSACK_SIZE;
      totalOverKnapsackSize = totalOverKnapsackSize + max( 0, itemsWeight );
   }

   return ( totalRewards - double( C * totalOverKnapsackSize / n ) );
}


double f_scenario( const int x[], double W[][SET], int scenarioNum )
{
   double totalRewards = 0;
   double totalOverKnapsackSize = 0;

   for ( int ii = 0; ii < SET; ii++ )
   {
      totalRewards = rewards[ii] * x[ii] + totalRewards;
   }

   double itemsWeight = 0;
   for ( int ii = 0; ii < SET; ii ++ )
      itemsWeight =  itemsWeight + W[scenarioNum][ii] * x[ii];
   itemsWeight = itemsWeight - KNAPSACK_SIZE;
   totalOverKnapsackSize = totalOverKnapsackSize + max( 0, itemsWeight );

   return ( totalRewards - double( C * totalOverKnapsackSize ) );
}


void Set_Aktiv( int Aktiv[][POWSET], double f[], int s )
{

   for( int ii = 1; ii < pow( 2, s + 1 ); ii = ii + 1 )
   {
      if ( !Aktiv[s][ii] )
      {
         Aktiv[s + 1][2 * ii] = 0;
         Aktiv[s + 1][2 * ii + 1] = 0;
      }
      else
      {
         if( ii % 2 == 1 )
         {
            if ( f[ii] < f[ii - 1] )
            {
               Aktiv[s + 1][2 * ii] = 0;
               Aktiv[s + 1][2 * ii + 1] = 0;
            }
         }
      }
   }
}


// CHANGED BY DAVID LOVE -- I was running into segmentation faults caused by
// delta[][] being too big to be allocated.  In order to save memory, I have
// changed the allocation to a static one, and I have changed the indices.
// Delta only needs o be as large as [n_k] in the first index, whereas before
// it was [start + n_k] and before that [COUNT] for a very large number. 
// However, only the final n_k rows were used.  Below, I first changed all 
// instances of "start" to "int(0)", to make the loop indices correct.
// Then I changed all instances of "W[jj]" to "W[start+jj]".  
// Finally, the declaration of delta[][] was changed to be dynamic.
// You can get the old code back by reversing these changes.
void Solve_SIP( int n_k, int start, double W[][SET], int x_sol[], double& z )
{
   // CHANGED BY DAVID LOVE -- made delta dynamically allocated.
   double** delta;     // Cf. Dynamic Programming Function (Solve_SIP)
   //double delta[int(0) + n_k][POWSET];     // Cf. Dynamic Programming Function (Solve_SIP)
   double f[POWSET];                // For use in dynamic programming algorithm
   int Active[SET][POWSET];         // Idem

   // CHANGED BY DAVID LOVE -- Dynamically allocating the array delta
   delta = new double*[int(0) + n_k];
   for( size_t ii = 0; ii < int(0) + n_k; ii++ )
      delta[ii] = new double[POWSET];

   /*Initializations*/
   for( int ii = 0; ii < pow( 2, SET ); ii++ )
   {
      f[ii] = 0;
      for ( int k = 0; k < SET; k++ )
         Active[k][ii] = 1;
      for( int jj = int(0); jj < int(0) + n_k; jj++ )
         delta[jj][ii] = - KNAPSACK_SIZE;
   }

   /*Stage Zero Computations*/
   f[0] = 0;
   f[1] = rewards[0];
   for( int jj = int(0); jj < int(0) + n_k; jj++ )
   {
      delta[jj][1] = delta[jj][1] + W[start+jj][0];
      if ( W[start+jj][0] > KNAPSACK_SIZE )
         f[1] = f[1] - ( double( C ) / double( n_k ) ) * ( W[start+jj][0] - KNAPSACK_SIZE );
   }

   /*Inactivate "inferior" combinations*/
   Set_Aktiv( Active, f, 0 );

   /*Continue with procedure up to stage SET - 2*/
   int s;
   for( s = 1; s < SET - 1; s++ )
   {
      for( int ii = int( pow( 2, s + 1 ) ) - 2; ii >= 0; ii = ii - 2 )
      {
         if( Active[s][ii] )
         {
            f[ii] = f[ii / 2];
            f[ii + 1] = f[ii] + rewards[s];
            for( int jj = int(0); jj < int(0) + n_k; jj++ )
            {
               delta[jj][ii] = delta[jj][ii / 2];
               delta[jj][ii + 1] = delta[jj][ii] + W[start+jj][s];
               if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] <= 0 ) )
                  f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * delta[jj][ii + 1];
               if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] > 0 ) )
                  f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * W[start+jj][s];
            }
         }
      }
      Set_Aktiv( Active, f, s );
   }

   /*Last Stage -- And finding optimal solution in decimal*/
   int sol = 0;
   z = 0;
   for( int ii = int( pow( 2, SET ) ) - 2; ii >= 0; ii = ii - 2 )
   {
      if( Active[SET - 1][ii] )
      {
         f[ii] = f[ii / 2];
         if( f[ii] > z )
         {
            sol = ii;
            z = f[ii];
         }
      }
      if( Active[SET - 1][ii + 1] )
      {
         f[ii + 1] = f[ii] + rewards[s];
         for( int jj = int(0); jj < int(0) + n_k; jj++ )
         {
            delta[jj][ii] = delta[jj][ii / 2];
            delta[jj][ii + 1] = delta[jj][ii] + W[start+jj][s];
            if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] <= 0 ) )
               f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * delta[jj][ii + 1];
            if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] > 0 ) )
               f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * W[start+jj][s];
         }
         if( f[ii + 1] > z )
         {
            sol = ii + 1;
            z = f[ii + 1];
         }
      }
   }

   /*Converting optimal solution to 0-1 format*/
   for( int ii = 0; ii < SET; ii++ )
      x_sol[ii] = 0;

   int k = SET - 1;
   while( sol > 1 )
   {
      x_sol[k] = sol % 2;
      sol = int( floor( double( sol / 2 ) ) );
      k--;
   }
   x_sol[k] = 1;

   // CHANGED BY DAVID LOVE -- Deleting the array delta
   for( size_t ii = 0; ii < int(0) + n_k; ii++ )
      delete[] delta[ii];
   delete[] delta;
}


void Solve_SIP_Orig( int n_k, int start, double W[][SET], int x_sol[], double& z )
{
   double delta[start + n_k][POWSET];     // Cf. Dynamic Programming Function (Solve_SIP)
   double f[POWSET];                // For use in dynamic programming algorithm
   int Active[SET][POWSET];         // Idem

   /*Initializations*/
   for( int ii = 0; ii < pow( 2, SET ); ii++ )
   {
      f[ii] = 0;
      for ( int k = 0; k < SET; k++ )
         Active[k][ii] = 1;
      for( int jj = start; jj < start + n_k; jj++ )
         delta[jj][ii] = - KNAPSACK_SIZE;
   }

   /*Stage Zero Computations*/
   f[0] = 0;
   f[1] = rewards[0];
   for( int jj = start; jj < start + n_k; jj++ )
   {
      delta[jj][1] = delta[jj][1] + W[jj][0];
      if ( W[jj][0] > KNAPSACK_SIZE )
         f[1] = f[1] - ( double( C ) / double( n_k ) ) * ( W[jj][0] - KNAPSACK_SIZE );
   }

   /*Inactivate "inferior" combinations*/
   Set_Aktiv( Active, f, 0 );

   /*Continue with procedure up to stage SET - 2*/
   int s;
   for( s = 1; s < SET - 1; s++ )
   {
      for( int ii = int( pow( 2, s + 1 ) ) - 2; ii >= 0; ii = ii - 2 )
      {
         if( Active[s][ii] )
         {
            f[ii] = f[ii / 2];
            f[ii + 1] = f[ii] + rewards[s];
            for( int jj = start; jj < start + n_k; jj++ )
            {
               delta[jj][ii] = delta[jj][ii / 2];
               delta[jj][ii + 1] = delta[jj][ii] + W[jj][s];
               if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] <= 0 ) )
                  f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * delta[jj][ii + 1];
               if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] > 0 ) )
                  f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * W[jj][s];
            }
         }
      }
      Set_Aktiv( Active, f, s );
   }

   /*Last Stage -- And finding optimal solution in decimal*/
   int sol = 0;
   z = 0;
   for( int ii = int( pow( 2, SET ) ) - 2; ii >= 0; ii = ii - 2 )
   {
      if( Active[SET - 1][ii] )
      {
         f[ii] = f[ii / 2];
         if( f[ii] > z )
         {
            sol = ii;
            z = f[ii];
         }
      }
      if( Active[SET - 1][ii + 1] )
      {
         f[ii + 1] = f[ii] + rewards[s];
         for( int jj = start; jj < start + n_k; jj++ )
         {
            delta[jj][ii] = delta[jj][ii / 2];
            delta[jj][ii + 1] = delta[jj][ii] + W[jj][s];
            if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] <= 0 ) )
               f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * delta[jj][ii + 1];
            if( ( delta[jj][ii + 1] > 0 ) && ( delta[jj][ii] > 0 ) )
               f[ii + 1] = f[ii + 1] - ( double( C ) / double( n_k ) ) * W[jj][s];
         }
         if( f[ii + 1] > z )
         {
            sol = ii + 1;
            z = f[ii + 1];
         }
      }
   }

   /*Converting optimal solution to 0-1 format*/
   for( int ii = 0; ii < SET; ii++ )
      x_sol[ii] = 0;

   int k = SET - 1;
   while( sol > 1 )
   {
      x_sol[k] = sol % 2;
      sol = int( floor( double( sol / 2 ) ) );
      k--;
   }
   x_sol[k] = 1;

}


double obj_func( const int x[] )
{
   double mu_x = 0;
   double sigma_x = 0;
   double zScore;
   double z_temp;

   for ( int ii = 0; ii < SET; ii++ )
      mu_x = mu_x + mu[ii] * x[ii];
   mu_x = mu_x - KNAPSACK_SIZE;

   double Sum = 0;
   for ( int ii = 0; ii < SET; ii++ )
      Sum = Sum + pow( sigma[ii], 2 ) * pow( x[ii], 2 );
   sigma_x = sqrt( Sum );

   double totalReward = 0;
   for ( int ii = 0; ii < SET; ii++ )
      totalReward = totalReward + rewards[ii] * x[ii];
   
   if ( sigma_x == 0 )
   {
      //zScore = double( mu_x / 0.00001 );
      z_temp = totalReward - C*( max( 0, mu_x ) );
   }
   else
   {
      zScore = double( mu_x / sigma_x );

      double P1, P2, Q;
      NORMP( zScore, &P1, &P2, &Q );
      z_temp = totalReward  - C*( mu_x*P1 + ( sigma_x/( sqrt( 2*PI )*exp( 0.5*pow( zScore, 2 ) ) ) ) );
   }

   return z_temp;

}


void Find_Next( int x[], int last )
{
   if ( x[last] == 0 )
      x[last] = 1;
   else
   {
      x[last] = 0;
      last = last - 1;
      if( !( last < 0 ) )
         Find_Next( x, last );
   }
}


void Optimal( int x_star[], double& z_star )
{
   int x[SET];
   double temp;
   z_star = -100000;


   for ( int ii = 0; ii < SET; ii++ )
      x[ii] = 0;

   for ( int jj = 0; jj < pow( 2, SET ); jj++ )
   {
      temp = obj_func( x );
      if ( temp > z_star )
      {
         z_star = temp;
         for ( int ii = 0; ii < SET; ii++ )
            x_star[ii] = x[ii];
      }
      Find_Next( x, SET - 1 );

   }

}


double S_variance( const int x_hat[], int x_gap[], const int m_k, int start, double W2[][SET] )
{
   double sum = 0;
   for ( int jj = start; jj < start + m_k; jj++ )
   {
      sum = sum + pow( f_scenario( x_hat, W2, jj ) - f_scenario( x_gap, W2, jj ) - f_bar_n( x_hat, W2, m_k, start ) + f_bar_n( x_gap, W2, m_k, start ) , 2 );
   }

   return double( sum / ( m_k - 1 ) );

}


double Gap_estimate( const int x_hat[], const int n_k, double& s_std, double W2[][SET] )
{
   static int x_gap[SET];                  // Stores solution of gap problem
   static double z_gap;                    // objective value of gap problem

   int start = 0;
   int m = COUNT / RP;
   double mu_hat = 0.0;
   s_std = 0.0;

   for( int ii = 0; ii < RP; ii++ )
   {
      start = ii * m;
      Solve_SIP( n_k, start, W2, x_gap, z_gap );
      mu_hat = mu_hat + z_gap - f_bar_n( x_hat, W2, n_k, start );
      s_std = s_std + S_variance( x_hat, x_gap, n_k, start, W2 );
   }

   mu_hat = double( mu_hat / RP );
   s_std = sqrt( double( s_std / RP ) );

   cout << "Test Output: Gap Estimate = " << mu_hat << endl;

   return mu_hat;
}
