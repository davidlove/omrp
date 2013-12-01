#include "defs.h"
#include "sip_solver.h"
#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes

/*============================================================================*/
/*============================================================================*/
/*============================================================================*/
/*=========================== Main Function ==================================*/
int main()
{
   /* Defining/Declaring Necessary Parameters*/

   time_t seconds1, seconds2, total_time;     // Keep Track of Time



   int k;                          // Iterations
   int n_k;                        // SAA sample size
   int m_k;                        // Gap Problem sample size
   double mu_hat;                  // Gap for each replication
   double gap;
   double CI[REPLICAT];
   int ref = 0;                    // Number of times we refresh
   // int FQ_Estimate1 = -1;           // First quadratic estimate
   double s_std;

   int x_star[SET];                 // True Solution
   double z_star;                   // True optimal objective

   double W[COUNT][SET];            // Stores random weights
   double W2[COUNT][SET];           // Stores random weights

   int x_hat[SET];                 // Cand. Solutions
   double z_hat;                   // Obj. Values

   int x_temp[SET];                // Temp Variables
   double z_temp;

   int Jumps;
   double IT1_GAPESTIMATE;         // Gap estimate at the first iteration
   double IT1SN;                   // sample standard deviation at the first iteration
   double ITSNE_1;                   // sample standard deviation at iteration before last
   double MUHATKE_1;                 // Gap estimate at iteration before last



   /*Save Output to File*/
   ofstream OutFile;
   OutFile.open( "fsp_amrp_log.txt" );

   /* OutFile << setw( 3 ) << "REP" << setw( 6 ) << "PROCED." << setw( 9 ) << "STOP_RULE" << setw( 5 ) << "N_0" << setw( 5 ) << "M_0" << setw( 5 ) << "d" << setw( 6 ) << "NE"
           << setw( 6 ) << "ME" << setw( 5 ) << "KE" << setw( 15 ) << "GAP" << setw( 15 ) << "GAP_ESTIMATE" << setw( 11 ) << "SOLUTION" << setw( 5 ) << "JUMPS"
           << setw( 8 ) << "DURATION" << setw( 5 ) << "KFN" << setw( 5 ) << "KFM" << setw( 5 ) << "TOTAL_REPLIC" << setw( 15 ) << "ZSTAR" << setw( 15 ) << "ZNSTAR"
           << setw( 15 ) << "RUN_ID" << setw( 7 ) << "REFRESH" << setw( 8 ) << "FQ_Estimate1" << setw( 15 ) << "IT1_GAPESTIMATE" << setw( 15 ) << "IT1SN"
           << setw( 15 ) << "ITSNE_1" << setw( 15 ) << "MUHATKE_1" << setw( 15 ) << "ITSNE_1" << endl << endl; */


   Optimal( x_star, z_star );                // True Optimals
   cout << "Test OutPut ; z* = " << z_star << endl;
   cout << "X Star = ";
   for ( int ii = 0 ; ii < SET; ii++ )
      cout << x_star[ii];
   cout << endl;

   int r_seed = 1;

   for( int ii = 0; ii < Nk; ii++ )
   {
      for( int uu = 0; uu < REPLICAT; uu++ )
      {
         ref = 0;                               // # of times resampling takes palce
         seconds1 = time ( NULL );              // Start Time of replication
         int seed = 1234567890 + 12345678 * uu; //int32 seed = time(0);                // random seed
         StochasticLib1 sto( seed );            // make instance of random library

         for ( int ii2 = 0; ii2 < SET; ii2++ )        // Generating Random Variates
         {
            for ( int jj = 0; jj < COUNT; jj++ )
            {
               W[jj][ii2] = sto.Normal( mu[ii2], sigma[ii2] );
               W2[jj][ii2] = sto.Normal( mu[ii2], sigma[ii2] );
            }
         }

         n_k = n[ii];
         m_k = m[ii];
         k = 1;
         Jumps = 0;

         Solve_SIP( m_k, 0, W, x_hat, z_temp );
         mu_hat = Gap_estimate( x_hat, n_k, s_std, W2 );

         int jj = 0;
         while ( jj < Dk )
         {
            if ( mu_hat + ( 1 + Z_ALPHA * s_std ) / sqrt( double( n_k ) ) <= d[jj]*z_star )
            {
               gap = z_star - obj_func( x_hat );
               seconds2 = time ( NULL );                     // Stopping Time
               total_time = seconds2 - seconds1;             // Total Time

               if( k == 1 )
               {
                  IT1_GAPESTIMATE = mu_hat;
                  IT1SN = s_std;
                  ITSNE_1 = s_std;
                  MUHATKE_1 = mu_hat;
               }

               /* OutFile << setw( 3 ) << uu << setw( 6 ) << "FSP" << setw( 9 ) << "A2RP" << setw( 5 ) << n[ii] << setw( 5 ) << m[ii] << setw( 5 ) << d[jj]
                       << setw( 6 ) << n_k << setw( 6 ) << m_k << setw( 5 ) << k << setw( 15 ) << gap << setw( 15 ) << mu_hat << setw( 11 );
               for ( int t = 0; t < SET; t++ )
                  OutFile << x_hat[t];
               OutFile << setw( 5 ) << Jumps << setw( 8 ) << total_time << setw( 5 ) << KFN << setw( 5 ) << KFM << setw( 5 ) << REPLICAT << setw( 15 ) << z_star
                       << setw( 15 ) << z_gap << setw( 15 ) << "FSP_A2RP_082608" << setw( 7 ) << ref << setw( 8 ) << FQ_Estimate1 << setw( 15 ) << IT1_GAPESTIMATE
                       << setw( 15 ) << IT1SN << setw( 15 ) << ITSNE_1 << setw( 15 ) << MUHATKE_1 << setw( 15 ) << s_std << endl; */

               jj++;
            }
            else
            {
               n_k = n_k + 1;  // Sample size increase
               m_k = m_k + 1;
               k++;

               ITSNE_1 = s_std;
               MUHATKE_1 = mu_hat;

               if( k == 2 )
               {
                  IT1_GAPESTIMATE = mu_hat;
                  IT1SN = s_std;
               }

               /*Refreshing*/
               if ( ( k % KFN ) == 0 )
               {
                  ref++;
                  seed = 1234567890 + 1313 * k + 12345678 * uu;
                  StochasticLib1 sto( seed );
                  for ( int ii = 0; ii < SET; ii++ )
                  {
                     for ( int jj = 0; jj < COUNT; jj++ )
                     {
                        W2[jj][ii] = sto.Normal( mu[ii], sigma[ii] );
                     }
                  }
               }

               Solve_SIP( m_k, 0, W, x_hat, z_temp );
               mu_hat = Gap_estimate( x_hat, n_k, s_std, W2 );
            }


         }

      }
   }

   OutFile.close();

   cout << "===================The End===========================" << endl;

   getchar();

   return 0;
}
/*============================= The End ======================================*/


