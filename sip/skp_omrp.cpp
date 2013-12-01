#include "defs.h"
#include "sip_solver.h"
#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes

#define DEBUG_DEMAND 0  // True to debug demand generation
#define DEBUG_GAP    0  // True to debug gap-estimation 
#define DEBUG_COVER  0  // True to debug coverage estimation
#define DEBUG_OBJ    0  // True to debug optimality gap determination
#define DEBUG        DEBUG_DEMAND || DEBUG_GAP || DEBUG_COVER || DEBUG_OBJ

// Find the coverage probability of the gap estimates
double FindCoverageProbability( const double* gapest, const double* gbar, const int rep, const double ng_nonoverlap, const int ng, const double za, const double OptGap, const int degreesFreedom, ofstream& debugFile, double* var = NULL, double* gapbar = NULL );

//calculates plain old sample variance 
double CalculateVariance ( const double* y1, double ybar1, int size, int repCounter, int degreesFreedom ); 

int main( int argc, char** argv )
{
   //general variables

   int ii,repCounter,batchCounter,gammaCounter;           
   clock_t start, finish; 
   double runtime; 
   ofstream outFile,
            debugFile;

   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // Get important parameter values from command line input
   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   char           optionCharacter;
   //extern char*   optarg;
   extern int     optind;
                  //optopt;

   int            inputN = 0,
                  inputM = 0,
                  inputK = 0,
                  inputG = 0,
                  inputR = 1000,
                  nkFlag = 0;
   string         outputFileName;

   while( ( optionCharacter = getopt( argc, argv, "n:m:k:g:r:o:" ) ) != -1 )
   {
      switch( optionCharacter )
      {
         case 'n':
            inputN = atoi( argv[optind-1] );
            nkFlag++;
            break;
         case 'm':
            inputM = atoi( argv[optind-1] );
            break;
         case 'k':
            inputK = atoi( argv[optind-1] );
            nkFlag++;
            break;
         case 'g':
            inputG = atoi( argv[optind-1] );
            break;
         case 'r':
            inputR = atoi( argv[optind-1] );
            break;
         case 'o':
            outputFileName = argv[optind-1];
            break;
         case ':':
         case '?':
            exit(1);
      }
   }

   // If g is not specified, use nonoverlapping batches
   if( !inputG )
   {
      inputG = inputM;
   }

   // If string is not specified, use a default
   if( outputFileName.empty() )
   {
      outputFileName = "10D_OMRP_soln.txt";
   }
   // Only one of n, k can be specified.  Enforce this
   if( !nkFlag )
   {
      inputK = 30;
      inputN = inputM * inputK;
   }
   else if( nkFlag > 1 )
   {
      cerr << "Can only specify one of n, k" << endl;
      exit(1);
   }
   else if( !inputN )
   {
      inputN = inputM * inputK;
   }
   else
   {
      inputK = inputN / inputM;
   }

   // m must be specified
   if( !inputM )
   {
      cerr << "Batch size m must be specified" << endl;
      exit(1);
   }
   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------
   // End of parsing the input
   // ------------------------------------------------------------------------
   // ------------------------------------------------------------------------

   //parameters for multiple replications procedure
   
   // Need new parameters: Batch Size (m), Sample Size (n), Overlap Parameter (gamma).
   const int      batchSize = inputM,  //sample size within each batch, originally 1000
                  gamma = inputG,      //number of entries until overlap
                  numBatches = floor( ( (double) inputN - batchSize ) / gamma + 1 ),
                                       //replication size 30
                  sampleSize = batchSize + (numBatches - 1) * gamma, 
                                       //total number of sample drawn per confidence interval, originally 30000
                  degreesFreedom = numBatches,
                                       // Old one
                                       //numBatches - ( mPrime / gamma ),
                                       // degrees of freedom of variance operator
                  numCIs = inputR;     //10000  -- formerly rep
   const double   numNonOverBatches = (double) sampleSize / (double) batchSize - 1.0;
                                       // number of nonoverlapping batches
   
   // This needs to change for overlapping batch means.  Stores the total unumber of gap estimates
   int            totssize = numCIs * numBatches;  
                                       //        -- just to store everything
   
   static double  za = 1.282;           //za=1.645 -> 95% quantile of a st normal
                                        //za=1.282 -> 90% quantile of a st normal

   //parameters related to solution 

   int            xhat[SET] = {0,0,1,1,1,0,1,1,0,1},
                                       //candidate solution
                  xnstar[SET] = {0};   //solution to n-scenario problem   
   double         gbar[numCIs],        //final gap estimate
                  W[sampleSize][SET];  //Stores random weights


   double         OptGap;              //stores the true optimality gap

   //parameters for calculation of gap and sample variance

   double         coverage = 0.0; 

   //other arrays

   double *       gapest,              //gap estimate of single numCIs procedure in each replication
          *       indivGaps;           //gap esimates for individual samples
   int    *       indivCounts;

   // Variables for calculating MSE

   double         var[numCIs];

   for( ii = 0; ii < numCIs; ii++ )
   {
      var[ii] = 0.0;
   }
   
   //=============================================================================================
   
#if DEBUG
   debugFile.open( "debug_10D.dat", ios::out );
   if ( !debugFile )
   {
      cout<<"error opening debugging file"<<endl;
      abort(); 
   }
#endif
   outFile.open( outputFileName.c_str(), ios::out );
   if ( !outFile )
   {
      cout<<"error opening output file"<<endl;
      abort(); 
   }

   start = clock();
   
   //first, some initializations

   // Find true optimality gap
   int xstar[SET];
   Optimal( xstar, OptGap );
   OptGap = OptGap - obj_func( xhat );

   int seed = 1234567890 + 12345678; //int seed = time(0);                // random seed
   StochasticLib1 sto( seed );            // make instance of random library

   gapest  = new double [totssize]; 
   indivGaps = new double[sampleSize];
   indivCounts = new int[sampleSize];

   for( repCounter = 0; repCounter < numCIs; repCounter++ )
   {
      for( ii = 0; ii < sampleSize; ii++ )
      {
         indivGaps[ii] = 0.0;
         indivCounts[ii] = 0;
      }
      // Create sample for this confidence interval
      for( ii = 0; ii < SET; ii++ )
         for( int jj = 0; jj < sampleSize; jj++ )
         {
            W[jj][ii] = sto.Normal( mu[ii], sigma[ii] );
         }

      for( batchCounter = 0; batchCounter < numBatches; batchCounter++ )
      { 
         // MODIFY -- Code doesn't work because gapest only accounts for xnstar, not xhat.
         Solve_SIP( batchSize, batchCounter*gamma, W, xnstar, gapest[(repCounter*numBatches)+batchCounter] );
         gapest[(repCounter*numBatches)+batchCounter] -= f_bar_n( xhat, W, batchSize, batchCounter*gamma );
         // Loop for calculating Gbar
         for( gammaCounter = 0; gammaCounter < batchSize; gammaCounter++ )
         {
            indivGaps[batchCounter*gamma + gammaCounter] += f_scenario( xnstar, W, batchCounter*gamma + gammaCounter ) - f_scenario( xhat, W, batchCounter*gamma+ gammaCounter );
            indivCounts[batchCounter*gamma + gammaCounter] += 1;
         }

      }  //end of loop with batchCounter (numBatches) 

      // Divide by the "vertical" component to get the true gap estimate
      // and finish by averaging
      gbar[repCounter] = 0.0;
      for( ii = 0; ii < sampleSize; ii++ )
      {
         indivGaps[ii] /= indivCounts[ii];
         gbar[repCounter] += indivGaps[ii];
      }
      gbar[repCounter] /= sampleSize;

   }  //end of for loop with numCIs

   double   avgGap = 0.0,
            bias = 0.0,
            avgVar = 0.0,
            varVar = 0.0;
   
   coverage = FindCoverageProbability( gapest, gbar, numCIs, numNonOverBatches, numBatches, za, OptGap, degreesFreedom, debugFile, var, &avgGap );
   
   bias = avgGap - OptGap;
   for( ii = 0, avgVar = 0; ii < numCIs; ii++ )
      avgVar += var[ii];
   avgVar /= numCIs;
   for( ii = 0, varVar = 0.0; ii < numCIs; ii++ )
      varVar += pow( var[ii] - avgVar, 2 );
   varVar /= (numCIs - 1);

   double widthCI;

   widthCI = za * sqrt( avgVar ) / sqrt( numNonOverBatches );
   
   finish = clock(); 

   //OUTPUT results:

   //outFile << "Candidate Solution: XHAT  =  " << xhat << endl; 
   outFile << "Candidate Solution: xhat  =  ";
   for( int ii = 0; ii < SET; ii++ )
      outFile << xhat[ii] << ' ';
   outFile << endl;
   outFile << "Optimal Solution: xstar   =  ";
   for( int ii = 0; ii < SET; ii++ )
      outFile << xstar[ii] << ' ';
   outFile << endl;
   outFile << "True Optimality Gap:      =  " << OptGap << endl;
   outFile << "Sample Size       :    n  =  " << sampleSize << endl;
   outFile << "Batch Size        :    m  =  " << batchSize << endl; 
   outFile << "Batch Non-overlap : gamma =  " << gamma << endl;
   outFile << "Number Batches    :       =  " << numBatches << endl; 
   outFile << "Number non-overlapping k  =  " << numNonOverBatches << endl;
   outFile << "degreesFreedom            =  " << degreesFreedom << endl;
   outFile << "Number of replications    =  " << numCIs << endl; 
   outFile << "Width of CI               =  " << widthCI << endl;
   outFile << "Bias                      =  " << bias << endl;
   outFile << "Average Variance          =  " << avgVar << endl;
   outFile << "Variance of VG            =  " << varVar << endl;
   //outFile << "Prob of Solution Overlap  =  " << (double) noChange / (double) (numBatches - 1) / (double) (numCIs) << endl;
   // outFile << "Bias = " << bias << ", avgVar = " << avgVar << endl;
   // outFile << "MSE = " << pow(bias,2) + avgVar << endl;

   outFile << endl<< coverage << "  +-  " << sqrt(coverage*(1-coverage)/numCIs)* za << endl; 

   runtime = double (finish - start)/double (CLOCKS_PER_SEC); 
   outFile << endl << "Time to run the program is :" << runtime << " seconds." << endl; 

   //outFile << endl << "Gap Estimates" << endl;
   //for( ii = 0; ii < totssize; ii++ )
      //outFile << gapest[ii] << endl;

   outFile.close();
   
   //it a good idea to clean up

   //if (demand)  delete [] demand; 
   if (gapest)  delete [] gapest; 
   if( indivGaps ) delete [] indivGaps; 
   if( indivCounts ) delete [] indivCounts; 

   //demand  = NULL; 
   gapest  = NULL; 

   return 0;
 } //end of main

 
 double FindCoverageProbability( const double* gapest, const double* gbar, const int numCIs, const double numNonOverBatches, const int numBatches, const double za, const double OptGap, const int degreesFreedom, ofstream& debugFile, double* var, double* gapbar )
{
   int      repCounter;
            //ii;
   double   //gapsum = 0,
            //gbar = 0,
            svar = 0,
            coverage = 0;
   double //* var,
          * CI;
   int    * check;
   bool     varPassed;     // true if var was passed in as an argument

   // gapbar  = new double [numCIs]; 
   if( !var )
   {
      varPassed = false;
      var  = new double [numCIs]; 
   }
   else
      varPassed = true;
   CI      = new double [numCIs]; 
   check   = new int    [numCIs]; 

   coverage = 0;
   for( repCounter = 0; repCounter < numCIs; repCounter++ )
   {
      //calculate gapbar (average gap of numBatches replications)

      // This is the average of the overlapping gaps, which doesn't put enough emphasis on samples near either end of the distribution.
      // This needs to be changed somehow...  Maybe an average of non-overlapping batches?  Enforce n = km?
      //gapsum = 0.0; 
      //for( ii = (repCounter*numBatches); ii < (repCounter*numBatches)+numBatches; ii++ )
         //gapsum += gapest[ii]; 
      
      //gbar = (double) gapsum / (double) numBatches; 
      //cout << "gbar = " << gbar << endl;
      //cout << gbar << endl;
      // gapbar[repCounter]  = gbar; 

      //calculate variance of the gap in numBatches replications 

      svar   = CalculateVariance( gapest, gbar[repCounter], numBatches, repCounter, degreesFreedom); 
      var[repCounter] = svar; 


      //calculate CI and fill in check, update coverage
      
      CI[repCounter] = gbar[repCounter] + za * sqrt(svar) / sqrt(numNonOverBatches); 

      
      if( CI[repCounter] < OptGap )
         check[repCounter] = 0; 
      else 
         check[repCounter] = 1; 

      coverage += check[repCounter]; 
   }
   coverage /= numCIs; 

   // return the average gap estimate, if it was passed in
   //if( gapbar )
      //*gapbar = gbar;


#if DEBUG_COVER
   int totssize = numBatches * numCIs;
   debugFile<<endl<<"Gap Estimate"<<endl; 
   debugFile<<"------------"<<endl; 
   for(repCounter=0; repCounter<totssize; repCounter++)
      debugFile<<gapest[repCounter]<<endl; 

   // debugFile<<endl<<"ZnStar"<<endl; 
   // debugFile<<"---------"<<endl; 
   // for(repCounter=0; repCounter<totssize; repCounter++)
      // debugFile<<zn[repCounter]<<endl; 

   // debugFile<<endl<<"ExpCost Estimate"<<endl; 
   // debugFile<<"------------------"<<endl; 
   // for(repCounter=0; repCounter<totssize; repCounter++)
      // debugFile<<expCost[repCounter]<<endl; 


   debugFile<<endl; 


   // debugFile<<endl<<"Gap_Bar"<<endl; 
   // debugFile<<"------------"<<endl; 
   // for(repCounter=0; repCounter<numCIs; repCounter++)
      // debugFile<<gapbar[repCounter]<<endl; 

   debugFile<<endl<<"Var estimate"<<endl; 
   debugFile<<"---------------"<<endl; 
   for(repCounter=0; repCounter<numCIs; repCounter++)
      debugFile<<var[repCounter]<<endl; 

   debugFile<<endl<<"CI Estimate"<<endl; 
   debugFile<<"------------------"<<endl; 
   for(repCounter=0; repCounter<numCIs; repCounter++)
      debugFile<<CI[repCounter]<<endl; 
   
   debugFile<<endl<<"Checks of CI"<<endl; 
   debugFile<<"------------------"<<endl; 
   for(repCounter=0; repCounter<numCIs; repCounter++)
      debugFile<<check[repCounter]<<endl; 
   debugFile.close();
#endif

   // if (gapbar)  delete [] gapbar; 
   if( !varPassed )
      if (var)     delete [] var; 
   if (CI)      delete [] CI; 
   if (check)   delete [] check; 
   // gapbar  = NULL; 
   if( !varPassed )
      var     = NULL; 
   CI      = NULL;
   check   = NULL;

   return coverage;
} // end of FindCoverageProbability


double CalculateVariance ( const double* y1, double ybar1, int size, int currRep, int degreesFreedom )
{
   int ii; 
   double sv1, sumsquare; 
   
   sumsquare = 0.0; 

   for( ii = 0; ii < size; ii++ )
   {
      // Change currRed*size to reflect overlapping batches.  I think currRep*gamma + ii should do it.
      sv1 = y1[(currRep*size)+ii] - ybar1; 
      sumsquare += pow(sv1,2);
   }

   sumsquare *= ( 1.0 / (double) degreesFreedom);

   return sumsquare;    
} // end of CalculateVariance 


