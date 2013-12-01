/***********************************************************************
   
    Multiple Replications procedure applied to newsvendor problem 

************************************************************************/

#include "multiplerep.h"

int main ( int argc, char** argv )
{

#if DEBUG
   srand(49873); 
#else
   srand(16326);
   //srand( time(NULL) );
   //srand( 49276 );
#endif
   
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
                  inputR = 10000,
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
      outputFileName = "Newsvendor_MRP_soln.txt";
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

   //parameters for the newsvendor problem

   const int      r = 15,               //price
                  c = 5;                //cost
   const double   a = 0,                //Uniform(a,b) parameters
                  b = 10;
// const double   mean = 10,            //mean of normal distribution
//                stdev = 5,            //stdev of normal distribution
   
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

   double         xhat,                //candidate solution
                  xnstar[numBatches],  //solution to n-scenario problem
                  gbar[numCIs],        //final gap estimate
                  quantile;            //quantile where optimum resides
   int            discreteQuantile;    //location of optimum solution
   double         here;                //to find location of optimum solution


   double         OptGap;              //stores the true optimality gap

   //parameters for calculation of gap and sample variance

   double         coverage = 0.0; 

   //other arrays

   double *       gapest,              //gap estimate of single numCIs procedure in each replication
          *       indivGaps;           //gap esimates for individual samples
   int    *       indivCounts;

   circ_sort      demand;

   // Variables for calculating MSE

   double         var[numCIs];

   for( ii = 0; ii < numCIs; ii++ )
   {
      var[ii] = 0.0;
   }
   
   //=============================================================================================
   
#if DEBUG
   debugFile.open( "debug_Newsvendor.dat", ios::out );
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

   //OptGap = 3.333;   
   //OptGap = 47.619;   
   //xhat   = 6.5121; 
   xhat   = 8.7749; 
   //xhat = b - (double) c * (b-a)/ r;
   OptGap = FindObjValue( xhat, r, c, a, b, debugFile );

   // prob = (double) (1.00/batchSize); 
   
   demand.Initialize( batchSize );
   gapest  = new double [totssize]; 
   indivGaps = new double[sampleSize];
   indivCounts = new int[sampleSize];

   //find the location of the optimum solution 
   //in the sorted demand of the sampling problem

   quantile = (double)(r - c)/ (double) r; 
   for( ii = 1; ii <= batchSize; ii++ )
   {
      here = (double) ii/batchSize; 
      if( here >= quantile )
      {
         break;   
      }
   }
   discreteQuantile = ii;   //this is where optimum resides


   // -------------------------------
   // This code recordes the number of time the optimal solutions overlap as the random numbers are changed.
   double oldOptSoln = 0.0;
   int noChange = 0;
   // -------------------------------


   // batchCounter and repCounter do nothing indpendent in this function.  They may need to be combined together for overlapping means
   for( repCounter = 0; repCounter < numCIs; repCounter++ )
   {
      for( ii = 0; ii < sampleSize; ii++ )
      {
         indivGaps[ii] = 0.0;
         indivCounts[ii] = 0;
      }
      // We want the demands for each CI to be completely independent
      GenerateDemand( demand, batchSize, a, b, debugFile );
      for( batchCounter = 0; batchCounter < numBatches; batchCounter++ )
      { 
         // Was: GenerateDemand( demand, batchSize, a, b, debugFile );
         xnstar[batchCounter] = demand[discreteQuantile-1]; 
         gapest[(repCounter*numBatches)+batchCounter] = FindGapEstimate( demand, batchSize, discreteQuantile, xnstar[batchCounter], xhat, r, c, debugFile );
         // Loop for the new Gbar setting
         for( gammaCounter = 0; gammaCounter < batchSize; gammaCounter++ )
         {
            indivGaps[batchCounter*gamma + gammaCounter] += FindSampleProfit(xnstar[batchCounter],demand(gammaCounter),r,c,debugFile) - FindSampleProfit(xhat,demand(gammaCounter),r,c,debugFile);
            indivCounts[batchCounter*gamma + gammaCounter] += 1;
         }
         GenerateDemand( demand, gamma, a, b, debugFile );

         // --------------------------
         // More code to test how often optimal solutions are changing as we change random numbers
         if( xnstar[batchCounter] == oldOptSoln )
            noChange++;
         oldOptSoln = xnstar[batchCounter];
         // --------------------------

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

      // GenerateDemand( demand, batchSize, a, b, debugFile );
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

   outFile << "Candidate Solution: XHAT  =  " << xhat << endl; 
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
   outFile << "Prob of Solution Overlap  =  " << (double) noChange / (double) (numBatches - 1) / (double) (numCIs) << endl;
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

// void GenerateDemand( double* demand, const int batchSize, const double lowerLimit, const double upperLimit, ofstream& debugFile )
void GenerateDemand( circ_sort& demand, const int numNewData, const double lowerLimit, const double upperLimit, ofstream& debugFile )
{
   int ii;

   for( ii = 0; ii < numNewData; ii++ )
   {
      //demand[ii] = TruncNormal(mean, stdev); 
      //demand[ii] = Uniform( lowerLimit, upperLimit ); 
      demand.AddEntry( Uniform( lowerLimit, upperLimit ) );
   } 

#if DEBUG_DEMAND
   debugFile << endl << "Unsorted Demand" << endl << "-----------" << endl;
   demand.PrintQueue( debugFile );
   // for( ii = 0; ii < batchSize; ii++ )
   // {
      // debugFile<<"d["<<ii<<"] = "<<demand[ii]<<endl; 
   // }
#endif

   //sort demand and find optimum solution
   //qsort(demand, (size_t) batchSize, sizeof(double), Compare); 
   demand.Sort();

#if DEBUG_DEMAND
   debugFile<<endl<<"Sorted Demand"; 
   debugFile<<endl<<"---------------"<<endl; 
   demand.PrintSorted( debugFile );
   //for( ii=0; ii < numNewData; ii++ )
   //{
      //debugFile<<"d["<<ii<<"] = "<<demand[ii]<<endl; 
   //}
   //debugFile << endl;
#endif

   return;
} // end of GenerateDemand

double FindGapEstimate( circ_sort& demand, const int batchSize, const int discreteQuantile, const double xnstar, const double xhat, const int r, const int c, ofstream& debugFile )
{
   int    ii;
   double optimalSum = 0,
          candidateSum = 0,
          znstar,
          obj;

   for( ii = 0; ii < discreteQuantile; ii++ )
   {
      optimalSum += demand[ii];
   }
   for( ii = discreteQuantile; ii < batchSize; ii++ )
   {
      optimalSum += demand[discreteQuantile-1];
   }

   znstar = -c*xnstar + r*optimalSum/batchSize; 
   // zn[(repCounter*numBatches)+batchCounter]  = znstar;

   for( ii = 0; ii <batchSize; ii++ )
   {
      if( demand[ii] <= xhat )
         candidateSum += demand[ii];
      else 
         candidateSum += xhat; 
   }
                     
   obj        = -c * xhat  +  r * candidateSum / batchSize;
   // expCost[(repCounter*numBatches)+batchCounter] = obj; 

#if DEBUG_GAP
   debugFile << endl << "discreteQuantile  = "        << discreteQuantile << endl; 
   debugFile << endl << "xnstar = "                   << xnstar << endl;
   debugFile << endl << "Min of Demand and xhat : "   << endl; 
   debugFile << endl << "optimalSum = "               << optimalSum; 
   debugFile << endl << "candidateSum = "             << candidateSum << endl;
   // debugFile << endl << "xhat objective value = "     << gap; 
   debugFile << endl << "              znstar = "     << znstar << endl;
   debugFile << endl << "Gap Estimate =         "     << znstar - obj << endl;
#endif

   return znstar - obj; 
} // end of FindGapEstimate

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


double FindObjValue( double xhat, double price, double cost, double a, double b, ofstream& debugFile )
{
   double OptimalProfit,
          SubOptimalProfit;

   OptimalProfit = FindExpectedProfit( b - cost * (b-a)/price, price, cost, a, b, debugFile );
   SubOptimalProfit = FindExpectedProfit( xhat, price, cost, a, b, debugFile );

#if DEBUG_OBJ
   debugFile << "Optimal Solution = " << b - cost * (b-a)/price << endl;
   debugFile << "Optimal Profit = " << OptimalProfit << endl;
   debugFile << "SubOptimal Profilt = " << SubOptimalProfit << endl;
#endif

   return OptimalProfit - SubOptimalProfit;
} // end of FindObjValue

double FindExpectedProfit( double xhat, double price, double cost, double a, double b, ofstream& debugfile )
{
   return price / (b-a) * (b*xhat - a*a - xhat*xhat/2) - cost * xhat;
} // end of FindExpectedProfit

double FindSampleProfit( double xhat, double xi, double price, double cost, ofstream& debugfile )
{
   if( xhat >= xi )
      return price * xi - cost * xhat;
   else
      return (price - cost) * xhat;
}

int Compare( const void *arg1, const void *arg2 )
{
   if(*(double *) arg1 < *(double *) arg2){
      return -1; 
   }
   else if ((*(double *) arg1 > *(double *) arg2)){
      return +1; 
   }
   else {
      return 0; 
   }
} // end of Compare

double TruncNormal ( double mean1, double stdev1)
{
   //normal variate generation using polar method
   //pp.491 of Law & Kelton, 2ed
   //generate truncated normal by rejecting the ones less than 0
   

   //Note: Not the best efficient way to use this as we are
   //      generating two truncated normal variates but using
   //      only one of them. 

   double u1, u2, v1, v2, w, y; 
   double norm1, norm2; 
   bool flag1, flag2; 

   norm1 = 0;
   norm2 = 0;

   flag2 = true; 
   
   while(flag2){

      flag1 = true; 

      while (flag1){ 
         u1 = (double) (rand()/ (double) RAND_MAX); 
         u2 = (double) (rand()/ (double) RAND_MAX);
         v1 = 2*u1 - 1;
         v2 = 2*u2 - 1; 
         w = v1*v1 + v2*v2; 
         if (w<=1){
            y = sqrt (-2*log(w)/w); 
            norm1 = v1 * y; 
            norm2 = v2 * y;
            norm1 = mean1 + stdev1 * norm1; 
            norm2 = mean1 + stdev1 * norm2; 
            flag1 = false; 
         }
      } //end of while with flag1

      if (norm1 >=0){
         return norm1; 
         flag2 = false; 
      }
      else 
         if (norm2 >=0){
            return norm2; 
            flag2 = false; 
         }

   }  //end of while with flag2

   return 0.0; // Function had no return value when I got it -- David Love, March 16, 2009
   
} // end of TruncNormal 

double Uniform ( double a, double b)
{
   double u1; 

   u1 = (double) (rand()/ (double) RAND_MAX); 
   u1 = a + (b-a)*u1; 
   return u1; 
} // end of Uniform 
   

double UnusedCalcVar ( double *y1, double *y2, double sum1, double sum2, int iss, int r )
{
   double sv1, sv2, ybar1, ybar2; 
   double sumsquare; 
   
   ybar1 = sum1 / (double) iss; 
   ybar2 = sum2 / (double) iss; 

   sumsquare = 0.0; 

   for(int ii = 0; ii<iss; ii++){
      sv1 = y1[ii] - ybar1; 
      sv2 = y2[ii] - ybar2; 
      sumsquare += pow((sv1-sv2),2);
   }

   sumsquare *= ( 1.0 / (double) (iss - 1));
   sumsquare *= pow((double) r,2); 

   return sumsquare; 
} // end of UnusedCalcVar 

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

