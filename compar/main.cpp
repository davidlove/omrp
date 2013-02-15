/*------------------------------------------------------------------------------
MODULE TYPE:		The main module.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski
--------------------------------------------------------------------------------

SOURCE FILE NAME:	main.cpp
CREATED:			1994.08.17
LAST MODIFIED:		1997.02.21

DEPENDENCIES:		error.h, vec_pool.h, smartptr.h, memblock.h, print.h,
					smartdcl.h, parsemps.h, stdtype.h, mps_lp.h, simplex.h, 
					sort_lab.h, std_tmpl.h, myalloc.h, smplx_lp.h, solv_lp.h, 
					compile.h, lp_codes.h, history.h, parsespc.h, inverse.h,
					scenario.h, linklist.h, subprobl.h, solver.h, solvcode.h,
					mast_sol.h, qdx_pub.h, lexer.h, work_vec.h, determlp.h,
					stochsol.h
					<stdio.h>, <stdlib.h>, <ctype.h>
--------------------------------------------------------------------------------
SOURCE FILE CONTENTS:
	x
--------------------------------------------------------------------------------
STATIC FUNCTIONS:
	x

STATIC DATA:
	x
--------------------------------------------------------------------------------
USED MACROS AND THEIR MEANING:
	XXXXX			- x
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Later... Modified by Guzin Bayraksan... 
The modifications (mostly) in between //@BEGIN--- and //@END----- 

Even Later... Modified by David Love...
Deleted most of the code.  Trying to make it all simpler.
--------------------------------------------------------------------------------*/
#include "main.h"

#define MAKE_TRD 0

// DESPERATION -- Need some new global variables that will hold scenarios and shit
vector<double> allScenarios,
	       allWeights;
double*	 batchScenarios;


/*------------------------------------------------------------------------------

	int main( int argc, char *argv[] )

PURPOSE:
	Calls the actual main procedure of the program. After it exits checks if
	all memory has been released properly.

PARAMETERS:
	int argc, char *argv[] 
		The command line arguments.

RETURN VALUE:
	Same as returned from the 'run()' function.

SIDE EFFECTS:
	None. No one would notice anyways, unless the program dumps core.

------------------------------------------------------------------------------*/

int main( int argc, char *argv[] )
{
	run( argc, argv );

	return 0;
}  //end of main


/*------------------------------------------------------------------------------

	void run( int argc, char *argv[] )

PURPOSE:
	The actual 'main()' function of the decomposition code. Calls argument
	parser, displays help if argument error is detected. Creates the necessary
	objects. Reads in the input files. Runs the solver. Reports results.

PARAMETERS:
	int argc, char *argv[] 
		The command line arguments.

RETURN VALUE:
	Zero on success, non-zero otherwise.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void run( int argc, char *argv[] )
{
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
                  inputR = 300,
                  nkFlag = 0;
   string	  scenFileName;

   while( ( optionCharacter = getopt( argc, argv, "n:m:k:g:r:f:" ) ) != -1 )
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
	 case 'f':
	    scenFileName = argv[optind-1];
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

   // If no scenario file name is given, quit
   if( scenFileName.empty() )
   {
      cerr << "Must specify a file for scenario information, -f filename" << endl;
      exit(1);
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

   // Get the distribution
   ReadScenarios( scenFileName );
   // Set up the vector to hold the batches
   batchScenarios = new double[inputM];

//****************************************************************************************

	//THIS SECTION IS THE MRP BY ITSELF...........
	
        // David Love -- Added number of nonoverlapping batches
        const int numNonOverBatches = inputK;  // Number of nonoverlapping batches
        const int batchSize = inputM;
        // David Love -- Delete gamma when I have a working gamma input
        const int gamma = inputG;

        // David Love -- Rewrote the number of batches to work with overlapping
	const int numBatches = floor( (numNonOverBatches - 1) * (double) batchSize / (double) gamma + 1 );      //this is the smaller loop for the MRP
        // David Love -- Degrees of freedom of variance operator
        const int sampleSize = batchSize + (numBatches-1)*gamma;
        // DESPERATION -- I'm chaning the degreesFreedom and ciDenom
        const double degreesFreedom = numBatches * ( (double) sampleSize / (double) batchSize - 1.0 );
        const double ciDenom = 1.0;

        // David Love -- Counters for calculating total average gap
        double *  indivGaps;     // Gap estimates for individual samples
        int    *  indivCounts;   // Counter for "vertical" part of total gap average
        double *  gapest;        // holds all of the gap estimators
	
	int kk, oo; 
	double *var, *gbar, *ci; 

	double *znstar, *znstarbar; 

	double zn; 
	
	double optg = 1282152.974;   //cep1, for xhat=(650,...,650, 150,...,150)  
		//2.833632377;   //pgp2, for xhat=(2.5, 6, 3.5, 4.5)
		//20.01728866212;  //pgp2, for xhat=(2, 7, 2, 5)
		//38129.093677509;   //cep1, for xhat=(0,125,875,2500,0,625,1375,3000)  
		//54.32289045997;  //pgp2, for xhat=(3, 4, 4,,  4) 
                //164.84;          //apl1p, for xhat = (1111.11, 2300)
		
		//11737.350015; //cep1, for xhat=(0,0,1166.67,2500,0,500,1666.67,3000) 
		
		
	double cov = 0.0; 

	double mns, vrs;		//mean and variance for the inner loop 

        // The following two variables are not defined in the code when I got it.  Had to add them myself.
	int numCIs = inputR;//1000;                  // Added by David Love
	double za = 1.282;              // Added by David Love


	var = new double[numCIs];
	gbar = new double[numCIs];
	ci = new double[numCIs];
	
	znstarbar = new double[numCIs];

	znstar = new double[numCIs*numBatches]; 

        // David Love -- Debug by listing all gap values
        //double gapValues[numCIs][numBatches];

        // David Love -- Testing that gamma, number of batches and degrees of freedom work correctly
        // printf( " gamma = %d \n num batches = %d \n degreesFreedom = %d\n", gamma, numBatches, degreesFreedom );
        // exit(0);
        
        // David Love -- Initializers for calculating total gap average
        indivGaps = new double[sampleSize];
        indivCounts = new int[sampleSize];
        gapest = new double[numBatches];

	for(kk=0; kk<numCIs; kk++){

		mns = 0.0; 
		vrs = 0.0; 
		zn  = 0.0; 

                // David Love -- Let the screen count where we are in numCIs
#if !MAKE_TRD
                printf( "Rep %d of %d, gamma = %d, m = %d, n = %d, nb = %d, degreesFreedom = %lf\n", kk, numCIs, gamma, batchSize, sampleSize, numBatches, degreesFreedom );
#endif

                // David Love -- Initialize gap counters
                for( int ii = 0; ii < sampleSize; ii++ )
                {
                   indivGaps[ii] = 0.0;
                   indivCounts[ii] = 0;
                }

		// David Love -- Generate a full selection of scenarios for each batch
                //printf( "j -- Beginning replication %d of %d\n", kk+1, numCIs );
		ReGenerateScenarios( batchSize, batchSize ); 
		for(oo = 0; oo<numBatches; oo++){
                        // David Love -- ReGenerate back at the end of the loop.  Initialized outside the looop
                        //printf( "j -- Batch %d of %d, overlap = %d\n", oo+1, numBatches, batchSize - gamma );
			ReGenerateScenarios( batchSize, gamma ); 


                        /* // @BEGIN Rebecca's code to print the scenarios
                        // Tested and Confirmed: Overlapping is working correctly!!!
                        printf( "SCEN Printing Scenarios\n" );
                        printf( "SCEN Batch %d of %d\n", oo+1, numBatches );
                        for (int scenCounter = 0; scenCounter < batchSize; scenCounter++)
                        {
                           printf( "SCEN " );
                           printf("%d ", scenCounter+1);
                           Scen->GetScenarioComponents(scenCounter);
                           printf("\n");
                        }
                        printf("SCEN\n");
                        */ // @END Rebecca's code

                        // DESPERATION -- Commenting out a couple things to see if it compiles
			//znstar[kk*numBatches+oo] = sol->result; 

                        // DESPERATION -- Commenting out a couple things to see if it compiles
			//zn += sol->result; 
			
			gapest[oo] = CalculateGap( batchSize );

                        // David Love -- Loop for the new "G Double Bar"
                        for( int gammaCounter = 0; gammaCounter < batchSize; gammaCounter++ )
                        {
                           //indivGaps[oo*gamma + gammaCounter] += ( master.GetIndivGap(gammaCounter) + master.GetIndivGap(batchSize) ) * pow(10.0, ObjScale);
                           indivGaps[oo*gamma + gammaCounter] += batchScenarios[gammaCounter];
#if MAKE_TRD
                           cout << batchScenarios[gammaCounter] << ' ';
#endif
                           indivCounts[oo*gamma + gammaCounter] += 1;
                           //printf( "indivGaps[%d] = %lf, indivCounts[%d] = %d\n", oo*gamma + gammaCounter, indivGaps[oo*gamma + gammaCounter], 
                                 //oo*gamma + gammaCounter, indivCounts[oo*gamma + gammaCounter] );
                           //printf( "ExpC[%d] = %lf, v[%d] = %lf\n", oo*gamma + gammaCounter, master.GetExpC(gammaCounter) + master.GetExpC(batchSize),
                                 //oo*gamma + gammaCounter, master.GetV(gammaCounter) + master.GetV(batchSize) );
                        }
#if MAKE_TRD
                        cout << endl;
#endif
		} //end of for with "oo" (0 == oo < numBatches). 


		//MRP calculations:

                // David Love -- New calculation of the gap estimate
                gbar[kk] = 0.0;
                for( int ii = 0; ii < sampleSize; ii++ )
                {
                   indivGaps[ii] /= indivCounts[ii];
                   gbar[kk] += indivGaps[ii];
                }
                gbar[kk] /= sampleSize;
                //printf( "gbar[%d]   = %lf\n", kk, gbar[kk] );
		//mns /= (double) numBatches; 
		//vrs /= (double) numBatches; 

                vrs = 0.0;
                for( int ii = 0; ii < numBatches; ii++ )
                {
                   vrs += pow( gapest[ii] - gbar[kk], 2 );
                }
                vrs /= (double) degreesFreedom;

		//vrs = vrs - pow(mns, 2) ;
                // David Love -- Had to change degrees of freedom
		// vrs = vrs * numBatches / (double) (numBatches-1); 
		// vrs = vrs * numBatches / (double) degreesFreedom; 

		//gbar[kk] = mns; 
		var[kk] = vrs; 

		ci[kk] = gbar[kk] + za*sqrt(var[kk]) / sqrt(ciDenom); 
		if(ci[kk] >= optg) cov++; 

		//also calculate znstarbar

		znstarbar[kk] = zn / (double) numBatches; 
		
	}  //end of loop with kk

        // David Love -- Statistics on variance calculations
        double avgVar = 0.0;
        for( kk = 0; kk < numCIs; kk++ )
                avgVar += var[kk];
        avgVar = avgVar / numCIs;
        double varVariance = 0.0;
        for( kk = 0; kk < numCIs; kk++ )
           varVariance += pow( var[kk] - avgVar, 2 );
        varVariance = varVariance / (numCIs-1);
        // David Love -- Statistics on confidence interval width
        double ciWidth = 0.0,
               varCIWidth = 0.0;
        for( kk = 0; kk < numCIs; kk++ )
                ciWidth += za * sqrt( var[kk] ) / sqrt( ciDenom );
        ciWidth = ciWidth / numCIs;
        for( kk = 0; kk < numCIs; kk++ )
           varCIWidth += pow( (za * sqrt( var[kk] ) / sqrt( ciDenom )) - ciWidth, 2 );
        varCIWidth = varCIWidth / (numCIs-1);
        // David Love -- Statistics on total CI size
        double ciSize = 0.0,
               varCISize = 0.0;
        for( kk = 0; kk < numCIs; kk++ )
                ciSize += znstarbar[kk] + za * sqrt( var[kk] ) / sqrt( ciDenom );
        ciSize = ciSize / numCIs;
        for( kk = 0; kk < numCIs; kk++ )
           varCISize += pow( (znstarbar[kk] + za * sqrt( var[kk] ) / sqrt( ciDenom )) - ciSize, 2 );
        varCISize = varCISize / (numCIs-1);

        //for( kk = 0; kk < numCIs; kk++ )
           //ciWidth += znstarbar[kk];
        //ciWidth /= numCIs;
        //ciWidth += za * sqrt(avgVar) / sqrt( ciDenom );

	//print results:
	cov = (double) cov / numCIs;  

        // David Love -- Clearing out printed information
	// printf ("Optimality Gap = %f\n\n", optg); 
	
        // David Love -- Clearing out printed information
	// printf ("\nMRP RESULTS:\n"); 
	// printf ("-----------------\n\n"); 
	// printf ("k:         GAP:			SVAR:		CI			znstarbar:\n"); 
	// for(kk=0; kk<numCIs; kk++){
		// printf ("%d\t%f\t%20.10E\t%f\t%f\n", kk+1, gbar[kk], var[kk], ci[kk], znstarbar[kk]);
	// }

        printf("BEGIN MRP\n");
        printf("Batch Size = %d\n", batchSize);
        printf("gamma = %d\n", gamma);
        printf("CI Width = %lf\n", ciWidth);
	printf ("Coverage, MRP: = %lf\n", cov); 
        printf ("Average Variance = %lf\n", avgVar );
        printf("Variance of Variance = %le\n", varVariance );
        printf("Number Batches = %d\n", numBatches);
        printf("Degrees Freedom = %lf\n", degreesFreedom);
        printf("CI Denom = %lf\n", (double) ciDenom);
        printf("Var CI Width = %lf\n", varCIWidth);
        printf("CI Sze = %lf\n", ciSize);
        printf("Var CI Size = %lf\n", varCISize);

	//print also the znbar
        // David Love -- Clearing out printed information
	// printf("\nznbar values:\n"); 
	// printf ("-----------------\n\n"); 
	// printf ("k:         znstar:\n"); 
	// for(kk=0; kk<numBatches*numCIs; kk++){
		// printf ("%d \t %f\n", kk+1, znstar[kk]);
	// }
        // David Love -- Some debugging code
        //printf( "Vales of mn\n" );
        //for( kk = 0; kk < numCIs; kk++ ) {
           //for( oo = 0; oo < numBatches; oo++ ) {
              //printf("%lf ", gapValues[kk][oo]);
           //}
           //printf("\n");
        //}
        //printf( "Variance Information\n" );
        //for( kk = 0; kk < numCIs; kk++ ) {
           //printf( "%lf\n", var[kk] );
        //}
	

	delete var; delete gbar; delete ci;
	var = NULL; gbar = NULL; ci = NULL;

        // David Love -- Clean up memory
        delete [] indivGaps;
        delete [] indivCounts;

        delete [] batchScenarios;

//****************************************************************************************
} // end of run

void ReGenerateScenarios( int scennum, int gamma )
{
	int s, i, j, l; 
	double prob; 
        int len = 1;
        // David Love ---  Counts the number of batches that we use.
        // static int batchNumber = 0;
        // David Love -- hold seeds

	// David Love -- Sets the same random seed for every set of scenarios
        // Random01::Seed( 7625, 3293, 41);

        static MTRand::uint32 batchState[ MTRand::SAVE ];
        static bool firstTime = true;
        if( firstTime ) 
        {
           Random01::GetSeed( batchState );
           firstTime = false;
        }
	//--------------------------------------------------------------------------
	//	If needed, resize. Also, Re-assign the scenario probabilities.
	//
	// DESPERATION -- I'm not sure what this does...
	// if( scennum < ScenNum ){
	// 	ScenNum = scennum; 
	// 	this ->ArrayOfScenarios.Resize( scennum );		
	// }

	// for(s = 0; s < scennum; s++ )
	// 		ArrayOfScenarios[s]->SetProbability ( 1.0/(double)ScenNum );

	//(Asumes that the probabilities of distributions have been filled by an  
	// earlier call of GenerateScenarios...)

	//--------------------------------------------------------------------------
	//	Loop on all scenarios. Then, loop on distributions and re-form scenarios
	//
 
        // David Love -- Set up the random number seed again
        // Print( "Seed Values: %d %d %d\n", seedx, seedy, seedz );
        Random01::ReSeed( batchState );

        // David Love -- Print out the distribution locations, check that overlapping is working properly
        // Print( "Batch random variables\n" );
	for(s = 0; s < scennum; s++ )
	{
		for(i = 0; i < len; i++ )
		{
			prob = Random01::Next();
                        // David Love -- quick debugging statement for probability
                        // printf( "prob = %lf\n", prob );

			l = allScenarios.size();
			for(j = 0; j < l; j++ )
                        {
                                // David Love -- more debugging check for probability break for large values
                                // printf( "j = %d, probability = %lf\n", j ,(*dist[i])[j].GetProbability() );
				if( allWeights[j] > prob ) break;
				//if( (*dist[i])[j].GetProbability() > prob ) break;
                        }

			assert( j < l );

			batchScenarios[s] = allScenarios[j];
			//ArrayOfScenarios[s]->SetAgain( i, &( (*dist[i])[j] ) );

                        // David Love -- Print out the distribution locations, check that overlapping is working properly
                        //Print( "j[%d][%d] = %3d\t", s, i, j );
		}
                // David Love -- Get the seed values to use for the overlapping
                if( s + 1 == gamma )
                {
                   Random01::GetSeed( batchState );
                }
                // David Love -- Print out the distribution locations, check that overlapping is working properly
                //Print("\n");
	}

        // David Love -- A quick ckech of how many batches we are using.
        // batchNumber++;
        // fprintf( stderr, "Batch # %d, len %d\n", batchNumber, len );


}  //end of ReGenerateScenarios

double CalculateGap( int batchSize )
{
   double sumGaps = 0.0;
   for( int ii = 0; ii < batchSize; ii++ )
   {
      sumGaps += batchScenarios[ii];
   }
   sumGaps /= batchSize;

   return sumGaps;
}

/*------------------------------------------------------------------------------
 * ReadScenarios( filename ) reads all scenario information from the file with 
 * the name given.
 * This file should consist of two columns of data, in the double format.  The
 * first column is the optimality gap of that scenario.  The second column is 
 * the probability of that solution appearing.
 * ---------------------------------------------------------------------------*/
void ReadScenarios( string filename )
{
   ifstream scenarioFile;
   double   scenHolder;

   scenarioFile.open( filename.c_str() );
   if( !scenarioFile.is_open() )
   {
      cerr << "Cannot open file " << filename << endl;
      exit(1);
   }

   while( scenarioFile >> scenHolder )
   {
      allScenarios.push_back( scenHolder );
      scenarioFile >> scenHolder;
      allWeights.push_back( scenHolder );
   }
   for( vector<double>::size_type ii = 1; ii < allWeights.size(); ii++ )
   {
      allWeights[ii] += allWeights[ii-1];
   }

   scenarioFile.close();
   return;
}
