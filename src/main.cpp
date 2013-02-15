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
--------------------------------------------------------------------------------*/



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __VEC_POOL_H__
#	include "vec_pool.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __TIME_CNT_H__
#	include "time_cnt.h"
#endif

#ifndef __PARSEMPS_H__
#	include "parsemps.h"
#endif
#ifndef __RD_SUBLP_H__
#	include "rd_sublp.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __SCENTREE_H__
#	include "scentree.h"
#endif
#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __MAST_SOL_H__
#	include "mast_sol.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif

#ifndef __READ_TIM_H__
#	include "read_tim.h"
#endif
#ifndef __MAIN_H__
#	include "main.h"
#endif
#ifndef __DETERMLP_H__
#	include "determlp.h"
#endif
#ifndef __STOCHSOL_H__
#	include "stochsol.h"
#endif

//==============================================================================
//	Static functions used in the module --- prototypes.
//
static DeterministicLP *ReadDeterministicLP( void );

static Bool_T ReadTimeFile( MPS_LP &LP, Int_T &Stage1Row, Int_T &Stage1Col,
	Int_T &Stage2Row, Int_T &Stage2Col );

static Scenarios *ReadAndGenerateScenarios( DeterministicLP &LP, TimeInfo &TI );

static Bool_T RD_Crash( const DeterministicLP &DetermLP, Spc &SPC,
	Array<Real_T> &x1, const Int_T x1n );

//
//	Static data.
//

static DecompOptions DecompOpt;

//
//	End of static data and function prototypes.
//==============================================================================


//@BEGIN--------------------------------------------------

//to get the value of the penalty (and value of "mg") 
//declaring it to be global for convenience...
Real_T PEN; 
Int_T MGG, NFIX2, M2;
Int_T InitScen;
//Bool_T addscen; 

//to calculate true optimality gap... 
Real_T Opt_Tol; 

//------------------------------
//For sadsam, scenarios are better generated w/o their 
//cumulative probabilities but actual probabilities,
//to assure that we use the following boolean in 
//GenerateScenarios... maybe, somewhere else in the future
bool sadsam = false; 
//------------------------------

Real_T whatistolcut_2;  //just to print out the eps-opt thing

//@END----------------------------------------------------




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
	if( argc <= 1 )
	{
		PrintCopyright();
		PrintHelpScreen( "decomp" );
		return -1;
	}

	run( argc, argv );
	

	WorkVectorPool::CleanUp();

#ifndef NDEBUG
	long cnt = MemoryBlock::Count();

	if( cnt > 0 ) FatalError( "%ld memory blocks were not freed!", cnt );
#endif
  	
	
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

enum
{
	TI_START = 0,
	TI_READ_LP,
	TI_READ_TIME,
	TI_RD_CRASH,
	TI_READ_SCEN,
	TI_GEN_SCEN,
	TI_LP_DIVIDE,
	TI_SOLVE,

	TI_END_MARKER
};



void run( int argc, char *argv[] )
{

	Int_T i;
	
	//--------------------------------------------------------------------------
	//	Set up and start the timer.
	//
	TimeInfo TI( TI_END_MARKER );

	TI.MarkTime( TI_START );

	//--------------------------------------------------------------------------
	//	Parse the arguments from the argument list.
	//
	if( !ParseArguments( --argc, ++argv, DecompOpt ) )
	{
		PrintCopyright();
		PrintHelpScreen( *(argv - 1) );
		FatalError( "Argument error(s). Cannot continue." );
	}

	if( DecompOpt.Verbosity > V_NONE ) PrintCopyright();

	//--------------------------------------------------------------------------
	//	Read in the deterministic LP.
	//
	DeterministicLP *DetermLP = ReadDeterministicLP();


	if( DetermLP == NULL ) FatalError( "Errors processing the core file." );
	TI.MarkTime( TI_READ_LP );

	//--------------------------------------------------------------------------
	//	Read in the time file.
	//
	Int_T Stage1Row, Stage1Col, Stage2Row, Stage2Col;

	if( !ReadTimeFile( *DetermLP, Stage1Row, Stage1Col, Stage2Row, Stage2Col ) )
		FatalError( "Errors processing the time file: %s\n",
			DecompOpt.TimeFile );
	TI.MarkTime( TI_READ_TIME );

	//--------------------------------------------------------------------------
	//	Check if there are any non-zeros in the part of the matrix, which should
	//	be empty. If there are, signal error.
	//
	if( !DetermLP->CheckStructure( Stage2Row, Stage2Col ) )
		FatalError( "Invalid structure of the problem in the core file." );
	
	//--------------------------------------------------------------------------
	//	Decide on a pricing scheme for subproblems. Possible values for
	//	'SPC.Pricing': PRS_RC, PRS_SE, PRS_ASE
	//

	Spc SPC;
	SPC.Pricing = DecompOpt.Pricing;

	//--------------------------------------------------------------------------
	//	Run the RD crash: solve a deterministic problem.
	//
	Int_T x1n = Stage2Col;
	Array<Real_T> x1( x1n, 0.0 );

	if( DecompOpt.DoCrash && !RD_Crash( *DetermLP, SPC, x1, x1n ) )
		FatalError( "RD crash failed." );


	Print ("\nAFTER RD_CRASH:\n");
	for (i = 0; i< x1n; i++){
		Print ("x1[%1d] = " 
			   "%10g\n", 
			   i+1, x1.start[i] );
	}

	
	TI.MarkTime( TI_RD_CRASH );

	//--------------------------------------------------------------------------
	//	Read the stochastic data file and call the scenario generating routine.
	//	Then zero the random data in the deterministic LP object.
	//
	Scenarios *Scen = ReadAndGenerateScenarios( *DetermLP, TI);
	
	if( Scen == NULL ) FatalError( "\nFailed to create the scenarios!" );

	
	//@BEGIN--------------------------------------------------------------------
	// For Two-Replication Procedures, creating a second set of scenarios.
	//
//	Scenarios *Scen2 = ReadAndGenerateScenarios( *DetermLP, TI);
//	if( Scen2 == NULL ) FatalError( "\nFailed to create the second set of scenarios!" );
	
	//@END----------------------------------------------------------------------
	
	
	DetermLP->ZeroRandomData( (*Scen)[0] );

	Int_T ObjScale = DetermLP->ScaleObjective( *Scen );

	if( ObjScale && DecompOpt.Verbosity >= V_LOW )
		Print( "\nSCALING THE OBJECTIVE BY 1e%d.\n", (int) ObjScale );

	//--------------------------------------------------------------------------
	//	Create two linear problems: a master problem ('A') and a Slave problem
	//	('W'). Also create technology matrix 'T'. Set up the RD subproblem 'W'.
	//
	if( DecompOpt.Verbosity >= V_LOW )
		Print( "\nCONSTRAINTS ARE BEING DIVIDED INTO STAGES..." );

	SolvableLP A;
	SolvableLP T;
	RD_SubproblemLP W( T );

	DetermLP->DivideIntoStages( A, T, W, Stage2Col, (*Scen)[0] );

	if( DecompOpt.Verbosity >= V_LOW )
	{
		Real_T St1Col	= Stage2Col - Stage1Col,
			St2Col		= DetermLP->GetN() - Stage2Col,
			St1Row		= DetermLP->GetM() - DetermLP->GetActualStage2Rows()+1,
			St2Row		= DetermLP->GetActualStage2Rows(),
			ScenNo		= Scen->GetNumberOfScenarios(),
			LP_Col		= St1Col + ScenNo * St2Col,
			LP_Row		= St1Row + ScenNo * St2Row,
			MastRow		= St1Row + 2 * ScenNo,
			MastCol		= St1Col + 2 * ScenNo;

		Print(
			"DONE\n"
			"\tStage 1 rows:             %10g\n"
			"\tStage 1 columns:          %10g\n"
			"\tStage 2 rows:             %10g\n"
			"\tStage 2 columns:          %10g\n"
			"\tDeterministic equivalent problem dimensions:\n"
			"\tRows:                     %10g\n"
			"\tColumns:                  %10g\n"
			"\tMaster problem maximum dimensions:\n"
			"\tRows:                     %10g\n"
			"\tColumns:                  %10g\n",

			St1Row, St1Col, St2Row, St2Col,
			LP_Row, LP_Col, MastRow, MastCol
		);
	}

	W.ToStandard( V_NONE );
	W.InitializeRD_Subproblem();

	//--------------------------------------------------------------------------
	//	Renumber the rows and columns in scenarios. After that the deterministic
	//	LP is no longer needed and may be disposed of.
	//
	DetermLP->RenumberIndiceInScenarios( *Scen, Stage2Col );

	//@BEGIN-----------------------------------------------------------
	//Need to renumber indices in the scenarios for the second sample 
//	DetermLP->RenumberIndiceInScenarios( *Scen2, Stage2Col );
	//@END-------------------------------------------------------------


	delete DetermLP; DetermLP = NULL;

	TI.MarkTime( TI_LP_DIVIDE );

	//--------------------------------------------------------------------------
	//	Initialize the data for subproblem solver. Create a simplex solver
	//	object (together with it's auxiliary objects).
	//
	RD_SubproblemManager SubMan( A, W, SPC, DecompOpt.Restart );

	SubMan.SetScenarios( *Scen );
	SubMan.SetVerbosity( DecompOpt.Verbosity );

	//--------------------------------------------------------------------------
	//	Create and initialize the master solver object. Solve the problem.
	//
	MasterSolver master( A, *Scen );

	master.SetSubproblemManager( SubMan );
	
	//@BEGIN--------
	//Just to try what happens when there are no cuts and initial point is optimal
	
/*	//this is for cep1
	x1.start[0]=0.00;
	x1.start[1]=0.0;
	x1.start[2]=1166.6666666667;
	x1.start[3]=2500.0;
	x1.start[4]=0.00;
	x1.start[5]=0.00;
	x1.start[6]=2333.3333333333;
	x1.start[7]=3000.00;

*/	//@END----------

	
	master.SetStartingPoint( x1, x1n );
	master.SetPenalty( DecompOpt.InitPen );

	if( DecompOpt.Verbosity >= V_LOW )
		Print(
			"\nINVOKING THE REGULARIZED DECOMPOSITION OPTMIZER.\n"
			"\tThe initial penalty  value is %g.\n\n", master.GetPenalty()
		);
		

	master.Solve();

	TI.MarkTime( TI_SOLVE );
	PrintTimings( TI, true );

	//--------------------------------------------------------------------------
	//	Prepare and output the solution.
	//
	StochSolution *sol = master.GetSolution();

	assert( sol != NULL );
	if( ObjScale )
		sol->result *= pow( 10.0, ObjScale );

	PrintSolution( DecompOpt.SolutionFile, sol );


	Print ("\nOPTIMAL DECISION VECTOR:\n");
	for (i = 0; i< x1n; i++){
		Print ("x*[%1d] = " 
			   "%10g\n", 
			   i+1, sol->x[i] );
	}



	
//****************************************************************************************

	//THIS SECTION IS THE MRP BY ITSELF...........
	
        // David Love -- Added number of nonoverlapping batches
        const int numNonOverBatches = DecompOpt.NumNonOverBatches;  // Number of nonoverlapping batches
        const int batchSize = DecompOpt.ScenNum;
        // David Love -- Delete gamma when I have a working gamma input
        const int gamma = DecompOpt.NonOverlap;

        // David Love -- Rewrote the number of batches to work with overlapping
	const int numBatches = floor( (numNonOverBatches - 1) * (double) batchSize / (double) gamma + 1 );      //this is the smaller loop for the MRP
        // David Love -- Degrees of freedom of variance operator
        const int sampleSize = batchSize + (numBatches-1)*gamma;
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
	int numCIs = DecompOpt.Replications;//1000;                  // Added by David Love
	double za = 1.282;              // Added by David Love

        // David Love -- Variables to measure run time
        time_t startTime = time( NULL ),
               endTime = 0;

	var = new double[numCIs];
	gbar = new double[numCIs];
	ci = new double[numCIs];
	
	znstarbar = new double[numCIs];

	znstar = new double[numCIs*numBatches]; 

        // David Love -- Debug by listing all gap values
        //double gapValues[numCIs][numBatches];

        // David Love -- Testing that gamma, number of batches and degrees of freedom work correctly
        // Print( " gamma = %d \n num batches = %d \n degreesFreedom = %d\n", gamma, numBatches, degreesFreedom );
        // exit(0);
        
        // David Love -- Initializers for calculating total gap average
        indivGaps = new double[sampleSize];
        indivCounts = new int[sampleSize];
        gapest = new double[numBatches];


	for(kk=0; kk<numCIs; kk++){
                time_t loopStartTime = time( NULL );

		mns = 0.0; 
		vrs = 0.0; 
		zn  = 0.0; 

                // David Love -- Let the screen count where we are in numCIs
                printf( "Rep %d of %d, gamma = %d, m = %d, n = %d, nb = %d, degreesFreedom = %lf\n", kk, numCIs, gamma, batchSize, sampleSize, numBatches, degreesFreedom );

                // David Love -- Initialize gap counters
                for( int ii = 0; ii < sampleSize; ii++ )
                {
                   indivGaps[ii] = 0.0;
                   indivCounts[ii] = 0;
                }

		// David Love -- Generate a full selection of scenarios for each batch
                //Print( "j -- Beginning replication %d of %d\n", kk+1, numCIs );
		Scen->ReGenerateScenarios( batchSize, batchSize ); 
		for(oo = 0; oo<numBatches; oo++){
                        // David Love -- ReGenerate back at the end of the loop.  Initialized outside the looop
                        //Print( "j -- Batch %d of %d, overlap = %d\n", oo+1, numBatches, batchSize - gamma );
			Scen->ReGenerateScenarios( batchSize, gamma ); 


                        /* // @BEGIN Rebecca's code to print the scenarios
                        // Tested and Confirmed: Overlapping is working correctly!!!
                        Print( "SCEN Printing Scenarios\n" );
                        Print( "SCEN Batch %d of %d\n", oo+1, numBatches );
                        for (int scenCounter = 0; scenCounter < batchSize; scenCounter++)
                        {
                           Print( "SCEN " );
                           Print("%d ", scenCounter+1);
                           Scen->GetScenarioComponents(scenCounter);
                           Print("\n");
                        }
                        Print("SCEN\n");
                        */ // @END Rebecca's code


			//solve sampling problem:
			
			//this is to accelerate sol of sampling problem (cep1)
			x1.start[0]= 0.00;
			x1.start[1]= 0.00;
			x1.start[2]= 1833.33333333;
			x1.start[3]= 2500.00;
			x1.start[4]= 0.00;
			x1.start[5]= 0.00;
			x1.start[6]= 2333.33333333;
			x1.start[7]= 3000.00;
			
			//this is to accelerate solution of the sampling problem (pgp2) 
			//x1.start[0]= 1.5;  
			//x1.start[1]= 5.5;
			//x1.start[2]= 5.0;
			//x1.start[3]= 5.5;

                        //this is to accelerate solution of the sampling problem (apl1p)
                        //x1.start[0]= 1800.0;
                        //x1.start[1]= 1571.43;

			master.SetStartingPoint( x1, x1n );  
			
			master.ReInit();
			master.SetiStop(0); 
			master.Solve();

			master.FillSolution( sol );
			if( ObjScale )	sol->result *= pow( 10.0, ObjScale );

                        // David Love -- Commented out print statements, save room in file
			// PrintSolution( DecompOpt.SolutionFile, sol );

			// Print ("\nOPTIMAL DECISION VECTOR:\n");
			// for (i = 0; i< x1n; i++){
				// Print ("x*[%1d] = " 
					   // "%10g\n", 
					   // i+1, sol->x[i] );
			// }

			//store znstar

			znstar[kk*numBatches+oo] = sol->result; 

			zn += sol->result; 
			
			//solve xhat solution:

			//pgp2, this point is 0.63% from optimal
			//x1.start[0]= 2.5;  
			//x1.start[1]= 6.0;
			//x1.start[2]= 3.5;
			//x1.start[3]= 4.5;

			///cep1, this point is 10.73% from the optimal, 
			//x1.start[0]= 0.00;
			//x1.start[1]= 125.0;
			//x1.start[2]= 875.0;
			//x1.start[3]= 2500.00;
			//x1.start[4]= 0.00;
			//x1.start[5]= 625.00;
			//x1.start[6]= 1375.0;
			//x1.start[7]= 3000.00;
		
			///cep1, this point is 361% from the optimal, 
			x1.start[0]= 650;
			x1.start[1]= 650;
			x1.start[2]= 650;
			x1.start[3]= 650;
			x1.start[4]= 150;
			x1.start[5]= 150;
			x1.start[6]= 150;
			x1.start[7]= 150;

                        //apl1p, this point is ?? from the optimal,
                        //x1.start[0]= 1111.11;
                        //x1.start[1]= 2300;
		

                        // David Love -- Clearling out printed Information
			// Print("\nSTARTED XHAT SOLUTION...\n");
				
			master.SetStartingPoint( x1, x1n );  
			master.SetiStop( 10 );	     
  			master.ReInit();
			master.Solve(); 

                        // David Love -- Clearing out printed information
			// Print("FINISHED XHAT SOLUTION...\n\n");
			
			gapest[oo] = (master.CalculateGap() * pow(10.0, ObjScale)); 
                        //Print( "gapest[%d] = %lf\n", oo, gapest[oo] );

                        // David Love -- Loop for the new "G Double Bar"
                        for( int gammaCounter = 0; gammaCounter < batchSize; gammaCounter++ )
                        {
                           indivGaps[oo*gamma + gammaCounter] += ( master.GetIndivGap(gammaCounter) + master.GetIndivGap(batchSize) ) * pow(10.0, ObjScale);
                           indivCounts[oo*gamma + gammaCounter] += 1;
                           //Print( "indivGaps[%d] = %lf, indivCounts[%d] = %d\n", oo*gamma + gammaCounter, indivGaps[oo*gamma + gammaCounter], 
                                 //oo*gamma + gammaCounter, indivCounts[oo*gamma + gammaCounter] );
                           //Print( "ExpC[%d] = %lf, v[%d] = %lf\n", oo*gamma + gammaCounter, master.GetExpC(gammaCounter) + master.GetExpC(batchSize),
                                 //oo*gamma + gammaCounter, master.GetV(gammaCounter) + master.GetV(batchSize) );
                        }
		} //end of for with "oo" (0 == oo < numBatches). 
                // David Love -- Print time required for each replication
                endTime = time( NULL );
                Print( "Rep Time: %0.0lf\n", difftime( endTime, loopStartTime ) );


		//MRP calculations:

                // David Love -- New calculation of the gap estimate
                gbar[kk] = 0.0;
                for( int ii = 0; ii < sampleSize; ii++ )
                {
                   indivGaps[ii] /= indivCounts[ii];
                   gbar[kk] += indivGaps[ii];
                }
                gbar[kk] /= sampleSize;
                Print( "gbar[%d]   = %lf\n", kk, gbar[kk] );
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
		
	}  //end of loop with kk, 0 <= kk < numCIs

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
        // David Love -- 4/26/12, I am pretty sure that the lines below tih znstarber[kk] should actually contain gbar[kk].  
        // But these things have already been calculated and stored above as ci[kk].
        // The original lines of core are commented below my changes in the next two for loops.
        double ciSize = 0.0,
               varCISize = 0.0;
        for( kk = 0; kk < numCIs; kk++ )
                ciSize += ci[kk];
                // ciSize += znstarbar[kk] + za * sqrt( var[kk] ) / sqrt( ciDenom );
        ciSize = ciSize / numCIs;
        for( kk = 0; kk < numCIs; kk++ )
           varCISize += pow( ci[kk] - ciSize, 2 );
           // varCISize += pow( (znstarbar[kk] + za * sqrt( var[kk] ) / sqrt( ciDenom )) - ciSize, 2 );
        varCISize = varCISize / (numCIs-1);

        //for( kk = 0; kk < numCIs; kk++ )
           //ciWidth += znstarbar[kk];
        //ciWidth /= numCIs;
        //ciWidth += za * sqrt(avgVar) / sqrt( ciDenom );

	//print results:
	cov = (double) cov / numCIs;  

        // David Love -- Clearing out printed information
	// Print ("Optimality Gap = %f\n\n", optg); 
	
        // David Love -- Clearing out printed information
	// Print ("\nMRP RESULTS:\n"); 
	// Print ("-----------------\n\n"); 
	// Print ("k:         GAP:			SVAR:		CI			znstarbar:\n"); 
	// for(kk=0; kk<numCIs; kk++){
		// Print ("%d\t%f\t%20.10E\t%f\t%f\n", kk+1, gbar[kk], var[kk], ci[kk], znstarbar[kk]);
	// }

        // David Love -- Calculating the solution time
        endTime = time( NULL );

        Print("BEGIN MRP\n");
        Print("Batch Size = %d\n", batchSize);
        Print("gamma = %d\n", gamma);
        Print("CI Width = %lf\n", ciWidth);
	Print ("Coverage, MRP: = %lf\n", cov); 
        Print ("Average Variance = %lf\n", avgVar );
        Print("Variance of Variance = %le\n", varVariance );
        Print("Number Batches = %d\n", numBatches);
        Print("Degrees Freedom = %lf\n", degreesFreedom);
        Print("CI Denom = %lf\n", (double) ciDenom);
        Print("Var CI Width = %lf\n", varCIWidth);
        Print("CI Sze = %lf\n", ciSize);
        Print("Var CI Size = %lf\n", varCISize);
        Print("Time = %0.0lf\n", difftime( endTime, startTime ) );

	//print also the znbar
        // David Love -- Clearing out printed information
	// Print("\nznbar values:\n"); 
	// Print ("-----------------\n\n"); 
	// Print ("k:         znstar:\n"); 
	// for(kk=0; kk<numBatches*numCIs; kk++){
		// Print ("%d \t %f\n", kk+1, znstar[kk]);
	// }
        // David Love -- Some debugging code
        //Print( "Vales of mn\n" );
        //for( kk = 0; kk < numCIs; kk++ ) {
           //for( oo = 0; oo < numBatches; oo++ ) {
              //Print("%lf ", gapValues[kk][oo]);
           //}
           //Print("\n");
        //}
        //Print( "Variance Information\n" );
        //for( kk = 0; kk < numCIs; kk++ ) {
           //Print( "%lf\n", var[kk] );
        //}
	

	delete var; delete gbar; delete ci;
	var = NULL; gbar = NULL; ci = NULL;

        // David Love -- Clean up memory
        delete [] indivGaps;
        delete [] indivCounts;



//****************************************************************************************

	delete Scen; Scen = NULL;
//	delete Scen2; Scen2 = NULL; 
	
	delete sol;
	
} // end of run





/*------------------------------------------------------------------------------

	void PrintCopyright( void )
	void PrintHelpScreen( const char *ProgName )
	void PrintSubproblemStatistics( const Solver &s )
	void PrintSolution( const char *file, const StochSolution *sol );

PURPOSE:
	These routines take care of some simple screen output. They print a copy-
right notice, a help message, a set of statistics from subproblem solution, a
solution (first stage variables only) and recourse costs for scenarios resp.

PARAMETERS:
	const char *ProgName
		The name of the program (to be printed in the help message).

	const Solver &s
		Subproblem solver object.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void PrintCopyright( void )
{
	Print(
"----------------------------------------------------------------------------\n"
"REGULARIZED DECOMPOSITION  solver for 2 stage  stochastic  programs, v. 2.04\n"
"Coordinator and Master problem solver    (c)  1985-1995  Andrzej Ruszczynski\n"
"Primal SIMPLEX solver for subproblems    (c)  1992-1996   Artur Swietanowski\n"
"\n"
//"         This application was developed at IIASA, Laxenburg, Austria\n"
//"    and in the Institute of Automatic Control & Computation Engineering,\n"
//"                 Warsaw Univ. of Technology, Warsaw, Poland\n"
"----------------------------------------------------------------------------\n"
	);
}


void PrintHelpScreen( const char *ProgName )
{
	Print(
"Call: \"%s [options]\" where options are:\n"
"  [-problem] <problem_name>    - common name of data files (assumed default\n"
"                                 extensions: .cor, .tim and .sto),\n"
"  -cor <core_file_name>        - core file name,\n"
"  -tim <time_file_name>        - time period data file name,\n"
"  -sto <stoch_file_name>       - stochastic data file name,\n"
"  -s[cen] {all|<number>}       - scenario sample size or directive\n"
"                                 to generate all possible scenarios,\n"
"  -restart {tree|random|self*} - mode for subproblem solution restarts,\n"
"  -v <none|low*|high>          - output verbosity level,\n"
"  -txt_sol <solution_file>     - optional name of solution file,\n"
"  -pric {rc|se|ase*}           - simplex optimizer pricing mode.\n"
"\n"
"Asterisk marks default values. You may specify the problem files using\n"
"either the first option or all of the next three options.\n",

ProgName
	);
}


void PrintTimings( TimeInfo &TI, bool first)
{
	if (first)
	{
		double ReadLP_Time	= TI.TimeDifference( TI_START,		TI_READ_LP ),
			ReadTimeTime	= TI.TimeDifference( TI_READ_LP,	TI_READ_TIME ),
			RD_CrashTime	= TI.TimeDifference( TI_READ_TIME,	TI_RD_CRASH ),
			ReadScenTime	= TI.TimeDifference( TI_RD_CRASH,	TI_READ_SCEN ),
			GenScenTime		= TI.TimeDifference( TI_READ_SCEN,	TI_GEN_SCEN ),
			DivideTime		= TI.TimeDifference( TI_GEN_SCEN,	TI_LP_DIVIDE ),
			SolveTime		= TI.TimeDifference( TI_LP_DIVIDE,	TI_SOLVE ),
			TotalTime		= TI.TimeDifference( TI_START,		TI_SOLVE );

		Print(
			"\nTimings:\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%s\n"
			"\t%-40s:%7.1f\n",

			"Read LP",						ReadLP_Time,
			"Read time file",				ReadTimeTime,
			"Read scenario file",			ReadScenTime,
			"Generate scenarios",			GenScenTime,
			"Determ. LP divide",			DivideTime,
			"RD crash",						RD_CrashTime,
			"Solve",						SolveTime,
			"----------------------------------------",
			"TOTAL",						TotalTime
		);
	}
	else
	{
		double AppendScen   = TI.TimeDifference( TI_START,		TI_GEN_SCEN),
			SolveTime		= TI.TimeDifference( TI_GEN_SCEN,	TI_SOLVE ),
			TotalTime		= TI.TimeDifference( TI_START,		TI_SOLVE );

		Print(
			"\nTimings:\n"
			"\t%-40s:%7.5f\n"
			"\t%-40s:%7.5f\n"
			"\t%s\n"
			"\t%-40s:%7.5f\n",

			"Append Scenario",				AppendScen,
			"Solve",						SolveTime,
			"----------------------------------------",
			"TOTAL",						TotalTime
		);
	}
}


void PrintSubproblemStatistics( const Solver &s )
{
	Print(
		"\nSUBPROBLEM STATISTICS:\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n",

		"Total number of iterations:",
			(long) s.ReadStatisticCounter( Solver::TotalIter ),
		"Number of residuals' checks:",
			(long) s.ReadStatisticCounter( Solver::ResidCheck ),
		"Number of steepest edge resets:",
			(long) s.ReadStatisticCounter( Solver::SE_Reset ),
		"Number of primal var's computations:",
			(long) s.ReadStatisticCounter( Solver::PrimVarCompute ),
		"Number of dual var's computations:",
			(long) s.ReadStatisticCounter( Solver::DualVarCompute ),
		"Number of split pricings:",
			(long) s.ReadStatisticCounter( Solver::AltPric ),
		"Number of penalty adjustments:",
			(long) s.ReadStatisticCounter( Solver::PenaltyAdjust ),
		"Number of infeasibility minimizations:",
			(long) s.ReadStatisticCounter( Solver::InfeasMin )
	);

	Print(
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n"
		"\t%-40s%10ld\n",

		"Number of refactorizations:",
			(long) s.ReadStatisticCounter( Inverse::Refact ),
		"Number of basis updates:",
			(long) s.ReadStatisticCounter( Inverse::Upd ),
		"Number of sparse FTRAN's:",
			(long) s.ReadStatisticCounter( Inverse::SpFTRAN ),
		"Number of dense FTRAN's:",
			(long) s.ReadStatisticCounter( Inverse::DenFTRAN ),
		"Number of sparse BTRAN's:",
			(long) s.ReadStatisticCounter( Inverse::SpBTRAN ),
		"Number of dense BTRAN's:",
			(long) s.ReadStatisticCounter( Inverse::DenBTRAN )
	);
}


void PrintSolution( const char *file, const StochSolution *sol )
{
	Print( "\nOPTIMAL VALUE: %20.10E\n", sol->result );

	if( file && *file )
		((Solution *)sol)->WriteText( file, StochSolution::PrimalWeights );
}


/*------------------------------------------------------------------------------

	static DeterministicLP *ReadDeterministicLP( void )

PURPOSE:
	Reads a deterministic LP from a named input file.

PARAMETERS:
	None.

RETURN VALUE:
	Deterministic LP object pointer on success, NULL pointer on failure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static DeterministicLP *ReadDeterministicLP( void )
{
	//--------------------------------------------------------------------------
	//	Allocate memory. Exit on failure.
	//
	DeterministicLP *LP = new DeterministicLP;

	if( !LP )
		FatalError( "Not enough memory to allocate the deterministic LP." );

	//--------------------------------------------------------------------------
	//	Read the input file.
	//
	if( DecompOpt.Verbosity == V_LOW )
		Print( "\nREADING THE CORE FILE..." );

	if( !LP->ReadLP( DecompOpt.CoreFile, DecompOpt.Verbosity ) )
	{
		delete LP;
		LP = NULL;
		return NULL;
	}

	if( DecompOpt.Verbosity == V_LOW )
		Print( "DONE.\n" );

	//----------------------------------------------------------------------
	//	Find the objective function row. Treat the first free row as the
	//	objective function.
	//
	if( LP->GetCostRow() < 0 )
	{
		Error( "Objective function not found in the core file." );
		delete LP;
		return NULL;
	}

	//--------------------------------------------------------------------------
	//	If reading successful, return pointer to the deterministic LP object.
	//	Otherwise return a NULL pointer.
	//
	assert( LP != NULL );
	return LP;
}


/*------------------------------------------------------------------------------

	static Bool_T ReadTimeFile( MPS_LP &LP, Int_T &Stage1Row, Int_T &Stage1Col, 
		Int_T &Stage2Row, Int_T &Stage2Col )

PURPOSE:
	This function calls the time file reader. It reads the file and stores
the information about deterministic problem (LP) division into stages in
the remaining reference arguments.

PARAMETERS:
	MPS_LP &LP
		Deterministic LP.

	Int_T &Stage1Row, Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col
		"Return values" - they define LP division into stages.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Bool_T ReadTimeFile( MPS_LP &LP, Int_T &Stage1Row, Int_T &Stage1Col, // )
	Int_T &Stage2Row, Int_T &Stage2Col )
{
	if( DecompOpt.Verbosity == V_LOW )
		Print( "\nREADING THE TIME FILE..." );

	FILE *fp = fopen( DecompOpt.TimeFile, "r" );
	if( !fp )
		FatalError( "Unable to open time file %s.", DecompOpt.TimeFile );

	Bool_T Success = ReadTimes( fp, LP, Stage1Row, Stage1Col, Stage2Row,
		Stage2Col, DecompOpt.Verbosity );
	fclose( fp );

	if( DecompOpt.Verbosity == V_LOW )
		Print( "DONE.\n" );

	if( Stage1Row < 0 || Stage1Col < 0 ||
		Stage2Row < Stage1Row || Stage2Col < Stage1Col )
		FatalError( "Inconsistent data in time file. Invalid labels or "
			"wrong data order." );

	return Success;
}


/*------------------------------------------------------------------------------

	static Scenarios *ReadAndGenerateScenarios( DeterministicLP &LP,
		TimeInfo &TI )

PURPOSE:
	This function opens the file (name stored in the "DecompOpt" structure.

PARAMETERS:
	DeterministicLP &LP
		Deterministic LP structure.

RETURN VALUE:
	A pointer to created scenario repository object, or NULL if the attempt to
create the object failed.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Scenarios *ReadAndGenerateScenarios( DeterministicLP &LP, TimeInfo &TI )
{
	//--------------------------------------------------------------------------
	//	Create an appropriate scenario repository object.
	//
	Scenarios *Scen = ( DecompOpt.Restart == RD_SubproblemManager::TREE ) ?
		new TreeOfScenarios : new Scenarios;

	if( Scen == NULL ) FatalError( "Out of memory." );

	//--------------------------------------------------------------------------
	//	Open, read and then close the file.
	//
	if( DecompOpt.Verbosity >= V_LOW )
		Print( "\nREADING THE STOCH FILE AND GENERATING SCENARIOS..." );

	FILE *fp = fopen( DecompOpt.StochFile, "r"  );
	if( !fp ) FatalError( "Unable to open stoch file %s.", DecompOpt.StochFile );

	Bool_T Success = Scen->ReadScenarioFile( DecompOpt.StochFile, fp, LP );
	TI.MarkTime( TI_READ_SCEN );

	if( Success )
		Success = Scen->GenerateScenarios( DecompOpt.ScenNum );
	TI.MarkTime( TI_GEN_SCEN );

	fclose( fp );

	if( Success == False )
	{
		delete Scen;
		Scen = NULL;
	}

	if( DecompOpt.Verbosity >= V_LOW )
		Print( "DONE.\n" );

	return Scen;
}


static Bool_T RD_Crash( const DeterministicLP &DetermLP, Spc &SPC, // )
	Array<Real_T> &x1, const Int_T x1n )
{
	if( DecompOpt.Verbosity >= V_LOW )
		Print( "\nRUNNING RD CRASH (a deterministic LP solution)..." );
	
	//--------------------------------------------------------------------------
	//	Now create a SimplexLP object representing an optimization problem for
	//	one scenario. Solve the problem and retain the solution.
	//
	SimplexLP OneScenLP;

	if( ! OneScenLP.CreateAsSubmatrix( DetermLP, 0, DetermLP.GetM(), 0,
		DetermLP.GetN(), DetermLP.GetCostRow() ) )
	{
		Error( "Problems creating the crash LP." );
		return False;
	}

	//--------------------------------------------------------------------------
	//	Some more processing prior to solution.
	//
	OneScenLP.ScaleLP( V_NONE );
	OneScenLP.ToStandard( V_NONE );

	//--------------------------------------------------------------------------
	//	Create a simplex solver object and solve the problem.
	//
	Solver Solv( OneScenLP, SPC );

	Solv.InitializeSolver( V_NONE );
	SOLVE_RESULT SolverStatus = Solv.Solve( V_NONE );

	if( SolverStatus != SR_OPTIMUM && SolverStatus != SR_INFEASIBLE &&
		SolverStatus != SR_UNBOUNDED )
	{
		Error( "Problems solving the crash LP." );
		return False;
	}

	//--------------------------------------------------------------------------
	//	Get the solution. Store it in the "x1" vector of length "x1n".
	//
	Solution *sol = Solv.GetSolution( Solution::Primal, False );
	OneScenLP.ProcessSolution( *sol );

	x1.Copy( sol->x, sol->GetAllocN(), x1n, x1n );
	delete sol;

	//--------------------------------------------------------------------------
	//	Return the success status.
	//
	if( DecompOpt.Verbosity >= V_LOW )
		Print( "DONE!\n" );
	return True;
}

