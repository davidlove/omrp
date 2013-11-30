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
	--------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>

//@BEGIN----------

#include <time.h>

#ifndef __RAND01_H__
#	include "rand01.h"
#endif

//@END------------

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

/* PPL ADDENDUM -- INCLUDES */
#include <vector>
/* PPL ADDENDUM -- TERMINUS */

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
//==============================================================================
//to get the value of the penalty (and value of "mg") 
//declaring it to be global for convenience...
//

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

//==============================================================================
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
	
	------------------------------------------------------------------------------
*/

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
	
	//if( cnt > 0 ) FatalError( "%ld memory blocks were not freed!", cnt );
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
	
	------------------------------------------------------------------------------
*/

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
	
	/* RS addendum - fix the seed */
	srand(DecompOpt.seed_value);
	//Print("*****SEED********** %d\n", DecompOpt.seed_value);
	
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
	
	//Print ("\nAFTER RD_CRASH:\n");
	for (int i = 0; i< x1n; i++){
		//Print ("x1[%1d] = " 
					// "%10g\n", 
					// i+1, x1.start[i] );
	}
	
	TI.MarkTime( TI_RD_CRASH );
	
	//@BEGIN-------------------------------------
	
	int seed1 = 1*DecompOpt.seed_value; 
	int seed2 = 2*DecompOpt.seed_value; 
	int seed3 = 3*DecompOpt.seed_value; 
	
	Random01::Seed( seed1, seed2, seed3 ); 

	//@END---------------------------------------
	
	//--------------------------------------------------------------------------
	//	Read the stochastic data file and call the scenario generating routine.
	//	Then zero the random data in the deterministic LP object.
	//
	Scenarios *Scen = ReadAndGenerateScenarios( *DetermLP, TI);
	
	if( Scen == NULL ) FatalError( "\nFailed to create the scenarios!" );
	
	//@BEGIN--------------------------------------------------------------------
	// For Two-Replication Procedures, creating a second set of scenarios.
	//
	DecompOpt.ScenNum /= 2;  // /= 4;  //again, because m_k = 2*n_k... and now, we
	// are now doing 2-Rep procedures... 
	
	Scenarios *Scen2 = ReadAndGenerateScenarios( *DetermLP, TI);
	if( Scen2 == NULL ) FatalError( "\nFailed to create the second set of scenarios!" );
	
	Scenarios *Scen3 = ReadAndGenerateScenarios( *DetermLP, TI);
	if( Scen3 == NULL ) FatalError( "\nFailed to create the second set of scenarios!" );
	
	DecompOpt.ScenNum *= 2;   //change back to initial value
	
	//@END----------------------------------------------------------------------
	
	DetermLP->ZeroRandomData( (*Scen)[0] );
	
	Int_T ObjScale = DetermLP->ScaleObjective( *Scen );
	
	//if( ObjScale && DecompOpt.Verbosity >= V_LOW )
		//Print( "\nSCALING THE OBJECTIVE BY 1e%d.\n", (int) ObjScale );
	
	//--------------------------------------------------------------------------
	//	Create two linear problems: a master problem ('A') and a Slave problem
	//	('W'). Also create technology matrix 'T'. Set up the RD subproblem 'W'.
	//
	//if( DecompOpt.Verbosity >= V_LOW )
		//Print( "\nCONSTRAINTS ARE BEING DIVIDED INTO STAGES..." );
	
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
			
			/*Print(
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
						); */
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
	DetermLP->RenumberIndiceInScenarios( *Scen2, Stage2Col );
	DetermLP->RenumberIndiceInScenarios( *Scen3, Stage2Col );
	
	//also get the number of scenarios to be used later on
	Int_T ScGiven = DecompOpt.ScenNum;
	
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
	MasterSolver master( A, *Scen, true );
	
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
	
	//if( DecompOpt.Verbosity >= V_LOW )
		/*Print(
					"\nINVOKING THE REGULARIZED DECOMPOSITION OPTMIZER.\n"
					"\tThe initial penalty  value is %g.\n\n", master.GetPenalty() 
					); */
	
	//@BEGIN----------------------------
	//   just to get timings
	const int rep = 10;
	
	clock_t start, finish;		//to get the whole running time
	double *runtime; 
	runtime = new double[rep];
	
	clock_t start1, finish1;	//to get the Assesing Sol Qual part time
	double *assesstime;
	assesstime = new double[rep];
	
	for (int a=0; a < rep; a++){
		runtime[a] = 0;
		assesstime[a] = 0;		
	}
	
	start = clock(); 
	
	//@END------------------------------
	
	//@BEGIN--------------------------------------------------------------------
	// Now, let's have a second set of master and a SubMan; 
	//      called master2 and SubMan2, obviously!  :o)
	// 
	
	// define variables to store mastersolve parameters of MGG, NFIX2 and M2
	Int_T  mg_store,  nfix_store,  m_store,		//to store parameters for master
		mg_store2, nfix_store2, m_store2,	//"    "		"	   "  master2
		mg_store3, nfix_store3, m_store3;	//"    "		"	   "  master3
	
	// store mastersolve parameters before switching to the second set
	mg_store   = MGG;
	nfix_store = NFIX2;
	m_store    = M2;
	
	// change Scen Number for 2 rep
	DecompOpt.ScenNum /= 2; 
	
	// now, create SubMan2 and master2
	RD_SubproblemManager SubMan2( A, W, SPC, DecompOpt.Restart );
	SubMan2.SetScenarios( *Scen2 );
	SubMan2.SetVerbosity( DecompOpt.Verbosity );
	
	MasterSolver master2( A, *Scen2, false );
	master2.SetSubproblemManager( SubMan2 );
	master2.SetStartingPoint( x1, x1n );
	master2.SetPenalty( DecompOpt.InitPen );
	
	// and, solve master2 with Scen2
	//Print("\nSolving the Second Sampling Problem with Scen2 and master2:\n"); 
	
	start1 = clock(); 
	
	master2.Solve( true, 1E-08 );
	
	finish1 = clock(); 
	
	assesstime[0] =  double (finish1 - start1) / double (CLOCKS_PER_SEC) ; 
	
	StochSolution *sol_2 = master2.GetSolution();
	
	assert( sol_2 != NULL );
	if( ObjScale )	sol_2->result *= pow( 10.0, ObjScale );
	//PrintSolution( DecompOpt.SolutionFile, sol_2 );
	//Print ("\nOPTIMAL DECISION VECTOR:\n");
	for (int i = 0; i< x1n; i++){
		/*Print ("x*[%1d] = " 
					 "%10g\n", 
					 i+1, sol_2->x[i] ); */	}
	
	//now, store the solve parameters for the second set
	mg_store2   = MGG;
	nfix_store2 = NFIX2;
	m_store2    = M2;
	
	//88888888888888888888888888888888888888888888888888888888888888888888888888
	//The same for the third set: Create SubMan3, master3, solve, store
	
	RD_SubproblemManager SubMan3( A, W, SPC, DecompOpt.Restart );
	SubMan3.SetScenarios( *Scen3 );
	SubMan3.SetVerbosity( DecompOpt.Verbosity );
	
	MasterSolver master3( A, *Scen3, false );
	master3.SetSubproblemManager( SubMan3 );
	master3.SetStartingPoint( x1, x1n );
	master3.SetPenalty( DecompOpt.InitPen );	
	
	//Print("\nSolving the Third Sampling Problem with Scen3 and master3:\n"); 
	
	start1 = clock(); 
	
	master3.Solve( true, 1E-08 );
	
	finish1 = clock(); 
	
	assesstime[0] +=  double (finish1 - start1) / double (CLOCKS_PER_SEC) ; 
	
	StochSolution *sol_3 = master3.GetSolution();
	
	assert( sol_3 != NULL );
	if( ObjScale )	sol_3->result *= pow( 10.0, ObjScale );
	//PrintSolution( DecompOpt.SolutionFile, sol_3 );
	//Print ("\nOPTIMAL DECISION VECTOR:\n");
	for (int i = 0; i< x1n; i++){
	/*	Print ("x*[%1d] = " 
					 "%10g\n", 
					 i+1, sol_3->x[i] ); */
	}
	
	//now, store the solve parameters for the second set
	mg_store3   = MGG;
	nfix_store3 = NFIX2;
	m_store3    = M2;
	
	//88888888888888888888888888888888888888888888888888888888888888888888888888
	//change back Scennum (original number of scenarios)
	DecompOpt.ScenNum *= 2; 
	
	//@END----------------------------------------------------------------------
	
	int r;					//counter for the outer loop
	
	//  --procedure parameters:
	const double Z_star = 24642.32;
	const double   za =  1.282; 
	
	//  --for running the procedure:
	//	Int_T LastScen,			//stores the last scenario number for master
	//	NewScen;			//stores the new scenario number for master
	
	//  --for generating xhat:
	int nk,				                 	//k = counter for iteration k
		s,				                        //s = start of the new xhat location in the array
		count, count2, count3;			     	//count = to count the number of times we solve with warm start. 
	
	//  --for assessing xhat:	
	double *gap,			//gap estimate
		*svar;			//variance estimate 
	
	double gap1, gap2,		//used in calculation of gap & var in 2 Replications
		svar1, svar2;		
		
	/* RS addendum */
#ifdef halfgap
	double gap_first_half_1, gap_first_half_2,
		gap_second_half_1, gap_second_half_2;
			
	double *gap_first_half, *gap_second_half;
#endif	
		
	//  --for output of the procedure:
	double *xhat_T; 		//stores the stopping candidate solution
	double *xhat;			                //stores the xhat_k for iteration k; 
	double *bound;			//stores the bound (with sk term in it) 
	double  ak; 			//the inflation term (we take it here as 1/sqrt(nk))
	
	FILE *fdetail2;
	
	/* RS addendum - this makes the output file names change dynamically with sample size, test problem, and opt/sub candidate solution */			
	char filename_iter[30];
	char filename_rep[30];
	char temp[10];
	
	sprintf(filename_iter, "%d", ScGiven);
	sprintf(filename_rep, "%d", ScGiven);
	sprintf(temp, "%d", DecompOpt.seed_value);
	strcat(filename_iter, "_");
	strcat(filename_rep, "_");
	strcat(filename_iter, temp);
	strcat(filename_rep, temp);
	
#ifdef match 
	strcat(filename_iter, "B_iter.txt");
	strcat(filename_rep, "B");
#else
	strcat(filename_iter, "_iter.txt");
#endif

#ifdef apl1p_opt 
	strcat(filename_rep, "_apl1p.txt");
#endif
#ifdef apl1p_sub 
	strcat(filename_rep, "_apl1p.txt");
#endif
#ifdef apl1p_sub2
	strcat(filename_rep, "_apl1p.txt");
#endif
#ifdef cep1_opt	
	strcat(filename_rep, "_cep1.txt");
#endif
#ifdef cep1_sub
	strcat(filename_rep, "_cep1.txt");
#endif
#ifdef pgp2_opt	
	strcat(filename_rep, "_pgp2.txt");
#endif
#ifdef pgp2_sub 
	strcat(filename_rep, "_pgp2.txt");
#endif
#ifdef pgp2_sub2 
	strcat(filename_rep, "_pgp2.txt");
#endif
#ifdef gbd_opt
	strcat(filename_rep, "_gbd.txt");
#endif
#ifdef gbd_sub
	strcat(filename_rep, "_gbd.txt");
#endif



	/* RS addendum */	
#ifdef match
	/* RS addendum */
	int *PerfectMatch;      // Stores the perfect matchings sequentially	
	/* end */
#endif
	
	fdetail2 = fopen(filename_iter, "w");
	if ( !fdetail2 )
		exit (1); 
	fprintf(fdetail2, "RUN_ID; Problem; n_0; Z_star; REP; xhat; xgap1; gap1; svar1; xgap2; gap2; svar2; bound; RunTime; AssessTime");
	
	//--------------------------------------------------------------------------
	//Initializations
	//--------------------------------------------------------------------------
	
	//to store the xhat I need the size of the x for the problem
	//but that is already stored in the x1n
	xhat = new double[x1n];
	xhat_T = new double[rep * x1n];

	/* RS addendum */	
#ifdef apl1p_opt 
	xhat[0] = 1800;         // optimal xhat for APL1P
	xhat[1] = 1571.43; 
#endif

	/* RS addendum */
#ifdef apl1p_sub 
	xhat[0] = 1111.11;      // suboptimal xhat for APL1P
	xhat[1] = 2300;
#endif

	/* RS addendum */		
#ifdef apl1p_sub2
	xhat[0] = 2000;         // second suboptimal xhat for APL1P
	xhat[1] = 500;			// optimality gap is 1242.69
#endif

	/* RS addendum */			
#ifdef cep1_opt 
	xhat[0] = 0;            //optimal xhat for CEP1
	xhat[1] = 0;				
	xhat[2] = 1833.33;
	xhat[3] = 2500;
	xhat[4] = 0;
	xhat[5] = 0;
	xhat[6] = 2333.33;
	xhat[7] = 3000;
#endif

	/* RS addendum */	
#ifdef cep1_sub 
	xhat[0] = 0;            //suboptimal xhat for CEP1
	xhat[1] = 125;				
	xhat[2] = 875;
	xhat[3] = 2500;
	xhat[4] = 0;
	xhat[5] = 625;
	xhat[6] = 1375;
	xhat[7] = 3000;
#endif			

	/* RS addendum */	
#ifdef pgp2_opt
	xhat[0] = 1.5;          //optimal xhat for PGP2
	xhat[1] = 5.5;
	xhat[2] = 5;
	xhat[3] = 5.5;
#endif
	
	/* RS addendum */
#ifdef pgp2_sub 
	xhat[0] = 1.5;          //suboptimal xhat for PGP2
	xhat[1] = 5.5;
	xhat[2] = 5;
	xhat[3] = 4.5;
#endif
	
	/* RS addendum */
#ifdef pgp2_sub2
	xhat[0] = 2;            //second suboptimal xhat for PGP2
	xhat[1] = 6;
	xhat[2] = 3.5;
	xhat[3] = 4.5;			//optimality gap is 15.338	
#endif

	/* RS addendum */
#ifdef gbd_opt 
	xhat[0] = 10;           //optimal xhat for GBD
	xhat[1] = 0;
	xhat[2] = 0;
	xhat[3] = 0;
	xhat[4] = 0;
	xhat[5] = 12.4793;
	xhat[6] = 1.18736;
	xhat[7] = 5.33333;
	xhat[8] = 0;
	xhat[9] = 4.24138;
	xhat[10] = 0;
	xhat[11] = 20.7586;
	xhat[12] = 7.80104;
	xhat[13] = 0;
	xhat[14] = 7.19896;
	xhat[15] = 0;
	xhat[16] = 0;
#endif

	/* RS addendum */
#ifdef gbd_sub 
	xhat[0] = 10;           //suboptimal xhat for GBD
	xhat[1] = 0;
	xhat[2] = 0;
	xhat[3] = 0;
	xhat[4] = 0;
	xhat[5] = 12.2724;
	xhat[6] = 2.92759;
	xhat[7] = 3.8;
	xhat[8] = 0;
	xhat[9] = 4.65517;
	xhat[10] = 0;
	xhat[11] = 20.3448;
	xhat[12] = 9.13574;
	xhat[13] = 0;
	xhat[14] = 5.86426;
	xhat[15] = 0;
	xhat[16] = 0;
#endif
	
	gap   = new double[rep]; 
	
	#ifdef halfgap
		gap_first_half = new double[rep];
		gap_second_half = new double[rep];
	#endif
	
	svar  = new double[rep]; 
	bound = new double[rep]; 

	/* RS addendum */
#ifdef match
	PerfectMatch = new int[1500];
#endif	

	//--------------------------------------------------------------------------
	//The outer loop with "rep"
	//--------------------------------------------------------------------------
	
		/* RS addendum */
#ifdef match 	
	int flag;	
	FILE *edges;
	FILE *matches;
#endif	
	
	for(r = 0; r < rep; r++){
		
		for(int i = 0; i<x1n; i++)
			xhat_T[ r*x1n + i ] = xhat[ i ]; 
		
		//Print("\nOUTER LOOP - REP = %d:\n", r+1); 	
		
		//--------------------------------------------------------------------------
		//The sequential procedure
		//--------------------------------------------------------------------------		
		// start = clock();
		
		//initialize stuff before entering sequential procedure
		nk = ScGiven;
		count = 0;
		count2 = 0;
		count3 = 0;
		
		fprintf(fdetail2, "\n120309;AP1LP;%d;%f;", nk, Z_star);
		fprintf(fdetail2, "%d;", r+1);
		
		fprintf(fdetail2, "(");
		for (int i = 0; i< x1n; i++){
			if (i ==0)
				fprintf(fdetail2, "%g, ", xhat[i]);
			else
				fprintf(fdetail2, "%g); ", xhat[i]);
		}

		/* RS addendum */
#ifdef match
		
		/* PPL BIAS ADDENDUM */
		
/* Calculating the distances and create edge file */
		edges = fopen("edge_file.txt", "w");
		if (!edges){
			printf("error opening edge file");
			exit (1);
		}
		
		fprintf(edges, "%d %d\n", nk, nk*(nk - 1)/2);
		for (int i = 0; i < nk - 1; i++ ){
			for (int j = i + 1; j < nk; j++ ){
				int dist = (int) 1000*(Scen->GetDistance(i, j));
				fprintf(edges, "%d %d %d\n", i, j, dist );
			}
		}
		fclose( edges );
		
		system("./blossom4 -e edge_file.txt -w match_file.txt");
		
		/* Store Matching in Array */
		matches = fopen("match_file.txt", "r");
		
		if (!matches){
			printf("error opening match file");
			exit (1); 
		}
		
		int NumIn = 0; 
		fscanf(matches, "%d", &NumIn);
		fscanf(matches, "%d", &NumIn);
		
		int i = 0;
		while(!feof(matches) && i < nk){
			fscanf(matches, "%d", &PerfectMatch[i]); 
			i++;
			fscanf(matches, "%d", &PerfectMatch[i]); 
			i++;	
			fscanf(matches, "%d", &NumIn); 			
		}
		
		fclose( matches ); 
		
		/* Partitioning Scen */
		for (int i = 0; i < nk; i++){
			Scen2->ScenCopy(i/2, Scen->GetScen(PerfectMatch[i]));  
			i = i+1;
			Scen3->ScenCopy(i/2, Scen->GetScen(PerfectMatch[i]));  
		}	
		
		for (int i = 0; i < nk/2; i++){
		//	Print("%f\n", Scen2->GetValue(i));
		}		
		for (int i = 0; i < nk/2; i++){
		//	Print("%f\n", Scen3->GetValue(i));
		}	
				

#else

		/* Partitioning Scen -- Alternating */
		for (int i = 0; i < nk; i++){
			Scen2->ScenCopy(i/2, Scen->GetScen(i));
			i = i+1;
			Scen3->ScenCopy(i/2, Scen->GetScen(i)); 
		}
		
		for (int i = 0; i < nk/2; i++){
			//Print("%f\n", Scen2->GetValue(i));
		}		
		for (int i = 0; i < nk/2; i++){
			//Print("%f\n", Scen3->GetValue(i));
		}	

#endif

/*		Print("**********Printing Scenarios***********\n");
		for (int i = 0; i < nk/2; i++)
		{
			Print("%d ", i+1);
			Scen2->GetScenarioComponents(i);
			Print("\n");
		}
		Print("\n");
		for (int i = 0; i < nk/2; i++)
		{
			Print("%d ", i+1);
			Scen3->GetScenarioComponents(i);
			Print("\n");
		}	*/			

		
		
		/* PPL BIAS END */
		
		//STEP 2: 
		//----------------------------------------------------------------------
		//Assess the quality of xhat
		//----------------------------------------------------------------------
		
		// change Scen Number  /* RS ADDENDUM ************* */
		//DecompOpt.ScenNum /= 2; 
		
		//prepare to solve from scratch.
		//Print("Solving Assessing Prob from Scratch..\n");
		SubMan2.ReSizeSolverState( DecompOpt.ScenNum, InitScen); 
		SubMan2.ChangePreviousBlockNumber(-2); 
		master2.SetWeights(); 
		master2.SetStartingPoint(x1, x1n);   //?? what is x1 right now ??
		master2.SetiStop(0); 
		master2.ReInit();
		
		//solve, get and print solution
		start1 = clock(); 
		master2.Solve( true, 1E-08 );
		finish1 = clock(); 
		assesstime[r] += (  double (finish1 - start1) / double (CLOCKS_PER_SEC) ); 
		
		master2.FillSolution( sol_2);
		if( ObjScale )  sol_2->result *= pow( 10.0, ObjScale );
		PrintSolution( DecompOpt.SolutionFile, sol_2 );
		//Print ("\nOPTIMAL DECISION VECTOR:\n");
		for (int i = 0; i< x1n; i++){
		//	Print ("x*[%1d] =  %10g\n", i+1, sol_2->x[i] );
		}
		
		//Store the solve parameters
		mg_store2   = MGG;
		nfix_store2 = NFIX2;
		m_store2    = M2;				
		
		//prepare to solve from scratch.
		//Print("Solving Assessing Prob from Scratch..\n");
		SubMan3.ReSizeSolverState( DecompOpt.ScenNum, InitScen); 
		SubMan3.ChangePreviousBlockNumber(-2); 
		master3.SetWeights(); 
		master3.SetStartingPoint(x1, x1n);   //?? what is x1 right now ??
		master3.SetiStop(0); 
		master3.ReInit();
		
		count3 = 0;   //THIS should be done at SCEN3... 
		
		//solve, get and print solution
		start1 = clock(); 
		master3.Solve( true, 1E-08 );
		finish1 = clock(); 
		assesstime[r] += (  double (finish1 - start1) / double (CLOCKS_PER_SEC) ); 
		
		master3.FillSolution( sol_3);
		if( ObjScale )  sol_3->result *= pow( 10.0, ObjScale );
		PrintSolution( DecompOpt.SolutionFile, sol_3 );
		//Print ("\nOPTIMAL DECISION VECTOR:\n");
		
		for (int i = 0; i< x1n; i++){
		//	Print ("x*[%1d] =  %10g\n", i+1, sol_3->x[i] );
		}
		
		//Store the solve parameters
		mg_store3   = MGG;
		nfix_store3 = NFIX2;
		m_store3    = M2;
		
		//-------------------------------------------------------------------
 		//STEP 2.2 (1): Assess the quality of the current candidate solution:
		//-------------------------------------------------------------------
		
		//Print("\n ...at rep = %d\n", r+1);
		
		//get xhat function values to form the gap and variance estimates
		//Print("\nASSESSING XHAT QUALITY...\n");
		
		start1 = clock(); 	
		
		/* RS addendum */		
#ifdef halfgap
	//	Print("First Half\n");
		master2.AssessSolQual(xhat, gap1, gap_first_half_1, gap_second_half_1, svar1, false);
#else								
		master2.AssessSolQual(xhat, gap1, svar1, false);
#endif
		
		finish1 = clock(); 
		assesstime[r] += (double (finish1 - start1) / double (CLOCKS_PER_SEC) ); 
		
		gap1   *= pow( 10.0, ObjScale );
		svar1  *= pow( 100.0, ObjScale ); 
	
		/* RS addendum */		
#ifdef halfgap
		gap_first_half_1 *= pow( 10.0, ObjScale );
		//Print("First half 1 %f\n ", gap_first_half_1);
		gap_second_half_1 *= pow( 10.0, ObjScale );
#endif	
		
		fprintf(fdetail2, "(");
		for (int i = 0; i< x1n; i++){
			if (i ==0)
				fprintf(fdetail2, "%g, ", sol_2->x[i]);
			else
				fprintf(fdetail2, "%g); ", sol_2->x[i]);
		}
		fprintf(fdetail2, "%g;%g;", gap1, svar1);
		

		
		//-------------------------------------------------------------------
		//STEP 2.2 (2): Assess the quality of the current candidate solution:
		//-------------------------------------------------------------------
		
		//Print("\n ...at rep = %d\n", r+1);
		//get xhat function values to form the gap and variance estimates
		//Print("\nASSESSING XHAT QUALITY...\n");
		
		start1 = clock(); 	
		
		/* RS addendum */				
#ifdef halfgap
		//Print("Second Half\n");
		master3.AssessSolQual(xhat, gap2, gap_first_half_2, gap_second_half_2, svar2, false);
#else								
		master3.AssessSolQual(xhat, gap2, svar2, false);
#endif
		
		finish1 = clock(); 
		assesstime[r] += (  double (finish1 - start1) / double (CLOCKS_PER_SEC) ); 
		
		gap2   *= pow( 10.0, ObjScale );
		svar2  *= pow( 100.0, ObjScale ); 

		/* RS addendum */
#ifdef halfgap
		gap_first_half_2 *= pow( 10.0, ObjScale );
		//Print("First half 2 %f\n ", gap_first_half_2);
		gap_second_half_2 *= pow( 10.0, ObjScale );
#endif
		
		fprintf(fdetail2, "(");
		for (int i = 0; i< x1n; i++){
			if (i ==0)
				fprintf(fdetail2, "%g, ", sol_3->x[i]);
			else
				fprintf(fdetail2, "%g); ", sol_3->x[i]);
		}
		fprintf(fdetail2, "%g;%g;", gap2, svar2);
		
		//888888888888888888888888888888888888888888888888888888888888888888
		
		//Now, calculate the "pooled" gap & var estimates for A2RP
		gap[r]  = 0.5 * (gap1 + gap2); 
		svar[r] = 0.5 * (svar1 + svar2); 
		
		/* RS addendum */
#ifdef halfgap
		gap_first_half[r]  = 0.5 * (gap_first_half_1 + gap_first_half_2); 
		gap_second_half[r] = 0.5 * (gap_second_half_1 + gap_second_half_2);
#endif
		
		//---------------------------------------
		//STEP 2.3: Check termination condition:
		//---------------------------------------
		
		ak = 1/ sqrt( (double) nk); 
		
		finish = clock();
		runtime[r] +=  double (finish - start) / double (CLOCKS_PER_SEC) ;
		
		bound[r] = gap[r] + za*sqrt(svar[r])/sqrt( (double) nk);
		fprintf(fdetail2, "%g;", bound[r]);
		fprintf(fdetail2, "%g;%g", runtime[r], assesstime[r]);
		
		//Print("The bound on mu_T is = %g\n", bound[r]); 
		
		//--------------------------------------------------------------------------
		//Regenerate Scenarios to start another run
		//--------------------------------------------------------------------------
		
		if( r != rep-1){
			
			start = clock();
			
			seed1 += 1;  
			seed2 += 2; 
			seed3 += 3; 
			
			Random01::Seed( seed1, seed2, seed3 ); 			
			
			DecompOpt.ScenNum = ScGiven;  //*2;  //because m_k=2*n_k 
			
			master.SetNewl( DecompOpt.ScenNum +1 ); 
			Scen->ReGenerateScenarios( DecompOpt.ScenNum ); 
			SubMan.SetNumOfScenarios( DecompOpt.ScenNum ); 
			//SubMan.ReSizeSolverState( DecompOpt.ScenNum, LastScen);
			SubMan.ChangePreviousBlockNumber( -2 );
			
		}	//end of if with the "r != rep-1" check
		
	}	//end of for loop with "rep"
	
	//--------------------------------------------------------------------------
	//Print summary of results
	//--------------------------------------------------------------------------
	
	/*Print("gap1 = %f\n", gap1);
	Print("gap2 = %f\n", gap2);
	Print("svar1 = %f\n", svar1);
	Print("svar2 = %f\n", svar2); */
	
	FILE *fdetail;
	
	fdetail = fopen(filename_rep, "w");
	if ( !fdetail )
		exit (1); 
	
	//Print("\nSUMMARY OF RESULTS:"); 
	//Print("\n-------------------\n"); 
	
	//here, use count as counter to how many times the procedure stopped
	count = 0; 
	
		/* RS addendum */
#ifdef halfgap 
		/*Print("\nREP: ;    GAP: ;      SVAR: ;   xhat:        ;  bound: "); 
		Print("\n-----;---------;------------;------;----------------------;---------"); */

		for(r = 0; r < rep; r++){
		
			if(r%20 == 0)
			{
				if(r == 0)
				{
					fprintf(fdetail, "REP; n_0; EST_GAP; SVAR; xhat; bound; Proced; Z_star; TotalRep; RUN_ID; AssessTime");
				}
				else
				{
					fprintf(fdetail, "\nREP; n_0; EST_GAP; SVAR; xhat; bound; Proced; Z_star; TotalRep; RUN_ID; AssessTime");

				}
			}
			int i;
			s = r*x1n;
		
			//Print("\n%d;%g;%g;", r+1, gap[r], svar[r]); 
			fprintf(fdetail, "\n%d;%d;%g;%g;%g;%g;", r+1, nk, gap[r], gap_first_half[r], gap_second_half[r], svar[r]);
		
			//Print("("); 
			fprintf(fdetail, "(");
			for (i = 0; i< x1n-1; i++){
				//Print ("%g, ", xhat_T[s+i]);
				fprintf(fdetail, "%g, ", xhat_T[s+i]);
		}
		//i++; 
		//Print("%g);", xhat_T[s+i]);
		fprintf(fdetail, "%g);", xhat_T[s+i]);
		
		//Print("%g", bound[r]); 
		fprintf(fdetail, "%g;", bound[r]);
		
		fprintf(fdetail, "%f;%d;120309;APL1P;%g", Z_star, rep, assesstime[r]);
			
#else
		//Print("\nREP: ;    GAP: ;      SVAR: ;   xhat:        ;  bound: "); 
		//Print("\n-----;---------;------------;------;----------------------;---------");

		for(r = 0; r < rep; r++){
		
			if(r%20 == 0)
			{
				if(r == 0)
				{
					fprintf(fdetail, "REP; n_0; EST_GAP; SVAR; xhat; bound; Proced; Z_star; TotalRep; RUN_ID; AssessTime");
				}
				else
				{
					fprintf(fdetail, "\nREP; n_0; EST_GAP; SVAR; xhat; bound; Proced; Z_star; TotalRep; RUN_ID; AssessTime");

				}
			}
			int i;
			s = r*x1n;
		
			//Print("\n%d;%g;%g;", r+1, gap[r], svar[r]); 
			fprintf(fdetail, "\n%d;%d;%g;%g;", r+1, nk, gap[r], svar[r]);
		
			//Print("("); 
			fprintf(fdetail, "(");
			for (i = 0; i< x1n-1; i++){
				//Print ("%g, ", xhat_T[s+i]);
				fprintf(fdetail, "%g, ", xhat_T[s+i]);
		}
		//i++; 
		//Print("%g);", xhat_T[s+i]);
		fprintf(fdetail, "%g);", xhat_T[s+i]);
		
		//Print("%g", bound[r]); 
		fprintf(fdetail, "%g;", bound[r]);
		
		fprintf(fdetail, "%f;%d;120309;APL1P;%g", Z_star, rep, assesstime[r]);
		
#endif	
		
	}
	
	//--------------------------------------------------------------------------
	//Clean up
	//--------------------------------------------------------------------------
	
	delete xhat;     xhat   = NULL; 
	delete xhat_T;   xhat_T = NULL;
	delete gap;      gap    = NULL; 
	
	/* RS addendum */
#ifdef halfgap
	delete gap_first_half;
	delete gap_second_half;
#endif
	
	delete svar;     svar   = NULL;
	delete bound;    bound  = NULL; 
	delete assesstime; assesstime = NULL;
	
	delete sol_2;    sol_2  = NULL;
	delete sol_3;    sol_3  = NULL;
	delete Scen; Scen = NULL;
	delete Scen2; Scen2 = NULL;
	delete Scen3; Scen3 = NULL;	
	
} 

//****************************************************************************************

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
	
	------------------------------------------------------------------------------
*/

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
	//Print( "\nOPTIMAL VALUE: %20.10E\n", sol->result );
	
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
	
	------------------------------------------------------------------------------
*/

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
	
	------------------------------------------------------------------------------
*/

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
	
	------------------------------------------------------------------------------
*/

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

