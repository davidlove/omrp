/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solver2.cpp
CREATED:            1992.09.29
LAST MODIFIED:		1996.09.21

DEPENDENCIES:       stdtype.h, std_tmpl.h, error.h, solver.h, history.h,
					solvcode.h, std_math.h, smartptr.h, print.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/


#include <math.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __HISTORY_H__
#	include "history.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif
#ifndef __SOLVTOL_H__
#	include "solvtol.h"
#endif


#define RC_UPDATE_ACCURACY			(2.0)
#define MAX_ARTIF_INCREASE			(3)
#define PENALTY_INCREASE_FACTOR		(10.0)
#define PENALTY_INCREASE_LIMIT		(1.0e+3)


/*------------------------------------------------------------------------------

	SOLVE_RESULT Solver::Solve( VerbLevel Verbosity, Long_T IterLimit )

PURPOSE:
	The most important procedure of the Solver class. Solves the given LP
problem using an advanced version of the revised simplex algorithm.
	Solver status and iteration information is output using "Print"
interface. Verbosity of the output depends on global variable "Verbosity".

PARAMETERS:
	VerbLevel Verbosity
		Report verbosity level.

	Long_T IterLimit
		Maximum number of iterations to perform. Default value of -1 is
		transformed into 3 * ( M + N );

RETURN VALUE:
	Function returns the solve status, which is one of:
1.	SR_OPTIMUM				when an optimal solution was found,
2.	SR_INFEASIBLE			when the primal problem is infeasible,
3.	SR_UNBOUNDED			when the primal is unbounded,
4.	SR_CANNOT_FIND_INIT_SOL	when initial solution was not found,
5.	SR_CANNOT_BACKTRACK		when program tried to backtrack and couldn't,
6.	SR_STALLED				when stalling occured,
7.	SR_RUNNING				when iteration limit has been exceeded and the
							problem was not solved.

SIDE EFFECTS:
	None. But obviuosly all 'Solver' class data is used, changed, moved around,
reallocated etc. After all this is what the 'Solve' function is all about.

------------------------------------------------------------------------------*/

SOLVE_RESULT Solver::Solve( VerbLevel Verbosity, Long_T IterLimit )
{
	assert( IterLimit == -1 || IterLimit >= 1 );

	if( !Initialized && !InitializeSolver( Verbosity ) )
		return Status = SR_CANNOT_FIND_INIT_SOL;

	if( IterLimit == -1 )
	{
		IterLimit = 5 * ( M + N );
		assert( IterLimit > 0 );
	}

	Int_T p, pp, q, bound;	// These six data will hold the basic
	Real_T theta, pivot;	// information during the simplex iteration.

	Int_T i;				// Various loop iterator.

	int CheckResidMode;		// Used when checking residuals.
	int UpdateStatus;		// Used for reading the basis update status.
	Int_T DegenIter;
	Int_T RC_ResetCounter;
	Long_T LastRC_ResetIter = -1;
	const Int_T RC_ResetFreq = 200;
	SOLVER_ERROR ErrorCode = SE_OK;

	//--------------------------------------------------------------------------
	//	Each iteration will output it's most important characteristics when
	//	in HIGH verbosity mode.
	//
	//	Output has to fit in 80 columns on a screen. The format will be:
	//		Iteration number:		6 digits + 2 spc				(8)
	//		Column leaving basis:	6 digits + 2 spc				(16)
	//		Column entering basis:	6 digits + 2 spc				(24)
	//		Reduced cost:			+x.xxExxx + 2 spc				(37)
	//		Step length:			+x.xxExxx + spc + [R ] + 2 spc	(48)
	//		Residuals:				+x.xxExxx + spc + [IPD] + 2 spc	(61)
	//		Objective function:		+x.xxxxExxx						(74)
	//
	if( Verbosity >= V_HIGH )
		Print(
			"\n\nSIMPLEX OPTIMIZER:\n"
			"%6s  %6s  %6s  %10s  %10s    %10s    %12s\n",
			"ITER", "IN", "OUT", "RC", "STEP", "RESID", "OBJ"
			);
	else if( Verbosity >= V_LOW )
		Print( "\nSIMPLEX OPTIMIZER INVOKED.\n" );

	ArtifIncrease = 0;

//	H.ResetHistory();
	for( IterCnt = 1, DegenIter = 0, RC_ResetCounter = 0;
		IterCnt <= IterLimit && DegenIter < CYCLE_CNT; IterCnt++ )
	{
#ifndef NDEBUG
//		CheckA2B_Consistency();
#endif

		//----------------------------------------------------------------------
		//	Call pricing and optimum search procedure. The pricing will either
		//	determine which column is to be inserted into the basis, or will
		//	declare the current solution to be optimal.
		InGoingColumnNumber( q, Verbosity );
		if( q == SLV_OPTIMUM )
			{ Status = SR_OPTIMUM; goto Epilogue; }
		else if( q == SLV_PROBLEM_INFEASIBLE )
			{ Status = SR_INFEASIBLE; goto Epilogue; }

#ifndef NDEBUG
		const Real_T ftol	= 10 * FEASIBILITY_TOL_DEF,
			utol			= u[q] - ftol;
#endif
		assert(
			( A2B[q] == A2B_LO && z[q] < -OPTIMALITY_TOL && x[q] <= ftol ) ||
			( A2B[q] == A2B_UP && z[q] > OPTIMALITY_TOL  && x[q] >= utol ) ||
			( A2B[q] == A2B_IN && fabs( z[q] ) > OPTIMALITY_TOL  &&
				( !( VarType[q] & VT_LO ) || x[q] > ftol ) && x[q] < utol ) ||
			( VarType[q] & VT_FIXED )
			);

		//----------------------------------------------------------------------
		//	Compute direction vector:
		//	1.	copy q-th column of the constraint matrix into 'w1' work vector,
		//	2.	multiply the resulting column by basis inverse.
		//
		{
			Ptr<Real_T> A;			// These two data members are used for
			Ptr<Int_T> Ind;			// reading columns of the constraint matrix.

			LP.GetColumn( q, A, Ind, w1Len );

			for( i = 0; i < w1Len;  ++i, ++A, ++Ind )
			{
				w1[i]		= *A;
				w1Ind[i]	= *Ind;
			}
			B->SparseFTRAN( w1, w1Ind, w1Len );
		}

		//----------------------------------------------------------------------
		//	Compute reduced cost and weight of variable 'q' and store it in
		//	local variables. If the freshly computed reduced cost indicates,
		//	that introducing this column into basis will not be profitable,
		//	we will recompute reduced costs and go to the beginnig of the loop.
		//
		//	If we are moving the in-coming variable towards its lower bound
		//	we multiply the direction vertor by (-1).
		//
		Real_T z_q = LP.GetC( q ),
			gamma_q = ( Pricing >= PRS_SE ) ? 1.0 : 0.0;

		for( i = 0; i < w1Len; i++ )
		{
			Real_T ww = w1[i];
			if( Pricing >= PRS_SE )
				gamma_q += ww * ww;

			Real_T cc = LP.GetC( B2A[w1Ind[i]] );
			z_q -= cc * ww;
		}

		assert( Pricing < PRS_SE || gamma_q > 0.0 );

		if( z_q > 0 )
			for( i = 0; i < w1Len; i++ )
				w1[i] = -w1[i];

		Bool_T UpdateRC		= False,
			UpdateResult	= False;

		if( ++RC_ResetCounter > RC_ResetFreq )
			UpdateRC = True;
		else if( VarType[q] & ( VT_FIXED | VT_ARTIF ) )
		{
			if( ( A2B[q] == A2B_LO && z_q > OPTIMALITY_TOL ) ||
				( A2B[q] == A2B_UP && z_q < -OPTIMALITY_TOL ) )
				UpdateRC = True;
		}
		else if( ( A2B[q] == A2B_LO && z_q > -OPTIMALITY_TOL ) ||
				( A2B[q] == A2B_UP && z_q < OPTIMALITY_TOL ) ||
				( A2B[q] == A2B_IN && fabs( z_q ) < OPTIMALITY_TOL ) )
			UpdateRC = True;
		else if( ( fabs( z_q ) < fabs( z[q] / RC_UPDATE_ACCURACY ) ) ||
			( fabs( z_q ) > fabs( z[q] * RC_UPDATE_ACCURACY ) ) )
			UpdateRC = True;
		else if( ArtifIncrease >= MAX_ARTIF_INCREASE )
		{
			if( LP.GetPenalty() < PENALTY_INCREASE_LIMIT )
			{
				LP.SetPenalty( Min( LP.GetPenalty() *
					PENALTY_INCREASE_FACTOR, PENALTY_INCREASE_LIMIT ) );
				PenaltyAdjustCnt++;
				UpdateResult = True;
				if( Verbosity >= V_LOW )
					Print( "Penalty increased to %12.4E.\n", LP.GetPenalty() );
			}

			ArtifIncrease = 0;
			UpdateRC = True;
		}

		if( UpdateRC )
		{

			if( LastRC_ResetIter == IterCnt )
				FatalError( "Pricing loop!!!" );

			ComputeDualVariables();
			CheckResidMode = CHK_DUAL;
			if( CheckResiduals( CheckResidMode ) > ALARM_RESID &&
				!HandleNumericalDifficulties( SE_DUAL_RESID_ALARM, Verbosity ) )
				return Status = SR_CANNOT_BACKTRACK;
			ComputeReducedCosts();
			if( UpdateResult )
				ComputeResult();

			RC_ResetCounter = 0;
			IterCnt--; RC_FaultCnt++;
			continue;
		}

		//----------------------------------------------------------------------
		//	Pivoting procedure will determine which basic variable (if any)
		//	is to be removed from the basis. Other possible results are:
		//	1.	problem is unbounded,
		//	2.	non-basic variable number 'q' will reach one of its bounds and
		//		the basis will remain unchanged.
		//
		pivot = OutGoingColumnNumber( q, z_q, p, theta, bound );
		if( p == SLV_UNBOUNDED )
			if( VerifyUnboundedness( q ) )
			{
				Status = SR_UNBOUNDED;
				goto Epilogue;
			}
			else
				continue;
		Result -= fabs( z_q ) * theta;

		if( !HandleNumericalDifficulties( UpdatePrimalVars( q, z_q, theta ),
			Verbosity ) )
			return Status = SR_CANNOT_BACKTRACK;

		//----------------------------------------------------------------------
		//	See if we have made a non-zero step. If not - increment degenerate
		//	step counter 'DegenIter'. If appropriate, try to take action: switch
		//	to steepest edge pricing, or even terminate.
		//
		if( theta < MIN_STEP_LENGTH )
		{
			if( ++DegenIter > DEGEN_CNT && Pricing == PRS_RC )
			{
				Pricing = PRS_ASE;
				if( Verbosity >= V_LOW )
					Print( "\nDegeneration: Forced switch to approx. steepest "
						"edge pricing.\n" );
			}
			else if( DegenIter > CYCLE_CNT )
				return Status = SR_STALLED;
		}
		else
			DegenIter = 0;

		if( Verbosity >= V_HIGH )
			Print( "%6ld  %6ld  ", (long) IterCnt, (long) q );
		else if( Verbosity == V_LOW && IterCnt % 25 == 0 )
			Print( "Iter: %8d\tResult: %10.2E\n",
				(int) IterCnt, (double) Result );

		if( p == SLV_TO_BND )
		{	
			if( Verbosity == V_HIGH )
				Print( "%6s  %10.2E  %10.2E %15s  %12.4E\n",
					"", (double) z_q, (double) theta, "", (double) Result );

			A2B[q] = ( z_q < 0.0 ) ? A2B_UP : A2B_LO;
			z[q] = z_q;
			if( Pricing >= PRS_SE )
				gamma[q] = gamma_q;

			if( ( VarType[q] & VT_ARTIF ) && ( A2B[q] == VT_LO ) )
			{
				LP.FixLambda( q );
				VarType[q] |= VT_FIXED;
				x[q] = u[q] = 0.0;
			}
			continue;
		}
		//
		//	If a non-basic variable just moved between its bounds, this is the
		//	end of the current iteration. Otherwise there's still some work to
		//	do.
		//----------------------------------------------------------------------

		assert( IsNonZero( pivot ) );
		assert( bound != A2B_UNDEF );

		pp = B2A[p];

		if( Verbosity >= V_HIGH )
		{
			if( theta == 0.0 )
				Print( "%6ld  %10.2E  %10s ",
					(long) pp, (double) z_q, "" );
			else
				Print( "%6ld  %10.2E  %10.2E ",
					(long) pp, (double) z_q, (double) theta );
		}

		//----------------------------------------------------------------------
		//	The basis will be changed! Continue with the standard iteration.
		//	Store basis change on the history list. Then update hash tables.
		//
		//	When appropriate fix artificial out-going variable at zero.
		//
//		H.PushStep( pp, q );
		A2B[pp] = bound;
		B2A[p] = q;
		A2B[q] = p;

		if( ( VarType[pp] & VT_ARTIF ) && A2B[pp] == A2B_LO )
		{
			LP.FixLambda( pp );
			VarType[pp] |= VT_FIXED;
			x[pp] = u[pp] = 0.0;
		}

		//----------------------------------------------------------------------
		//	If we are performing steepest edge pricing, we will need to multiply
		//	the direction vector by basis transpose inverse.
		//	We unpack 'w1' using 'w3' and create 'mark' vector in 'w2Ind'.
		//	Compute a new work vector: multiply 'w1' by basis transpose inverse.
		//
		if( Pricing == PRS_SE )
		{
			Array<Int_T> &mark = w2Ind;

			for( i = 0; i < w1Len; i++ ) w3[i] = w1[i];
			w1.Fill( 0.0, M );
			mark.Fill( 0, M );
			for( i = 0; i < w1Len; i++ )
			{
				w1[ w1Ind[i] ]		= w3[i];
				mark[ w1Ind[i] ]	= 1;
			}

			B->SparseBTRAN( w1, mark );
		}

		//----------------------------------------------------------------------
		//	Update basis representation.
		//
		UpdateStatus = UpdateBasis( p );
		if( Verbosity >= V_HIGH )
		{
			char c = ' ';
			switch( UpdateStatus )
			{
			case SLV_BASIS_REFACTORIZED:	c = 'R'; break;
			case SLV_BASIS_UPDATED:			c = ' '; break;
			case SLV_UPDATE_FAILED:			c = 'F'; break;
			}
			Print( "%c  ", c );
		}

		if( UpdateStatus == SLV_UPDATE_FAILED &&
			!HandleNumericalDifficulties( SE_SINGULAR_BASIS, Verbosity ) )
		{
			if( Verbosity == V_LOW )
				Print( "\n" );
			return Status = SR_CANNOT_BACKTRACK;
		}

		//----------------------------------------------------------------------
		//	Once every "RESID_CHECK_FREQ" check residuals.
		//
		if( !PeriodicalResidCheck( ErrorCode, Verbosity ) &&
			!HandleNumericalDifficulties( ErrorCode, Verbosity ) )
		{
			if( Verbosity == V_LOW )
				Print( "\n" );
			return Status = SR_CANNOT_BACKTRACK;
		}
		if( Verbosity >= V_HIGH )
			Print( "%12.4E\n", (double) Result );

		//----------------------------------------------------------------------
		//	Update reduced costs and their weights (according to Goldfarb's
		//	"steepest edge" formulas.
		//
		UpdateReducedCosts( p, pp, z_q, gamma_q, pivot );
		continue;
	}
	//
	//	End of simplex algorithm loop.
	//
	//--------------------------------------------------------------------------

	Status = ( IterCnt > IterLimit ) ? SR_RUNNING : SR_STALLED;

	//--------------------------------------------------------------------------
	//	Report the result from the subproblem solver.
	//
Epilogue:
	IterCnt--;
	TotalIterCnt += IterCnt;

	if( Verbosity == V_LOW ) Print( "\n" );
	ComputeResult();

	PrimalResiduals			= 0.0,
	BoxConstraintViolation	= 0.0;

	switch( Status )
	{
	case SR_OPTIMUM:
	case SR_INFEASIBLE:
	case SR_UNBOUNDED:
		{
			int Resid;

			PrimalResiduals			= CheckResiduals( Resid = CHK_PRIM );
			BoxConstraintViolation	= CheckResiduals( Resid = CHK_INF );
		}
		break;

	case SR_CANNOT_FIND_INIT_SOL:
	case SR_CANNOT_BACKTRACK:
	case SR_STALLED:
	case SR_RUNNING:
	case SR_OK:
		break;

	case SR_UNKNOWN:
	case SR_UNINITIALIZED:
#ifndef NDEBUG
		abort();
#endif
		break;
	}

	return Status;
}


/*------------------------------------------------------------------------------

	Bool_T Solver::InitializeSolver( void )

PURPOSE:
	This function initializes the simplex solver. It::
	1.	finds the initial basis (calling the crash procedure),
	2.	factorizes the basis matrix,
	3.	computes the initial solution 'x' (which is likely to be infeasible),
	4.	checks the primal residuals (if the residuals are bigger than
		"SATISF_RESID" retries basis construction with tighter numerical
		tolerance),
	5.	computes the dual variables 'y',
	6.	checks the dual residuals (if the residuals are bigger than
		"SATISF_RESID" retries basis construction with tighter numerical
		tolerance),
	7.	computes the infeasibility as 'lambda = b - B x_B'
	8.	adds a vector to the constraint matrix, so that the solution becomes
		feasible,
	9.	computes reduced costs and fills the steepest edge weight vector with
		ones.

	Thus we obtain a feasible basic solution and an initial set of reduced costs
	and weights for the steepest edge algorithm.

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	As residuals are computed by means of using 'Solver::CheckResiduals()'
	function, work vector 'w3' is overwritten. Function
	'Solver::MakeSolutionFeasible()' uses (and overwrites) vectors 'w1' and
	'w1Ind' and variable 'w1Len'. It may also result in problem dimension change
	(when new columns are added to the problem).

------------------------------------------------------------------------------*/

Bool_T Solver::InitializeSolver( VerbLevel Verbosity, // )
	FSP FindStartingPoint, FIB FindInitialBasis, Bool_T BasisFirst )
{
	FeasibilityRestoreCount = 0;
	ResetTolerances( Verbosity );
	N = Int_T( LP.GetStructN() + LP.GetSlackN() );

	//--------------------------------------------------------------------------
	//	When NULL pointers are given (e.g. by default parameter values) ---
	//	the default actions are taken.
	//
	if( FindStartingPoint == NULL )
		if( UseExternalSolution )
			FindStartingPoint	= &Solver::EXTERNAL_STARTING_POINT;
		else
			FindStartingPoint	= &Solver::DEFAULT_STARTING_POINT;

	if( FindInitialBasis == NULL )
//		if( UseExternalBasis )
//			FindInitialBasis	= &Solver::EXTERNAL_INITIAL_BASIS;
//		else
			FindInitialBasis	= &Solver::DEFAULT_INITIAL_BASIS;

	Int_T j;

	//--------------------------------------------------------------------------
	//	Initialize vector of variable types 'VarType' and upper bounds.
	//
	for( j = 0; j < N; j++ ) VarType[j] = LP.GetVarType( j );
	for( j = 0; j < N; j++ ) u[j] = LP.GetU( j );

	VarType.Fill( VT_NORM,	AllocN, N );
	u.Fill( +INFINITY,		AllocN, N );
	A2B.Fill( A2B_LO,		AllocN, N );

	//--------------------------------------------------------------------------
	//	Look for a stable initial basis in loop. Exit the loop if you cannot
	//	obtain satisfactory numerical stability even when pivot tolerance is
	//	nearly equal to 1.
	//
	for( LU_PIVOT_TOL = 0.1; LU_PIVOT_TOL < 0.99;
		LU_PIVOT_TOL = Min( 0.9999, 3 * LU_PIVOT_TOL ) )
	{
		if( BasisFirst )
		{
			if( (this->*FindInitialBasis)( LU_PIVOT_TOL, Verbosity ) &&
				(this->*FindStartingPoint)( Verbosity ) )
				break;
		}
		else
		{
			if( (this->*FindStartingPoint)( Verbosity ) &&
				(this->*FindInitialBasis)( LU_PIVOT_TOL, Verbosity ) )
				break;
		}
	}

	//--------------------------------------------------------------------------
	//	Now is the time to reduce infeasibility by means of adding super
	//	artificial columns. After that we may compute the objective function.
	//
	MakeSolutionFeasible( FEASIBILITY_TOL, Verbosity, True );

	Real_T MinM, MaxM;

	LP.GetPenaltyEstimates( MinM, MaxM );

	if( LP.GetPenalty() < MaxM * 2.0 )
		LP.SetPenalty( MaxM * 2.0 );

	//--------------------------------------------------------------------------
	//	Initialize data for steepest edge pricing: weights 'gamma' and reduced
	//	cost vector 'z'.
	//
	ComputeResult();
	ComputeDualVariables();
	ComputeReducedCosts();

	//
	//	Now we have a feasible starting solution. All solver data is properly
	//	initialized.
	//
	return Initialized = True;
}


/*------------------------------------------------------------------------------

	SOLVER_ERROR Solver::UpdatePrimalVars( Int_T q, Real_T z_q, Real_T theta )

PURPOSE:
	This function updates the values of the primal variables 'x' after making a
step in direction given in sparse vector ('w1', 'w1Ind', 'w1Len') of length
'theta'. Depending on the sign of reduced cost 'z_q' of the non-basic variable
'x[q]' changes value by plus or minus 'theta'.
	In order to ensure that the variable simple bounds are not violated,
residuals are checked.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	Error code 'SE_INF_RESID' if the infeasibility has occured, or 'SE_OK'
success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SOLVER_ERROR Solver::UpdatePrimalVars( Int_T q, Real_T z_q, Real_T theta )
{
	Int_T i;

	//--------------------------------------------------------------------------
	//	Solution update: the basic variables.
	//
	for( i = 0; i < w1Len; i++ )
	{
		const Int_T j = w1Ind[i];
		Real_T &xx = x[ B2A[j] ];
		
		//----------------------------------------------------------------------
		//	Fix artificial variables that were reduced to zero.
		//
		xx -= theta * w1[i];
		if( VarType[j] & VT_ARTIF )
		{
			if( IsNegative( w1[i] ) )	// Increase an artificial variable.
				ArtifIncrease++;
				
			if( IsZero( xx ) )
			{
				LP.FixLambda( j );
				VarType[j] |= VT_FIXED;
				xx = u[j] = 0.0;
			}
		}
		else if( IsZero( xx ) )
			xx = 0.0;
	}

	//--------------------------------------------------------------------------
	//	Solution update: the non-basic variable.
	//
	if( IsZero( x[q] += ( z_q < 0.0 ) ? theta : -theta ) ) x[q] = 0.0;
	if( VarType[q] & VT_ARTIF && ( z_q < 0.0 ) )
		ArtifIncrease++;

	//--------------------------------------------------------------------------
	//	In DEBUG mode only!
	//
	//	Check if in feasibility restoring mode (i.e. when
	//	FeasibilityRestoreCount > 0 ) whether any artificial variables are below
	//	zero.
	//
/*
#ifndef NDEBUG
	if( FeasibilityRestoreCount > 0 )
	{
		for( Int_T j = LP.GetStructN() + LP.GetSlackN(); j < N; j++ )
			if( VarType[j] & VT_ARTIF )
			{
				if( x[j] < -1.1*FEASIBILITY_TOL )
					Print( "Negative artificial: x[%d] = %14.10e\n",
						(int)j, x[j] );
			}
	}
#endif
*/

	return SE_OK;
}


/*------------------------------------------------------------------------------

	void Solver::UpdateReducedCosts( Int_T p, Int_T pp, Real_T z_q,
		Real_T gamma_q, Real_T pivot )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

void Solver::UpdateReducedCosts( Int_T p, Int_T pp, Real_T z_q, // )
	Real_T gamma_q, Real_T pivot )
{
	//--------------------------------------------------------------------------
	//	Compute dense work vector 'w2' as a product of 'e_p' and basis inverse.
	//
	Array<Int_T> &mark = w1Ind;
	w2.Fill( 0.0, M );		w2[p]	= 1.0;
	mark.Fill( 0, M );		mark[p]	= 1;
	B->SparseBTRAN( w2, mark );

	//--------------------------------------------------------------------------
	//	Update all reduced costs and weights. Use work vectors 'w1' and
	//	'w2' (both are stored in dense form by now).
	//

	//
	//	Pass 1.	Compute work vectors 'alpha' and 'beta'.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T slCol, laCol, len;
	Real_T sl, la, x;

	alpha.Fill( 0.0,	AllocN );
	beta.Fill( 0.0,		AllocN );

	for( Int_T i = 0; i < M; i++ )
		if( mark[i] )
		{
			x = w2[i];
			LP.GetRow( i, a, col, len, slCol, sl, laCol, la );
			for( ; len; --len, ++a, ++col )
				if( A2B[ *col ] < 0 )
					alpha[ *col ] += x * *a;

			if( slCol >= 0 ) alpha[ slCol ] += x * sl;
			alpha[ laCol ] += x * la;
		}

	if( Pricing == PRS_SE )
	{
		beta.Fill( 0.0,		AllocN );

		for( Int_T i = 0; i < M; i++ )
			if( IsNonZero( x = w1[i] ) )
			{
				LP.GetRow( i, a, col, len, slCol, sl, laCol, la );
				for( ; len; --len, ++a, ++col )
					if( A2B[ *col ] < 0 )
						beta[ *col ] += x * *a;

				if( slCol >= 0 ) beta[ slCol ] += x * sl;
				beta[ laCol ] += x * la;
			}
	}

	//
	//	Pass 2.	Update all reduced costs (and perhaps SE weights).
	//
	for( Int_T j = 0; j < N; j++ )
		if( A2B[j] >= 0 )
		{
			z[j] = 0;
			if( Pricing >= PRS_SE )
				gamma[j] = 1.0;				// The value here is not important.
		}
		else if( j == pp )
		{
			z[pp] = fabs( z_q ) / pivot;
			if( Pricing >= PRS_SE )
				gamma[pp] = gamma_q / ( pivot * pivot );
		}
		else
		{
			z[j] -= alpha[j] * z_q;

			if( IsNonZero( alpha[j] ) && Pricing >= PRS_SE )
			{
				Real_T sqr_alpha_j = alpha[j] * alpha[j];

				switch( Pricing )
				{
				case PRS_SE:
					gamma[j] = Max( gamma[j] - 2*alpha[j]*beta[j] +
						sqr_alpha_j * gamma_q, 1 + sqr_alpha_j );
					break;

				case PRS_ASE:
					{
						Real_T old_alpha_sqr_alpha_j = sqr_alpha_j *
							pivot * pivot;

						gamma[j] = Max( gamma[j], old_alpha_sqr_alpha_j + 1.0 )
							- 2.0 * old_alpha_sqr_alpha_j
							+ sqr_alpha_j * gamma_q;
					}
					break;

				case PRS_RC:
					assert( Pricing >= PRS_SE );
				}

				assert( gamma[j] + 1e-6 >= sqr_alpha_j  + 1.0 );
			}
		}

/*
	//
	//	Now: a rigorous check of all reduced costs.
	//
	{
		WorkVector<Real_T> z_copy( AllocN );

		z_copy.Copy( z, AllocN, AllocN );

		ComputeDualVariables();
		ComputeReducedCosts();

		Real_T MaxDiff = 0.0;

		for( Int_T j = 0; j < N; j++ )
			MaxDiff = Max( MaxDiff, fabs( z[j] - z_copy[j] ) );

		if( IsNonZero( MaxDiff ) )
			Print( "Maximum reduced cost error after update: %12.2E.\n",
				MaxDiff );
	}
*/
}


/*------------------------------------------------------------------------------

	Bool_T Solver::PeriodicalResidCheck( SOLVER_ERROR &ErrorCode,
		VerbLevel Verbosity )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Bool_T Solver::PeriodicalResidCheck( SOLVER_ERROR &ErrorCode, // )
	VerbLevel Verbosity )
{
	//--------------------------------------------------------------------------
	//	Act only every 'RESID_CHECK_FREQ' calls. During the remaining calls just
	//	return.
	//
	if( ++ResidCheckClock < RESID_CHECK_FREQ )
	{
		ErrorCode = SE_OK;
		if( Verbosity >= V_HIGH ) Print( "%10s    ", "" );
		return True;
	}

	ResidCheckClock = 0;

	Real_T Resid = 0.0;
	int CheckResidMode;

	//--------------------------------------------------------------------------
	//	Compute primal residuals.
	//
	CheckResidMode = CHK_PRIM | CHK_INF;
	Resid = CheckResiduals( CheckResidMode );

	if( Verbosity >= V_HIGH )
	{
		char c = ' ';
		switch( CheckResidMode )
		{
		case CHK_PRIM:	c = 'P'; break;
		case CHK_INF:	c = 'I'; break;
		}

		//----------------------------------------------------------------------
		//	Print the computed residuals, or space instead.
		//
		if( IsPositive( Resid ) )
			Print( "%10.2E %c  ", Resid, c );
		else
			Print( "%10s    ", "" );
	}
/*
	if( Resid > 10 * FEASIBILITY_TOL )
		Print( "Iter %4d: resid %10.2E\n", (int)IterCnt, Resid );
*/

	//--------------------------------------------------------------------------
	//	In case of extremely large residuals set "ErrorCode" to appropriate
	//	value. Otherwise set it to 'SE_OK'.
	//
	if( Resid > ALARM_RESID )
	{
		switch( CheckResidMode )
		{
		case CHK_PRIM:	ErrorCode = SE_PRIM_RESID_ALARM;	break;
		case CHK_INF:	ErrorCode = SE_INF_RESID_ALARM;		break;
		}
		return False;
	}
	else if( Resid > POOR_RESID )
	{
		switch( CheckResidMode )
		{
		case CHK_PRIM:	ErrorCode = SE_PRIM_RESID_POOR;		break;
		case CHK_INF:	ErrorCode = SE_INF_RESID_POOR;		break;
		}
		return False;
	}

	ErrorCode = SE_OK;
	return True;
}


/*------------------------------------------------------------------------------

	SOLVE_RESULT Solver::RestartAndSolve( VerbLevel Verbosity,
		Long_T IterLimit, Bool_T ComputeDuals )
	SOLVE_RESULT Solver::RestartAndSolve( VerbLevel Verbosity,
		const SolverStateDump *dump, Long_T IterLimit )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

SOLVE_RESULT Solver::RestartAndSolve( VerbLevel Verbosity, Long_T IterLimit, //)
	Bool_T ComputeDuals )
{
	//--------------------------------------------------------------------------
	// Solver ought to have been previously initialized.
	//
	assert( Initialized );

	Int_T StructN = N = Int_T( LP.GetStructN() + LP.GetSlackN() );

	FeasibilityRestoreCount = 0;
	ResetTolerances( Verbosity );

	//--------------------------------------------------------------------------
	//	Re-read vectors of variable types 'VarType' and upper bounds, which may
	//	have changed since the previous solution..
	//
	Int_T i;
	for( i = 0; i < N; i++ ) VarType[i] = LP.GetVarType( i );
	VarType.Fill( VT_NORM, AllocN, N );
	for( i = 0; i < N; i++ ) u[i] = LP.GetU( i );
	u.Fill( +INFINITY, AllocN, N );

	//--------------------------------------------------------------------------
	//	Now is the time to reduce infeasibility by means of adding super
	//	artificial columns. After that we may compute the objective function.
	//
	//	Initialize data for steepest edge pricing: weights 'gamma' and reduced
	//	cost vector 'z'.
	//
	MakeSolutionFeasible( FEASIBILITY_TOL, Verbosity, False );
	ComputeResult();
	if( ComputeDuals )
	{
		ComputeDualVariables();
		ComputeReducedCosts();
	}
	else
		ComputeReducedCosts( StructN );

	//
	//	Now we have a feasible starting solution. All solver data is properly
	//	initialized.
	//
	return Solve( Verbosity, IterLimit );
}


SOLVE_RESULT Solver::RestartAndSolve( VerbLevel Verbosity, // )
	const SolverStateDump *dump, Long_T IterLimit )
{
	assert( dump != NULL );

	//--------------------------------------------------------------------------
	//	Read data from the "dump" structure.
	//
	x.Copy( dump->X, AllocN, AllocN );
	A2B.Copy( dump->A2B, AllocN, AllocN );

	//--------------------------------------------------------------------------
	//	Restore the B2A array. Factorize the basis.
	//
	B2A.Fill( B2A_UNDEF, M );
	Int_T i, j;
	for( j = 0, i = 0; j < N; j++ )
		if( A2B[j] >= 0 )
		{
			assert( B2A[A2B[j]] == B2A_UNDEF );
			B2A[A2B[j]] = j;
			i++;
		}

	assert( i == M );
	UpdateBasis();

	//--------------------------------------------------------------------------
	//	We have restored a solution stored in the "dump" structure. Now we
	//	can run the solver.
	//
	return RestartAndSolve( Verbosity, IterLimit, True );
}
