/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solver3.cpp
CREATED:			1993.11.01
LAST MODIFIED:		1996.09.16

DEPENDENCIES:		error.h, std_tmpl.h, stdtype.h, solver.h, solv_lp.h,
					lp_codes.h, mps_lp.h, inverse.h, solvcode.h, simplex.h,
					std_math.h, smartptr.h, print.h

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

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
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
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T Solver::FactorizeBasis( void )

PURPOSE:
	Performs an unconditional factorisation of the current simplex basis.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status ('True' on success).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Solver::FactorizeBasis( void )
{
	int code;
	Ptr<Real_T> A;
	Ptr<Int_T> Row;
	Int_T Len;

	//--------------------------------------------------------------------------
	//	Empty the inverse data structure. Fill it with the columns of the
	//	current basis. Factorize.
	//
	B->Clear();
	for( Int_T i = 0; i < M; i++ )
	{
		LP.GetColumn( B2A[i], A, Row, Len );
		B->AddCol( i, A, Row, Len );
	}
	code = B->Factor();

	if( code != 1 )	return False;

	//--------------------------------------------------------------------------
	//	Remember the number of non-zeros and the largest non-zero absolute
	//	value.
	//
	Alen = Ulen = B->GetFactorLen();
	Amax = Umax = B->GetUMax();

	return True;
}


/*------------------------------------------------------------------------------

	void Solver::FixArtificialVariables( void )

PURPOSE:
	Fix (at zero, of course) any artificial variables that have been reduced to
zero (with respect to feasibility tolerance, of course).

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::FixArtificialVariables( void )
{
	for( Int_T j = Int_T( LP.GetStructN() + LP.GetSlackN() ); j < N; j++ )
		if( ( VarType[j] & VT_ARTIF ) && fabs( x[j] ) < FEASIBILITY_TOL )
		{
			VarType[j] |= VT_FIXED;
			x[j] = u[j] = 0.0;
			LP.FixLambda( j );
		}
}


/*------------------------------------------------------------------------------

	void Solver::ComputePrimalVariables( void )

PURPOSE:
	Compute afresh the vector of primal variables 'x' in two steps (shown in
comments included in the function body).

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ComputePrimalVariables( void )
{
	Int_T i;
	Ptr<Real_T> A;
	Ptr<Int_T> Row;
	Int_T Len;

	//----------------------------------------------------------------------
	//	Calculate b' = b - N * xN. Use 'w3' work vector.
	//
	LP.GetRHS( w3 );

	for( i = 0; i < N; i++ )
		if( A2B[i] < 0 && IsNonZero( x[i] ) )
		{
			for( LP.GetColumn( i, A, Row, Len ); Len; --Len, ++A, ++Row )
				w3[ *Row ] -= x[i] * *A;
		}

	for( i = 0; i < M; i++ )
		if( IsZero( w3[i] ) ) w3[i] = 0.0;

	//----------------------------------------------------------------------
	//	Calculate x = B^(-1)b'.
	//
	B->DenseFTRAN( w3 );
	for( i = 0; i < M; i++ ) x[ B2A[i] ] = w3[i];

	PrimVarComputeCnt++;
}


/*------------------------------------------------------------------------------

	void Solver::ComputeDualVariables( Bool_T CountArtif )
	void Solver::ComputeDualVariables( Real_T penalty );

PURPOSE:
	The first version computes the vector of dual variables by storing the
basic costs in the vector 'y' and then performing a dense BTRAN (i.e.
y := B^(-T) * c_B).
	The second version uses the (previously computed) split dual variables
to calculate the standard dual var. vector.

PARAMETERS:
	Bool_T CountArtif
		If "True", the prices by the artificial variables will be counted.

	Real_T penalty
		The current value of penalty.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ComputeDualVariables( Bool_T CountArtif )
{
	if( CountArtif )
	{
		for( Int_T i = 0; i < M; i++ )
			y[i] = LP.GetC( B2A[i] );
	}
	else
	{
		for( Int_T i = 0; i < M; i++ )
			y[i] = ( VarType[B2A[i]] & VT_ARTIF ) ? 0.0 : LP.GetC( B2A[i] );
	}

	B->DenseBTRAN( y );

	DualVarComputeCnt++;
}


void Solver::ComputeDualVariables( Real_T penalty )
{
	for( Int_T i = 0; i < M; i++ )
		y[i] = y_x[i] + penalty * y_t[i];
}


/*------------------------------------------------------------------------------

	void Solver::ComputeResult( void )

PURPOSE:
	Computes the current value of the objective function as: obj := c^T * x

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ComputeResult( void )
{
	Result = LP.GetF();
	for( Int_T i = 0; i < N; i++ )
		Result += LP.GetC( i ) * x[i];
}


/*------------------------------------------------------------------------------

	void Solver::ComputeReducedCosts( void )
	void Solver::ComputeReducedCosts( Int_T Start )

PURPOSE:
	Assuming that the dual variables ('y') are already computed, we compute in
loop the values of reduced costs of the primal variables. The two versions of
the function are given for functionality and speed.
	One always calculates ALL reduced costs. It uses the row form of constraint
matrix and scans only rows to which non-zero dual variables correspond.
	The other takes one argument (see below) and computes reduced costs only for
some of the variables. It is used only when very few reduced costs need to be
recomputed. In fact, it is called during the restarts of the algorithm, when
only the reduced costs corresponding to the artificial variables need
refreshing.

PARAMETERS:
	Int_T Start
		The starting column for RC calculation.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ComputeReducedCosts( void )
{
	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Real_T sl, la;
	Int_T slCol, laCol, len, i;

	for( i = 0; i < N; i++ )
		z[i] = LP.GetC( i );

	for( i = 0; i < M; i++ )
	{
		Real_T yy = y[i];

		if( IsNonZero( yy ) )
		{
			LP.GetRow( i, a, col, len,  slCol, sl, laCol, la );
			for( ; len; --len, ++a, ++col )
				z[ *col ] -= *a * yy;
			if( slCol >= 0 ) z[ slCol ] -= sl * yy;
			z[ laCol ] -= la * yy;
		}
	}

	for( i = 0; i < N; i++ )
		if( IsZero( z[i] ) )
			z[i] = 0.0;

	ResetWeights();

	SE_ResetCnt++;
}


void Solver::ComputeReducedCosts( Int_T Start )
{
	Ptr<Real_T> A;
	Ptr<Int_T> Row;
	Int_T Len;

	for( Int_T j = Start; j < N; j++ )
	{
		Real_T rc = LP.GetC( j );

		for( LP.GetColumn( j, A, Row, Len ); Len;  --Len, ++A, ++Row )
			rc -= *A * y[ *Row ];
		z[j] = ( IsNonZero( rc ) ) ? rc : 0.0;
	}

	ResetWeights();

	SE_ResetCnt++;
}


/*------------------------------------------------------------------------------

	void Solver::ResetWeights( void )

PURPOSE:
	Reset steepest edge or approximate steepest edge weights.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ResetWeights( void )
{
	if( Pricing >= PRS_SE )
	{
		Ptr<Real_T> a;
		Ptr<Int_T> row;
		Int_T len;

		for( Int_T j = 0; j < N; j++ )
		{
			LP.GetColumn( j, a, row, len );
			gamma[j] = (double)len + 1.0;
		}
	}
}


/*------------------------------------------------------------------------------

	Real_T Solver::CheckResiduals( int &Mode )

PURPOSE:
	Calculates infeasibility, primal and dual residuals (depending on the
parameter 'Mode', which is a bit mask). In case of overly large residuals
outputs the report (this is what the 'Iter' and 'Result' may be needed for).
	The first stage of the check that reveals large inaccuracy (larger than or
equal to 'ALARM_RESID') breaks the process. On the other hand individual
residuals smaller than 'GOOD_RESID' are considered to be zero.
	Residual estimator 'R' is calculated as the maximum of the appropriate
residual or infeasibility vector. If none of the three categories of residuals
should trigger alert (and return from the procedure) the maximum of their values
is returned.
	Value of mode on return indicates which type of residual was the largest
(with respect to second norm) or which one trigerred a premature return.

PARAMETERS:
	int &Mode
		On entry:	a bit mask marking what needs to be checked.
		On exit:	a bit mask saying which residuals were biggest (no more than
					one bit shall be set to one). If all residuals are found
					to small to bother, zero is returned.

RETURN VALUE:
	x

SIDE EFFECTS:
	Overwrites 'w3' if primal or dual residuals are computed.

------------------------------------------------------------------------------*/

Real_T Solver::CheckResiduals( int &Mode )
{
	if( Mode & ( CHK_PRIM | CHK_DUAL ) ) ResidCheckCnt++;

	Real_T InfR	= 0.0,						// Infeasibility.
		pR		= 0.0,						// Primal residuals.
		dR		= 0.0,						// dual residuals.
		R		= 0.0;
	Int_T i;
	int RetMode	= 0;

	//--------------------------------------------------------------------------
	//	Phase I: Infeasibility check.
	//		If variables' infeasibilities are greater than ALARM_RESID
	//		other phases are not performed.
	//
	if( Mode & CHK_INF )
	{
		for( i = 0; i < N; i++ )
			if( ( VarType[i] & VT_HAS_LO_BND ) && x[i] < 0.0 )
				InfR = Max( fabs( x[i] ), InfR );
			else if( ( VarType[i] & VT_HAS_UP_BND ) && x[i] > u[i] )
				InfR = Max( fabs( x[i] - u[i] ), InfR );

		if( InfR > ALARM_RESID )
		{
			Mode = CHK_INF;
			return InfR;
		}
	}
	if( IsNonZero( InfR ) )
	{
		R		= InfR;
		RetMode	= CHK_INF;
	}

	//--------------------------------------------------------------------------
	//	Phase II: Primal residuals check.
	//
	Ptr<Real_T> A;
	Ptr<Int_T> Row;
	Int_T Len;

	if( Mode & CHK_PRIM )
	{
		LP.GetRHS( w3 );

		for( i = 0; i < N; i++ )
			for( LP.GetColumn( i, A, Row, Len ); Len; --Len, ++A, ++Row )
				w3[ *Row ] -= x[i] * *A;

		for( i = 0; i < M; i++ )
			pR = Max( pR, fabs( w3[i] ) );

		if( pR > ALARM_RESID )
		{
			Mode = CHK_PRIM;
			return pR;
		}
	}
	if( pR > R )
	{
		R		= pR;
		RetMode	= CHK_PRIM;
	}

	//--------------------------------------------------------------------------
	//	Phase III: Dual residuals check.
	//
	if( Mode & CHK_DUAL )
	{
		for( i = 0; i < M; i++ )
		{
			w3[i] = LP.GetC( B2A[i] );
			for( LP.GetColumn( B2A[i], A, Row, Len ); Len; --Len, ++A, ++Row)
				w3[i] -= y[ *Row ] * *A;
		}

		for( i = 0; i < M; i++ )
			dR = Max( dR, fabs( w3[i] ) );

		if( dR >= ALARM_RESID )
		{
			Mode = CHK_DUAL;
			return dR;
		}
	}
	if( dR > R )
	{
		R		= dR;
		RetMode	= CHK_PRIM;
	}

	Mode = RetMode;
	return R;
}


/*------------------------------------------------------------------------------

	void Solver::MakeSolutionFeasible( Real_T ResidLevel, VerbLevel Verbosity,
		Bool_T MinimizeInfeasibility )

PURPOSE:
	Used for restoring the feasibility of the solution of the original problem
(i.e. a problem without the artificial variables, which are disregarded).
Projects the structural and slack variables on the intervals delimited by their
simple bounds. Then computes the primal residuals and, if necessary, adds (or
replaces) a vector of artificial variables that will contain the infeasibility.
If necessary, the simplex basis is factorized.

PARAMETERS:
	Real_T ResidLevel
		Greatest absolute value of primal residuals accepted for an initial
		solution.

	VerbLevel Verbosity
		Verbosity level.

	Bool_T MinimizeInfeasibility
		If 'True', additional heuristic for infeasibility minimization is run.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Overwrites work vectors 'w1', 'w1Ind', 'w3' (indirectly - by calling
'Solver::CheckResiduals()' and variable 'w1Len'. Possibly changes the dimension
(number of columns) of the constraint matrix.

------------------------------------------------------------------------------*/

void Solver::MakeSolutionFeasible( Real_T ResidLevel, VerbLevel Verbosity, // )
	Bool_T MinimizeInfeasibility )
{
	assert( ResidLevel > 1.0e-16 );

	//--------------------------------------------------------------------------
	//	We do not need to take into account the previous 'Lambda' vector.
	//
	Int_T OrigN = Int_T( LP.GetStructN() + LP.GetSlackN() );

	//--------------------------------------------------------------------------
	//	Move variables which violate their simple bounds to the bound that is
	//	closest to their current position. Adjust 'A2B'. Remember to leave the
	//	basis unchanged!
	//
	Int_T j;
	for( j = 0; j < OrigN; j++ )
		if( A2B[j] < 0 )
			if( VarType[j] & VT_HAS_LO_BND && x[j] < 0.0 )
			{
				x[j] = 0.0;
				A2B[j] = A2B_LO;
			}
			else if( VarType[j] & VT_HAS_UP_BND && x[j] > u[j] )
			{
				x[j] = u[j];
				A2B[j] = A2B_UP;
			}
//	ComputePrimalVariables();
	Int_T i;
	for( i = 0; i < M; i++ )
	{
		Int_T j = B2A[i];

		if( VarType[j] & VT_HAS_LO_BND && x[j] < 0.0 )
			x[j] = 0.0;
		else if( VarType[j] & VT_HAS_UP_BND && x[j] > u[j] )
			x[j] = u[j];
	}

	N = LP.GetN();

	Bool_T ArtifInBasis = False;

	for( j = OrigN; j < N; j++ )
		if( A2B[j] >= 0 )
		{
			ArtifInBasis = True;
			break;
		}

	N = OrigN;

	//--------------------------------------------------------------------------
	//	Compute residual vector (using 'CheckResiduals') in 'w3'. If residuals
	//	are small, exit.
	//
	int CheckResidMode = CHK_PRIM;
	if( CheckResiduals( CheckResidMode ) < ResidLevel )
		goto FeasibleSolution;

	//--------------------------------------------------------------------------
	//	Minimize the lambda vector (in dense form in 'w3' work vector) by means
	//	of shifting primal variables between their simple bounds..
	//
	CheckResidMode = CHK_PRIM;
	if( MinimizeInfeasibility )
		if( MinimizeLambda( w3, Verbosity ) )
			if( CheckResiduals( CheckResidMode ) < ResidLevel )
				goto FeasibleSolution;

	//--------------------------------------------------------------------------
	//	If the residuals were larger than 'ResidLevel' threshold, we must
	//	reduce infeasibilities. We update vector 'x' at the same time.
	//	Ready 'Lambda' vector is then passed to 'SimplexLP' class function.
	//	Finally the dimensions of the problem and the problem data ('x', 'u',
	//	'VarType', 'A2B') are updated.
	//
	for( i = 0, w1Len = 0; i < M; i++ )
		if( IsNonZero( w3[i] ) )
		{
			x[N+i] = fabs( w3[i] );
			w1Len++;
		}
		else
			x[N+i] = 0.0;

	if( w1Len == 0 )
	{
		if( Verbosity >= V_LOW )
			Print( "\nSolution feasible.\n" );
		goto FeasibleSolution;
	}
	else
	{
		if( Verbosity >= V_LOW )
			Print( "\nReducing infeasibilities:\n"
				"\tNo. of infeasibilities: %d\n", (int)w1Len );
	}

	LP.CreateLambda( w3 );

	//--------------------------------------------------------------------------
	//	Scan the artificial variables. Adjust A2B[j], u[j], VarType[j].
	//
	N = LP.GetN();
	for( j = OrigN; j < N; j++ )
	{
		Short_T vt = VarType[j] = LP.GetVarType( j );

		assert( vt & VT_ARTIF );

		u[j] = LP.GetU( j );
		if( vt & VT_FX )
		{
			assert( IsZero( x[j] ) );
			assert( IsZero( u[j] ) );

			if( A2B[j] < 0 ) A2B[j]	= A2B_LO;
		}
		else
		{
			assert( IsNonZero( x[j] ) );
			assert( IsNonZero( u[j] ) );

			if( A2B[j] < 0 ) A2B[j] = ( vt & VT_UP ) ? A2B_UP : A2B_IN;
		}
	}

	if( ArtifInBasis )
		UpdateBasis();

	//--------------------------------------------------------------------------
	//	Now some simple checking.
	//
	assert( CheckResiduals( CheckResidMode = CHK_PRIM | CHK_INF ) <=
		1.1 * ResidLevel );

	if( Verbosity >= V_LOW )
		Print( "Reduction successful.\n" );

	return;

	//--------------------------------------------------------------------------
	//	If the initial solution was found feasible, fix all artificial variables
	//	at zero and proceed.
	//
FeasibleSolution:
	w3.Fill( 0.0, M );
	LP.CreateLambda( w3 );

	u.Fill( 0.0, AllocN, OrigN );
	for( j = OrigN; j < AllocN; j++ )
		if( A2B[j] < 0 )
			A2B[j] = A2B_LO;

	x.Fill( 0.0, AllocN, OrigN );
	VarType.Fill( VT_FIXED | VT_ARTIF, AllocN, OrigN );
	N = LP.GetN();

	if( ArtifInBasis )
		UpdateBasis();
	
	assert( CheckResiduals( CheckResidMode = CHK_PRIM | CHK_INF ) <=
		1.1 * ResidLevel );
}


/*------------------------------------------------------------------------------

	int Solver::UpdateBasis( Int_T p = -1 )

PURPOSE:
	Performs an elementary basis column exchange. The basis is assumed to be
already stored in 'A2B' and 'B2A' index tables.
	If 'p' is equal to -1 a full factorization is forced. This may be trigerred
by various events (from initial basis constructon to numerical difficulties).
If however the calling procedure does not force refactorization by any of those
methods, the procedure itself is free to decide whether to update the basis
representation, or to perform a new LU factorization.
	On return the procedure returns the success / failure flag. In case of
failure the possible cause is reported. 'CB_CHANGED_COLS' and 'CB_REFACTORIZED'
are success codes, while other codes indicate some kind of difficulties (of
course numerical difficulties during factori- zation or update).

PARAMETERS:
	Int_T p
		A number of basic coumn to be removed from the basis and replaces by
		the pivot column (stored during the last FTRAN). If negative, the basis
		is forcibly factorized (see above).

RETURN VALUE:
	A return code, one of:
	SLV_BASIS_REFACTORIZED	- successfull update by full factorization,
	SLV_BASIS_UPDATED		- successfull update by column exchange,
	SLV_UPDATE_FAILED		- update failure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int Solver::UpdateBasis( Int_T p )
{
	//--------------------------------------------------------------------------
	//	There are multiple criteria of choosing if we should update or
	//	refactorize. We refactorize if:
	//	-	both 'p' and 'q' are equal zero (forced refactorization), or
	//	-	any of growth factors exceeds limit, or
	//	-	refactorization frequency is exceeded, or
	//	-	residuals are checked and they exceed preset limit (indicated
	//		by forced refatorization).
	//
	Ulen = B->GetFactorLen();		Umax = B->GetUMax();

	Bool_T Refact = ( ( p < 0 ) || ++UpdateCnt >= REFACT_FREQ ||
		( UpdateCnt >= 5 &&
		( Ulen > Alen * LENGTH_FACTOR || Umax > Amax * GROWTH_FACTOR) ) ) ?
		True : False;

	int Retval = ( Refact ) ? SLV_BASIS_REFACTORIZED : SLV_BASIS_UPDATED;

	if( Refact ) UpdateCnt = 0;

	//--------------------------------------------------------------------------
	//	First attempt whatever action seems suitable.
	//
	if( ( Refact && FactorizeBasis() ) || ( B->Update( p ) == 1 ) )
	{
		goto End;						// Successful update / refactorization.
	}
	//--------------------------------------------------------------------------
	//	In case of update failure, we may still be able to refactorize and thus
	//	recover from the numerical difficulties.
	//
	else if( !Refact && FactorizeBasis() )
	{
		Retval = SLV_BASIS_REFACTORIZED;
		goto End;
	}
	//--------------------------------------------------------------------------
	//	When refactorization failed, there is nothing more we can do in this
	//	procedure. Numerical difficulties status is returned to the calling
	//	function. In some verbosity modes a suitable message is output.
	//
	else
	{
		Retval = SLV_UPDATE_FAILED;
		goto End;
	}

End:
	return Retval;
}


/*------------------------------------------------------------------------------

	Bool_T Solver::HandleNumericalDifficulties( SOLVER_ERROR Problem,
		VerbLevel Verbosity )

PURPOSE:

	!! SOON TO BE OBSOLETED BY AUTOMATIC BASIS RECOVERY DURING FACTORIZATION !!

	This function is supposed to work as a numerical error handler for the main
	function of the class: "Solver::Solve". Whenever "Solver::Solve" encounters
	numerical problems, like large residuals, infeasibility, basis singularity
	etc., it calls this function with appropriate parameter setting.

	This function has two sets of methods for error handling:
	1.	factorization followed by computing the variables from the scratch, when
		the problem does not seem to be extremely difficult and
	2.	repeated backtracking and basis factorization, followed by computation
		of variables, as soon as a stable basis is found.

	Repeated use of "goto" is perhaps unnecessary, but convenient. It allows us
	to recover from numerical difficulties that arise during the function
	without recursion. Any doubts about the meaning of the "goto's" may be
	resolved by drawing a simple flow chart.

PARAMETERS:
	NUM_PROBLEM Problem
		Defines the kind of numerical problem encountered. For possible values:
		see "solvcode.h".

RETURN VALUE:
	Boolean success flag.

SIDE EFFECTS:
	Possible recomputations of primal and dual variables and the objective
	function. Possible change of the basis. Repeated use (and overwriting) of
	work vectors.

------------------------------------------------------------------------------*/

#define BACK_STEPS		( Min( Int_T( M/10 ), Int_T( 15 ) ) )

Bool_T Solver::HandleNumericalDifficulties( SOLVER_ERROR Problem, // )
	VerbLevel Verbosity )
{
	Bool_T ComputeDuals	= False,
		FactorizeBasis	= True;
	int CheckResidMode;

	switch( Problem )
	{
	case SE_OK:					return True;

	case SE_DUAL_RESID_ALARM:	ComputeDuals = True;
	case SE_SINGULAR_BASIS:
	case SE_PRIM_RESID_ALARM:
	case SE_INF_RESID_ALARM:	goto SevereProblems;

	case SE_DUAL_RESID_POOR:	ComputeDuals = True;
	case SE_PRIM_RESID_POOR:
	case SE_INF_RESID_POOR:		goto MinorProblems;
	}

SevereProblems:

	if( Verbosity >= V_LOW )
		Print( "\n\nSevere numerical difficulties.\n" );

	do
	{
		if( BackTrack( BACK_STEPS ) <= 0 ) return False;
		HardLPTolerances( Verbosity );
	} while( UpdateBasis() != SLV_BASIS_REFACTORIZED );

	FactorizeBasis = False;

MinorProblems:

	if( Verbosity >= V_LOW )
		Print( "\n\nMinor numerical difficulties.\n" );

	if( FactorizeBasis && UpdateBasis() != SLV_BASIS_REFACTORIZED )
		goto SevereProblems;

	if( ComputeDuals )
	{
		ComputeDualVariables();
		CheckResidMode = CHK_DUAL;
		if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
			goto SevereProblems;
	}

	ComputePrimalVariables();
	CheckResidMode = CHK_PRIM | CHK_INF;
 	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
 		goto SevereProblems;

	ComputeResult();

	return True;
}
