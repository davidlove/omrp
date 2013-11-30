/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	init_bas.cpp
CREATED:            1994.04.16
LAST MODIFIED:		1995.08.12

DEPENDENCIES:       stdtype.h, std_tmpl.h, error.h, solver.h, solvcode.h,
					std_math.h, smartptr.h, print.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

STATIC FUNCTIONS:
	None.

STATIC DATA:
	None.

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __SOLVER_H__
#	include "solver.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T Solver::FindInitialBasis_Bixby( const Real_T PIV_TOL,
		VerbLevel Verbosity )
	Bool_T Solver::FindInitialBasis_External( const Real_T LU_PIVOT_TOL,
		VerbLevel Verbosity )

PURPOSE:
	Create an initial basis for the simplex method. The basis will be an upper
	triangular matrix created (in as much as possible) from the columns of the
	constraint matrix of the linear problem (LP). The algorithm used is a
	variant of Bob Bixby's basis creation method.

PARAMETERS:
	const Real_T PIV_TOL
		Relative tolerance for initial basis pivot size.

	VerbLevel Verbosity
		Verbosity level.

RETURN VALUE:
	Success / failure status.

SIDE EFFECTS:
	Fills A2B and B2A with new values. In case of failure, A2B and B2A are not
	guaranteed to be consistent.

------------------------------------------------------------------------------*/


Bool_T Solver::FindInitialBasis_Bixby( const Real_T PIV_TOL, // )
	VerbLevel Verbosity )
{
	//--------------------------------------------------------------------------
	//	Run crash procedure to obtain a basis proposal. Then reconstruct 'A2B'
	//	and 'B2A' hash tables. WARNING: Constraint matrix dimension (number of
	//	columns) may change.
	//
	A2B.Fill( A2B_UNDEF, LP.GetN() );
	B2A.Fill( B2A_UNDEF, LP.GetM() );

	LP.InitialBasis( PIV_TOL, A2B, Verbosity );
	N = LP.GetN();

	Int_T i, j;
	for( j = i = 0; i < N; i++ )
		if( A2B[i] >= 0 )
		{
			A2B[i] = j;
			assert( B2A[j] == B2A_UNDEF );
			B2A[ j++ ] = i;
		}
		else
			assert( A2B[i] == A2B_UNDEF );

	assert( j == M );

	//--------------------------------------------------------------------------
	//	Factorize basis. Compute dual variables. Return failure status if
	//	residuals are too big.
	//
	if( !FactorizeBasis() )
	{
		if( Verbosity >= V_LOW )
			Print( "\tCannot factorize the initial basis.\n"
				"\tRetrying with tighter pivot tolerance.\n" );
		return False;
	}

	ComputeDualVariables();

	int CheckResidMode = CHK_DUAL;

	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
	{
		if( Verbosity >= V_LOW )
			Print( "\tInitial dual residuals too large.\n"
				"\tRetrying with tighter pivot tolerances\n" );
		return False;
	}

	return True;
}

/*
Bool_T Solver::FindInitialBasis_External( const Real_T , // )
	VerbLevel Verbosity )
{
	assert( False );

	//--------------------------------------------------------------------------
	//	Copy 'A2B' from the 'ExternA2B' table. Then remake the 'B2A' hash table.
	//
	B2A.Fill( B2A_UNDEF, LP.GetM() );

	N = LP.GetN();

	Int_T i, j;
	for( j = i = 0; j < N; j++ )
	{
		assert( ExternalA2B[j] == A2B_LO || ExternalA2B[j] == A2B_UP ||
			ExternalA2B[j] == A2B_IN || ExternalA2B[j] >= 0 );
	
		A2B[j] = ExternalA2B[j];

		if( A2B[j] >= 0 )
		{
			A2B[j] = i;
			assert( B2A[i] == B2A_UNDEF );
			B2A[ i++ ] = j;
		}
		else
			assert( A2B[j] == A2B_UNDEF );
	}

	assert( i == M );

	//--------------------------------------------------------------------------
	//	Factorize basis. Compute dual variables. Return failure status if
	//	residuals are too big.
	//
	if( !FactorizeBasis() )
	{
		if( Verbosity >= V_LOW )
			Print( "\tCannot factorize the initial basis.\n" );
		return False;
	}

	ComputeDualVariables();

	int CheckResidMode = CHK_DUAL;

	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
	{
		if( Verbosity >= V_LOW )
			Print( "\tInitial dual residuals too large.\n" );
		return False;
	}

	return True;
}
*/
