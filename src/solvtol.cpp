/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solvtol.cpp
CREATED:			1993.10.24
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		solver.h, stdtype.h, solv_lp.h, std_tmpl.h, error.h,
					lp_codes.h, mps_lp.h, inverse.h, solvcode.h, simplex.h,
					solvtol.h, std_math.h, print.h

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

------------------------------------------------------------------------------*/

#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __SOLVTOL_H__
#	include "solvtol.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif


/*------------------------------------------------------------------------------

	void Solver::ResetTolerances( VerbLevel Verbosity )


PURPOSE:
	x

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ResetTolerances( VerbLevel Verbosity )
{
	if( Initialized )
	{
		if( Verbosity >= V_HIGH  )
			Print( "\nReseting the tolerances.\n\n" );
		else if( Verbosity >= V_LOW )
			Print( "Reseting the tolerances.\n" );
	}

	SetZeroTolerance( SMALL_ELEM_DEF );
	FEASIBILITY_TOL = FEASIBILITY_TOL_DEF;
	if( Pricing >= PRS_SE )
		OPTIMALITY_TOL = OPTIMALITY_TOL_DEF * OPTIMALITY_TOL_DEF;
	else
		OPTIMALITY_TOL = OPTIMALITY_TOL_DEF;
	LU_PIVOT_TOL		= S.LU_PIVOT_TOL;
	PIVOT_TOL			= S.PIVOT_TOL;
	REFACT_FREQ			= S.REFACT_FREQ;
	RESID_CHECK_FREQ	= S.RESID_CHECK_FREQ;
	B->SetLU_Tol( LU_PIVOT_TOL );
}


/*------------------------------------------------------------------------------

	void Solver::HardLPTolerances( VerbLevel Verbosity )

PURPOSE:
	When the problem seems to be too difficult, we adjust the tolerances. We
expect the numerical errors to be larger (and therefore increase zero tolerance
and 'PIVOT_TOL'). To improve factorization stability we also increase
'LU_UPDATE_TOL'. In anticipation of more numerical difficulties we reduce the
maximum number of iterations between refactorizations and computation of
residuals.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::HardLPTolerances( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH  )
		Print( "\nTolerances for harder LP.\n\n" );
	else if( Verbosity >= V_LOW )
		Print( "Iter: %8d\tResult: %10.2E\t"
			"Tolerances for harder LP.\n", (int) IterCnt, Result );

	SetZeroTolerance( Min( GetZeroTolerance() * 5.0e0, SMALL_ELEM_HI ) );
	LU_PIVOT_TOL = Min( LU_PIVOT_TOL * 2.0e0, LU_PIVOT_TOL_HI );
	PIVOT_TOL = Min( PIVOT_TOL * 2.0e0, PIVOT_TOL_HI );
	REFACT_FREQ = Max( Int_T( 2 * REFACT_FREQ / 3 ), (Int_T) REFACT_FREQ_LO );
	RESID_CHECK_FREQ = Max( Int_T( 2 * RESID_CHECK_FREQ / 3 ),
		(Int_T) RESID_CHECK_FREQ_LO );

	B->SetLU_Tol( LU_PIVOT_TOL );
}


/*------------------------------------------------------------------------------

	void Solver::EasyLPTolerances( VerbLevel Verbosity )

PURPOSE:
	When the problem seems to be easy, we adjust the tolerances. Since the
problem seems to be easy we may allow smaller pivots (reduced 'PIVOT_TOL' and
'LU_PIVOT_TOL'). Numerical errors are not expected to grow very large, therefore
zero tolerance is reduced in size in order to treat more small values as
non-zeros and not just round-off errors.  We may also allow the solver to
refactorize and check residuals less frequently (increased values of
'REFACT_FREQ' and 'RESID_CHECK_FREQ').

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::EasyLPTolerances( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH  )
		Print( "\nTolerances for easier LP.\n\n" );
	else if( Verbosity >= V_LOW )
		Print( "Iter: %8d\tResult: %10.2E\t"
			"Tolerances for easier LP.\n", (int) IterCnt, Result );

	SetZeroTolerance( Max( GetZeroTolerance() / 5.0e0, SMALL_ELEM_LO ) );
	LU_PIVOT_TOL = Max( LU_PIVOT_TOL / 2.0e0, LU_PIVOT_TOL_LO );
	PIVOT_TOL = Max( PIVOT_TOL / 2.0e0, PIVOT_TOL_LO );
	REFACT_FREQ = Min( Int_T( 3 * REFACT_FREQ / 2 ), (Int_T) REFACT_FREQ_HI );
	RESID_CHECK_FREQ = Min( Int_T( 3 * RESID_CHECK_FREQ / 2 ),
		(Int_T) RESID_CHECK_FREQ_HI );

	B->SetLU_Tol( LU_PIVOT_TOL );
}


/*------------------------------------------------------------------------------

	void Solver::FinalTolerances( VerbLevel Verbosity )

PURPOSE:
	When we are getting close to the solution, we need to have an exact
solution. For this reason we enforce highest possible factorization accuracy
(through improved stability - increased 'LU_PIVOT_TOL').  We tend to treat
more numbers seriously (reduced zero tolerance and 'PIVOT_TOL'), beacause the
factorization (and solves) will be more accurate and the conditioning of the
optimal basis is usually good.  To be sure the solution is feasible we
significantly reduce 'FEASIBILITY_TOL'.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::FinalTolerances( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH  )
		Print( "\nTaking final tolerances.\n\n" );
	else if( Verbosity >= V_LOW )
		Print( "Iter: %8d\tResult: %10.2E\t"
			"Taking final tolerances.\n", (int) IterCnt, Result );

	SetZeroTolerance( SMALL_ELEM_DEF );
	FEASIBILITY_TOL = FEASIBILITY_TOL_LO;
	if( Pricing >= PRS_SE )
		OPTIMALITY_TOL = OPTIMALITY_TOL_LO * OPTIMALITY_TOL_LO;
	else
		OPTIMALITY_TOL = OPTIMALITY_TOL_LO;
	PIVOT_TOL = PIVOT_TOL_DEF;
	LU_PIVOT_TOL = LU_PIVOT_TOL_HI;

	B->SetLU_Tol( LU_PIVOT_TOL );
}
