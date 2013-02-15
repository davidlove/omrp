/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	start_pt.cpp
CREATED:            1994.04.16
LAST MODIFIED:		1995.09.12

DEPENDENCIES:       stdtype.h, std_tmpl.h, error.h, solver.h, solvcode.h,
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

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif



/*------------------------------------------------------------------------------

	Bool_T Solver::FindStartingPoint_Basic( VerbLevel Verbosity )
	Bool_T Solver::FindStartingPoint_AllOnes( VerbLevel Verbosity )
	Bool_T Solver::FindStartingPoint_AllNonBasicOnes( VerbLevel Verbosity )
	Bool_T Solver::FindStartingPoint_External( VerbLevel Verbosity )

PURPOSE:
	These functions are meant to provide the initial starting point for the
simplex algorithm. They all share the same functionality:
-	they set the values of all structural variables (according to some rule),
-	they initialize the 'A2B' hash table.

	The selection of the starting point for the simplex method has been made
independent from the selection of the initial basis. However some of those
functions may depend on previous construction (and successful factorization) of
the basis. If that should be the case, the corresponding basis building function
must be called first.

	FindStartingPoint_Basic
		assumes, that the non-basic variables are at their bounds, the basic
		ones are computed accordingly (implicitly - this will be done by the
		calling function); requires that the basis is constructed first.

	FindStartingPoint_AllOnes
		all variables are initially set to 1; independent of basis construction.

	FindStartingPoint_AllNonBasicOnes
		all non-basic variables are set to one, the basic ones are computed
		accordingly; requires that the basis is constructed first.

	FindStartingPoint_External
		we accept a solution that has been passed by "SetExternalSolution"
		member function call.

PARAMETERS:
	VerbLevel Verbosity
		Controls report output verbosity level.

RETURN VALUE:
	Boolean success status ('True' on success, 'False' on failure).

SIDE EFFECTS:
	Some of the functions check the primal residuals, and thus overwrite some
work vectors.

------------------------------------------------------------------------------*/

Bool_T Solver::FindStartingPoint_Basic( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH )
		Print( "\tGenerating a standard simplex starting point.\n" );

	//--------------------------------------------------------------------------
	//	Construct a standard simplex basic solution - put all the non-basic
	//	variables on their finite bounds (if they have them).
	//
	for( Int_T j = 0; j < N; j++ )
		if( A2B[j] == A2B_UNDEF )
		{
			if( VarType[j] & ( VT_FX | VT_HAS_LO_BND | VT_ARTIF ) )
			{
				A2B[j]	= A2B_LO;
				x[j]	= 0.0;
			}
			else if( VarType[j] & VT_HAS_UP_BND )
			{
				A2B[j]	= A2B_UP;
				x[j]	= LP.GetU( j );
			}
			else
			{
				A2B[j]	= A2B_IN;
				x[j]	= 0.0;
			}
		}
		else
			x[j] = 0.0;

	//--------------------------------------------------------------------------
	//	Function "FindInitialBasis()" has already factorized the basis and
	//	computed the dual variables. Thus we may now compute the primal
	//	variables.
	//
	ComputePrimalVariables();

	int CheckResidMode = CHK_PRIM;

	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
	{
		if( Verbosity >= V_LOW )
			Print( "\tInitial primal residuals too large.\n"
				"\tRetrying with tighter pivot tolerance.\n" );
		return False;
	}

	return True;
}


/*----------------------------------------------------------------------------*/

Bool_T Solver::FindStartingPoint_AllOnes( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH )
		Print( "\tGenerating an interior starting point:\n"
			"\tall variables set to one.\n" );

	//--------------------------------------------------------------------------
	//	Construct an interior starting solution. All variables with a finite
	//	lower bound (which is guaranteed to have been shifted to zero) will be
	//	placed at one. Wariables without lower bound and with a finite upper
	//	bound 'u[j]' will be placed at 'u[j] - 1'. Free and fixed variables will
	//	be placed at zero.
	//
	for( Int_T j = 0; j < N; j++ )
	{
		if( ( VarType[j] & VT_HAS_LO_BND ) && !( VarType[j] & VT_FX ) )
			x[j] = 1.0;
		else if( VarType[j] & VT_HAS_UP_BND )
			x[j] = LP.GetU( j ) - 1.0;
		else
			x[j] = 0.0;
		if( A2B[j] == A2B_UNDEF ) A2B[j] = A2B_IN;
	}

	return True;
}


/*----------------------------------------------------------------------------*/

Bool_T Solver::FindStartingPoint_AllNonBasicOnes( VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH )
		Print( "\tGenerating an interior starting point:\n"
			"\tAll non-basic variables set to one.\n" );

	if( !FindStartingPoint_AllOnes( V_NONE ) ) return False;

	//--------------------------------------------------------------------------
	//	Function "FindInitialBasis()" has already factorized the basis and
	//	computed the dual variables. Thus we may now compute the primal
	//	variables.
	//
	ComputePrimalVariables();

	int CheckResidMode = CHK_PRIM;

	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
	{
		if( Verbosity >= V_LOW )
			Print( "\tInitial primal residuals too large.\n"
				"\tRetrying with tighter pivot tolerance.\n" );
		return False;
	}

	return True;
}


/*----------------------------------------------------------------------------*/

Bool_T Solver::FindStartingPoint_External( VerbLevel Verbosity )
{
	if( !UseExternalSolution )
		return False;

	if( Verbosity >= V_HIGH )
		Print( "\tUsing user-selected starting point.\n" );

	//--------------------------------------------------------------------------
	//	Copy the user-defined starting point from the array "ExternalX".
	//
	Int_T StructN = LP.GetStructN();

	Int_T j;
	for( j = 0; j < StructN; j++ )
	{
		x[j]	= ExternalX[j];
		if( A2B[j] < 0 )
			A2B[j]	= A2B_IN;

		if( VarType[j] & VT_HAS_LO_BND && x[j] <= 0.0 )
		{
			if( A2B[j] < 0 )
				A2B[j]	= A2B_LO;
			x[j]	= 0.0;
		}

		if( VarType[j] & VT_HAS_UP_BND && x[j] >= u[j] )
		{
			if( A2B[j] < 0 )
				A2B[j]	= A2B_UP;
			x[j]	= u[j];
		}
	}

	for( ; j < N; j++ )
	{
		x[j] = 0.0;
		if( A2B[j] < 0 )
		{
			if( VarType[j] & VT_HAS_LO_BND )
				A2B[j]	= A2B_LO;
			else if( VarType[j] & VT_HAS_UP_BND && u[j] < 0.0 )
				A2B[j]	= A2B_UP;
			else
				A2B[j]	= A2B_IN;
		}
	}

	//--------------------------------------------------------------------------
	//	Function "FindInitialBasis()" has already factorized the basis and
	//	computed the dual variables. Thus we may now compute the primal
	//	variables.
	//
	ComputePrimalVariables();

	int CheckResidMode = CHK_PRIM;

	if( CheckResiduals( CheckResidMode ) > SATISF_RESID )
	{
		if( Verbosity >= V_LOW )
			Print( "\tInitial primal residuals too large.\n"
				"\tRetrying with tighter pivot tolerance.\n" );
		return False;
	}

	//--------------------------------------------------------------------------
	//	Unmark the external solution flag, so that it is not used more than
	//	once.
	//
	UseExternalSolution = False;

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Solver::SetExternalSolution( const Array<Real_T> &ExternX,
		Int_T Len )

PURPOSE:

PARAMETERS:
	const Array<Real_T> &ExternX
		External solution vector.

	Int_T Len
		Vector length. Should be equal to LP.GetStructN().

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	The next time solver is initialized the external solution will become the
	default initial solution.

------------------------------------------------------------------------------*/

Bool_T Solver::SetExternalSolution( const Array<Real_T> &ExternX, Int_T Len )
{
	if( Len != LP.GetStructN() )
		return False;

	ExternalX.Resize( Len );
	ExternalX.Copy( ExternX, Len, Len );

	UseExternalSolution = True;

	return True;
}
