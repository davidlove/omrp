/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	pp_simpl.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1995.10.24

DEPENDENCIES:		std_tmpl.h, presolve.h, simplex.h, error.h, smartptr.h,
					memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h, stdtype.h,
					lp_codes.h, solv_lp.h, compile.h, mps_lp.h, sort_lab.h,
					myalloc.h, std_math.h, cl_list.h
					<math.h>, assert.h

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


#include <math.h>
#include <assert.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif


/*------------------------------------------------------------------------------

	Int_T Presolver::EliminateEmptyRows( void )

PURPOSE:
	Finds empty rows in the problem and removes them. Uses linked lists of rows
ordered by row length.

PARAMETERS:
	None.

RETURN VALUE:
	Number of removed rows.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::EliminateEmptyRows( void )
{
	Int_T Elim = 0, i;

	for( i = RowList.GetFirst( 0 ); i != -1; i = RowList.GetNext(), Elim++ )
	{
		assert( i >= 0 && i < m && !ExcludeRows[i] );

		//----------------------------------------------------------------------
		//	See if the problem is feasible (i.e. if 0 \in < bl[i], bu[i] >.
		//
		if( bl[i] > 0.0 || bu[i] < 0.0 )
		{
			PrimalInf = Max( Max( -bl[i], bu[i] ), PrimalInf );
			if( bl[i] > ftol || bu[i] < -ftol )
			{
				Status = LPS_INFEASIBLE;
				break;
			}
		}


		RemoveRow( i );
	}

	return Elim;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::EliminateEmptyColumns( void )

PURPOSE:
	Finds empty columns in the problem and removes them. Uses linked lists of
columns ordered by column length.

PARAMETERS:
	None.

RETURN VALUE:
	Number of removed columns.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::EliminateEmptyColumns( void )
{
	Int_T Elim = 0, j;

	for( j = ColList.GetFirst( 0 ); j != -1; j = ColList.GetNext(), Elim++ )
	{
		assert( !ExcludeCols[j] && j >= 0 && j < n );

		//----------------------------------------------------------------------
		//	Attempt to fix a variable with an empty column (we may fix it on any
		//	bound, or between the bounds - it all depends on the value of its
		//	objective function coefficient).
		//
		if( LP.GetC( j ) < 0.0 )
		{
			if( !( LP.GetVarType( j ) & VT_UP ) )
			{
				DualInf = Max( DualInf, -LP.GetC( j ) );
				if( LP.GetC( j ) < -ftol )
				{
					Status = LPS_UNBOUNDED;
					break;
				}
				FixVariable( j, 0.0 );
			}
			else
				FixVariable( j, VT_UP );
		}
		else if( LP.GetC( j ) > 0.0 )
		{
			if( !LP.GetVarType( j ) & VT_LO )
			{
				DualInf = Max( DualInf, LP.GetC( j ) );
				if( LP.GetC( j ) > ftol )
				{
					Status = LPS_UNBOUNDED;
					break;
				}
				else
					FixVariable( j, 0.0 );
			}
			else
				FixVariable( j, VT_LO );
		}
		else
		{
			if( LP.GetVarType( j ) & VT_LO )
				FixVariable( j, VT_LO );
			else if( LP.GetVarType( j ) & VT_UP )
				FixVariable( j, VT_UP );
			else
				FixVariable( j, 0.0 );
		}
	}

	return Elim;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::EliminateFixedVariables( void )

PURPOSE:
	Eliminates the fixed variables from the LP (only those that were found to be
fixed in the original problem formulation). Is called only once at the beginning
of presolve analysis.

PARAMETERS:
	None.

RETURN VALUE:
	Number of removed columns.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::EliminateFixedVariables( void )
{
	Int_T Elim = 0;
	for( Int_T j = 0; j < n; j++ )
		if( !ExcludeCols[j] && ( LP.GetVarType( j ) & VT_FX ) )
		{
			FixVariable( j, VT_FX );
			Elim++;
		}

	return Elim;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::EliminateSingletonRows( void )

PURPOSE:
	Finds singleton rows in the problem and converts them into simple bounds on
variables. Uses linked lists of rows ordered by row length.

PARAMETERS:
	None.

RETURN VALUE:
	Number of removed (converted) rows.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::EliminateSingletonRows( void )
{
	//--------------------------------------------------------------------------
	//	Loop on singleton rows.
	//
	Int_T Elim = 0;
	for( Int_T i = RowList.GetFirst( 1 ); i != -1;
		i = RowList.GetNext(), Elim++ )
	{
		//----------------------------------------------------------------------
		//	See if this really is a singleton and if it hasn't been removed
		//	before. See if 'i' is in range 0,...,M-1.
		//
		assert( RowLen[i] == 1 && !ExcludeRows[i] && i >= 0 && i < m );

		//----------------------------------------------------------------------
		//	Find the non-zero column and value.
		//
		Ptr<Real_T> a;
		Ptr<Int_T> col;
		Int_T len;

		for( LP.GetRow( i, a, col, len ); len; ++col, ++a, --len )
			if( !ExcludeCols[ *col ] ) break;

		//----------------------------------------------------------------------
		//	Remove row as redundant.
		//
		RemoveRow( i );

		//----------------------------------------------------------------------
		//	Change singleton row data to bounds on the variable. Remember to
		//	change direction of inequalities if a_ij is negative.
		//
		if( rt[i] & VT_FX )
			FixVariable( *col, bl[i] / *a );
		else
		{
			if( rt[i] & VT_LO )
				if( IsPositive( *a ) )
					SetL( *col, bl[i] / *a );
				else
					SetU( *col, bl[i] / *a );

			//
			//	The previous call may have fixed the variable.
			//
			if( !ExcludeCols[ *col ] && rt[i] & VT_UP )
				if( IsPositive( *a ) )
					SetU( *col, bu[i] / *a );
				else
					SetL( *col, bu[i] / *a );
		}

		if( Status != LPS_UNKNOWN ) break;
	}
	//
	//	End of loop on singleton rows.
	//--------------------------------------------------------------------------

	return Elim;
}
