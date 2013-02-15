/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	pp_primi.cpp
CREATED:			1993.10.07
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		std_tmpl.h, std_math.h, presolve.h, simplex.h, error.h,
					smartptr.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h,
					stdtype.h, lp_codes.h, solv_lp.h, compile.h, mps_lp.h,
					myalloc.h, cl_list.h, postsolv.h, my_defs.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Presolver primitives - to be used as the low level functions for removing
columns and rows of the constraint matrix.

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
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif

#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif
#ifndef __POSTSOLV_H__
#	include "postsolv.h"
#endif


/*------------------------------------------------------------------------------

	void Presolver::RemoveColumn( Int_T i )

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

void Presolver::RemoveColumn( Int_T j )
{
	assert( j >= 0 && j < n && !ExcludeCols[j] );

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	for( LP.GetColumn( j, a, row, len ); len; ++a, ++row, --len )
		if( !ExcludeRows[ *row ]  )
			RowList.TruncateLength( *row, 1 );
	
	ColList.RmFromList( j );
	ExcludeCols[j] = True;
	EliminatedCols++;
}


/*------------------------------------------------------------------------------

	void Presolver::FixVariable( Int_T j, int vt )

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

void Presolver::FixVariable( Int_T j, int vt )
{
	switch( vt )
	{
	//--------------------------------------------------------------------------
	//	Attempt to fix (and remove) a variable on its lower bound.
	//
	case VT_LO:
		if( LP.GetVarType( j ) & VT_LO )
			FixVariable( j, LP.GetL( j ) );
		else
			Status = LPS_UNBOUNDED;
		break;

	//--------------------------------------------------------------------------
	//	Attempt to fix (and remove) a variable on its upper bound.
	//
	case VT_UP:
		if( LP.GetVarType( j ) & VT_UP )
			FixVariable( j, LP.GetU( j ) );
		else
			Status = LPS_UNBOUNDED;
		break;

	//--------------------------------------------------------------------------
	//	Removal of a fixed variable.
	//
	case VT_FX:
		assert( LP.GetVarType( j ) & VT_FX );

		FixVariable( j, LP.GetL( j ) );
		break;

	//--------------------------------------------------------------------------
	//	Obviously an error.
	//
#ifndef NDEBUG
	default:
		abort();
#endif
	}
}


/*------------------------------------------------------------------------------

	void Presolver::FixVariable( Int_T j, Real_T val )

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


void Presolver::FixVariable( Int_T j, Real_T val )
{
	assert( j >= 0 && j < n && !ExcludeCols[j] );

	LP.SetL( j, val );
	LP.SetU( j, val );

	assert( LP.GetVarType( j ) & VT_FX );

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	if( IsNonZero( val ) )
	{
		for( LP.GetColumn( j, a, row, len ); len; ++a, ++row, --len )
			if( !ExcludeRows[ *row ]  )
			{
				Real_T shift = *a * val;

				if( rt[ *row ] & VT_LO && IsZero( bl[ *row ] -= shift ) )
					bl[ *row ] = 0.0;

				if( rt[ *row ] & VT_UP && IsZero( bu[ *row ] -= shift ) )
					bu[ *row ] = 0.0;
			}

		f+= val * LP.GetC( j );
	}

	RemoveColumn( j );

	if( PostSolve )
		PostSolve->VariableFixing( j, val );
}


/*------------------------------------------------------------------------------

	void Presolver::RemoveRow( Int_T i )

PURPOSE:
	Remove a row from lists of rows and truncate all columns, in which the
	row in question has non-zeros.

PARAMETERS:
	Int_T i
		Row number.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/


void Presolver::RemoveRow( Int_T i )
{
	assert( i >= 0 && i < m && !ExcludeRows[i] );

	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T len;

	for( LP.GetRow( i, a, col, len ); len; ++col, --len )
		if( !ExcludeCols[ *col ] )
			ColList.TruncateLength( *col, 1 );

	RowList.RmFromList( i );
	ExcludeRows[i] = True;
	EliminatedRows++;
}


/*------------------------------------------------------------------------------

	void Presolver::EliminateFreeSingletonColumn( Int_T i, Int_T j,
		Real_T a_ij )

PURPOSE:
	Remove the free singleton variable and the row in which it has its only non-
zero. Adjust cost coefficients for all other variables that have non-zeros in
this row.

PARAMETERS:
	Int_T i, Int_T j, Real_T a_ij
		Row, column and non-zero value of a column singleton that was found to
		be a free or implied free variable, and may thus be eliminated.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/


void Presolver::EliminateFreeSingletonColumn( Int_T i, Int_T j, Real_T a_ij )
{
	assert( !ExcludeCols[j] && j >= 0 && j < n );

	Real_T cf_af	= LP.GetC( j ) / a_ij;
	const Short_T RT	= rt[i];
	Real_T SlOpt		= 0.0,
		rhs				= bu[i];

	//--------------------------------------------------------------------------
	//	Update the cost vector and the fixed adjustment to the objective
	//	function.
	//
	if( IsNonZero( cf_af ) )
	{
		if( RT == VT_FIXED )
			f += cf_af * bl[i];
		else if( RT == ( VT_LO | VT_UP ) )
		{
			if( cf_af > 0 )
			{
				f += cf_af * bl[i];
				SlOpt = bu[i] - bl[i];
			}
			else
				f += cf_af * bu[i];
		}
		else if( RT == VT_UP )
		{
			if( cf_af > 0 )
			{
				DualInf = Max( DualInf, cf_af );
				if( cf_af > ftol )
				{
					Status = LPS_UNBOUNDED;
					return;
				}
				else
					cf_af = 0.0;
			}
			else
				f += cf_af * bu[i];
		}
		else if( RT == VT_LO )
		{
			rhs = bl[i];

			if( cf_af > 0 )
				f += cf_af * bl[i];
			else
			{
				DualInf = Max( DualInf, -cf_af );
				if( cf_af < -ftol )
				{
					Status = LPS_UNBOUNDED;
					return;
				}
				else
					cf_af = 0.0;
			}
		}
#ifndef NDEBUG
		else
			abort();
#endif
	}

	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T len;

	RemoveColumn( j );

	if( IsNonZero( cf_af ) )
		for( LP.GetRow( i, a, col, len ); len; ++a, ++col, --len )
			if( !ExcludeCols[ *col ] )
				LP.SetC( *col, LP.GetC( *col ) - cf_af * *a );

	LP.GetRow( i, a, col, len );
	if( PostSolve )
		PostSolve->FreeSingletonColumnRemoval( j, a, col, len, a_ij, rhs,
			SlOpt );

	//--------------------------------------------------------------------------
	//	Finally remove the row, which becomes free after the free variable is
	//	removed.
	//
	RemoveRow( i );
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::SetU( Int_T j, Real_T val )

PURPOSE:
	Set the upper bound on variable 'j' to new value 'val' if and only if the
new bound is tighter than the one that was set before.

PARAMETERS:
	Int_T j
		Number of variable.
	
	Real_T val
		The new (proposed) upper bound.

RETURN VALUE:
	'True' if the upper bound was tightened. 'False' otherwise.

SIDE EFFECTS:
	Possible detection of infeasibility.

------------------------------------------------------------------------------*/

Bool_T Presolver::SetU( Int_T j, Real_T val )
{
	assert( j >= 0 && j < n && !ExcludeCols[j] );

	if( val > LP.GetU( j ) ) return False;
	if( val < LP.GetL( j ) )
	{
		PrimalInf = Max( LP.GetL( j ) - val, PrimalInf );
		if( val < LP.GetL( j ) - ftol )
		{
			Status = LPS_INFEASIBLE;
			return False;
		}
		else
			val = LP.GetL( j );
	}

	LP.SetU( j, val );

	if( LP.GetVarType( j ) & VT_FX )
		FixVariable( j, VT_FX );

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::SetL( Int_T j, Real_T val )

PURPOSE:
	Set the lower bound on variable 'j' to new value 'val' if and only if the
new bound is tighter than the one that was set before.

PARAMETERS:
	Int_T j
		Number of variable.
	
	Real_T val
		The new (proposed) lower bound.

RETURN VALUE:
	'True' if the lower bound was tightened. 'False' otherwise.

SIDE EFFECTS:
	Possible detection of infeasibility.

------------------------------------------------------------------------------*/

Bool_T Presolver::SetL( Int_T j, Real_T val )
{
	assert( j >= 0 && j < n && !ExcludeCols[j] );

	if( val < LP.GetL( j ) ) return False;
	if( val > LP.GetU( j ) )
	{
		PrimalInf = Max( val - LP.GetU( j ), PrimalInf );
		if( val > LP.GetU( j ) + ftol )
		{
			Status = LPS_INFEASIBLE;
			return False;
		}
		else
			val = LP.GetU( j );
	}

	LP.SetL( j, val );

	if( LP.GetVarType( j ) & VT_FX )
		FixVariable( j, VT_FX );

	return True;
}


/*------------------------------------------------------------------------------

	void Presolver::RHS2BLU( void )

PURPOSE:
	Converts the right hand side and range vectors (with appropriate row type
	information) into the format more suitable for presolving: lower and upper
	bounds on row activity plus a flag denoting the row type.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::RHS2BLU( void )
{
	for( Int_T i = 0; i < m; i++ )
	{
		LP.GetConstraintRange( i, bl[i], bu[i] );

		Short_T RowType = LP.GetRowType( i );

		if( RowType & RT_RNG ) 
			rt[i] = VT_BOUNDED;
		else
			switch( RowType & RT_TYPE )
			{
			case RT_LE:	rt[i] = VT_MI;		break;
			case RT_GE:	rt[i] = VT_NORM;	break;
			case RT_EQ:	rt[i] = VT_FIXED;	break;
			}
	}
}


/*------------------------------------------------------------------------------

	void Presolver::BLU2RHS( void )

PURPOSE:
	Converts the lower/upper bounds on row activity back to the original form
	of right hand side and range vectors.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::BLU2RHS( void )
{
	for( Int_T i = 0; i < m; i++ )
		LP.SetConstraintRange( i, bl[i], (rt[i] & VT_LO) ? True : False, bu[i],
			(rt[i] & VT_UP) ? True : False );
}
