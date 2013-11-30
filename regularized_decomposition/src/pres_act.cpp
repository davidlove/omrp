/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	pres_act.cpp
CREATED:			1994.05.08
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		std_math.h, std_tmpl.h, smartptr.h, error.h, memblock.h,
					smartdcl.h, sptr_deb.h, sptr_ndb.h, solution.h, compile.h,
					stdtype.h, lp_codes.h, postsolv.h, my_defs.h,
					simplex.h, myalloc.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif

#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __POSTSOLV_H__
#	include "postsolv.h"
#endif


//==============================================================================
//
//	Class "FixedAdjustment" member functions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	Bool_T FixedAdjustment::Write( FILE *fp )

PURPOSE:
	Write a line of text to a text file.

PARAMETERS:
	FILE *fp
		Output file.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T FixedAdjustment::Write( FILE *fp )
{
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	First line - marks the beginning of a new description.
	//
	if( fprintf( fp, "%-8s\n", "ADJUST" ) == EOF ||
		fprintf( fp, "  %-8s  %.12G\n", "VALUE", Val ) == EOF )
		return False;

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T FixedAdjustment::Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols )

PURPOSE:
	x

PARAMETERS:
	Solution &solution
		Solution structure.

	Real_T FTOL
		Ignored.

	Array<Bool_T> &ExcludeCols
		Ignored.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T FixedAdjustment::Undo( Solution &solution, Real_T, Array<Bool_T> & )
{
	solution.result += Val;
	return True;
}


//==============================================================================
//
//	Class "VariableFixing" member functions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	Bool_T VariableFixing::Write( FILE *fp )

PURPOSE:
	Write a line of text to a text file.

PARAMETERS:
	FILE *fp
		Output file.
	
RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T VariableFixing::Write( FILE *fp )
{
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	First line - marks the beginning of a new description. Second line gives
	//	the variable value.
	//
	if( fprintf( fp, "%-8s  %-8d  %-8s\n", "ACTION", (int)j, "FIX" ) == EOF ||
		fprintf( fp, "  %-8s  %.12G\n", "VALUE", Val ) == EOF )
		return False;

	return True;
}


//==============================================================================
//
//	Class "ExplicitSlackRemoval" member functions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	ExplicitSlackRemoval::ExplicitSlackRemoval( Int_T jj,
		const Array<Bool_T> *ExcludeCols,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
		Short_T vartype, Real_T ls, Real_T us,
		Short_T rowtype, Real_T bl, Real_T bu )

PURPOSE:
	Constructor.

PARAMETERS:
	Int_T jj
		Explicit slack variable's index.

	const Array<Bool_T> *ExcludeCols
		Used for writing down only meaningful part of the row.

	const Ptr<Real_T> &a, const Ptr<Label> &col, Int_T len
		Matrix row.

	Real_T as
		The coefficient in front of the explicit slack.

	Short_T vartype, Real_T ls, Real_T us
		Variable type, lower and upper bounds on variable.

	Short_T rowtype, Real_T bl, Real_T bu
		Row type, lower and upper bounds on row activity.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

ExplicitSlackRemoval::ExplicitSlackRemoval( Int_T jj, // )
	const Array<Bool_T> *ExcludeCols,
	const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
	Short_T vartype, Real_T ls, Real_T us,
	Short_T rowtype, Real_T bl, Real_T bu )
	: j( jj ), A(), Col(), Len( 0 ), As( as ),
	Ls( ( vartype & VT_LO ) ? ls : -INFINITY ),
	Us( ( vartype & VT_UP ) ? us : INFINITY ),
	bL( ( rowtype & VT_LO ) ? bl : -INFINITY ),
	bU( ( rowtype & VT_UP ) ? bu : INFINITY ),
	VarType( vartype ), RowType( rowtype )
{
	assert( jj >= 0 );
	assert( ( vartype & ( VT_LO | VT_UP | VT_FX ) ) == vartype );
	assert( ( rowtype & ( VT_LO | VT_UP | VT_FX ) ) == rowtype );

	Int_T i, RowLen;

	for( i = RowLen = 0; i < len; i++ )
		if( col[i] != j && ( !ExcludeCols || !(*ExcludeCols)[ col[i] ] ) )
			RowLen++;

	assert( RowLen > 0 );

	Len = RowLen;
	A.Resize( RowLen );
	Col.Resize( RowLen );
	for( i = 0; i < len; i++ )
		if( col[i] != j && ( !ExcludeCols || !(*ExcludeCols)[ col[i] ] ) )
		{
			RowLen--;
			A[RowLen]	= a[i];
			Col[RowLen]	= col[i];
		}
}


/*------------------------------------------------------------------------------

	Bool_T ExplicitSlackRemoval::Write( FILE *fp )

PURPOSE:
	Writes the information stored during the explicit slack removal to the
	output file.

PARAMETERS:
	FILE *fp
		Output file.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ExplicitSlackRemoval::Write( FILE *fp )
{
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	First line - marks the beginning of a new description.
	//
	if( fprintf( fp, "%-8s  %-8d  %12s\n", "ACTION", (int)j, "EXPL_SLACK" )
		== EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Section dependent on the row type: row type and row activity limits.
	//
	switch( RowType & VT_TYPE )
	{
	case VT_LO:
		if( fprintf( fp, "  %-8s  %.12G\n", "GE", bL ) == EOF )
			return False;
		break;

	case VT_UP:
		if( fprintf( fp, "  %-8s  %.12G\n", "LE", bU ) == EOF )
			return False;
		break;

	case VT_LO | VT_UP:
		if( fprintf( fp, "  %-8s  %.12G  %.12G\n", "RG", bL, bU ) == EOF )
			return False;
		break;

	case VT_LO | VT_UP | VT_FX:
		if( fprintf( fp, "  %-8s  %.12G\n", "EQ", bL ) == EOF )
			return False;
		break;

#ifndef NDEBUG
	default:
		abort();
#endif
	}

	//--------------------------------------------------------------------------
	//	Section dependent on the variable type: variable simple bounds.
	//
	if( VarType & VT_FX )
	{
		if( fprintf( fp, "  %-8s  %.12G\n", "FX", Ls ) == EOF )
			return False;
	}
	else if( ( VarType & VT_LO ) && ( VarType & VT_UP ) )
	{
		if( fprintf( fp, "  %-8s  %.12G  %.12G\n", "UP", Ls, Us ) == EOF )
			return False;
	}
	else if( VarType & VT_LO )
	{
		if( fprintf( fp, "  %-8s  %.12G\n", "PL", Ls ) == EOF )
			return False;
	}
	else if( VarType & VT_UP )
	{
		if( fprintf( fp, "  %-8s  %.12G\n", "MI", Us ) == EOF )
			return False;
	}

	//--------------------------------------------------------------------------
	//	Explicit slack's non-zero value.
	//
	if( fprintf( fp, "  %-8s  %.12G\n", "VALUE", As ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	The row (all non-zeros coefficients).
	//
	for( Int_T i = 0; i < Len; i++ )
		if( fprintf( fp, "  %-8s  %-8d  %.12G\n", "COEFF", Col[i], A[i] )
			== EOF )
			return False;

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T ExplicitSlackRemoval::Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols )

PURPOSE:
	x

PARAMETERS:
	Solution &solution
		Solution structure.

	Real_T FTOL
		Feasibility tolerance.

	Array<Bool_T> &ExcludeCols
		"ExcludeCols" is used to check if we're not trying to write the same
		value twice. it is updated on exit to indicate that the variable that
		was just restored is no longer "excluded".


RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ExplicitSlackRemoval::Undo( Solution &solution, Real_T FTOL, // )
	Array<Bool_T> &ExcludeCols )
{
	assert( ExcludeCols[j] );
	ExcludeCols[j] = False;

	Real_T xl	= -INFINITY,		// Lower bound on 'x'.
		xu		= INFINITY,			// Upper bound on 'x'.
		ax		= 0.0,				// Dot product of row and vector 'x'.
		&x		= solution.x[j];	// Reference to the expl. slack's value.
	x = 0.0;

	//--------------------------------------------------------------------------
	//	Compute ax = a^T * x
	//
	for( Int_T i = 0; i < Len; i++ )
	{
		assert( j != Col[i] );
		ax += A[i] * solution.x[Col[i]];
	}

	assert( IsNonZero( As ) );

	//--------------------------------------------------------------------------
	//	Compute lower and upper bounds on 'x' (explicit slack) implied by
	//	row activity 'ax' and bounds on row constraints.
	//
	if( IsPositive( As ) )
	{
		if( RowType & VT_LO )	xl = ( bL - ax ) / As;
		if( RowType & VT_UP )	xu = ( bU - ax ) / As;

		assert( xu >= xl );
	}
	else
	{
		if( RowType & VT_LO )	xu = ( bL - ax ) / As;
		if( RowType & VT_UP )	xl = ( bU - ax ) / As;

		assert( xu >= xl );
	}

	//--------------------------------------------------------------------------
	//	Variable simple bounds (Ls, Us) and bounds implied by row activity
	//	(xl, xu) must not contradict each other. This is checked. Afterwards
	//	value of 'x' (expl. slack) is chosen so as to be on one of the simple
	//	bounds.
	//
	if( ( VarType & VT_LO ) &&		// Finite lower simple bound.
		!IsNegative( Ls - xl ) &&	// Ls >= xl
		!IsPositive( Ls - xu ) )	// Ls <= xu
		x = Ls;
	else if( ( VarType & VT_UP ) &&	// Finite upper sime bound.
		!IsNegative( Us - xl ) &&	// Us >= xl
		!IsPositive( Us - xu ) )	// Us <= xu
		x = Us;
	else if( ( RowType & VT_LO ) &&	// Finite lower bound on row activity.
		!IsNegative( xl - Ls ) )	// xl >= Ls
		x = xl;
	else if( ( RowType & VT_UP ) &&	// Finite upper bound on row activity.
		!IsPositive( xu - Us ) )	// xu <= Us
		x = xu;
	else
	{
		if( ( VarType & VT_LO ) &&		// Finite lower simple bound.
			xl > Us )
			x = Us;
		else if( ( VarType & VT_UP ) &&	// Finite upper simple bound.
			xu < Ls )
			x = Ls;

		if( x > Us + FTOL )
			Warning( "Explicit slack x[%d] infeasibile: u=%g, x-u=%g.",
				j, Us, x - Us );
		else if( x < Ls - FTOL )
			Warning( "Explicit slack x[%d] infeasibile: l=%g, l-x=%g.",
				j, Ls, Ls - x );
	}

	return True;
}



/*------------------------------------------------------------------------------

	FreeSingletonColumnRemoval::FreeSingletonColumnRemoval( Int_T jj,
		const Array<Bool_T> *ExcludeCols,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
		Real_T af, Real_T b, Real_T slOpt )

PURPOSE:
	Constructor.

PARAMETERS:
	Int_T jj
		Variable index.

	const Array<Bool_T> *ExcludeCols
		Array of column previously excluded.

	const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len
		The row in which a free singleton variable is found.

	Real_T af
		The coefficient by the free singleton.

	Real_T b
		The right hand side value for the row.

	Real_T slOpt
		Slack's value.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

FreeSingletonColumnRemoval::FreeSingletonColumnRemoval( Int_T jj, // )
	const Array<Bool_T> *ExcludeCols,
	const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
	Real_T af, Real_T b, Real_T slOpt )
	: j( jj ), A(), Col(), Len( 0 ), Af( af ), B( b ), SlOpt( slOpt )
{
	assert( jj >= 0 );

	Int_T i, RowLen;

	for( i = RowLen = 0; i < len; i++ )
		if( col[i] != j && ( !ExcludeCols || !(*ExcludeCols)[ col[i] ] ) )
			RowLen++;

	Len = RowLen;
	A.Resize( RowLen );
	Col.Resize( RowLen );
	for( i = 0; i < len; i++ )
		if( col[i] != j && ( !ExcludeCols || !(*ExcludeCols)[ col[i] ] ) )
		{
			RowLen--;
			A[RowLen]	= a[i];
			Col[RowLen]	= col[i];
		}
}


/*------------------------------------------------------------------------------

	Bool_T FreeSingletonColumnRemoval::Write( FILE *fp )

PURPOSE:
	Writes the information stored during the free singleton variable removal to
	the output file.

PARAMETERS:
	FILE *fp
		Output file.
	
RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T FreeSingletonColumnRemoval::Write( FILE *fp )
{
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	First line - marks the beginning of a new description.
	//
	if( fprintf( fp, "%-8s  %-8d  %12s\n", "ACTION", (int)j, "FREE_SINGL" )
		== EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Slack's optimal value (if non-zero).
	//
	if( IsNonZero( SlOpt ) && fprintf( fp, "  %-8s  %.12G\n", "SLACK", SlOpt )
		== EOF )
		return False;

	//--------------------------------------------------------------------------
	//	The RHS value (always one - regardless of the row type).
	//
	if( fprintf( fp, "  %-8s  %.12G\n", "RHS", B ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Free singleton's non-zero value.
	//
	if( fprintf( fp, "  %-8s  %.12G\n", "VALUE", Af ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	The row (all non-zeros coefficients).
	//
	for( Int_T i = 0; i < Len; i++ )
		if( fprintf( fp, "  %-8s  %-8d  %.12G\n", "COEFF", int(Col[i]), A[i] )
			== EOF )
			return False;

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T FreeSingletonColumnRemoval::Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols )
		

PURPOSE:
	x

PARAMETERS:
	Solution &solution
		Solution structure.

	Real_T FTOL
		Ignored.

	Array<Bool_T> &ExcludeCols
		"ExcludeCols" is used to check if we're not trying to write the same
		value twice. it is updated on exit to indicate that the variable that
		was just restored is no longer "excluded".


RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T FreeSingletonColumnRemoval::Undo( Solution &solution, Real_T, // )
	Array<Bool_T> & )
{
	Real_T &x = solution.x[j];

	x = B - SlOpt;

	for( Int_T i = 0; i < Len; i++ )
	{
		assert( j != Col[i] );
		x -= A[i] * solution.x[Col[i]];
	}

	assert( IsNonZero( Af ) );

	x /= Af;

	return True;
}
