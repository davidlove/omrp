/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	sing_col.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1995.10.29

DEPENDENCIES:		std_tmpl.h, std_math.h, presolve.h, simplex.h, error.h,
					smartptr.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h,
					stdtype.h, lp_codes.h, solv_lp.h, compile.h, mps_lp.h,
					sort_lab.h, myalloc.h, cl_list.h, postsolv.h, my_defs.h
					<assert.h>

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

	Int_T Presolver::DealWithSigletonColumns( const int Mode )

PURPOSE:
	Loops on all singleton columns and sees whether each one of them is:
a)	a free singleton,
b)	an impled free singleton or
c)	an explict slack (this is optional and depends on the argument ).

For each of the cases one of the elimination procedures is invoked.

PARAMETERS:
	Bool_T ElimExplSlacks
		Decides whether explicit slacks are to be eliminated or not.

RETURN VALUE:
	The number of variables removed from the problem.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::DealWithSigletonColumns( const int Mode )
{
	if( !( Mode & ( LPR_SINGL_COLS | LPR_EXPLICIT_SLACKS ) ) )
		return 0;

	Int_T Elim = 0, j;

	//--------------------------------------------------------------------------
	//	Loop on singletom columns only.
	//
	for( j = ColList.GetFirst( 1 ); j != -1; j = ColList.GetNext() )
	{
		assert( !ExcludeCols[j] && j >= 0 && j < n );

		Ptr<Real_T> a;
		Ptr<Int_T> row;
		Int_T len;

		//----------------------------------------------------------------------
		//	Find the non-zero row and its coefficient.
		//
		for( LP.GetColumn( j, a, row, len ); len; ++a, ++row, --len )
			if( !ExcludeRows[ *row ] ) break;

		// Free singleton column.
		//
		if( ( Mode & LPR_SINGL_COLS ) &&
			! ( LP.GetVarType( j ) & ( VT_LO | VT_UP ) ) )
		{
			EliminateFreeSingletonColumn( *row, j, *a );
			OrigFreeSingletonCols++;
			Elim++;
		}
		else if( ( Mode & LPR_SINGL_COLS ) &&
			IsImpliedFreeSingleton( *row, j, *a ) )
		{
			EliminateFreeSingletonColumn( *row, j, *a );
			ImpliedFreeSingletonCols++;
			Elim++;
		}
		//	"IsImpliedFreeSingleton()" may fix the variable and thus remove it
		//	from the problem.
		//
		else if( ExcludeCols[j] )
			Elim++;
		//	Possible explicit slack removal.
		//
		else if( ( Mode & LPR_EXPLICIT_SLACKS ) && !ExplSlackRemoved[*row] &&
			( rt[*row] == VT_FIXED || IsZero( LP.GetC( j ) ) ) )
		{
			if( EliminateExplicitSlack( *row, j, *a ) )
			{
				RelaxedConstraints++;
				Elim++;
			}
			ExplSlackRemoved[*row] = True;
		}

		if( Status != LPS_UNKNOWN ) break;
	}

	return Elim;
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::EliminateExplicitSlack( Int_T i, Int_T j, Real_T a_ij )

PURPOSE:
	Eliminates explicit slack variables. Adjusts bounds on row activity as well
as cost function coefficients and fixed adjustment.

PARAMETERS:
	Int_T i, Int_T j, Real_T a_ij
		Coordinates and non-zero value of a suspected explicit slack.

RETURN VALUE:
	'True' if the slack was converted, 'False' otherwise.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Presolver::EliminateExplicitSlack( Int_T i, Int_T j, Real_T a_ij )
{
	assert( !ExcludeRows[i] );
	assert( !ExplSlackRemoved[i] );

	Short_T	&cur_rt	= rt[i];
	Real_T	&cur_bl	= bl[i],
		&cur_bu		= bu[i];
	Short_T	old_rt	= cur_rt,
		vt			= LP.GetVarType( j );
	Real_T	old_bl	= cur_bl,
		old_bu		= cur_bu,
		l			= LP.GetL( j ),
		u			= LP.GetU( j ),
		c			= LP.GetC( j );

	//--------------------------------------------------------------------------
	//	Since 'a_ij * l or a_ij * u' will be used for updating the lower
	//	and upper bounds on row activity, the appropriate values have to be
	//	of the same (or similar) order of magnitude. Otherwise numerical errors
	//	will occur.
	//
	{
		Real_T v_bnd	= Max( ( vt & VT_LO ) ? fabs( l ) : 0.0,
						( vt & VT_UP ) ? fabs( u ) : 0.0 ),
			r_bnd		= Max( ( old_rt & VT_LO ) ? fabs( old_bl ) : 0.0,
						( old_rt & VT_UP ) ? fabs( old_bu ) : 0.0 );

		if( fabs( a_ij ) * v_bnd > 100.0 * r_bnd )
			return False;
	}

	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T len;
	
	LP.GetRow( i, a, col, len );

	if( old_rt == VT_FIXED && IsNonZero( c ) )
	{
		Real_T cf_af = c / a_ij;

		f += cf_af * cur_bl;
		for( ; len; ++a, ++col, --len )
			if( !ExcludeCols[ *col ] )
				LP.SetC( *col, LP.GetC( *col ) - cf_af * *a );
	}

	if( IsPositive( a_ij ) )
	{
		if( cur_rt & VT_LO )
			if( vt & VT_UP )
				cur_bl -= a_ij * u;
			else
			{
				cur_bl = -INFINITY;
				cur_rt &= ~VT_LO;
			}

		if( cur_rt & VT_UP )
			if( vt & VT_LO )
				cur_bu -= a_ij * l;
			else
			{
				cur_bu = +INFINITY;
				cur_rt &= ~VT_UP;
			}
	}
	else if( IsNegative( a_ij ) )
	{
		if( cur_rt & VT_LO )
			if( vt & VT_LO )
				cur_bl -= a_ij * l;
			else
			{
				cur_bl = -INFINITY;
				cur_rt &= ~VT_LO;
			}

		if( cur_rt & VT_UP )
			if( vt & VT_UP )
				cur_bu -= a_ij * u;
			else
			{
				cur_bu = +INFINITY;
				cur_rt &= ~VT_UP;
			}
	}

	if( ( cur_rt & VT_LO ) && ( cur_rt & VT_UP ) &&
		IsEqual( cur_bl, cur_bu ) )
		cur_rt |= VT_FX;
	else
		cur_rt &= ~VT_FX;

	if( cur_rt & ( VT_LO | VT_UP ) )
	{
		LP.GetRow( i, a, col, len );
		if( PostSolve )
			PostSolve->ExplicitSlackRemoval( j, a, col, len, a_ij, vt, l, u,
				old_rt, old_bl, old_bu );

		RemoveColumn( j );
	}
	else
		EliminateFreeSingletonColumn( i, j, a_ij );

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::IsImpliedFreeSingleton( Int_T i, Int_T j, Real_T a_ij )

PURPOSE:
	Checks whether the variable 'j' is an implied free singleton (with its only
non-zero in row 'i' having value of 'a_ij').

PARAMETERS:
	Int_T i, Int_T j, Real_T a_ij
		Row and column number, non-zero value of the suspected implied free
		singleton column (see above).

RETURN VALUE:
	'True' if the answer is yes.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Presolver::IsImpliedFreeSingleton( Int_T i, Int_T j, Real_T a_ij )
{
	assert( ColLen[j] == 1 && j >= 0 && j < n );

	Real_T il, iu;
	Short_T vt;

	ComputeImpliedBounds( i, j, a_ij, il, iu, vt );

	//--------------------------------------------------------------------------
	//	We compare implied and explicit bounds on 'x_j' to determine, if the
	//	explicit ones may be dropped. If implied bounds are at least as tight as
	//	the explicit ones, the latter are superfluous. Logical value is
	//	returned.
	//
	if( ( !( LP.GetVarType( j ) & VT_LO ) || il >= LP.GetL( j ) ) &&
		( !( LP.GetVarType( j ) & VT_UP ) || iu <= LP.GetU( j ) ) )
		return True;
	else
	{
		if( vt & VT_LO ) SetL( j, il );
		if( vt & VT_UP ) SetU( j, iu );
	}

	return False;
}


/*------------------------------------------------------------------------------

	void Presolver::ComputeImpliedBounds( Int_T i, Int_T j, Real_T a_ij,
		Real_T &il, Real_T &iu, Short_T &vt )

PURPOSE:
	Compute lower and upper bounds on variable number 'j' implied by the 'i'-th
row. Put the computed values into 'il' and 'iu' reference arguments and their
respective infinity marker 'vt'.

PARAMETERS:
	Int_T i, Int_T j, Real_T a_ij
		Row and column number, non-zero value.

	Real_T &il, Real_T &iu, Short_T &vt
		Implied bounds together with bit mask indicating which values are
		finite.

RETURN VALUE:
	None. Results are returned in the last three reference arguments.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::ComputeImpliedBounds( Int_T i, Int_T j, Real_T a_ij, // )
	Real_T &il, Real_T &iu, Short_T &vt )
{
	assert( j >= 0 && j < n && !ExcludeCols[j] &&
		i >= 0 && i < m && !ExcludeRows[i] );

	//--------------------------------------------------------------------------
	//	First we compute lower and upper bounds on (a_ij * x_j) implied by the
	//	matrix row.
	//	Use the following formulas:
	//
	//	                  T                               T
	//	infAX = b - sup (a x)           supAX = b - inf (a x)
	//	         l   x    i                      u   x    i
	//
	//	and finally obtain:
	//
	//	          T                       T
	//	b - sup (a x) ó a  x  ó b - inf (a x)
	//	 l   x    i      ij j    u   x    i
	//
	//	where a  denotes 'i'-th row of the constraint matrix, b  and b  denote
	//	       i                                               l      u
	//
	//	lower and upper bounds on 'i'-th row activity respectively.
	//
	//	Note opposite order of infAX and supAX arguments. As a result of the
	//	call supAX will be assigned LOWER bound on row activity (excluding
	//	column 'j'), and infAX - UPPER bound. Then the substractions will be
	//	performed. The 'ifs' are there to handle infinite bound values.
	//
	Real_T infAX, supAX;
	Short_T RT;
	Bool_T FiniteLo = True, FiniteUp = True;

	ComputeRowActivityLimits( i, supAX, infAX, RT, j );
	if( ( RT & VT_UP ) && ( rt[ i ] & VT_LO ) )
		infAX		= bl[ i ] - infAX;
	else
	{
		infAX		= -INFINITY;
		FiniteLo	= False;
	}

	if( ( RT & VT_LO ) && ( rt[ i ] & VT_UP ) )
		supAX		= bu[ i ] - supAX;
	else
	{
		supAX		= +INFINITY;
		FiniteUp	= False;
	}

	//--------------------------------------------------------------------------
	//	When we know implied bounds on (a_ij * x_j), we may compute implied
	//	bounds on 'x_j' ('il', 'iu'). We assume, that a_ij is non-zero.
	//
	if( a_ij > 0.0 )
	{
		il = ( FiniteLo ) ? infAX / a_ij : -INFINITY;
		iu = ( FiniteUp ) ? supAX / a_ij : +INFINITY;
		vt = Short_T( ( FiniteLo ? (Short_T) VT_LO : 0 ) |
			( FiniteUp ? (Short_T) VT_UP : 0 ) );
	}
	else
	{
		il = ( FiniteUp ) ? supAX / a_ij : -INFINITY;
		iu = ( FiniteLo ) ? infAX / a_ij : +INFINITY;
		vt = Short_T( ( FiniteLo ? (Short_T) VT_UP : 0 ) |
			( FiniteUp ? (Short_T) VT_LO : 0 ) );
	}

	if( ( vt & VT_LO ) && ( vt & VT_UP ) && IsEqual( il, iu ) )
		vt = VT_FIXED;
}
