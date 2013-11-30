/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	row_anal.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1995.10.24

DEPENDENCIES:		std_math.h, std_tmpl.h, presolve.h, simplex.h, error.h,
					smartptr.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h,
					stdtype.h, lp_codes.h, solv_lp.h, compile.h, mps_lp.h,
					sort_lab.h, myalloc.h, cl_list.h

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

#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif


/*------------------------------------------------------------------------------

	Int_T Presolver::ForcingAndDominatedRows( void )

PURPOSE:
	x

PARAMETERS:
	None.

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Int_T Presolver::ForcingAndDominatedRows( void )
{
	Int_T ElNum = 0;

	for( Int_T rlen = 2; rlen < n; rlen++ )
		for( Int_T i = RowList.GetFirst( rlen ); i >= 0; i = RowList.GetNext() )
		{
			Real_T AL, AU;
			Short_T RT;

			ComputeRowActivityLimits( i, AL, AU, RT );

			if( AL > bu[i] )
			{
				PrimalInf = Max( AL - bu[i], PrimalInf );
				if( AL > bu[i] + ftol )
				{
					Status = LPS_INFEASIBLE;
					return ElNum;
				}
				else
					AL = bu[i];
			}

			if( AU < bl[i] )
			{
				PrimalInf = Max( bl[i] - AU, PrimalInf );
				if( AU < bl[i] - ftol )
				{
					Status = LPS_INFEASIBLE;
					return ElNum;
				}
				else
					AU = bl[i];
			}

			if( IsDominatedRow( i, AL, AU, RT ) )
			{
				DominatedRows++;
				ElNum++;
			}
			else if( IsForcingRow( i, AL, AU, RT ) )
			{
				ForcingRows++;
				ElNum++;
			}
			else
				VariableBoundsTightened += TightenBounds( i, AL, AU, RT );

			if( Status != LPS_UNKNOWN ) return ElNum;
		}

	return ElNum;
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::IsDominatedRow( Int_T row, Real_T AL, Real_T AU,
		Short_T RT )

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

Bool_T Presolver::IsDominatedRow( Int_T row, Real_T AL, Real_T AU, Short_T RT )
{
	if( ( ! ( rt[ row ] & VT_LO ) || ( ( RT & VT_LO ) && bl[ row ] <= AL ) ) &&
		( ! ( rt[ row ] & VT_UP ) || ( ( RT & VT_UP ) && bu[ row ] >= AU ) ) )
	{
		RemoveRow( row );
		return True;
	}
	else
		return False;
}


/*------------------------------------------------------------------------------

	Bool_T Presolver::IsForcingRow( Int_T row, Real_T AL, Real_T AU,
		Short_T RT )

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

Bool_T Presolver::IsForcingRow( Int_T row, Real_T AL, Real_T AU, Short_T RT )
{
	Ptr<Int_T> col;
	Ptr<Real_T> a;
	Int_T len;

	if( ( rt[ row ] & VT_LO ) && ( RT & VT_UP ) && bl[ row ] == AU )
		for( LP.GetRow( row, a, col, len ); len; ++a, ++col, --len )
		{
			if( ExcludeCols[ *col ] ) continue;
			
			FixVariable( *col, ( *a > 0.0 ) ? VT_UP : VT_LO );
			if( Status != LPS_UNKNOWN ) return True;
		}
	else if( ( rt[ row ] & VT_UP ) && ( RT & VT_LO ) && bu[ row ] == AL )
		for( LP.GetRow( row, a, col, len ); len; ++a, ++col, --len )
		{
			if( ExcludeCols[ *col ] ) continue;

			FixVariable( *col, ( *a > 0.0 ) ? VT_LO : VT_UP );
			if( Status != LPS_UNKNOWN ) return True;
		}
	else
		return False;

	return True;
}


/*------------------------------------------------------------------------------

	Int_T TightenBounds( Int_T row, Real_T AL, Real_T AU, Short_T RT )

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

Int_T Presolver::TightenBounds( Int_T row, Real_T AL, Real_T AU, Short_T RT )
{
	assert( !ExcludeRows[row] );

	//--------------------------------------------------------------------------
	//	At least one finite bound needed to go on.
	//
	if( !( RT & ( VT_LO | VT_UP ) ) )
		return 0;

	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Short_T _rt	= rt[row];
	Int_T len, Count	= 0;

	//--------------------------------------------------------------------------
	//	Loop on row non-zeros (those that have not been eliminated yet) and try
	//	to compute implied bounds on all the variables.
	//
	for( LP.GetRow( row, a, col, len ); len; ++a, ++col, --len )
		if( !ExcludeCols[ *col ] )
		{
			Short_T vt = LP.GetVarType( *col );

			//------------------------------------------------------------------
			//	At least one finite bound needed to go on.
			//
			if( !( vt & ( VT_LO | VT_UP ) ) )
				continue;

			Short_T ivt	= 0;
			Real_T il	= -INFINITY,
				iu		= +INFINITY;

			assert( IsNonZero( *a ) );

			//------------------------------------------------------------------
			//	Use one of the formulas to compute the implied bounds (depending
			//	on the sign of the non-zero coefficient).
			//
			if( *a > 0 )
			{
				if( ( _rt & VT_LO ) && ( RT & VT_UP ) && ( vt & VT_UP ) )
				{
					il = ( bl[row] - AU ) / *a + LP.GetU( *col );
					ivt |= VT_LO;
				}
				if( ( _rt & VT_UP ) && ( RT & VT_LO ) && ( vt & VT_LO ) )
				{
					iu = ( bu[row] - AL ) / *a + LP.GetL( *col );
					ivt |= VT_UP;
				}
			}
			else
			{
				if( ( _rt & VT_UP ) && ( RT & VT_LO ) && ( vt & VT_UP ) )
				{
					il = ( bu[row] - AL ) / *a + LP.GetU( *col );
					ivt |= VT_LO;
				}
				if( ( _rt & VT_LO ) && ( RT & VT_UP ) && ( vt & VT_LO ) )
				{
					iu = ( bl[row] - AU ) / *a + LP.GetL( *col );
					ivt |= VT_UP;
				}
			}

			//------------------------------------------------------------------
			//	Compare the implied bounds with the original bounds.
			//
			if( ( VT_LO & vt & ivt ) && SetL( *col, il ) )
				Count++;

			if( ( VT_UP & vt & ivt ) && SetU( *col, iu ) )
				Count++;
		}

	return Count;
}


/*------------------------------------------------------------------------------

	void Presolver::ComputeRowActivityLimits( Int_T row, Real_T &infAX,
		Real_T &supAX, Short_T &rt, Int_T xcol )

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

void Presolver::ComputeRowActivityLimits( Int_T row, Real_T &infAX, // )
	Real_T &supAX, Short_T &RT, Int_T xcol )
{
	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T len;

	infAX = supAX = 0.0;
	RT = VT_LO | VT_UP;

	//--------------------------------------------------------------------------
	//	Compute upper and lower bound on row activity (excluding column 'xcol').
	//	Use the following formulas:
	//
	//			  T
	//	a = sup (a x) =  ä   a u  +  ä   a l  ó +ì
	//	 u   x         	a >0  i i   a <0  i i
	//                   i           i
	//	          T
	//	a = inf (a x) =  ä   a l  +  ä   a u  ò -ì
	//	 l   x		   	a >0  i i   a <0  i i
	//				     i           i
	//
	//
	for( LP.GetRow( row, a, col, len ); len; ++a, ++col, --len )
		if( !ExcludeCols[ *col ] && *col != xcol )
			if( *a > 0.0 )
				if( LP.GetVarType( *col ) & VT_UP )
					supAX	+= *a * LP.GetU( *col );
				else
				{
					supAX	= +INFINITY;
					RT		&= ~VT_UP;
					break;
				}
			else
				if( LP.GetVarType( *col ) & VT_LO )
					supAX	+= *a * LP.GetL( *col );
				else
				{
					supAX	= +INFINITY;
					RT		&= ~VT_UP;
					break;
				}

	for( LP.GetRow( row, a, col, len ); len; ++a, ++col, --len )
		if( !ExcludeCols[ *col ] && *col != xcol )
			if( *a > 0.0 )
				if( LP.GetVarType( *col ) & VT_LO )
					infAX	+= *a * LP.GetL( *col );
				else
				{
					infAX	= -INFINITY;
					RT		&= ~VT_LO;
					break;
				}
			else
				if( LP.GetVarType( *col ) & VT_UP )
					infAX	+= *a * LP.GetU( *col );
				else
				{
					infAX	= -INFINITY;
					RT		&= ~VT_LO;
					break;
				}
}
