/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solvpiv.cpp
CREATED:			1993.10.24
LAST MODIFIED:		1995.11.24

DEPENDENCIES:		solver.h, stdtype.h, solv_lp.h, std_tmpl.h, error.h,
					lp_codes.h, mps_lp.h, inverse.h, solvcode.h 
					<math.h>

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
	COMP_XXXXX			- x

------------------------------------------------------------------------------*/


#include <math.h>

#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif


/*------------------------------------------------------------------------------

	void Solver::OutGoingColumnNumber( Int_T q, Int_T &p, Real_T &theta,
		Int_T &pBound )

PURPOSE:
	Employs Harris rule to find a variable that will be removed from the basis
and replaced by column 'q'. There are some additional criteria for variable
removal (in case of ties).
	Artificial variables are given preference for ousting.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Real_T Solver::OutGoingColumnNumber( Int_T q, Real_T z_q, Int_T &p, // )
	Real_T &ThetaMin, Int_T &pBound )
{
	Int_T i, j, b2a;
	Real_T Theta, a, Pivot = 0.0;
	Short_T bm;
	pBound = A2B_UNDEF;

	//--------------------------------------------------------------------------
	// First run: estabilishing maximum steplength 'ThetaMin'.
	//
	for( j = 0, ThetaMin = INFINITY, p = SLV_UNBOUNDED; j < w1Len; j++ )
	{
		i	= w1Ind[j];
		b2a	= B2A[i];
		bm	= VarType[b2a];
		a	= w1[j];

		if( ( bm & VT_FX ) && FeasibilityRestoreCount > 0 )
		{
			//------------------------------------------------------------------
			//	In the feasibility restoring phase do not move artificial
			//	variables beyond their simple bounds.
			//
			if( ( bm & VT_HAS_LO_BND ) && a >= PIVOT_TOL )
				Theta = x[b2a] / a;
			else if( ( bm & VT_HAS_UP_BND ) && a <= -PIVOT_TOL )
				Theta = ( x[b2a] - u[b2a] ) / a;
			else
				continue;
		}
		else
		{
			if( ( bm & VT_ARTIF ) && ( bm & VT_FX ) && fabs( a ) > PIVOT_TOL )
				// Do not move fixed artificial variables!
				Theta = 0.0;
			else if( bm & VT_HAS_LO_BND && a >= PIVOT_TOL )
				Theta = ( x[b2a] + FEASIBILITY_TOL ) / a;
			else if( bm & VT_HAS_UP_BND && a <= -PIVOT_TOL )
				Theta = ( x[b2a] - FEASIBILITY_TOL - u[b2a] ) / a;
			else
				continue;
		}
		if( Theta > ThetaMin ) continue;

		//----------------------------------------------------------------------
		//	Due to numerical roundoff error it is possible that 'Theta' is < 0
		//	or very small. Both those tendencies need to be eliminated.
		//	If steplength reaches zero, there's no point in continuing the loop.
		//
		if( Theta < 0.0 ) Theta = 0.0;
		pBound = ( a > 0.0 ) ? A2B_LO: A2B_UP;
		p = i;
		Pivot = a;
		ThetaMin = Theta;
		if( ThetaMin == 0.0 ) break;
	}

	//--------------------------------------------------------------------------
	// Here we handle the situations in which the in-coming variable will sooner
	// reach its opposite bound, than move one of the basic variables to their
	// bounds.
	//
	if( ( A2B[q] == A2B_IN && ThetaMin > 0.0 ) ||
		( ThetaMin > 0.0 && ( VarType[q] & VT_HAS_UP_BND ) &&
		( VarType[q] & VT_HAS_LO_BND ) ) )
	{
		Real_T uq = u[q], xq = x[q];
		Short_T vt = VarType[q];

		switch( A2B[q] )
		{
		case A2B_UP:
		case A2B_LO:
			if( uq <= ThetaMin )
			{
				p = SLV_TO_BND;
				ThetaMin = uq;
			}
			break;

		case A2B_IN:
			if( vt & VT_HAS_LO_BND && z_q > 0.0 && xq < ThetaMin )
			{
				ThetaMin = xq;
				p = SLV_TO_BND;
			}
			else if( vt & VT_HAS_UP_BND && z_q < 0.0 && uq - xq < ThetaMin )
			{
				ThetaMin = uq - xq;
				p = SLV_TO_BND;
			}
			break;

		default:
#			ifndef NDEBUG
				FatalError( "solvpiv.cpp: OutGoingColumnNumber: "
					"Invalid variable enters the basis." );
#			endif
			break;
		}
	}

	if( p < 0 ) return 0.0;

	//--------------------------------------------------------------------------
	//	To have a warranty that:
	//	1.	There will be a basic column candidate on the list.
	//	2.	A non-zero step will be preferred if possible (regardless of
	//		Harris' rule).
	//	We first set the minimum step length and its corresponding row number.
	//
	Real_T MaxA = 0.0;			// The best candidate's pivot absolute value.

	//--------------------------------------------------------------------------
	//	Second run: we find biggest alpha element among those that would not
	//	violate the maximum steplength calculated in the previous run.
	//
	for( j = 0; j < w1Len; j++ )
	{
		i	= w1Ind[j];
		b2a	= B2A[i];
		bm	= VarType[b2a];
		a	= w1[j];

		if( bm & VT_HAS_LO_BND && a >= PIVOT_TOL )
			Theta = x[b2a] / a;
		else if( bm & VT_HAS_UP_BND && a <= -PIVOT_TOL )
			Theta = ( x[b2a] - u[b2a] ) / a;
		else
			continue;

		if( Theta > ThetaMin ) continue;			// Observe steplength.

		//----------------------------------------------------------------------
		//	We check if this candidate is any better than the previous ones.
		//	Artificial and fixed columns are preferred as candidates for leaving
		//	the basis (a weight factor of 10 is used).
		//
		if( fabs( ( bm & ( VT_ARTIF | VT_FX ) ) ? a : a * 10.0 ) > MaxA )
		{
			MaxA = fabs( a );
			ThetaMin = Theta;
			p = i;
			pBound = ( a > 0.0 ) ? A2B_LO: A2B_UP;
			Pivot = a;
		}
	}

	return Pivot;
}
