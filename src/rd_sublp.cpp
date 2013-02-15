/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	rd_sublp.cpp
CREATED:			1995.07.27
LAST MODIFIED:		1996.04.05

DEPENDENCIES:		rd_sublp.h, smplx_lp.h, solv_lp.h, mps_lp.h, scenario.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __RD_SUBLP_H__
#	include "rd_sublp.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif


void RD_SubproblemLP::InitializeRD_Subproblem( void )
{
	if( State == UNINIT )
	{
		n1st	= T_base.GetN();
		m2st	= T_base.GetM();

		assert( m2st == m );

		d_base.Resize( m2st );		d_base.Copy( b, m2st, m2st );
		q_base.Resize( n );			q_base.Copy( c, n, n );
		h.Resize( m2st );			h.Fill( 0.0, m2st );
		q.Resize( n );				q.Fill( 0.0, n );
		T_base_y.Resize( m2st );	T_base_y.Fill( 0.0, m2st );

		State = INIT;
	}
}


/*------------------------------------------------------------------------------

	void RD_SubproblemLP::ApplyScenario( const Scenario &Sc,
		Bool_T NewTrialPoint, Int_T nn, const Real_T *TrialPoint )

PURPOSE:
	Modify the underlying "SimplexLP" linear subproblem according to the
scenario.

PARAMETERS:
	const Scenario &Sc
		Scenario.

	Bool_T NewTrialPoint
		"True" if the trial point has been changed.

	Int_T nn, const Real_T *TrialPoint
		First stage variables (as a vector).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void RD_SubproblemLP::ApplyScenario( const Scenario &Sc,  // )
	Bool_T NewTrialPoint,
#ifndef NDEBUG
	Int_T nn,
#else
	Int_T,
#endif
	const Real_T *TrialPoint )
{
	assert( nn == n1st && TrialPoint != NULL );

	//--------------------------------------------------------------------------
	//	Compute the new right hand side and cost to match the current scenario.
	//
	//	Algorithm:
	//	1.	h <- d_base; q <- q_base
	//	2.	h += Delta_d; h -= Delta_T * y; q += Delta_q
	//	3.	if( y changed )
	//			T_base_y <- T_base * y
	//	4.	h -= T_base_y
	//
	//	Step 2 is performed as one because Delta_d, Delta_T and Delta_q are
	//	stored together on the same list corresponding to one scenario.
	//
	h.Copy( d_base, m2st, m2st );
	q.Copy( q_base, n, n );

	//
	//	Now add the scenario-specific modifications.
	//
	for( Int_T l = Sc.GetLength(), j = 0; j < l; j++ )
	{
		for( Int_T i = 0, bl = Sc[j].Len(); i < bl; i++ )
		{
			const Delta &d = Sc[j][i];

			switch( d.type )
			{
			case Delta::RHS:
				h[d.row] += d.value;
				break;

			case Delta::MATRIX:
				assert( d.col >= 0 && d.col <= nn );
				h[d.row] -= d.value * TrialPoint[d.col];
				break;

			case Delta::COST:
				q[d.col] = d.value;
				break;

			default:
#ifndef NDEBUG
				abort();
#endif
				break;
			}
		}
	}

	//--------------------------------------------------------------------------
	//	Check if 'TrialPoint' has changed. If so - recompute 'T_base_y'.
	//
	if( NewTrialPoint )
	{
		T_base_y.Fill( 0.0, m2st );

		Ptr<Real_T> a;
		Ptr<Int_T> row;
		Int_T len;

		for( Int_T j = 0; j < n1st; j++ )
			for( T_base.GetColumn( j, a, row, len ); len; --len, ++a, ++row )
				T_base_y[*row] += *a * TrialPoint[j];
	}

	for( Int_T i = 0; i < m2st; i++ )
		h[i] -= T_base_y[i];

	//--------------------------------------------------------------------------
	//	Replace the right hand side and cost of the LP with computed ones.
	//
	SetRHS( h );
	c.Copy( q, n, n );
}
