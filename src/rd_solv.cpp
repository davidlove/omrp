/*------------------------------------------------------------------------------
MODULE TYPE:		Subproblem solution routines.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	rd_solv.cpp
CREATED:			1994.12.07
LAST MODIFIED:		1996.09.16

DEPENDENCIES:		rd_solv.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/


#ifndef __RD_SOLV_H__
#	include "rd_solv.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif


/*------------------------------------------------------------------------------

	void RD_SubproblemSolver::GetFeasibilityCut( Real_T &value, Real_T *grad,
		Int_T n, const Scenario &Scen )
	void RD_SubproblemSolver::GetOptimalityCut( Real_T &value, Real_T *grad,
		Int_T n, const Scenario &Scen )

PURPOSE:
	These functions are used to generate an outer approximation of subproblem's
objective function's epigraph. A feasibility cut is produced when a subproblems
is found to be infeasible. It cuts off a portion of the feasible region of the
first stage problem. An optimality cut (calculated when a subproblem has an
optimal solution) is a support of the objective function at the current trial
point.

PARAMETERS:
	Real_T &value, Real_T *grad, Int_T n
		Optimality cut: objective's value and support's gradient.
		Feasibility cut: distance to the feasible region and cutting plane's
			normal.
		"n" is the dimension of the first stage variables.

	const Scenario &Scen
		A reference to the scenario repository (needed in calculations).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void RD_SubproblemSolver::GetFeasibilityCut( Real_T &value, Real_T *grad, // )
	Int_T n, const Scenario &Scen )
{
	assert( grad != NULL );

	//--------------------------------------------------------------------------
	//	Here we generate a feasibility cut. For that purpose we need to know
	//	which row of the constraint matrix was most infeasible and what was the
	//	infeasibility.
	//

	//
	//	Find the most violated constraint corresponding to a non-zero basic
	//	artificial variable. Put the index in 'infrow'.
	//
	Int_T infrow = -1;
	value = 0.0;

	for( Int_T j = Int_T( LP.GetStructN() + LP.GetSlackN() ); j < N; j++ )
		if( VarType[j] & VT_ARTIF && x[j] > FEASIBILITY_TOL )
			if( x[j] > value && A2B[j] >= 0 )
			{
				infrow = A2B[j];
				value = x[j];
			}

	//
	//	It is possible that the problem is infeasible, but all artificial
	//	variables (some of them have non-zero values) are non-basic. In such
	//	case a dual optimal solution to the first stage problem ('y_t') is
	//	used to generate the cut.
	//
	if( infrow < 0 )
	{
#ifndef NDEBUG
		Bool_T found = False;

		for( Int_T i = 0; i < M; i++ )
			if( IsNonZero( y_t[i] ) )
			{
				found = True;
				break;
			}
		if( !found )
			FatalError( "Internal error when generating a feasibility cut." );
#endif
		for( Int_T j = Int_T( LP.GetStructN() + LP.GetSlackN() ); j < N; j++ )
			value += fabs( x[j] );
		CalculateGradient( grad, n, y_t, Scen );
	}
	else
	{
		//----------------------------------------------------------------------
		//	Extract the appropriate row of the basis inverse into "pi" vector.
		//	Use "mark" work vector as a work space for sparsity pattern.
		//
		WorkVector<Real_T> pi( M );
		WorkVector<Int_T> mark( M );

		pi.Fill( 0.0, M );	mark.Fill( 0, M );
		pi[ infrow ] = 1.0;	mark[ infrow ] = 1;
		B->SparseBTRAN( pi, mark );

		CalculateGradient( grad, n, pi, Scen );
	}
}


//
//	Calculate a cut as
//		g_omega = -Transp( T_omega ) * pi
//	where
//		g_omega				is the cutting plane,
//		Transp( T_omega )	is technology matrix transpose,
//		pi					is either a row of basis inverse, or a dual
//							vector of the first stage problem.
//
void RD_SubproblemSolver::CalculateGradient( Real_T *grad, Int_T n, // )
	Array<Real_T> &pi, const Scenario &Scen )
const
{
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T j, len;

	//--------------------------------------------------------------------------
	//	First calculate the basis part of the gradient vector.
	//
	for( j = 0; j < n; j++ )
	{
		Real_T &g = grad[j];
		
		g = 0.0;
		for( LP.T_base.GetColumn( j, a, row, len ); len; --len, ++a, ++row )
			g -= *a * pi[*row];
	}

	//--------------------------------------------------------------------------
	//	Then calculate the scenario-dependent part.
	//
	Int_T l = Scen.GetLength();
	for( j = 0; j < l; j++ )
		for( Int_T i = 0, bl = Scen[j].Len(); i < bl; i++ )
		{
			const Delta &d = Scen[j][i];

			if( d.type == Delta::MATRIX )
			{
				assert( d.col >= 0 && d.col < n );
				grad[ d.col ] -= pi[ d.row ] * d.value;
			}
		}
}


void RD_SubproblemSolver::GetOptimalityCut( Real_T &value, Real_T *grad, // )
	Int_T n, const Scenario &Scen )
{
	assert( grad != NULL );

	//----------------------------------------------------------------------
	//	Time to generate an optimality cut g_omega such that:
	//		g_omega = - Transp( T_omega ) * pi_opt
	//	where
	//		g_omega				is the gradient,
	//		Transp( T_omega )	is technology matrix transpose and
	//		pi_opt				is the vector of optimal dual variables.
	//
	CalculateGradient( grad, n, y, Scen );

	//--------------------------------------------------------------------------
	//	Compute the duality gap - verify the objective value in this manner.
	//
	assert( DualityGap() < 1e-6 );
	value = Result;
}
