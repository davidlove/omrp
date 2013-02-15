/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solvpric.cpp
CREATED:			1993.10.24
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		smartptr.h, solver.h, std_math.h, inverse.h
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
	XXXXX			- x

------------------------------------------------------------------------------*/

#include <math.h>

#ifndef __SOLVER_H__
#	include "solver.h"
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
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif


/*------------------------------------------------------------------------------

	void Solver::InGoingColumnNumber( Int_T &q, VerbLevel Verbosity )

PURPOSE:
	x

PARAMETERS:
	Int_T &q
		This parameter is actually a return value. Number of column found most
		atractive is returned. If no candidates are found, 'q' is set to
		'SLV_OPTIMUM'.

RETURN VALUE:
	None.

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

void Solver::InGoingColumnNumber( Int_T &q, VerbLevel Verbosity )
{
	//--------------------------------------------------------------------------
	//	Try it the easy way. Should work most of the time!
	//
	if( FindColumnCandidates( q ) != SLV_OPTIMUM )
		return;

	//--------------------------------------------------------------------------
	//	Tighten the tolerances and (if necessary) restore feasibility of the
	//	solution by adding artificial variables.
	//
	int ResidCheckMode = CHK_PRIM;
	PrimalResiduals = CheckResiduals( ResidCheckMode );
	if( PrimalResiduals > FEASIBILITY_TOL )
		UpdateBasis();
	FinalTolerances( Verbosity );
	if( FeasibilityRestoreCount++ < MAX_FEAS_RESTORE )
		MakeSolutionFeasible( FEASIBILITY_TOL, Verbosity, False );
	FixArtificialVariables();

	ComputeDualVariables();
	ComputeReducedCosts();

	if( FindColumnCandidates( q ) != SLV_OPTIMUM )
		return;

	//--------------------------------------------------------------------------
	//	If it still was no good - tighten tolerances and repeat the previous
	//	process.
	//
	//	Before computing the new reduced costs, see if there are any non-zero
	//	artificial non-basic variables. If so, increase the Big M value (using
	//	'SimplexLP' member function, to which the dual variables will be passed.
	//
	ToBounds();
	ComputePrimalVariables();

	Bool_T NZ_Art	= NonZeroArtificials(),
		BasArt		= BasicArtificials();

	if( NZ_Art && !BasArt )
		FatalError( "Internal error in pricing routine." );

	if( NZ_Art || BasArt )
	{
		AltPricCnt++;

		ComputeSplitDualVariables();

		Real_T MinM = MinimumPenaltyToMove( BasArt );
    	if( MinM <= LP.GetPenalty() )
		{
        	q = ( NZ_Art ) ? SLV_PROBLEM_INFEASIBLE : SLV_OPTIMUM;
			return;
		}

		PenaltyAdjustCnt++;
		LP.SetPenalty( 10.0 * MinM );
		ComputeDualVariables( LP.GetPenalty() );

		ComputeReducedCosts();
		ComputeResult();
		FindColumnCandidates( q );
		assert( q >= 0 );
	}
	else
	{
		ResidCheckMode = CHK_DUAL;
		DualResiduals = CheckResiduals( ResidCheckMode );
		q = SLV_OPTIMUM;
	}
}


/*------------------------------------------------------------------------------

	int Solver::FindColumnCandidates( Int_T &qMin )

PURPOSE:
	xxx 
	Assumes that 'pi' contains current dual variables.

PARAMETERS:
	x

RETURN VALUE:
	If column cadidate is found 'SLV_OPTIMUM' is returned. Otherwise
	'SLV_CANDIDATE_FOUND' is returned and 'q' and 'rc' are set correspondingly.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SLV_FCC Solver::FindColumnCandidates( Int_T &qMin )
{
	Bool_T FoundMarked;
	Real_T w = 0.0, wMin;

	do
	{
		FoundMarked = False;
		qMin = SLV_OPTIMUM;
//		wMin = -1.0;
		wMin = OPTIMALITY_TOL;

		//----------------------------------------------------------------------
		//	Find best weighted reduced cost (optional steepest edge or
		//	approx. steepest edge pricing).
		//
		for( Int_T j = 0; j < N; j++ )
		{
			Int_T a2b = A2B[j];

			//------------------------------------------------------------------
			//	Skip basic and fixed variables. Skip marked variables and note
			//	that you have encountered any.
			//
			if( a2b >= 0 || VarType[j] & VT_FX ) continue;

			if( VarType[j] & VT_MARKED ) { FoundMarked = True; continue; }

			switch( a2b )
			{
			case A2B_UP:	w = -z[j];			break;
			case A2B_LO:	w = z[j];			break;
			case A2B_IN:	w = -fabs( z[j] );	break;
			}

			if( w >= 0.0 ) continue;

			switch( Pricing )
			{
			case PRS_SE:			// Steepest edge or approximate
			case PRS_ASE:			// steepest edge.
				w = w * w / gamma[j];	
				break;

			case PRS_RC:
				w = -w;
				break;
			}

			//------------------------------------------------------------------
			//	Introduce into basis variables with a profitable reduced cost
			//	or non-zero artificial variables with a reduced cost of zero.
			//
//			if( ( w > OPTIMALITY_TOL ||
//				( ( VarType[j] & VT_ARTIF ) && A2B[j] != A2B_LO ) ) &&
//				w > wMin )
			if( w > wMin )
			{
				wMin = w;
				qMin = j;
			}
		}

		//----------------------------------------------------------------------
		//	If no more candidates to enter the basis, but some columns were
		//	marked we need to unmark them and search again.
		//
		if( qMin == SLV_OPTIMUM && FoundMarked ) UnMark();						

	} while( FoundMarked );

	return ( qMin == SLV_OPTIMUM ) ? SLV_OPTIMUM: SLV_CANDIDATE_FOUND;
}


/*------------------------------------------------------------------------------

	void Solver::UnMark( void )

PURPOSE:
	Dangerous variables may have been marked while solving the problem. If this
	should happen, they will also have to be unmarked before the solver finally
	reports optimal solution.

	This function unmarks any variables that may have been previously marked as
	dangerous and thus skipped during a search for column candidates to enter
	the basis. Otherwise the 'VarType' vector remains unchanged.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	The 'VarType' table is updated, bits corresponding to 'VT_MARKED' are
	cleared.

------------------------------------------------------------------------------*/

void Solver::UnMark( void )
{
	for( Int_T i = 0; i < N; i++ ) VarType[i] &= VT_UNMARK;
}



/*------------------------------------------------------------------------------

	Bool_T Solver::NonZeroArtificials( void )
	Bool_T Solver::BasicArtificials( void )

PURPOSE:
	The first function checks if all artificial variables are already reduced
to zero. Since they were assigned positive values only in order to reflect
infeasibility, only positive values are to be taken into account. On the other
hand significant negative values are unacceptable.
	The other one checks whether there are any artificial variables in the
basis.

PARAMETERS:
	None.

RETURN VALUE:
	Returns 'True' if some artificial variables are non-zero, 'False' otherwise.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Solver::NonZeroArtificials( void )
{
	for( Int_T j = 0; j < N; j++ )
		if( VarType[j] & VT_ARTIF )
		{
//			assert( x[j] > -FEASIBILITY_TOL );
			if( x[j] > FEASIBILITY_TOL )
				return True;
		}

	return False;
}


Bool_T Solver::BasicArtificials( void )
{
	for( Int_T i = 0; i < M; i++ )
		if( VarType[B2A[i]] & VT_ARTIF )
			return True;

	return False;
}


/*------------------------------------------------------------------------------

	void Solver::ToBounds( void )

PURPOSE:
	Moves the non-basic variables which should be between their bounds and are
	not (they violate their simple bounds by ~FEASIBILITY_TOL) back into their
	range.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ToBounds( void )
{
	for( Int_T j = 0; j < N; j++ )
		switch( A2B[j] )
		{
		case A2B_IN:
			//------------------------------------------------------------------
			//	This complicated expression is here because we demand that
			//	non-basic variables which are between their bounds sholudn't
			//	even reach the bounds without changing their type, let alone
			//	exceed them.
			//
			assert( 
					!( ( VarType[j] & VT_HAS_LO_BND ) && x[j] < 0.0 ) &&
					!( ( VarType[j] & VT_HAS_UP_BND ) && x[j] > u[j] )
				);
			break;

		case A2B_UP:
			assert( VarType[j] & VT_HAS_UP_BND );
			x[j]	= u[j];
			break;

		case A2B_LO:
			assert( VarType[j] & VT_HAS_LO_BND );
			x[j]	= 0.0;
			break;
		}
}


/*------------------------------------------------------------------------------

	void Solver::ComputeSplitDualVariables( void );

PURPOSE:
	Calculates the split dual variables (the results are stored in 'y_t' and
'y_x' Solver class data members).   

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::ComputeSplitDualVariables( void )
{
	y_x.Fill( 0.0, M );
	y_t.Fill( 0.0, M );

	for( Int_T i = 0; i < M; i++ )
	{
		Int_T j = B2A[i];

		if( VarType[j] & VT_ARTIF )
			y_t[i] = IsNonZero( LP.GetC( j ) ) ? 1.0 : 0.0;
		else
			y_x[i] = LP.GetC( j );
	}

	B->DenseBTRAN( y_x );
	B->DenseBTRAN( y_t );

#ifndef NDEBUG
	{
		const Real_T p = LP.GetPenalty();

		for( Int_T i = 0; i < M; i++ )
			assert( IsZero( y[i] - p * y_t[i] - y_x[i] ) );
	}
#endif
}


/*------------------------------------------------------------------------------

	Real_T Solver::MinimumPenaltyToMove( const Bool_T DependsOnArtif )
	Real_T Solver::MinimumPenaltyToStayBounded( Int_T j,
		const Bool_T DependsOnArtif )

PURPOSE:
	The first function: determine the value of penalty parameter 'M' which
would allow to move the values of as many artificial variables as possible.
	The second one: if possible increase the value of the penalty so that
too small a value of the penalty would not cause unboundedness.

PARAMETERS:
	Int_T j
		Variable which generated the unbounded objective decrease direction.

	const Array<Real_T> &y_x
	const Array<Real_T> &y_t
		Arrays of dual variables computed separately for the original and the
		artificial basic variables.

	const Bool_T DependsOnArtif
		"True" if there are artificial variables in the basis.

RETURN VALUE:
	The minimum value of the penalty parameter required or (-1) if no the
problems is infeasible / unbounded.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T Solver::MinimumPenaltyToMove( const Bool_T DependsOnArtif )
{
	Real_T MinPenalty = -1.0;

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	//--------------------------------------------------------------------------
	//	Loop on all columns; skip basic variables.
	//
	for( Int_T j = 0; j < N; j++ )
	{
		//
		//	Skip basic and fixed variables.
		//
		if( A2B[j] >= 0 || ( VarType[j] & VT_FX ) ) continue;

		if( ( VarType[j] & VT_ARTIF ) && IsPositive( x[j] ) )
		//
		//	See if 'j'-th artificial variable's reduced cost may be made
		//	non-zero.
		//
		{
			LP.GetColumn( j, a, row, len );
			assert( len == 1 );
			assert( *a == 1.0 || *a == -1.0 );

			Int_T i = *row;

			if( *a > 0.0 )
			{
				if( IsNotEqual( y_t[i], 1.0 ) )
				{
					Real_T Penalty = y_x[i] / ( 1.0 - y_t[i] );

					if( !( A2B[j] == A2B_UP && y_t[i] > 1.0 ) )
					{
						assert( 1.01 * Penalty > LP.GetPenalty() );
						MinPenalty = Max( MinPenalty, Penalty );
					}
				}
			}
			else
			{
				if( IsNotEqual( y_t[i], -1.0 ) )
				{
					Real_T Penalty = - y_x[i] / ( 1.0 + y_t[i] );

					if( !( A2B[j] == A2B_UP && y_t[i] < -1.0 ) )
					{
						assert( 1.01 * Penalty > LP.GetPenalty() );
						MinPenalty = Max( MinPenalty, Penalty );
					}
				}
			}
		}
		else if( DependsOnArtif )
		//
		//	See if 'j'-th structural variable may be moved to basis.
		//	If 'y_t^T a_j == 0' then no structural variable can be moved into
		//	the basis by changes to penalty 'M' (then it's reduced cost does
		//	not depend on 'M').
		//
		{
			Real_T y_t_a = 0.0, y_x_a = LP.GetC( j );

			for( LP.GetColumn( j, a, row, len ); len; --len, ++a, ++row )
			{
				y_t_a += *a * y_t[*row];
				y_x_a -= *a * y_x[*row];
			}

			if( IsZero( y_t_a ) ) continue;
			y_x_a /= y_t_a;

			//------------------------------------------------------------------
			//	This is the proposed penalty value (we still have to check, if
			//	'y_t_a' has appropriate value). All formulas for lower bound on
			//	penalty eventually come to this formula.
			//
			switch( A2B[j] )
			{
			case A2B_LO:
				if( IsPositive( y_t_a ) )
				{
					assert( 1.01 * y_x_a > LP.GetPenalty() );
					MinPenalty = Max( MinPenalty, y_x_a );
				}
				break;

			case A2B_UP:
				if( IsNegative( y_t_a ) )
				{
					assert( 1.01 * y_x_a > LP.GetPenalty() );
					MinPenalty = Max( MinPenalty, y_x_a );
				}
				break;

			case A2B_IN:
				if( IsNonZero( y_t_a ) )
				{
					assert( 1.01 * y_x_a > LP.GetPenalty() );
					MinPenalty = Max( MinPenalty, y_x_a );
				}
				break;
			}
		}
	}

	return MinPenalty;
}


Real_T Solver::MinimumPenaltyToStayBounded( Int_T j, // )
	const Bool_T DependsOnArtif )
{
	Real_T MinPenalty = -1.0;

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	assert( A2B[j] < 0 && !( VarType[j] & VT_FX ) );

	if( ( VarType[j] & VT_ARTIF ) && IsPositive( x[j] ) )
	//
	//	See if 'j'-th artificial variable's reduced cost may be made
	//	non-zero.
	//
	{
		LP.GetColumn( j, a, row, len );
		assert( len == 1 );
		assert( *a == 1.0 || *a == -1.0 );

		Int_T i = *row;

		if( *a > 0.0 )
		{
			if( IsNotEqual( y_t[i], 1.0 ) )
			{
				Real_T Penalty = y_x[i] / ( 1.0 - y_t[i] );

				if( !( A2B[j] == A2B_UP && y_t[i] < 1.0 ) )
				{
					assert( 1.01 * Penalty > LP.GetPenalty() );
					MinPenalty = Max( MinPenalty, Penalty );
				}
			}
		}
		else
		{
			if( IsNotEqual( y_t[i], -1.0 ) )
			{
				Real_T Penalty = - y_x[i] / ( 1.0 + y_t[i] );

				if( !( A2B[j] == A2B_UP && y_t[i] > -1.0 ) )
				{
					assert( 1.01 * Penalty > LP.GetPenalty() );
					MinPenalty = Max( MinPenalty, Penalty );
				}
			}
		}
	}
	else if( DependsOnArtif )
	//
	//	See if 'j'-th structural variable may be moved to basis.
	//	If 'y_t^T a_j == 0' then no structural variable can be moved into
	//	the basis by changes to penalty 'M' (then it's reduced cost does
	//	not depend on 'M').
	//
	{
		Real_T y_t_a = 0.0, y_x_a = LP.GetC( j );

		for( LP.GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			y_t_a += *a * y_t[*row];
			y_x_a -= *a * y_x[*row];
		}

		y_x_a /= y_t_a;

		//------------------------------------------------------------------
		//	This is the proposed penalty value (we still have to check, if
		//	'y_t_a' has appropriate value). All formulas for lower bound on
		//	penalty eventually come to this formula.
		//
		switch( A2B[j] )
		{
		case A2B_UP:
			if( IsPositive( y_t_a ) )
			{
				assert( 1.01 * y_x_a > LP.GetPenalty() );
				MinPenalty = Max( MinPenalty, y_x_a );
			}
			break;

		case A2B_LO:
			if( IsNegative( y_t_a ) )
			{
				assert( 1.01 * y_x_a > LP.GetPenalty() );
				MinPenalty = Max( MinPenalty, y_x_a );
			}
			break;

		case A2B_IN:
			if( IsNonZero( y_t_a ) )
			{
				assert( 1.01 * y_x_a > LP.GetPenalty() );
				MinPenalty = Max( MinPenalty, y_x_a );
			}
			break;
		}
	}

	return MinPenalty;
}


Bool_T Solver::VerifyUnboundedness( Int_T j )
{
	WorkVector<Real_T> y_x( M ), y_t( M );

#ifndef NDEBUG
	ComputeDualVariables();
#endif

	ComputeSplitDualVariables();

	Real_T MinM = MinimumPenaltyToStayBounded( j, BasicArtificials() );
   	if( MinM <= LP.GetPenalty() )
		return True;

	PenaltyAdjustCnt++;
	LP.SetPenalty( 10.0 * MinM );
	ComputeDualVariables();

	ComputeReducedCosts();
	ComputeResult();

	return False;
}
