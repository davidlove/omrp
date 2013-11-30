/*------------------------------------------------------------------------------
MODULE TYPE:		Pseudo - simplex type optimizer core code.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of the revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio
ADD. SUPPORT:		prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	rm_infeas.cpp
CREATED:			1994.06.11
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		solver.h, std_math.h, smplx_lp.h, std_tmpl.h, work_vec.h
					<assert.h>, <math.h>

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

#include <assert.h>
#include <math.h>

#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T Solver::MinimizeLambda( Array<Real_T> &lambda, VerbLevel Verbosity )

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

Bool_T Solver::MinimizeLambda( Array<Real_T> &lambda, VerbLevel Verbosity )
{
	InfeasMinCnt++;

	if( Verbosity >= V_LOW )
		Print( "\nINFEASIBILITY MINIMIZATION ROUTINE INVOKED:\n" );

	//--------------------------------------------------------------------------
	//	'NonZero' work vector will store the current non-zero pattern of the
	//	'lambda' vector.
	//
	WorkVector<Bool_T> NonZero( M );
	NonZero.Fill( False, M );

	//--------------------------------------------------------------------------
	//	'InfCnt' will store the number of non-zero entries of the lambda vector.
	//	'InfNorm' will store a square of the Euclidean norm of the lambda
	//	vector. Thus updates of 'InfNorm' will not be very expensive. The root
	//	will only have to be computed for outputting the current value of
	//	'InfNorm'.
	//
	Int_T InfCnt	= 0;
	Real_T InfNorm	= 0.0;

	//--------------------------------------------------------------------------
	//	Scan lambda vector to determine the non-zero pattern and remove small
	//	non-zero entries.
	//
	Int_T i;
	for( i = 0; i < M; i++ )
		if( IsNonZero( lambda[i] ) )
		{
			NonZero[i]	= True;
			InfNorm		+= lambda[i] * lambda[i];
			InfCnt++;
		}
		else
			lambda[i] = 0.0;

	if( Verbosity >= V_LOW )
		Print(
			"\t%-10s  %-10s  %-10s\n"
			"\t%-10s  %-10s  %-10s\n"
			"\t%-10d  %-10d  %-10.3E\n",

			"Pass",			"Inf. num.",	"Inf. norm",
			"----------",	"----------",	"----------",
			0,				int(InfCnt),	double( sqrt( InfNorm ) )
		);

	//--------------------------------------------------------------------------
	//	No infeasibility -- we may as well return to the caller now.
	//
	if( InfCnt == 0 )
		return False;

	//--------------------------------------------------------------------------
	//	Now is the time to sort the problem's columns from the shortest to the
	//	longest one.
	//	ATTENTION: We do not take artificial variables into account!
	//
	Int_T OrigN = Int_T( LP.GetStructN() + LP.GetSlackN() );
	WorkVector<Int_T> ColList( OrigN ),
		Start( M + 1 );

	Start.Fill( 0, M + 1 );
	ColList.Fill( -1, OrigN );

	Ptr<Int_T> row;
	Ptr<Real_T> a;
	Int_T len;

	//
	//	Count the numbers of columns of identical length.
	//
	for( i = 0; i < OrigN; i++ )
	{
		LP.GetColumn( i, a, row, len );
		Start[ len ] ++;
	}

	//
	//	Make the 'Start' table entries point at the ends of the (future)
	//	lists of columns of the same length.
	//
	for( i = 1; i <= M; i++ )
		Start[i] += Start[i-1];

	//
	//	Fill the column list with actual numbers.
	//
	for( i = 0; i < OrigN; i++ )
	{
		LP.GetColumn( i, a, row, len );
		ColList[ --Start[ len ] ] = i;
	}

	//--------------------------------------------------------------------------
	//	Here begins the part of the procedure, where the reductions may be
	//	performed.
	//
	Int_T MaxPasses			= 2,
		MaxFillIn			= 0;
	Bool_T Improved			= True,
		OverallImproved		= False;

	//--------------------------------------------------------------------------
	//	Loop on pass counter. Unlimited number of passes if 'MaxPasses == 0'.
	//
	for( Int_T pass = 0;
		( MaxPasses == 0 || pass < MaxPasses ) && Improved && InfCnt > 0;
		pass++ )
	{
		Improved = False;

		//----------------------------------------------------------------------
		//	Loop on columns of lengths in the range <1,InfCnt+MaxFill>.
		//
		Int_T j;
		for( len = 1; len < Min( Int_T(InfCnt + MaxFillIn), M ); len++ )
			for( j = Start[len]; j < Start[len+1]; j++ )
			{
				if( MinimizeInfOneVariable( ColList[j], lambda, NonZero,
					MaxFillIn, InfCnt, InfNorm ) )
					Improved = True;
			}

		//----------------------------------------------------------------------
		//	Now report the current pass results: number of non-zeros and
		//	norm of infeasibility vector.
		//
		if( Verbosity >= V_LOW )
			Print(
				"\t%-10d  %-10d  %-10.3E\n",
				int( pass + 1 ), int(InfCnt), double( sqrt( InfNorm ) )
			);

		if( Improved )
			OverallImproved = True;
	}

	return OverallImproved;
}

/*------------------------------------------------------------------------------

	Bool_T Solver::MinimizeInfOneVariable( Int_T col, Array<Real_T> &lambda,
		Array<Bool_T> &NonZero, Int_T MaxFillIn, Int_T &InfCnt,
		Real_T &InfNorm )

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

Bool_T Solver::MinimizeInfOneVariable( Int_T col, Array<Real_T> &lambda,
	Array<Bool_T> &NonZero, Int_T MaxFillIn, Int_T &InfCnt, Real_T &InfNorm )
{
	//--------------------------------------------------------------------------
	//	Check the resulting fill-in. Try next column if fill in
	//	too big.
	//
	Int_T FillIn	= 0,
		NonZeroHits	= 0;
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	for( LP.GetColumn( col, a, row, len ); len; --len, ++row )
		if( !NonZero[ *row ] )
		{
			if(  ++FillIn > MaxFillIn )
				return False;
		}
		else
			NonZeroHits++;

	if( NonZeroHits == 0 )
		return False;

	//--------------------------------------------------------------------------
	//	Compute optimal value of x_j.
	//
	Real_T a_l	= 0.0,
		a_a		= 0.0;

	for( LP.GetColumn( col, a, row, len ); len; --len, ++a, ++row )
	{
		if( NonZero[ *row ] )
			a_l	+= *a * lambda[ *row ];
		a_a	+= *a * *a;
	}

	Real_T xOpt = a_l / a_a;

	//--------------------------------------------------------------------------
	//	If no (or almost no) change, try next column from the list.
	//
	if( fabs( xOpt ) < 1.0e-3 )
		return False;

	//--------------------------------------------------------------------------
	//	Find closest feasible value of 'x[col]' and set it to this
	//	value.
	//
	if( VarType[col] & VT_HAS_LO_BND && x[col] + xOpt < 0.0 )
	{
		xOpt	= -x[col];
		x[col]	= 0.0;
		if( A2B[col] < 0 )
			A2B[col] = A2B_LO;
	}
	else if( VarType[col] & VT_HAS_UP_BND && x[col] + xOpt > u[col] )
	{
		xOpt	= u[col] - x[col];
		x[col]	= u[col];
		if( A2B[col] < 0 )
			A2B[col] = A2B_UP;
	}
	else
	{
		x[col]	+= xOpt;
		if( A2B[col] < 0 )
			A2B[col] = A2B_IN;
	}

	if( fabs( xOpt ) < 1.0e-3 )
		return False;

#ifndef NDEBUG
	Real_T OldInfNorm = InfNorm;
#endif

	//--------------------------------------------------------------------------
	//	Update 'lambda', its non-zero pattern and norm.
	//
	for( LP.GetColumn( col, a, row, len ); len; --len, ++a, ++row )
	{
		if( NonZero[ *row ] )
		{
			if( ( InfNorm -= lambda[ *row ] * lambda[ *row ] ) < 0.0 )
				InfNorm = 0.0;

			if( IsZero( lambda[ *row ] -= xOpt * *a ) )
			{
				lambda[ *row ]	= 0.0;
				NonZero[ *row ]	= False;
				InfCnt--;
			}
			else if( ( InfNorm += lambda[ *row ] * lambda[ *row ] ) < 0.0 )
				InfNorm = 0.0;
		}
		else
		{
			lambda[ *row ] = - xOpt * *a;
			InfNorm += lambda[ *row ] * lambda[ *row ];
			NonZero[ *row ]	= True;
			InfCnt++;
		}
	}

#ifndef NDEBUG
	OldInfNorm -= InfNorm;
#endif

	return True;
}
