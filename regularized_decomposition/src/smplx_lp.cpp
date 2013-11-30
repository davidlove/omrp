/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	smplx_lp.cpp
CREATED:			1994.03.11
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		smplx_lp.h, solv_lp.h, mps_lp.h, smartptr.h, std_math.h,
					print.h
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

#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif


/*------------------------------------------------------------------------------

	SimplexLP::SimplexLP( void )

PURPOSE:
	Object constructor. Creates an empty object.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SimplexLP::SimplexLP( void )
	: Standard( False ), M( 1.0 ), SlackLen( 0 )
{}


/*------------------------------------------------------------------------------

	Real_T SimplexLP::GetL( Int_T j ) const
	Real_T SimplexLP::GetU( Int_T j ) const
	Real_T SimplexLP::GetC( Int_T j ) const

PURPOSE:
	Return lower bound, upper bound and cost respectively for a variable "j".

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T SimplexLP::GetL( Int_T j )
const
{
	assert( j >= 0 && j < n + SlackLen + m );

	return ( j < n ) ?
		MPS_LP::GetL( j ) :		// Structural variable.
		0.0;					// Slack or artificial variable.
}


Real_T SimplexLP::GetU( Int_T j )
const
{
	assert( j >= 0 && j < n + SlackLen + m );

	if( j < n )						// Structural variable.
		return MPS_LP::GetU( j );
	j -= n;
	if( j < SlackLen )				// Slack variable.
		return SlackU[j];

	// Super-artificial variable: fixed or does not have an upper bound.
	//
	return ( LambdaVT[j-SlackLen] & VT_FX ) ? 0.0 : +INFINITY;
}


Real_T SimplexLP::GetC( Int_T j )
const
{
	assert( j >= 0 && j < n + SlackLen + m );

	if( j < n )						// Structural variable.
		return SolvableLP::GetC( j );
	j -= n;
	if( j < SlackLen )				// Slack variable.
		return 0.0;

	j -= SlackLen;
									// Artificial variable (of any kind).
	return /* (LambdaVT[j] & VT_FX) ? 0.0 : */ M;
}


/*------------------------------------------------------------------------------

	Short_T SimplexLP::GetVarType( Int_T j ) const

PURPOSE:
	Returns variable type for "j"-th variable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Short_T SimplexLP::GetVarType( Int_T j )
const
{
	assert( j >= 0 && j < n + SlackLen + m );

	if( j < n )						// Structural variable.
		return MPS_LP::GetVarType( j );
	j -= n;
	if( j < SlackLen )				// Slack variable.
		return SlackVT[j];
	return LambdaVT[j-SlackLen];	// Super-artificial variable.
}


/*------------------------------------------------------------------------------

	void SimplexLP::CreateLambda( const Array<Real_T> &v )

PURPOSE:
	Create a vector that will be capable of holding the infeasibility of the
problem. It (logically) corresponds to a diagonal matrix, from which (most
likely) some columns were taken out. The infeasibility that needs to be handled
is passed in vector 'v'.

PARAMETERS:
	const Array<Real_T> &v
		A dense residual vector.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::CreateLambda( const Array<Real_T> &v )
{
	Lambda.Resize( m );		Lambda.Fill( 1.0, m );
	LambdaVT.Resize( m );	LambdaVT.Fill( VT_FIXED | VT_ARTIF, m );

	for( Int_T i = 0; i < m; i++ )
		if( IsNonZero( v[i] ) )
		{
			if( v[i] < 0.0 ) Lambda[i] = -1.0;

			LambdaVT[i]	= VT_LO | VT_ARTIF;
		}
}


/*------------------------------------------------------------------------------

	void SimplexLP::ToStandard( void )

PURPOSE:
	Converts the problem to a standard form, in which all rows are equalities
(this is done by adding slacks) and variable lower bounds are at zero.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::ToStandard( VerbLevel Verbosity )
{
	//--------------------------------------------------------------------------
	//	Resize all the vectors to match the constraint matrix dimensions.
	//
	SlackU.Resize( m );				SlackU.Fill( +INFINITY, m );
	SlackVT.Resize( m );			SlackVT.Fill( VT_NORM, m );

	Slack.Resize( m );
	SlackRow.Resize( m );
	SlackPosByRows.Resize( m );		SlackPosByRows.Fill( -1, m );
	Lambda.Resize( m );				Lambda.Fill( 1.0, m );
	LambdaVT.Resize( m );			LambdaVT.Fill( VT_FIXED | VT_ARTIF, m );
	LambdaRow.Resize( m );

	Int_T i;
	for( i = 0; i < m; i++ ) LambdaRow[i] = i;

	ComputePenalty();
	ShiftLowerBoundsToZero();

	//--------------------------------------------------------------------------
	//	Add slacks to non-equality rows.
	//
	for( i = 0; i < m; i++ )
	{
		//----------------------------------------------------------------------
		//	Row with a range => upper bound on a slack variable has to be
		//	imposed.
		//
		if( RowType[i] & RT_RNG )
		{
			switch( RowType[i] & RT_TYPE )
			{
			case RT_LE:	Slack[ SlackLen ]	= 1.0e0;					break;
			case RT_GE:	Slack[ SlackLen ]	= -1.0e0;					break;
			case RT_EQ: Slack[ SlackLen ]	= r[i] > 0.0 ? -1.0 : 1.0;	break;
#ifndef NDEBUG
			default:	abort();
#endif
			}
		 	SlackU[ SlackLen ]		= fabs( r[i] );
			SlackVT[ SlackLen ]		= VT_BOUNDED;
			SlackRow[ SlackLen ]	= i;
			SlackPosByRows[ i ]		= SlackLen++;
		}
		else
		//----------------------------------------------------------------------
		//	Row without a range. Slack declared PL or MI depending on row type.
		//
		{
			Slack[ SlackLen ] = 1.0e0;
			switch( RowType[i] & RT_TYPE )
			{
			case RT_LE:
				SlackVT[ SlackLen ]	= VT_NORM;
				break;

			case RT_GE:
				SlackVT[ SlackLen ]	= VT_MI;
				SlackU[ SlackLen ]	= 0.0;
				break;

			case RT_EQ:	break;

#ifndef NDEBUG
			default:	abort();
#endif
			}

			if( !(RowType[i] & RT_EQ) )
			{
				SlackRow[ SlackLen ]	= i;
				SlackPosByRows[ i ]		= SlackLen++;
			}
		}
	}

	//--------------------------------------------------------------------------
	//	Report on current problem statistics.
	//
	if( Verbosity >= V_HIGH )
		Print(

			/* Stream output FORMAT (two lines at the top, then - a table). */

			"\nConversion to standard form completed.\n"
			"Constraint matrix data:\n"
			"\t%-18s%10d\n"
			"\t%-18s%10d\n"
			"\t%-18s%10ld\n"
			"\t%-18s%10d\n"
			"\t%-18s%10.4f %%\n"
			"\t%-18s%10.4f %%\n",

			/* DATA. */

			"No. of rows:",				(int) GetM(),
			"No. of vars:",				(int) GetN(),
			"No. of non-zeros:",		(long int) GetNZ(),
			"No. of slacks:",			(int) SlackLen,
			"Density (orig.):",			(double)( 100 * nz )
										/ (double) ( m * n ),
			"Density (total):",			double( 100 * GetNZ() )
										/ (double) ( GetM() * GetN() )
			);
	Standard = True;
}


void SimplexLP::UndoStandard( void )
{
	//--------------------------------------------------------------------------
	//	Restore the original positions of the lower bounds.
	//
	Int_T i, len;
	Ptr<Real_T> a;
	Ptr<Int_T> row;

	for( i = 0; i < n; i++ )
	{
		Real_T ll = GetL( i );
		Bool_T up = ( GetVarType( i ) & VT_HAS_LO_BND ) ? True : False;

		if( up && IsNonZero( ll ) )
		{
			for( MPS_LP::GetColumn( i, a, row, len ); len; --len, ++a, ++row )
				b[ *row ] += *a * ll;

			if( up ) u[i] += ll;
		}
	}
}


/*------------------------------------------------------------------------------

	void SimplexLP::GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row,
		Int_T &len )
		const

	void SimplexLP::GetRow( Int_T row, Ptr<Real_T> &a, Ptr<Int_T> &col,
		Int_T &len, Int_T &SlackCol, Real_T &theSlack,
		Int_T &LambdaCol, Real_T &theLambda )
		const

PURPOSE:
	This function gives a uniform way of accessing all columns of the standard
form (and perhaps presolved) LP problem. It returns pointers INTO the tables
which store the LP problem. That's why the "Col" and "Row" are pointers to
"const" data.

PARAMETERS:
	Int_T j/i
		Column number (in range 0 : "n + SlackLen + m" ) / row number
		(in range 0 : m) of column to be accessed.

	Ptr<Real_T> &a
		Array of column "j"/row "i" non-zeros. It's value on entry is
		irrelevant.

	Ptr<Int_T> &Row/&Col
		Array which holds column "j"/row "i" non-zero pattern. It's value
		on entry is irrelevant.

	Int_T &len
		Number of non zeros in column/row.

	Int_T &SlackCol, Real_T &theSlack
		Number of slack column for thr "i"-th row (or -1 if none) and
		the value of the corresp. non-zero.

	Int_T &LambdaCol, Real_T &theLambda
		Number of the lambda column corresp. to the "i"-th row and it's
		non-zero value.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row, // )
	Int_T &len )
	const
{
	assert( j >= 0 && j < n + SlackLen + m );

	if( j < n )										// Original column.
		MPS_LP::GetColumn( j, a, row, len );
	else if( ( j -= n ) < SlackLen )				// Slack column.
	{
		a.ExtractFragment(		Slack,		j, j+1 );
		row.ExtractFragment(	SlackRow,	j, j+1 );
		len = 1;
	}
	else
	{
		j	-= SlackLen;							// Super artificial column.

		a.ExtractFragment(		Lambda,		j, j+1 );
		row.ExtractFragment(	LambdaRow,	j, j+1 );
		len = 1;
	}
}


void SimplexLP::GetRow( Int_T row, Ptr<Real_T> &a, Ptr<Int_T> &col, Int_T &len,
	Int_T &SlackCol, Real_T &theSlack, Int_T &LambdaCol, Real_T &theLambda )
	const
{
	SolvableLP::GetRow( row, a, col, len );

	Int_T pos = SlackPosByRows[ row ];

	if( pos >= 0 )
	{
		theSlack = Slack[ pos ];
		SlackCol = Int_T( n + pos );
	}
	else
	{
		theSlack = 0.0;
		SlackCol = -1;
	}

	LambdaCol = Int_T( n + SlackLen + row );
	theLambda = Lambda[ row ];
}

/*------------------------------------------------------------------------------

	void SimplexLP::ComputePenalty( void )

PURPOSE:
	Computes an initial low value of the penalty "M" as twice the maximum cost
absolute value in the original objective.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::ComputePenalty( void )
{
	//--------------------------------------------------------------------------
	//	Compute the penalty 'M'.
	//
	Int_T i, nn = GetStructN();
	for( M = 1.0, i = 0; i < nn; i++ )
		M = Max( M, fabs( SolvableLP::GetC( i ) ) );
	M *= 2.0;
}


void SimplexLP::GetPenaltyEstimates( Real_T &MinM, Real_T &MaxM )
{
	//--------------------------------------------------------------------------
	//	Estimate the minimum penalty.
	//
	MaxM = -INFINITY;	// "Guarantees" that the problem will not be unbounded.
	MinM = +INFINITY;	// Guarantees that the problem will be unbounded.

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T j, len, nn = GetStructN();

	for( j = 0; j < nn; j++ )
	{
		Real_T x = 0.0;

		for( GetColumn( j, a, row, len ); len; --len, ++a, ++row )
			if( ! ( LambdaVT[ *row ] & VT_FX ) )
				x += Lambda[ *row ] * *a;

		if( IsNonZero( x ) )
		{
			x = GetC( j ) / x;

			if( MaxM < x ) MaxM = x;
			if( MinM > x ) MinM = x;
		}
	}
}


/*------------------------------------------------------------------------------

	void SimplexLP::ShiftLowerBoundsToZero( void )

PURPOSE:
	Shifts all finite lower bounds on variables to zero. Adjusts the right hand
side and the objective fixed adjustment.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::ShiftLowerBoundsToZero( void )
{
	//--------------------------------------------------------------------------
	//	Shift lower bounds to zero (where possible, i.e. where the lower bounds
	//	are finite. Adjust the finite upper bound on variable as well as the
	//	value of the objective function fixed adjustment accordingly.
	//
	Int_T i, nn = GetN(), len;
	Ptr<Real_T> a;
	Ptr<Int_T> row;

	for( i = 0; i < nn; i++ )
	{
		Real_T ll = GetL( i );
		Bool_T up = ( GetVarType( i ) & VT_HAS_LO_BND ) ? True : False;

		if( up && IsNonZero( ll ) )
		{
			for( MPS_LP::GetColumn( i, a, row, len ); len; --len, ++a, ++row )
				b[ *row ] -= *a * ll;

			if( up ) u[i] -= ll;

			f += GetC( i ) * ll;
		}
	}
}


/*------------------------------------------------------------------------------

	void SimplexLP::ProcessSolution( Solution &sol )

PURPOSE:
	This function takes the solution to a linear problem and (if the problem
was scaled before solving) retrieves the solution to the iriginal problem.

PARAMETERS:
	Solution &sol
		the solution structure, which holds:
		-	solution vector's dimensions,
		-	primal and dual variables,
		-	the objective function value.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::ProcessSolution( Solution &sol )
{
	Array<Real_T> &x	= sol.x;

	int contents = sol.GetContents();

	assert( !( contents & Solution::Primal ) || sol.GetN() == n );
	assert( !( contents & Solution::Dual ) || sol.GetM() == m );

	if( Standard && ( contents & Solution::Primal ) )
		for( Int_T j = 0; j < n; j++ )
			if( VarType[j] & VT_HAS_LO_BND && IsNonZero( l[j] ) )
				x[j] += l[j];

	SolvableLP::ProcessSolution( sol );
}
