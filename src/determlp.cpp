/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	determlp.cpp
CREATED:			1995.07.31
LAST MODIFIED:		1996.04.16

DEPENDENCIES:		detrmlp.h, mps_lp.h, lp_codes.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif
#ifndef __DETERMLP_H__
#	include "determlp.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif

#ifndef __PRINT_H__
#	include "print.h"
#endif


#define LN_10 (2.302585092994)


/*------------------------------------------------------------------------------

	Int_T DeterministicLP::GetCostRow( void ) const

PURPOSE:
	Find the first free row in the constraint matrux (of the underlying MPS_LP
object).

PARAMETERS:
	None.

RETURN VALUE:
	Row number of -1 if none found.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T DeterministicLP::GetCostRow( void )
const
{
	if( CostRow < 0 )
	{
		Int_T i;
		for( i = 0; i < m; i++ )
			if( GetRowType( i ) & RT_FR )
				break;

		//
		//	Explicit cast overrides "const".
		//
		*( (Int_T *) &CostRow ) = Int_T( ( i < m ) ? i : -1 );
	}

	return CostRow;
}


const char *DeterministicLP::RevealObjectiveName( void )
const
{
	if( ObjName == NULL )
	{
		GetCostRow();
		assert( CostRow >= 0 );

		*((char **)&ObjName) = (char *)RowLabels.FindLabel( CostRow );
	}

	assert( ObjName != NULL );

	return ObjName;
}


/*------------------------------------------------------------------------------

	void DeterministicLP::ScaleObjective( Scenarios &scen )

PURPOSE:
	If necessary, scales the objective function by a integer power of ten.

PARAMETERS:
	None.

RETURN VALUE:
	The scaling factor used (or zero if no scaling performed).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T DeterministicLP::ScaleObjective( Scenarios &scen )
{
	GetCostRow();
	assert( CostRow >= 0 );

	//--------------------------------------------------------------------------
	//	Optional objective scaling.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;
	Real_T aMin = +INFINITY, aMax = 0.0;
	Bool_T Found = False;

	//
	//	Find the largest and the smallest cost vector non-zero.
	//
	for( Int_T j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++row, ++a )
		{
			if( *row != CostRow || IsZero( *a ) ) continue;

			aMax = Max( aMax, fabs( *a ) );
			aMin = Min( aMin, fabs( *a ) );
			Found = True;
		}

	if( !Found ) return 0;

	assert( aMin <= aMax );
	assert( IsPositive( aMin * aMax ) );
	//
	//	Find the scaling factor as an integer power of 10.
	//
	int Div = (int)( log(aMin * aMax) / (2.0 * LN_10) + 0.5 );

	//
	//	If necessary scale the objective.
	//
	if( Div )
	{
		aMax = pow( 10.0, (double) -Div );

		//----------------------------------------------------------------------
		//	Scale the LP first.
		//
		for( Int_T j = 0; j < n; j++ )
			for( GetColumn( j, a, row, len ); len; --len, ++row, ++a )
				if( *row == CostRow )
					*a *= aMax;

		//----------------------------------------------------------------------
		//	Then scale the random costs in scenarios.
		//
		scen.ScaleObjective( aMax );

		return (Int_T)Div;
	}
	else
		return 0;
}


/*------------------------------------------------------------------------------

	Bool_T DeterministicLP::CheckStructure( Int_T Stage2Row, Int_T Stage2Col )

PURPOSE:
	Checks if the LP structure is correct for a two stage stochastic LP in
deterministic formulation. The structure should be such, that there are no
non-zeros in the constraint matrix in rows 0 to "Stage2Row" and columns
"Stage2Col" till the end.

PARAMETERS:
	Int_T Stage2Row, Int_T Stage2Col
		Define the area of the matrix which should be empty (see above).

RETURN VALUE:
	Boolean value; "True" for correct structure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T DeterministicLP::CheckStructure( Int_T Stage2Row, Int_T Stage2Col )
{
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	for( Int_T j = Stage2Col; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++row )
			if( *row < Stage2Row && *row != CostRow )
				return False;
	return True;
}


/*------------------------------------------------------------------------------

	Int_T DeterministicLP::DivideIntoStages( MPS_LP &A, MPS_LP &T, MPS_LP &W,
		Int_T Stage2Col, const Scenario& Scen )

PURPOSE:
	Finds out which rows and columns actually belong in the second stage
subproblem. Second stage columns are those that start with "Stage2Col". Only
rows which have non-zeros in those columns or which have random data in the
technology matrix or the right hand side need be included in the second stage
subproblem.

PARAMETERS:
	MPS_LP &A, MPS_LP &T, MPS_LP &W
		The matrices into which the deterministic LP matrix should be divided:

		           | A 0 |
		DetermLP = |     |
		           | T W |

	Int_T Stage2Col
		The number of columns of A.

	const Scenario& Scen
		One scenario which shall be used to determine which parts of the
		deterministic LP formulation contain random data in stochastic form.

RETURN VALUE:
	The actual number of the second stage constraints.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T DeterministicLP::DivideIntoStages( MPS_LP &A, MPS_LP &T, MPS_LP &W, // )
	Int_T Stage2Col, const Scenario& Scen )
{
	//--------------------------------------------------------------------------
	//	Now see which rows are actually 1st stage constraints and which are
	//	second stage constraints. Arrays 'IncludeRows' and 'IncludeCols'
	//	store boolean data:
	//		True	if row/column belongs to the second stage and
	//		False	otherwise.
	//	The 'Stage1Row'/'Stage2Row' variables are ignored.
	//
	IncludeRows.Resize( m );
	IncludeRows.Fill( False, m );

	IncludeCols.Resize( n );
	IncludeCols.Fill( False, n );
	IncludeCols.Fill( True, n, Stage2Col );

	//--------------------------------------------------------------------------
	//	First mark rows which depend on the second stage variables. Mark the
	//	cost row as well. All in all this should be all 2nd stage rows!
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	Int_T j;
	for( j = Stage2Col; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++row )
			IncludeRows[ *row ] = True;
	IncludeRows[ CostRow ] = True;

#ifndef NDEBUG
	{
		for( Int_T jj = 0, l = Scen.GetLength(); jj < l; jj++ )
			for( Int_T i = 0, bl = Scen[jj].Len(); i < bl; i++ )
			{
				const Delta &d = Scen[jj][i];
				if( !( d.type != Delta::COST || d.col >= Stage2Col ) ||
					!( d.type != Delta::MATRIX || d.col < Stage2Col ) )
				{
					Print( "d.type = %d, d.col = %d, Stage2Col = %d\n",
						d.type, d.col, (int)Stage2Col );
					abort();
				}
			}
	}
#endif

	//--------------------------------------------------------------------------
	//	Then scan one scenario to find out which rows depend on random data.
	//
	Int_T l = Scen.GetLength();
	for( j = 0; j < l; j++ )
		for( Int_T i = 0, bl = Scen[j].Len(); i < bl; i++ )
		{
			const Delta &d = Scen[j][i];

			switch( d.type )
			{
			case Delta::RHS:
			case Delta::MATRIX:
				IncludeRows[d.row] = True;
				break;

			case Delta::COST:
				assert( IncludeCols[d.col] );
				break;

			default:
#ifndef NDEBUG
				abort();
#endif
				break;
			}
		}

	//--------------------------------------------------------------------------
	//	Count how many rows are in the second stage now.
	//
	ActualStage2Rows = 0;

	Int_T i;
	for( i = 0; i < m; i++ )
		if( IncludeRows[i] )
			ActualStage2Rows++;

	//--------------------------------------------------------------------------
	//	Finally create all three matrices as submatrices of 'DetermLP'.
	//
	W.CreateAsSubmatrix( *this, IncludeRows, IncludeCols );

	for( j = 0; j < n; j++ )
		IncludeCols[j] = ( IncludeCols[j] ) ? False : True;

	T.CreateAsSubmatrix( *this, IncludeRows, IncludeCols );

	for( i = 0; i < m; i++ )
		IncludeRows[i] = ( IncludeRows[i] ) ? False : True;
	IncludeRows[ CostRow ] = True;

	A.CreateAsSubmatrix( *this, IncludeRows, IncludeCols );

	return ActualStage2Rows;
}


/*------------------------------------------------------------------------------

	void DeterministicLP::ZeroRandomData( const Scenario &Scen )

PURPOSE:
	Replace all random data by zero's.

PARAMETERS:
	const Scenario &Scen
		A scenario (any one).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void DeterministicLP::ZeroRandomData( const Scenario &Scen )
{
	for( Int_T l = Scen.GetLength(), j = 0; j < l; j++ )
		for( Int_T i = 0, bl = Scen[j].Len(); i < bl; i++ )
			switch( Scen[j][i].type )
			{
			case Delta::RHS:
				SetRHS_Elem( Scen[j][i].row, 0.0 );
				break;

			case Delta::MATRIX:
				SetMatrixElement( Scen[j][i].row, Scen[j][i].col, 0.0 );
				break;

			case Delta::COST:
				SetC( Scen[j][i].col, 0.0 );
				break;

			default:
#ifndef NDEBUG
				abort();
#endif
				break;
			}
}


/*------------------------------------------------------------------------------

	void DeterministicLP::RenumberIndiceInScenarios( Scenarios &Scen )

PURPOSE:
	Random data in scenarios is represented by row and column indice and some
values from the distribution. The indice are relative to the deterministic LP
formulation. This procedure recalculates the indice so that they will be correct
relative to a two stage formulation.

PARAMETERS:
	Scenarios &Scen
		All scenarios.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void DeterministicLP::RenumberIndiceInScenarios( Scenarios &Scen, Int_T // )
	Stage2Col )
{
	//--------------------------------------------------------------------------
	//	Free rows (including the objective row) are to be discarded in counting.
	//
	Int_T i;
	for( i = 0; i < m; i++ )
		if( GetRowType( i ) & RT_FR )
			IncludeRows[i] = True;

	//--------------------------------------------------------------------------
	//	Compute new row/column numbers (assuming that rows/columns for which
	//	'IncludeXxxs[i] == True' belong to the first stage constraints).
	//
	WorkVector<Int_T> NewRowNumber( m ),
		NewColNumber( n );
	Int_T i1, i2;

	//--------------------------------------------------------------------------
	//	First the constraints: 'i2' counts the second stage constraints only
	//	(nothing random may be present in the first stage rows).
	//
	NewRowNumber.Fill( -1, m );
	for( i = 0, i2 = 0; i < m; i++ )
		if( !IncludeRows[i] )
			NewRowNumber[i] = i2++;

	//--------------------------------------------------------------------------
	//	Then the variables (which appear in the technology matrix, both in its
	//	deterministic part and in the random one). 'i1' and 'i2' behave as
	//	before.
	//
	NewColNumber.Fill( -1, n );
	for( i = 0, i1 = 0, i2 = 0; i < n; i++ )
		if( IncludeCols[i] )
			NewColNumber[i] = i1++;
		else
			NewColNumber[i] = i2++;

	//--------------------------------------------------------------------------
	//	Renumber the rows and columns in scenarios.
	//
	Scen.RenumberIndiceInScenarios( NewRowNumber, m, NewColNumber, n,
		ActualStage2Rows, Stage2Col );
}
