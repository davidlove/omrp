/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, simplex type basis construction.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr. Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	crash.cpp
CREATED:			1993.05.17
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		smartptr.h, stdtype.h, std_tmpl.h, error.h, mps_lp.h,
					solv_lp.h, smplx_lp.h, solvcode.h, crash.h, print.h
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

------------------------------------------------------------------------------*/

#include <math.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __COMPILE_H__
#	include "compile.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __SOLVCODE_H__
#	include "solvcode.h"
#endif


//==============================================================================
//
//	Local type and symbolic constant definitions.
//

enum PREF_SET {
	PS_FREE	= 0,
	PS_NORM	= 1,
	PS_REST	= 2,
	PS_NONE	= 3
};

//------------------------------------------------------------------------------
//							Structure used in crash
//
struct PS			// Preference sets
{
	Int_T next,		// Next variable number.
		previous,	// Previous variable number.
		len,		// Constraint matrix column working length.
		origLen;	// Total column length.
	PREF_SET type;	// Variable type.
	Real_T price;	// Variable price (not cost!).
};

//
//	End of local type and symbolic constant definitions.
//
//==============================================================================

//------------------------------------------------------------------------------
//	Implicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )
	template class SmartPointerBase<PS> ;
	template class Array<PS>;
	template class Ptr<PS>;

	template class WorkVector<PS>;
	int WorkVector<PS>::VecType = WorkVectorPool::NO_TYPE;

	template PS *MALLOC( PS *& Table, size_t len );
	template PS *REALLOC( PS *& Table, size_t len );
	template PS *FREE( PS *& Table );
#endif
//
//------------------------------------------------------------------------------

//==============================================================================
//
//	Local type and symbolic constant definitions.
//
//==============================================================================

enum {
	EMPTY		= 0x0000,
	NON_ZERO	= 0x0001,
	PIVOT		= 0x0002
};

//==============================================================================
//
//	End of local type and symbolic cooonstant definitions.
//
//==============================================================================

//==============================================================================
//
//	Static function prototypes.
//
//==============================================================================

#ifndef NDEBUG
	static void CheckSets( Array<PS> &ps, Int_T n, Array<Int_T> &PsStart,
		Array<Int_T> &PsEnd, Array<Int_T> &PsLen );
#endif

static void RmColFromList( Array<PS> &ps, Int_T Col, Int_T &lStart,
	Int_T &lEnd );
static void RmRowFromList( const SimplexLP &LP, Array<PS> &ps, Int_T PivotRow,
	Int_T &ListStart, Int_T &ListEnd );
static Int_T FindPivotCol( Array<PS> &ps, Int_T m, Int_T ColNum );
static Int_T FindPivotRow( const SimplexLP *LP, Array<Short_T> &RowMark,
	Array<Real_T> &RowMax, Int_T PivotCol, Real_T PivTol );

//==============================================================================
//
//	End of static function prototypes.
//
//==============================================================================


/*------------------------------------------------------------------------------

	void SimplexLP::InitialBasis( Real_T PivotTol, Array<Int_T> A2B,
		VerbLevel Verbosity )

PURPOSE:
	Constructs an initial basis as a permutation of an upper triangular matrix.
Fills in the empty places with artificial variables' columns (added slacks for
equality rows).
	The columns of the original matrix are divided into four 'preference sets'
(a term suggested by Bixby) according to the level of freedom of associated
variables and probability of their presence in the optimal basis.
	The variables are then selected using the sparsity pattern of their columns
in order to form an upper triangular matrix. The row cliques are needed to
update lengths of columns after adding a column into basis and therefore
excluding some columns from further consideration.

PARAMETERS:
	None.

RETURN VALUE:
	Returns an initial basis layout expressed as a table of entries marking
basis / out-of-basis status of all variables (structural and slack). The "Int_T"
array is allocated inside this function and returned to the calling function.
Responsibility to deallocate the the array when it is no longer needed remains
with the user.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SimplexLP::InitialBasis( Real_T PivotTol, Array<Int_T> &A2B, // )
	VerbLevel Verbosity )
{
	//==========================================================================
	//
	//	ALLOCATE AND CONSTRUCT PREFERENCE SETS
	//
	//==========================================================================

	WorkVector<PS> ps( GetStructN() );

	Array<Int_T> PsStart( PS_NONE, -1 ),			// Mark empty lists.
		PsEnd( PS_NONE, -1 ),						// Mark end of lists.
		PsLen( PS_NONE, 0 );						// Sets' zero lengths.

	//--------------------------------------------------------------------------
	//	Allocate memory and initialize a table of maximum row entries.
	//
	WorkVector<Real_T> RowMax( m );
	RowMax.Fill( 0.0, m );

	//--------------------------------------------------------------------------
	//	Divide variables into three sets:
	//		1.	Free variables,
	//		2.	Normal variables,
	//		3.	Bounded and fixed.
	//	Seting up element lists (PsStart, PsEnd, PsLen) linked with fields
	//	"next" and "previous" of PS structure.
	//	Scan structural columns only.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len, j;

	for( j = 0; j < GetStructN(); j++ )
	{
		//----------------------------------------------------------------------
		//	Fill the table of maximum row entries.
		//
		MPS_LP::GetColumn( j, a, row, len );

		ps[j].origLen = ps[j].len = len;		// Store column length now.

		for( ; len; --len, ++a, ++row )
			RowMax[ *row ] = Max( fabs( *a ), RowMax[ *row ] );

		PREF_SET type = PS_NONE;

		//----------------------------------------------------------------------
		//	Price is a measure of freedom of variable (1 for PS_FREE, 0 for
		//	PS_NORM, '-c[j]+u[j]' for BOUNDED and '-INFINITY/2' for FIXED.
		//		
		//	Higher priced columns are preferred when considering candidacy for
		//	basis (in the same preference set).
		//
		switch( GetVarType( j ) & VT_TYPE )
		{
		case VT_FIXED:
			type			= PS_REST;
			ps[j].price		= -INFINITY/2;
			break;

		case VT_BOUNDED:
			type			= PS_REST;
			ps[j].price		= -GetC( j )/1.0e3 + GetU( j );
			break;

		case VT_MI:
		case VT_NORM:
			type			= PS_NORM;
			ps[j].price		= 0.0e0;
			break;

		case VT_FREE:
			type			= PS_FREE;
			ps[j].price		= 1.0e0;
			break;

		default:
			FatalError( "crash: InitialBasis: "
				"Invalid variable type detected: %d.",
				int( GetVarType( j ) & VT_TYPE ) );
		}

		ps[j].type		= type;				// Variable category.

		if( ps[j].origLen == 0 )			// Do not put empty columns on the
											// lists.
			ps[j].previous = ps[j].next = -1;
		else
		{
			PsLen[ type ]++;				// Increment suitable list's length.
			if( PsStart[ type ] == -1 )		// Mark list start if list empty.
				PsStart[ type ] = j;
			ps[j].previous = PsEnd[ type ];	// Link to previous element.
			if( PsEnd[ type ] != -1 )		// Link the previous el. with the
				ps[ PsEnd[ type ] ].next = j;// new one.
			PsEnd[ type ] = j;				// Mark new list end.
		}
	}

	//--------------------------------------------------------------------------
	// Make sure the doubly linked lists are properly terminated.
	//
	Int_T i;
	for( i = PS_FREE; i <= PS_REST; i++ )
	{
		if( PsStart[i] >= 0 )	ps[ PsStart[i] ].previous = -1;
		if( PsEnd[i] >= 0 )		ps[ PsEnd[i] ].next = -1;
	}

	//--------------------------------------------------------------------------
	//	Check the preference sets (lists).
	//
#ifndef NDEBUG
	CheckSets( ps, GetStructN(), PsStart, PsEnd, PsLen );
#endif

	//--------------------------------------------------------------------------
	//	Make one continuous list of the three sets in their order of preference:
	//	first PS_FREE, then PS_NORM, PS_REST.
	//
	Int_T ListStart = -1,
		ListEnd = -1;

	for( i = PS_FREE; i <= PS_REST; i++ )		// Find first non-empty set.
		if( PsStart[i] >= 0 )
		{
			ListStart	= PsStart[i];
			ListEnd		= PsEnd[i];
			break;
		}

	assert( ListStart != -1 );

	for( ++i; i <= PS_REST; i++ )				// Link the remaining sets.
		if( PsStart[i] != -1 )
		{
			ps[ ListEnd ].next			= PsStart[i];
			ps[ PsStart[i] ].previous	= ListEnd;
			ListEnd						= PsEnd[i];
		}

#ifndef NDEBUG
	{
		Int_T ColCnt, jj;

		for( ColCnt = 0, jj = ListStart; jj != -1 && ColCnt < GetStructN();
			jj = ps[jj].next, ColCnt++ )
			;

		assert( jj == -1 && ColCnt >= GetStructN() - 1 );
	}
#endif

	//==========================================================================
	//	
	//	Construct as much of the initial basis as possible from the original
	//	problem's columns.
	//	
	//==========================================================================

	WorkVector<Short_T> RowMark( m );	// Array of row marks for rows that have
	RowMark.Fill( EMPTY, m );			// already been pivotal, or have
										// non-pivotal non-zeros.

	//--------------------------------------------------------------------------
	//	Fill basis with columns.
	//	Construct a permutation of an upper triangular matrix.
	//
	//	First insert all slacks into basis, then scan each of the remaining sets
	//	of variables in order of preference.
	//
	Int_T PivotRow, PivotCol, ColNum, PivotCnt;

	//--------------------------------------------------------------------------
	//	Loop in which slacks are inserted into basis.
	//
	for( PivotCnt = j = 0; j < SlackLen; j++ )
	{
		//----------------------------------------------------------------------
		//	Find the pivot (that is the only non-zero in slack column).
		//
		PivotCol = Int_T( GetStructN() + j );
		PivotRow = SlackRow[j];

		//----------------------------------------------------------------------
		//	Store basis column number.
		//
		assert( A2B[ PivotCol ] == A2B_UNDEF );
		A2B[ PivotCol ] = A2B_BASIC;
		PivotCnt++;

		//----------------------------------------------------------------------
		//	Update 'ps', 'ListStart', 'ListEnd':
		//	-	shorten column lenghts of columns that have non-zeros in pivot
		//		row,
		//	-	mark sets that became empty in the process.
		//
		if( RowMark[ PivotRow ] == EMPTY )
		{
			RmRowFromList( *this, ps, PivotRow, ListStart, ListEnd );
			RowMark[ PivotRow ] = NON_ZERO | PIVOT;
		}
	}

	//--------------------------------------------------------------------------
	//	Loop in which original variables are inserted into basis.
	//
	do
	{
		ColNum = ListStart;

		do									// Single pivot search loop.
		{
			//------------------------------------------------------------------
			//	Find a pivot column. Break if no appropriate column found.
			//
			PivotRow = PivotCol = -1;
			PivotCol = FindPivotCol( ps, m, ColNum );
			if( PivotCol == -1 ) break;

			//------------------------------------------------------------------
			//	Try to find an acceptable pivot in column 'PivotCol'.
			//
			PivotRow = FindPivotRow( this, RowMark, RowMax, PivotCol,
				PivotTol );

			//------------------------------------------------------------------
			//	The column (whether accepted, or rejected) is not to be scanned
			//	any more. Therefore we remove it from the list.
			//
			RmColFromList( ps, PivotCol, ListStart, ListEnd );

			//------------------------------------------------------------------
			//	If the non-zeros of column "PivotCol" all failed the stability
			//	test, we move the "ColNum" to the next column. If there is
			//	no next column - we break out of the loop.
			//
			if( PivotRow == -1 && ( ColNum = ps[ PivotCol ].next ) == -1 )
				break;
		} while( PivotRow == -1 );
		// End of pivot search loop.
		//----------------------------------------------------------------------

		//----------------------------------------------------------------------
		//	If there is a good pivot at "[PivotRow,PivotCol]", we need to update
		//	"ps" and "RowMark" now.
		//
		if( PivotRow >= 0 )
		{
			assert( PivotRow < m );
			assert( A2B[ PivotCol ] == A2B_UNDEF );

			//------------------------------------------------------------------
			//	Store basic column number.
			//
			A2B[ PivotCol ] = A2B_BASIC;
			PivotCnt++;

			//------------------------------------------------------------------
			//	Update "ps", "ListStart", "ListEnd":
			//	-	shorten column lenghts of columns where that have non-zeros
			//		in rows in which NEW non-zeros appeared,
			//	-	mark sets that became empty in the process.
			//
			Ptr<Real_T> a;
			Ptr<Int_T> row;
			Int_T len;

			for( MPS_LP::GetColumn( PivotCol, a, row, len ); len; --len, ++row )
			{
				//--------------------------------------------------------------
				//	If a new non zero exists at position [row,PivotCol],
				//	update sets and 'RowMark'.
				//
				if( RowMark[ *row ] == EMPTY )
				{
					RmRowFromList( *this, ps, *row, ListStart, ListEnd );
					RowMark[ *row ] |= NON_ZERO;
				}
			}

			assert( RowMark[ PivotRow ] != EMPTY );

			RowMark[ PivotRow ] |= PIVOT;	// Mark pivot row in 'RowMark'

			assert( ps[ PivotCol ].len == -1 );
		}
	} while( PivotCol >= 0 && PivotCnt < m );
	//
	//	End search if the whole basis was constructed or if no more pivots
	//	can be found.
	//--------------------------------------------------------------------------
	
	//--------------------------------------------------------------------------
	//	Report the results for original columns.
	//	
	if( Verbosity >= V_HIGH )
		Print( "\nCRASH:\n\tCols: %6d original,", (int) PivotCnt );

	//--------------------------------------------------------------------------
	//	Fill the remainder of basis with artificial columns (if necessary).
	//	Put the artificial columns in the lambda vector. Don't worry about the
	//	sign of the lambda coefficient. It may be changed later.
	//
	if( PivotCnt < MPS_LP::GetM() )
	{
		for( PivotCnt = i = 0; i < MPS_LP::GetM(); i++ )
			if( !( RowMark[i] & PIVOT ) )
			{
				// Artificial var. coefficient.
				Lambda[i] = 1.0e0;

				// Variable goes into basis and is marked as artificial.
				A2B[ n + SlackLen + i ]	= A2B_BASIC;

				LambdaVT[i]		= Short_T( VT_NORM | VT_SLACK | VT_ARTIF );

				PivotCnt++;
			}
	}
	else
		PivotCnt = 0;

	if( Verbosity >= V_HIGH )
		Print( "%6d artificial.\n", (int) PivotCnt );
}


#ifndef NDEBUG
/*------------------------------------------------------------------------------

	static void CheckSets( Array<PS> &ps, Int_T n, Array<Int_T> &PsStart,
		Array<Int_T> &PsEnd, Array<Int_T> &PsLen )

PURPOSE:
	Checks integrity of preference sets after their creation and modification.
WARNING: Compiles only when not in NDEBUG mode (see top of the file).

PARAMETERS:
	Array<PS> &ps
		Preference sets (thre of them) are passed by this argument.

	Int_T n
		Number of original matrix columns variables.

	Int_T *PsStart
		This is the table of indices of the first elements of the sets.

	Int_T *PsEnd
		This is the table of indices of the last elements of the sets.

	Int_T *PsLen
		This is the table of sets' lengths.

RETURN VALUE:
	None.

SIDE EFFECTS:
	May abort the program if problems are detected in preference sets data
structures.

------------------------------------------------------------------------------*/

static void CheckSets( Array<PS> &ps, Int_T n, Array<Int_T> &PsStart, // )
	Array<Int_T> &PsEnd, Array<Int_T> &PsLen )
{
	Int_T Cnt[3] = { 0, 0, 0 },
		tp, tn ;

	for( Int_T i = PS_FREE; i <= PS_REST; i++ )
	{
		Int_T j;
		for( j = PsStart[i]; j != -1; j = ps[j].next )
		{
			Cnt[i]++;
			tp = ps[j].previous;
			assert( tp >= -1 && tp < n );

			if( tp == -1 ) continue;		// Skip first elements of the sets.

			tn = ps[ tp ].next;
			assert( tn >= -1 && tn < n );
			assert( tn == j );
		}
		assert( Cnt[i] == PsLen[i] );
		assert( !( j != -1 && ps[j].next == -1 && PsEnd[i] != j ) );
	}
}
#endif


/*------------------------------------------------------------------------------

	static Int_T FindPivotCol( Array<PS> &ps, Int_T m, Int_T ColNum )

PURPOSE:
	Scan the given list of columns for a column that ensures the best
triangularity - will introduce a smallest number of new non-zeros into the
basis. Two other criteria of column choice are:
-	belonging to a higher preference set,
-	having higher price (meaningful only for columns corresponding to bounded
	and fixed variables).

Current preference order is:
-	length (number of new non-zeros),
-	type (number of preference set) and
-	price (computed during preference sets' setup).

PARAMETERS:
	Array<PS> &ps
		Preference sets are passed with this argument.
	
	Int_T ColNum
		This is the number of column from which the search proceeds.

RETURN VALUE:
	Number of found candidate pivot column, or (if none is found) -1.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Int_T FindPivotCol( Array<PS> &ps, Int_T m, Int_T ColNum )
{
	//--------------------------------------------------------------------------
	//	Find shortest (but not EMPTY) column with the biggest price and lowest
	//	set number.
	//
	Int_T Len = Int_T( m + 1 ),
		OrigLen = Len,
		PivotCol;
	Short_T Type = PS_NONE;
	Real_T Price = -INFINITY;

	for( PivotCol = -1; ColNum != -1; ColNum = ps[ ColNum ].next )
	{
		PS &p = ps[ ColNum ];

		assert( p.len > 0 );

		//----------------------------------------------------------------------
		//	Preference sets are placed on the list in their order of
		//	preference. Therefore we may skip scanning other columns if
		//	we have already found an acceptable column in "better" set.
		//
		if( p.type > Type )
			break;
		else if( p.type == Type )
		{
			if( p.len > Len )
				continue;
			else if( p.len == Len )
			{
				if( ( p.origLen > OrigLen ) ||
					( p.origLen == OrigLen && p.price <= Price ) )
					continue;
			}
		}

		OrigLen		= p.origLen;
		Len			= p.len;
		Price		= p.price;
		PivotCol	= ColNum;
		if( Len == 1 ) break;
	} // End pivot search loop.

	return PivotCol;
}


/*------------------------------------------------------------------------------

	static Int_T FindPivotRow( const SimplexLP *LP, Array<Short_T> &RowMark,
		Array<Real_T> &RowMax, Int_T PivotCol, Real_T PivTol )

PURPOSE:
	Finds a pivot in given column. Accept almost any column singleton. In
non-singleton columns choose the biggest new (not marked in 'RowMark') non-zero
element. See if it passes the stability test:

	pivot >= PivTol * maximum element in pivot row

If no suitable element is found return -1.

PARAMETERS:
	Short_T *RowMark
		Table which holds information on which rows have been pivotal and which
		have not.

	Real_T *RowMax
		Table of maximum row non-zeros (needed to perform the stabili tests.

	Int_T PivotCol
		Pivot column number.

	Real_T PivTol
		Pivot tolerance (see equation above).

RETURN VALUE:
	Pivot row number - if one was found - or -1 otherwise.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Int_T FindPivotRow( const SimplexLP *LP, Array<Short_T> &RowMark, // )
	Array<Real_T> &RowMax, Int_T PivotCol, Real_T PivTol )
{
	Int_T PivotRow, len;
	Ptr<Real_T> a;
	Ptr<Int_T> row;

	LP->MPS_LP::GetColumn( PivotCol, a, row, len );

	//--------------------------------------------------------------------------
	//	Accept almost any column singleton.
	//
	if( len == 1 && fabs( *a ) > 1.0e-3 && RowMark[ *row ] == EMPTY )
		return *row;

	//--------------------------------------------------------------------------
	//	Compare non-zeros with "RowMax[]".
	//	Choose the biggest of those that pass the stability test.
	//
	Real_T aMax = 0.0e0, aa = 0.0;

	for( PivotRow = -1; len; --len, ++a, ++row )
		if( RowMark[ *row ] == EMPTY && ( aa = fabs( *a ) ) > aMax )
		{
			PivotRow = *row;
			aMax = aa;
		}

	if( aMax >= PivTol * RowMax[ PivotRow ] )
		return PivotRow;
	else
		return -1;
}


/*------------------------------------------------------------------------------

	static void RmColFromList( Array<PS> &ps, Int_T Col, Int_T &lStart,
		Int_T &lEnd )

PURPOSE:
	Removes a single column from the given list. "PS" structure is a single
list node (fields "next" and "previous" link the list). "lStart" and "lEnd"
(indices of the first and the last element of the list) may need to be updated
if "Col" is respectively the first or the last column on the list.

PARAMETERS:
	Array<PS> &ps
		Preference sets are passed here.

	Int_T Col
		Number of column to be removed.

	Int_T &lStart
		Reference to the number of the first element of the list.

	Int_T &lEnd
		Reference to the number of the last element of the list.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static void RmColFromList( Array<PS> &ps, Int_T Col, Int_T &lStart, // )
	Int_T &lEnd )
{
	PS &p = ps[ Col ];

	if( p.previous != -1 )
		ps[ p.previous ].next = p.next;
	else
		lStart = p.next;

	if( p.next != -1 )
		ps[ p.next ].previous = p.previous;
	else
		lEnd = p.previous;

	p.previous = p.next = p.len = -1;
}


/*------------------------------------------------------------------------------

	static void RmRowFromList( const SimplexLP &LP, Array<PS> &ps,
		Int_T PivotRow, Int_T &ListStart, Int_T &ListEnd )

PURPOSE:
	Updates preference sets by removing a row from the sets. This is done by
shortening the working lengths of the columns which have a non-zero in the row
being removed. Columns whose lengths are reduced to zero are removed from the
preference sets.
	Uses "SimplexLP::RmColFromList" to actually remove the columns.

PARAMETERS:
	Array<PS> &ps
		Preference sets.

	Int_T PivotRow
		Number of the row to be removed.

	Int_T &ListStart
		Reference to the number of the first element of the list (passed to
		"SimplexLP::RmColFromList".

	Int_T &ListEnd
		Reference to the number of the last element of the list (passed to
		"SimplexLP::RmColFromList"

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static void RmRowFromList( const SimplexLP &LP, Array<PS> &ps, // )
	Int_T PivotRow, Int_T &ListStart, Int_T &ListEnd )
{
	Ptr<Real_T> a;
	Ptr<Int_T> col;
	Int_T len;

	for( LP.SolvableLP::GetRow( PivotRow, a, col, len ); len; --len, ++col )
	{
		Int_T &l = ps[ *col ].len;

		if( l > 0 )
		{
			//------------------------------------------------------------------
			//	Remove columns that were reduced to zero length.
			//	
			if( --l == 0 )
				RmColFromList( ps, *col, ListStart, ListEnd );
		}
		else
			assert( l == -1 );
	}
}
