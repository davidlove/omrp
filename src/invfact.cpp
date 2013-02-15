/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines for sparse matrices.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	invfact.cpp
CREATED:			1991.12.29 by prof. A. Ruszczynski
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		smartptr.h, stdtype.h, std_tmpl.h, error.h, inverse.h,
					invaux.h, std_math.h, work_vec.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Inverse::lufactor()

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

------------------------------------------------------------------------------*/

#include <math.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif

#ifndef __INVERSE_H__
#	include "inverse.h"
#endif
#ifndef __INVAUX_H__
#	include "invaux.h"
#endif



#define MIN_RATIO		(1.0e-5)
#define MAX_COMPRESSES	(5)

/*------------------------------------------------------------------------------

	Short_T Inverse::Factor( void )

PURPOSE:
	Factorize a square basis matrix ('a') of dimension 'n' with row numbers in 
'irow[]' and column numbers in 'jcol[]'.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	Return code:
	 	1	:	successful,
		-1	:	'a' has an empty row or column,
		-2	:	'a' is singular,
		-3	:	global bound stabglo on factors is too small.

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Short_T Inverse::Factor( void )
{
	RefactCnt++;

	size_t				indent;
	WorkVector<Int_T>	adr( n );
	Real_T				M_entry;
	Short_T				code = 1;
	Int_T				ir, nz, ii, iter, ipiv, jpiv,
						i, j, k, kl, kp, klc, kc, kr, kc0;
	WorkVector<Int_T>	rpre( n ), rsuc( n ), cpre( n ), csuc( n );

	adr.Fill( 0L, n );

	//--------------------------------------------------------------------------
	//	Remove small entries, count elements in rows and columns.
	//
	CountRowsAndCols();
	M_len = 0;
	U_max = A_max;

	//--------------------------------------------------------------------------
	//	Check for empty row or column (if the diagonal matrix is not available).
	//
	for( i = 0; i < n; i++ )
		if( !rlen[i] || !clen[i] )
		{
			if( !DiagPresent )
				return -1;
			else
				ReplaceEmptyRowsAndColumns();
		}

	//--------------------------------------------------------------------------
	//	'cmprs' is the number of compresses since the last reallocation.
	//
	cmprs = 0;

	//--------------------------------------------------------------------------
	//	Reorder 'a' and 'jcol' by rows.
	//	Put the largest element to the front of each row. Construct row cliques.
	//	Set up bidirectional lists of rows/columns of equal lengths.
	//
	Reorder( n, U_len, a, jcol, rptr, irow );
	LargestInAllRowsToFront();
	ConstructColumnCliques();
	ConstructEqualLengthLists( rpre, rsuc, cpre, csuc );

	//-------------------------------------------------------------------------
	//
	//	S T A R T   O F   T H E   M A I N   E L I M I N A T I O N   L O O P .
	//
	//-------------------------------------------------------------------------

	for( iter = 0; iter < n; iter++ )
	{
#		ifdef FACTOR_DEBUG
			CheckIntegrity( rpre, cpre );
#		endif

		if( !FindPivot( ipiv, jpiv, rsuc, csuc ) )
		{
			if( ipiv < 0 || jpiv < 0 )
			{							//	No acceptable pivot has been found.
				if( !alert )			//	Singular matrix.
					return -2;
				else					//	Stabglo must be too big.
					return -3;
			}
		}

		assert( ipiv >= 0 && jpiv >= 0 && ipiv < n && jpiv < n );

		//----------------------------------------------------------------------
		//	Remove rows and columns involved in elimination from the doubly
		//	linked lists.
		//
		RemoveRowsFromList( jpiv, rsuc, rpre );
		RemoveColsFromList( ipiv, csuc, cpre );

		//----------------------------------------------------------------------
		//	Store pivot (and permutation data).
		//
		rpre[ ipiv ] = Int_T( -2 - iter );
	 	cpre[ jpiv ] = Int_T( -2 - iter );

		//----------------------------------------------------------------------
		//	Eliminate pivotal row from column file and find pivot in row file.
		//	'kp' is first in pivot row. 'kr' is pivot's position in the pivot
		//	row.
		//
		kr = ElimPivotRowFromColumnFile( ipiv, jpiv );

		//----------------------------------------------------------------------
		//	Bring pivot to front of pivotal row.
		//
		kp = rptr[ ipiv ];
		assert( kp >= 0 );

		if( kp != kr )
		{
		 	Swap( a[ kr ], a[ kp ] );
		  	Swap( jcol[ kr ], jcol[ kp ] );
		}

		//----------------------------------------------------------------------
		//	Perform elimination itself, looping on entries of pivot column
		//	clique.
		//
		klc = clen[ jpiv ];

		for( kc0 = 0; kc0 < klc; kc0++ )
		{
			kc = Int_T( cptr[ jpiv ] + kc0 );
			assert( kc >= 0 && kc < colend );

	 		ir = irow[ kc ];// The row for which elimination is be done.

			//
			//	Pack the row file to make room for new row and M-entry.
			//
			Int_T rMax = Int_T( 2 * rlen[ ir ] + 2 * rlen[ ipiv ] + M_len + 1 ),
				cMax = Int_T( n - iter + M_len + 3 ),
				maxLen = (Int_T) Max( rowend + rMax, colend + cMax );

			MakeRoomInTheRowFile( maxLen, Max( rMax, cMax ) );

			assert( rowend + rMax <= ia );

			//
			//	Place pivot row indices (excluding pivot) in 'adr'.
			//
			UnpackPivotRowIndice( ipiv, adr );

			M_entry = EliminateElementInPivotColumn( ipiv, jpiv, ir );

			//------------------------------------------------------------------
			//	Compute new row at positions where row ir has non-zeros.
			//	Put the updated row at the end of the row file if necessary.
			//	Set entries of the adr array to NULL as they are used.
			//
			EliminateNoFillIn( ir, M_entry, adr );

			//------------------------------------------------------------------
			//	Find fill-ins by scanning the work array for indice. Add
			//	fill-ins to the row file, at the end, and to column cliques.
			//
			EliminateWithFillIn( ipiv, ir, M_entry, adr );

			if( rlen[ir] == 0 )
				return -2;
			else
				assert( rlen[ir] >= 0 );

			//------------------------------------------------------------------
			//	Put the largest element to the front of the transformed row.
			//
			U_max = Max( U_max, LargestInRowToFront( ir ) );

			//------------------------------------------------------------------
			//	Store multiplier. Make an entry for l_inverse at the other end
			//	of 'a', 'irow' and 'jcol'.
			//
			CompressColumnFile( 1 );
			StoreMultiplier( ipiv, ir, M_entry );
		}

		//----------------------------------------------------------------------
		//	Insert rows and columns involved in elimination into linked lists of
		//	equal numbers of non-zeros.
		//
		indent	= cptr[ jpiv ];
		k		= clen[ jpiv ];
		while( k-- )
		{
			ir					= irow[ indent ];
			irow[ indent++ ]	= -1;
	 		nz					= rlen[ ir ];

			//------------------------------------------------------------------
			//	This is where we detect dependent rows. A row that was
			//	involved in the elimination has no more non-zeros in the
			//	remaining sub-matrix.
			//
	 		if( nz )
			{
				if( rlst[ nz - 1 ] == -1 )
				{
					rlst[ nz - 1 ] = ir;
					rpre[ ir ] = rsuc[ ir ] = ir;
				}
				else
				{
					rsuc[ ir ] = rlst[ nz - 1 ];
					rpre[ ir ] = rpre[ rsuc[ ir ] ];
					rpre[ rsuc[ ir ] ] = rsuc[ rpre[ ir ] ] = ir;
				}
	 		}
			else
				code = -2;
		}
	 	clen[ jpiv ]	= 0;
		indent			= rptr[ ipiv ] + 1;
		k				= Int_T( rlen[ ipiv ] - 1 );

		while( k-- )
		{
			ir = jcol[ indent++ ];
	 		nz = clen[ ir ];

			//------------------------------------------------------------------
			//	Detect dependent columns.
			//
	 		if( nz )
			{
				if( clst[ nz - 1 ] == -1 )
				{
					clst[ nz - 1 ] = ir;
					cpre[ ir ] = csuc[ ir ] = ir;
				}
				else
				{
					csuc[ ir ] = clst[ nz - 1 ];
					cpre[ ir ] = cpre[ csuc[ ir ] ];
					cpre[ csuc[ ir ] ] = csuc[ cpre[ ir ] ] = ir;
				}
			}
	 		else
				code = -2;
		}
	}
	//
	//	End of main elimination loop.
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	//	Reset column cliques for U and store row/column numbers in pivotal order
	//	in 'rlst' and 'clst'.
	//
	for( i = 0; i < n; i++ )
	{
		j	   	= Int_T( -2 - rpre[i] );
		rlst[j]	= i;
		j	   	= Int_T( -2 - cpre[i] );
		clst[j]	= i;
		clen[i]	= 0;
	}

	for( i = 0; i < n; i++ )
 		for( k = rptr[i], kl = Int_T( k + rlen[i] ); k < kl; k++ )
			clen[ jcol[ k ] ]++;

	for( k = i = 0; i < n; i++ )
		cptr[i] = k += clen[i];

	colend = k;
	for( ii = 0; ii < n; ii++ )
	{
 		i	= rlst[ ii ];
 		kp	= rptr[i];
 		kl	= Int_T( rlen[i] + kp );
		for( k = kp; k < kl; k++ )
			irow[ --cptr[ jcol[ k ] ] ] = i;
	}

//==============================================================================


//------------------------------------------------------------------------------
//	code == -1: row/column 'i' has no elements
//	code == -2: singular matrix detected by factorization
//	code == -3: numerical difficulties
//------------------------------------------------------------------------------

	return code;
}


/*------------------------------------------------------------------------------

	f()

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

void Inverse::CountRowsAndCols( void )
{
	if( clen[0] == 0 )
		for( Int_T l = 0; l < U_len; )
		{
			Real_T aa = fabs( a[l] );

			if( IsZero( aa ) )	// Small entry
			{
				U_len--;
				a[l]	= a[ U_len ];
				irow[l]	= irow[ U_len ];
				jcol[l]	= jcol[ U_len ];
			}
			else				// Count elements in rows and columns.
			{
				rlen[ irow[l] ]++;
				clen[ jcol[l] ]++;
				if( A_max < aa ) A_max = aa;
				l++;
			}
		}

	colend = rowend = U_len;
}


/*------------------------------------------------------------------------------

	void Inverse::LargestInAllRowsToFront( void )

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

void Inverse::LargestInAllRowsToFront( void )
{
	for( Int_T i = 0; i < n; i++ )
		LargestInRowToFront( i );
}


/*------------------------------------------------------------------------------

	void Inverse::ConstructColumnCliques( void )

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

void Inverse::ConstructColumnCliques( void )
{
	//--------------------------------------------------------------------------
	//	Initialize cptr[i] to point just beyond where the last component of
	//	column 'i' of 'a' will be stored.
	//
	for( Int_T k = 0, i = 0; i < n; i++ )
		cptr[i] = k += clen[i];

	//--------------------------------------------------------------------------
	//	Reconstruct the file of column cliques containing row numbers for each
	//	column. By putting the entries in backwards and decreasing 'cptr[j]'
	//	each time it is used we automatically leave it pointing to the first
	//	element of clique j.
	//
	for( Int_T row = Int_T( n - 1 ), kEnd = U_len; row >= 0; row-- )
	{
 		Int_T kStart = rptr[ row ];

 		for( Int_T kk = kStart; kk < kEnd; kk++ )
 			irow[ --cptr[ jcol[kk] ] ] = row;

		kEnd = kStart;
	}
}


/*------------------------------------------------------------------------------

	void Inverse::ConstructEqualLengthLists( Array<Int_T> &rpre,
		Array<Int_T> &rsuc, Array<Int_T> &cpre, Array<Int_T> &csuc )

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

void Inverse::ConstructEqualLengthLists( Array<Int_T> &rpre, // )
	Array<Int_T> &rsuc, Array<Int_T> &cpre, Array<Int_T> &csuc )
{
	//--------------------------------------------------------------------------
	//	Initialize the lists.
	//
	rpre.Fill( -1, n );
	rsuc.Fill( -1, n );
	cpre.Fill( -1, n );
	csuc.Fill( -1, n );

	rlst.Fill( -1, n );
	clst.Fill( -1, n );

	//--------------------------------------------------------------------------
	//	Set up bidirectional circular lists of rows and cols with equal numbers
	//	of non-zeros.
	//
	for( Int_T i = 0; i < n; i++ )
	{
		Int_T nz;	// Row/column length

		nz			= rlen[i];
		rsuc[i]		= rlst[--nz];
		rlst[nz]	= i;
		if( rsuc[ i ] >= 0 )
			rpre[ rsuc[i] ] = i;

		nz 			= clen[ i ];
		csuc[ i ]	= clst[ --nz ];
		clst[ nz ]	= i;
		if( csuc[ i ] >= 0 ) cpre[ csuc[ i ] ] = i;
	}

	//--------------------------------------------------------------------------
	//	Let the last list entry point to the first (bidirectionally). Loop on
	//	the row/column lengths.
	//
	for( Int_T nz = 0; nz < n; nz++ )
	{
		Int_T i;

		if( ( i = rlst[ nz ] ) != -1 )
		{
			Int_T ii = i;
			for( Int_T j; ( j = rsuc[ii] ) != -1; ii = j );

			rsuc[ii]	= i;
			rpre[i]		= ii;
		}

		if( ( i = clst[ nz ] ) != -1 )
		{
			Int_T ii = i;
			for( Int_T j; ( j = csuc[ii] ) != -1; ii = j );

			csuc[ii]	= i;
			cpre[i]		= ii;
		}
	}
}


/*------------------------------------------------------------------------------

	Bool_T FindPivot( Int_T &ipiv, Int_T &jpiv, const Ptr<Int_T> rsuc,
		const Ptr<Int_T> csuc )

PURPOSE:
	Finds a pivot element (ipiv,jpiv) which is the best according to Markowitz's
criterion and passes the stability test.

PARAMETERS:
	Int_T &ipiv, Int_T &jpiv
		The pivot coordinates. The values on entry is ignored. On exit it shall
		contain the coord's or (-1,-1) if no suitable pivot element was found.

	const Ptr<Int_T> rpre, const Ptr<Int_T> rsuc
	const Ptr<Int_T> cpre, const Ptr<Int_T> csuc
		Bidirectional lists of rows and columns, respectively, of the same
		length.

RETURN VALUE:
	'True' if pivot search is over; 'False' otherwise.
	Note: 'False;' does not have to mean that a pivot has not been found.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Inverse::FindPivot( Int_T &ipiv, Int_T &jpiv, const Ptr<Int_T> rsuc, // )
	const Ptr<Int_T> csuc )
{
	//--------------------------------------------------------------------------
	//	Find pivot. 'bestCost' is markowitz cost of the cheapest pivot found so
	//	far, which is in row 'ipiv' and column 'jpiv'.
	//
	Real_T bestCost	= Real_T( n ) * Real_T( n );
	Real_T ratio	= 0.0;
	Int_T ties		= 0;

	ipiv = jpiv = -1;

	//--------------------------------------------------------------------------
	//	Pivot selection loop is on length of column and row to be searched.
	//	'nz' is the number of non-zeros excluding the potential pivot.
	//
	for( Int_T nz = 0; nz < n; nz++ )
	{
		if( bestCost <= Real_T( nz ) * Real_T( nz ) ) break;

		//----------------------------------------------------------------------
		//	Search columns with 'nz + 1' non-zeros. Follow linked list.
		//
		Int_T start = clst[nz];

		if( start != -1 )
		{
			Int_T col = start;

			do
			{
				assert( col >= 0 && col < n );

				if( FindPivotInCol( col, nz, ipiv, jpiv, bestCost, ties,
					ratio ) )
					return True;

			} while( ( col = csuc[col] ) != start );
		}

		//----------------------------------------------------------------------
		//	Search rows with 'nz + 1' non-zeros. Follow linked list.
		//
		start = rlst[ nz ];

		if( start != -1 )
		{
			Int_T row = start;

			do
			{
				assert( row >= 0 && row < n );

				if( FindPivotInRow( row, nz, ipiv, jpiv, bestCost, ties,
					ratio ) )
					return True;

			} while ( ( row = rsuc[row] ) != start );
		}
	}

	return False;
}



/*------------------------------------------------------------------------------

	Bool_T Inverse::FindPivotInCol( const Int_T col, const Int_T nz,
		Int_T &ipiv, Int_T &jcol, Real_T &bestCost, Int_T &ties,
		Real_T &bestRatio )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	'True' if pivot search is over; 'False' otherwise.
	Note: 'False;' does not have to mean that a pivot has not been found.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Inverse::FindPivotInCol( const Int_T col, const Int_T nz, // )
	Int_T &ipiv, Int_T &jpiv, Real_T &bestCost, Int_T &ties, Real_T &bestRatio )
{
	//--------------------------------------------------------------------------
	//	Consider each row 'row' covered by column 'col'. Compute the Markowitz
	//	cost for the (row,col) pivot and see if it passes the stability test.
	//
	for( Int_T k = cptr[col], colEnd = Int_T( k + clen[col] ); k < colEnd; k++ )
	{
		Int_T row	= irow[k];
		Real_T cost	= (Real_T)nz * Real_T( rlen[row] - 1 ); // Markowitz cost.

		if( cost > bestCost )
			continue;

		//----------------------------------------------------------------------
		//	Find the non-zero value and perform a stability test. The largest
		//	element is first in its row.
		//
		Int_T pos = PositionInRowFile( row, col );

		assert( pos >= 0 );

		Real_T ratio = fabs( a[ pos ] / a[ rptr[row] ] );

		if( ( nz > 0 && ratio < stabrow ) ||
			( nz == 0 && ratio < MIN_RATIO ) ||
			( alert && fabs( a[ pos ] ) < A_max * stabglo) )
			continue;

		//----------------------------------------------------------------------
		//	Test for an acceptable element (Markowitz cost comparison, possible
		//	tiebreaking)
		//
		if( cost < bestCost || ( ++ties <= maxties && ratio > bestRatio ) )
		{
			bestCost	= cost;
			ipiv		= row;
			jpiv		= col;
			bestRatio	= ratio;
			ties		= 0;
		}
		else
			return True;

		if( bestCost <= Real_T( nz ) * Real_T( nz ) )
			return True;
	}

	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Inverse::FindPivotInRow( const Int_T row, const Int_T nz,
		Int_T &ipiv, Int_T &jpiv, Real_T &bestCost, Int_T &ties,
		Real_T &bestRatio )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	'True' if pivot search is over; 'False' otherwise.
	Note: 'False;' does not have to mean that a pivot has not been found.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Inverse::FindPivotInRow( const Int_T row, const Int_T nz, // )
	Int_T &ipiv, Int_T &jpiv, Real_T &bestCost, Int_T &ties, Real_T &bestRatio )
{
	Int_T rowStart	= rptr[row],
		rowEnd		= Int_T( rowStart + rlen[row] );

	assert( rowStart >= 0 && rowStart < rowend );
	assert( rowEnd >= 0 && rowEnd <= rowend );

	//--------------------------------------------------------------------------
	//	The largest element is in front of the row.
	//
	Real_T rowmax	= fabs( a[ rowStart ] ),
		trsh		= rowmax * stabrow;

	//--------------------------------------------------------------------------
	//	Find the column 'j' in row 'row' that passes the stability test and has
	//	the lowest Markowitz cost.
	//
	for( Int_T pos = rowStart; pos < rowEnd; pos++ )
	{
		Real_T cost = (Real_T)nz * Real_T( clen[ jcol[pos] ] - 1 );
											// Markowitz cost.

		if( cost > bestCost ) continue;

		Real_T x = fabs( a[pos] );

		//----------------------------------------------------------------------
		//	Perform stability test. Allow less stable column singletons.
		//
		if( ( nz > 0 && x < trsh ) ||
			( nz == 0 && x < MIN_RATIO ) ||
			( alert && x < A_max * stabglo ) )
			continue;

		Real_T ratio = x / rowmax;

		if( cost < bestCost || ( ++ties < maxties && ratio > bestRatio ) )
		{
			bestCost	= cost;
			ipiv		= row;
			jpiv		= jcol[pos];
			bestRatio	= ratio;
			ties		= 0;
		}
		else
			return True;

		if( bestCost <= Real_T( nz ) * Real_T( nz ) )
			return True;
	}

	return False;
}


/*------------------------------------------------------------------------------

	void Inverse::ReplaceEmptyRowsAndColumns( void )

PURPOSE:
	NOT IMPLEMENTED YET!

	Will (hopefully) replace empty rows/columns with entries from a diagonal
column passed to the Inverse object by a call to "RegisterDiagonal()" function.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::ReplaceEmptyRowsAndColumns( void )
{
	FatalError( "Empty row/column replacement not implemented yet." );
}


/*------------------------------------------------------------------------------

	void Inverse::RemoveRowsFromList( const Int_T col, Ptr<Int_T> rsuc,
		Ptr<Int_T> rpre )

PURPOSE:
	Updates the lists of rows of the same lengths as if the column "col" was
removed from the matrix. Each row entry in the list of rows is "shifted one
position": moved to the list of rows shorter by one element.

PARAMETERS:
	const Int_T col
		The number of column to be "removed".

	Ptr<Int_T> rsuc, Ptr<Int_T> rpre
		The lists of rows of the same lengths.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::RemoveRowsFromList( const Int_T col, Ptr<Int_T> rsuc, // )
	Ptr<Int_T> rpre )
{
	Int_T start	= cptr[ col ],
		len		= clen[ col ];

	while( len-- )
	{
		Int_T row	= irow[ start++ ],
			nz		= Int_T( rlen[row] - 1 );

		if( rsuc[row] == row )
			rlst[nz] = -1;
		else
		{
			rpre[ rsuc[row] ] = rpre[row];
			rsuc[ rpre[row] ] = rsuc[row];
			if( rlst[nz] == row )
				rlst[nz] = rsuc[row];
		}
	}
}


/*------------------------------------------------------------------------------

	void Inverse::RemoveColsFromList( const Int_T row, Ptr<Int_T> csuc,
		Ptr<Int_T> cpre )

PURPOSE:
	Updates the lists of columns of the same lengths as if the column "row" was
removed from the matrix. Each column entry in the list of columns is "shifted
one position": moved to the list of columns shorter by one element.

PARAMETERS:
	const Int_T row
		The number of row to be "removed".

	Ptr<Int_T> csuc, Ptr<Int_T> cpre
		The lists of rows of the same lengths.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::RemoveColsFromList( const Int_T row, Ptr<Int_T> csuc, // )
	Ptr<Int_T> cpre )
{
	Int_T start	= rptr[ row ],
	 	len		= rlen[ row ];

	while( len-- )
	{
		Int_T col	= jcol[ start++ ],
			nz		= Int_T( clen[col] - 1 );

		if( csuc[col]==col )
			clst[ nz ] = -1;
		else
		{
			cpre[ csuc[col] ] = cpre[col];
			csuc[ cpre[col] ] = csuc[col];
			if( clst[ nz ]==col )
				clst[ nz ] = csuc[col];
		}
	}
}


/*------------------------------------------------------------------------------

	void Inverse::MakeRoomInTheRowFile( Int_T maxLen, Int_T increase )

PURPOSE:
	Makes room in the row file by compression and possibly reallocation of data.
Reallocates if either compression would not suffice, or the number of compesses
is to high.

PARAMETERS:
	Int_T maxLen
		Required (maximum) length of the row file.

	Int_T increase
		Maximum increase of the current length of the file.

RETURN VALUE:
	None.

SIDE EFFECTS:
	The array 'rptr' may change. 'a', 'irow', 'jcol' may be reallocated (and
thus change their addresses.

------------------------------------------------------------------------------*/

void Inverse::MakeRoomInTheRowFile( Int_T maxLen, Int_T increase )
{
	if( maxLen + M_len > ia )
	{
		while( cmprs > MAX_COMPRESSES || U_len + U_len/5 + increase > ia )
		{
			cmprs = 0;

			Int_T OldM_End = Int_T( ia - 1 );

			ia = (Int_T) Max( Int_T( ia + ia/2 ), maxLen );

			Int_T NewM_End = Int_T( ia - 1 );

			a.Resize( ia );
			irow.Resize( ia );
			jcol.Resize( ia );

			for( Int_T i = 0; i < M_len; i++ )
			{
				a[ NewM_End ]		= a[ OldM_End ];
				irow[ NewM_End ]	= irow[ OldM_End ];
				jcol[ NewM_End-- ]	= jcol[ OldM_End-- ];
			}
		}

		Pack( a, jcol, rptr, n, rlen, True, rowend );
		cmprs++;
	}
}


/*------------------------------------------------------------------------------

	Int_T Inverse::ElimPivotRowFromColumnFile( Int_T ipiv, Int_T jpiv )

PURPOSE:
	Eliminate the pivot row 'ipiv' from the column file. Find the pivot
position in row file.

PARAMETERS:
	Int_T ipiv, Int_T jpiv
		Pivot row and column.

RETURN VALUE:
	Returns the pivot position in the row file as offset from the beginning of
the tables 'a' and 'jcol'.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Inverse::ElimPivotRowFromColumnFile( Int_T ipiv, Int_T jpiv )
{
	Int_T rowStart	= rptr[ ipiv ],
		rowEnd		= Int_T( rowStart + rlen[ ipiv ] );

	Int_T pivotPos = -1;
	for( Int_T i = rowStart; i < rowEnd; i++ )
	{
		Int_T col	= jcol[i],
			colEnd	= Int_T( cptr[col] + clen[col] ),
			pivPos	= PositionInColumnFile( ipiv, col );

		assert( pivPos >= 0 );

		irow[ pivPos ]	= irow[ --colEnd ];
		irow[ colEnd ]	= -1;
		clen[ col ]--;

		if( col == jpiv )
			pivotPos = i;
	}

	assert( pivotPos >= 0 && pivotPos < rowEnd );

	return pivotPos;
}


/*------------------------------------------------------------------------------

	void Inverse::CompressColumnFile( Int_T RequestedSpace )

PURPOSE:
	If necessary, the function compresses the column file. It is assumed (and
checked by a call to "assert()") that the compression will make enough space
available.

PARAMETERS:
	Int_T RequestedSpace
		The number of positions that may need to be added to the column file.

RETURN VALUE:
	None.

SIDE EFFECTS:
	The array 'cptr' may change.

------------------------------------------------------------------------------*/

void Inverse::CompressColumnFile( Int_T RequestedSpace )
{
	if( colend + RequestedSpace + M_len > ia )
	{
		Pack( a, irow, cptr, n, clen, False, colend );
		cmprs++;
	}
	assert( colend + M_len + RequestedSpace <= ia );
}


/*------------------------------------------------------------------------------

	void Inverse::StoreMultiplier( Int_T row, Int_T col, Real_T val )

PURPOSE:
	Make an entry for the Gaussian elimination factors at the other end of 'a',
'irow' and 'jcol'. Put the entries in backwards.

PARAMETERS:
	Int_T row, Int_T col, Real_T val
	The multiplier and its position in the matrix.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::StoreMultiplier( Int_T row, Int_T col, Real_T val )
{
	Int_T k	= Int_T( ia - ++M_len );

	a[k]	= val;
	irow[k]	= row;
	jcol[k]	= col;
}


/*------------------------------------------------------------------------------

	Real_T Inverse::LargestInRowToFront( Int_T row )

PURPOSE:
	This function rearanges a row so that the largest element is moved to the
front. This simplifies numerical stability checks.

PARAMETERS:
	Int_T row
		The row number.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T Inverse::LargestInRowToFront( Int_T row )
{
	Real_T rowMax	= 0.0;
	Int_T rowStart	= rptr[ row ],
		rowEnd		= Int_T( rowStart + rlen[ row ] ),
		pos			= -1;

	for( Int_T i = rowStart; i < rowEnd; i++ )
		if( fabs( a[i] ) > rowMax )
		{
			assert( jcol[i] >= 0 );	// One more place, where we can check if
									// the row still contains valid data.
			rowMax	= fabs( a[i] );
			pos		= i;
		}

	assert( pos >= 0 );

	if( pos != rowStart )
	{
		Swap( a[ rowStart ],	a[ pos ] );
		Swap( jcol[ rowStart ],	jcol[ pos ] );
	}

	return rowMax;
}


/*------------------------------------------------------------------------------

	void Inverse::MoveRowToEndOfFile( Int_T row )

PURPOSE:
	Moves the row to the end of the row file.

PARAMETERS:
	Int_T row
		The row number.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::MoveRowToEndOfFile( Int_T row )
{
	assert( row >= 0 && row < n );
	assert( rowend + M_len + rlen[row] <= ia );

	Int_T rowStart	= rptr[ row ],
		rowEnd		= Int_T( rowStart + rlen[ row ] );

	rptr[ row ] = rowend;

	for( Int_T i = rowStart; i < rowEnd; i++, rowend++ )
	{
		assert( jcol[i] >= 0 );

		a[ rowend ]		= a[i];
		jcol[ rowend ]	= jcol[i];
		jcol[i]			= -1;
	}
}


/*------------------------------------------------------------------------------

	void Inverse::MoveCliqueToEndOfFile( Int_T col )

PURPOSE:
	Moves column clique to the end of column file.

PARAMETERS:
	Int_T col
		Column number.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::MoveCliqueToEndOfFile( Int_T col )
{
	Int_T OldStart	= cptr[col],
		NewStart	= colend;

	cptr[col]		= colend;

	assert( colend + M_len + clen[col] <= ia );

	for( Int_T i = clen[col]; i; i-- )
	{
		irow[ NewStart++ ]	= irow[ OldStart ];
		irow[ OldStart++ ]	= -1;
	}

	colend += clen[col];
}


/*------------------------------------------------------------------------------

	void Inverse::UnpackPivotRowIndice( Int_T row, Array<Int_T> &adr )

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

void Inverse::UnpackPivotRowIndice( Int_T row, Array<Int_T> &adr )
{
	Int_T start = rptr[row],
		end = Int_T( start + rlen[row] );

	for( Int_T i = Int_T( start + 1 ); i < end; i++ )
		adr[ jcol[i] ] = i;
}


/*------------------------------------------------------------------------------

	Real_T Inverse::EliminateElementInPivotColumn( Int_T ipiv, Int_T jpiv,
		Int_T ielim )

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

Real_T Inverse::EliminateElementInPivotColumn( Int_T ipiv, Int_T jpiv, // )
	Int_T ielim )
{
	//--------------------------------------------------------------------------
	//	Search non-pivot row for element to be eliminated, that is,
	//	element in the pivot column.
	//
	Int_T elimEnd	= Int_T( rptr[ ielim ] + rlen[ ielim ] ),
		elimPos		= PositionInRowFile( ielim, jpiv );

	assert( elimPos >= 0 );

	//--------------------------------------------------------------------------
	//	Remove element to be eliminated from its row.
	//
	rlen[ ielim ]--;
	if( elimEnd >= rowend )
	{
		assert( elimEnd == rowend );

		rowend--;
	}
	elimEnd--;
	U_len--;

	Real_T M_entry	= - a[ elimPos ] / a[ rptr[ ipiv ] ];

	a[ elimPos ]	= a[ elimEnd ];
	jcol[ elimPos ]	= jcol[ elimEnd ];
	jcol[ elimEnd ]	= -1;

	return M_entry;
}


/*------------------------------------------------------------------------------

	void Inverse::EliminateNoFillIn( Int_T ielim, Real_T factor,
		Array<Int_T> &adr )

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

void Inverse::EliminateNoFillIn( Int_T ielim, Real_T factor, Array<Int_T> &adr )
{
	//--------------------------------------------------------------------------
	//	Compute new row at positions where row 'ielim' has non-zeros.
	//	Put the updated row at the end of the row file if necessary.
	//	Set entries of the adr array to ZERO as they are used.
	//
	for( Int_T i = rptr[ ielim ], rowEnd = Int_T( i + rlen[ ielim ] );
		i < rowEnd; )
	{
		Int_T col = jcol[i];

		if( adr[col] )
		{
			//
			//	Update element.
			//
			a[i]		+= factor * a[ adr[col] ];
			adr[col]	= 0;

			if( IsNonZero( a[i] ) )
				i++;
			else
			{
				//--------------------------------------------------------------
				// Remove small elements from U:
				//

				//	From the row file ...
				//
				rowEnd--;
				rlen[ ielim ]--;

				assert( jcol[ rowEnd ] >= 0 );

				a[i]			= a[ rowEnd ];
				jcol[i]			= jcol[ rowEnd ];
				jcol[ rowEnd ]	= -1;
				U_len--;

				//
				//	... and from the column file.
				//
				Int_T colPos	= PositionInColumnFile( ielim, col ),
					colEnd		= Int_T( cptr[col] + --clen[col] );

				assert( colPos >= 0 );

				irow[ colPos ]	= irow[ colEnd ];
				irow[ colEnd ]	= -1;
			}
		}
		else
			i++;
	}
}


/*------------------------------------------------------------------------------

	void Inverse::EliminateWithFillIn( Int_T ipiv, Int_T ielim, Real_T factor,
		Array<Int_T> &adr )

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

void Inverse::EliminateWithFillIn( Int_T ipiv, Int_T ielim, Real_T factor, // )
	Array<Int_T> &adr )
{
	Int_T pivStart	= rptr[ ipiv ],
		pivEnd		= Int_T( pivStart + rlen[ ipiv ] ),
		elimEnd		= Int_T( rptr[ ielim ] + rlen[ ielim ] );

	assert( pivStart >= 0 && pivStart < rowend );
	assert( pivEnd >= 0 && pivEnd <= rowend );
	assert( elimEnd >= 0 && elimEnd <= rowend );

	for( Int_T pivPos = Int_T( pivStart + 1 ); pivPos < pivEnd; pivPos++ )
	{
		Int_T col = jcol[ pivPos ];

		if( !adr[col] ) continue;

		Real_T U_entry = factor * a[ adr[col] ];

		//----------------------------------------------------------------------
		//	If necessary - move the row to end of file.
		//
		if( elimEnd < rowend && jcol[ elimEnd ] >= 0 )
		{
			assert( rowend + rlen[ ielim ] + M_len < ia );

			MoveRowToEndOfFile( ielim );
			elimEnd = rowend;
		}

		//----------------------------------------------------------------------
		//	Add new entry at end of the row.
		//
		assert( elimEnd + M_len < ia );

		a[ elimEnd ]	= U_entry;
		jcol[ elimEnd ]	= col;
		rlen[ ielim ]++;
		if( elimEnd >= rowend )
		{
			assert( elimEnd == rowend );
			rowend++;
		}
		elimEnd++;
		U_len++;

		//----------------------------------------------------------------------
		//	Create fill in column file.
		//
		Int_T colPos = Int_T( cptr[col] + clen[col] );

		assert( colPos <= colend );

		//----------------------------------------------------------------------
		//	Try to place new element at the end of the present entry.
		//
		if( ( colPos == colend && colend + M_len + 2 > ia ) ||
			( irow[ colPos ] >= 0 ) )
		{
			//
			//	Move the old clique to the end of column file before
			//	a new element is appended.
			//
			CompressColumnFile( Int_T( clen[col] + 2 ) );
			MoveCliqueToEndOfFile( col );

			colPos = colend++;
		}
		else if( colPos == colend )
			colend++;
		
		irow[ colPos ] = ielim;
		clen[col]++;
		adr[col] = 0;
	}
	//
	//	End of fill-in loop.
	//--------------------------------------------------------------------------

	rlen[ ielim ] = Int_T( elimEnd - rptr[ ielim ] );
}


#ifdef FACTOR_DEBUG
void Inverse::CheckIntegrity( const Ptr<Int_T> rpre, const Ptr<Int_T> cpre )
{
	//--------------------------------------------------------------------------
	//	Find all row file entries in the column file.
	//
	for( Int_T row = 0; row < n; row++ )
	{
		if( rpre[ row ] <= -2 ) continue;

		Int_T rowStart	= rptr[ row ],
			rowEnd		= rowStart + rlen[ row ];

		for( Int_T i = rowStart; i < rowEnd; i++ )
		{
			Int_T col = jcol[i];
			
			assert( col >= 0 && col < n );

			if( cpre[ col ] <= -2 ) continue;
			
			Int_T pos = PositionInColumnFile( row, col );
			
			assert( pos >= 0 );
		}
	}


	//--------------------------------------------------------------------------
	//	Find all column file entries in the row file.
	//
	for( Int_T col = 0; col < n; col++ )
	{
		if( cpre[ col ] <= -2 ) continue;

		Int_T colStart	= cptr[ col ],
			colEnd		= colStart + clen[ col ];

		for( Int_T i = colStart; i < colEnd; i++ )
		{
			Int_T row = irow[i];

			if( rpre[ row ] <= -2 ) continue;

			assert( row >= 0 && row < n );

			Int_T pos = PositionInRowFile( row, col );

			assert( pos >= 0 );
		}
	}
}
#endif
