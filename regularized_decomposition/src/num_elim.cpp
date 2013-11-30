/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	num_elim.cpp
CREATED:			1994.02.20
LAST MODIFIED:		1995.10.05

DEPENDENCIES:		std_math.h, std_tmpl.h, work_vec.h, smartptr.h, error.h,
					memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h, vec_pool.h,
					presolve.h, simplex.h, stdtype.h, lp_codes.h, solv_lp.h,
					compile.h, mps_lp.h, sort_lab.h, myalloc.h, cl_list.h
					<assert.h>

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

#include <assert.h>

#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif

#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif



/*------------------------------------------------------------------------------

	Int_T Presolver::NumericalEliminations( Real_T PivMax )

PURPOSE:
	This function is designed to make the sparse constraint matrix sparser by
	means of numerical eliminations performed on the rows of the constraint
	matrix. A row is used as a pivot row in elimination, if there exists another
	row, which has a non-zero pattern which is a superset of the non-zero
	pattern of the pivot row. Only equality rows may be pivot rows.

ALGORITHM DESCRIPTION:
	We search for pivot rows starting with the shortest ones, and then we pursue
	the search until all the rows are scanned. Since a singleton eqality row
	will be marked and changed to fixed bounds on a variable, we start the scan
	from doubleton rows.

	One row may be eliminated by another if the eliminating row is an equality
	row. Additionally we introduce a numerical stability criterion: we only
	allow pivots, for which absolute value of logarithm of absolute value

		|         |
		| log |p| |
		|         |

	is smaller than a prescribed tolerance

		log P

	where P > 1 is passed as argument 'PivMax'.

	When equality row 'a' is used as a pivot row in elimination on row 'b' with
	factor 'F', the type of the eliminated row 'b' does not change.

	Actual elimination is performed in blocks. One procedure finds a pivot row
	and numerical elimination factors for each row, that may be eliminated.
	Another procedure performs this block elimination on both row and column
	structures.

PARAMETERS:
	See above: ALGORITHM DESCRIPTION

RETURN VALUE:
	Function returns number of eliminated non-zeros.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Presolver::NumericalEliminations( Real_T PivMax )
{
	assert( PivMax > 1.0 );

	//--------------------------------------------------------------------------
	//	'Scan' is an 'm'-vector. Its entries are:
	//	-	'True' if the corresp. row is a potential pivot candidate or
	//	-	'False' if it has been scanned previously and has not changed (by
	//		elimination with some other row), or if it is not an equality row..
	//
	//	'ScanCnt' holds the number of equality rows still to be scanned.
	//	Initially it is set to the current number of rows still in the
	//	constraint matrix.
	//
	WorkVector<Bool_T> Scan( m );
	Int_T ScanCnt;
	Scan.Fill( False, m );

	//--------------------------------------------------------------------------
	//	'Factor' is an 'm'-vector of potential elimination factors. It is
	//	initialized to zeros. Later (for every candidate pivot row) only
	//	admissible values are stored at positions corresponding to the rows to
	//	be eliminated.
	//
	WorkVector<Real_T> Factor( m );
	Int_T i;
	Int_T ElimNZ = 0;

	for( ScanCnt = i = 0; i < m; i++ )
		if( !ExcludeRows[i] && rt[i] == VT_FIXED && RowLen[i] >= 2 )
		{
			Scan[i] = True;
			ScanCnt++;
		}

	//--------------------------------------------------------------------------
	//	Loop on row lengths until no more rows to be scanned remain.
	//
	for( Int_T RowLength = 0; ScanCnt > 0; RowLength++,
		RowLength %= Int_T( n + 1 ) )
	{

		//----------------------------------------------------------------------
		//	Check out all the rows of length 'RowLength'.
		//
		for( i = RowList.GetFirst( RowLength ); i >= 0; i = RowList.GetNext() )
		{
			if( Scan[i] && ( RowLen[i] < 2 || ExcludeRows[i] ) )
			{
				Scan[i] = False;
				ScanCnt--;
				continue;
			}
			else if( !Scan[i] || ExcludeRows[i] )
				continue;

			if( IsNonZeroSubset( i, Factor, PivMax ) )
			{
				Int_T elr	= NumElimRows( i, Factor );
#ifndef NDEBUG
				Int_T elc	=
#endif
							  NumElimCols( i, Factor, Scan, ScanCnt );

				assert( elc == elr );

				ElimNZ += elr;
			}

			Scan[i] = False;
			if( --ScanCnt <= 0 ) break;

		} // End of loop on rows of length 'RowLength'.

	} // End of loop on row legths.

	return ElimNZ;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::IsNonZeroSubset( Int_T piv, Array<Real_T> &Factor,
		const Real_T fMax )

	(auxiliary function)
		inline double Q( double f )

PURPOSE:
		This function determines if row number 'piv' may be used as a pivot row
	in numerical elimination on other rows. A row may be used as a pivot row if
	and only if:
	0.	row 'piv' is an equality row,
	1.	there are other rows, with a non-zero pattern, which is a superset of
		the non-zero pattern of the pivot row.
	2.	there exists at least one elimation factor which fulfills the numerical
		stability criterion described in function
		'Presolver::NumericalEliminations()' comments.

	ATTENTION:

	The numerical elimination factors are chosen by the following rules:
	-	factor 'a' is considered better than 'b' if absolute value of 'a' is
		logarithmically closer to 1, i.e.

			|         |   |         |
			| log |a| | < | log |b| |
			|         |   |         |

	-	absolute value of eligible factor 'a' must be in range

			1/fMax <= |a| <= fMax

		where 'fMax > 1'.

	Instead of cumbersome logarithmic computations we choose a combination
	of absolute value, comparison with one and comparison of resulting
	factors:

		fA = |F|
		fQ = ( fA < 1 ) ? 1/fA : fA

	then require that

		fQ < fMax

	and compare 'fQ1' and 'fQ2' and choose the smaller one.

	Computation of fQ in the above scheme is done by an auxiliary inline
	function 'Q( double )' defined just after this comment.

PARAMETERS:
	See above.

RETURN VALUE:
	Number of rows that may be eliminated wih the candidate pivot row.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
double Q( double f )
{
	if( IsZero( f ) ) return -INFINITY;

	double fA = fabs( f );

	return ( fA < 1.0 ) ? ( 1.0 / fA ) : fA;
}


Int_T Presolver::IsNonZeroSubset( Int_T piv, Array<Real_T> &Factor, // )
	const Real_T fMax )
{
	assert( piv >= 0 && piv < m && !ExcludeRows[ piv ] &&
		rt[ piv ] == VT_FIXED );

	//--------------------------------------------------------------------------
	//	'Elim' is an 'M'-vector. Its entries are the numbers of non-zeros common
	//	to potentially eliminated rows and the candidate pivot row 'piv'.
	//
	//	'ElimMax' stores the desired current value of the 'Elim' vector entry.
	//	If a certain 'Elim[i]' entry is different than 'ElimMax', then row 'i'
	//	is rejected, as its non-zero pattern is not a superset of pivot row's
	//	pattern.
	//
	WorkVector<Int_T> Elim( m );

	//--------------------------------------------------------------------------
	//
	//	'ElimCnt' stores the number of rows likely to be eligible for
	//	elimination. When it reaches zero, the function returns 0.
	//
	//	NOTE: 'ElimCnt' is only an upper bound on the number of possible row
	//	eliminations, it is not an exact value because it ignores the numerical
	//	stability criterion.
	//
	Int_T ElimCnt = 0;

	Factor.Fill( 0.0, m );
	Elim.Fill( 0, m );

	//--------------------------------------------------------------------------
	//	Loop on candidate pivot row's non-zeros. This is the loop in
	//	which it is estabilished whether this is a good pivot row.
	//
	//	ANALYSIS OF THE FIRST COLUMN OF THE PIVOT ROW.
	//
	Ptr<Real_T> AR, AC;
	Ptr<Int_T> col, row;
	Int_T rlen, clen;

	//==========================================================================
	//	First column is located and analysed.
	//

	//	Find the first column, that has not been eliminated before.
	//	(Since we search only non-zero length rows there has to be one).
	//
	for( LP.GetRow( piv, AR, col, rlen ); ExcludeCols[*col] && rlen;
		++AR, ++col, --rlen )
		;
	assert( !ExcludeCols[*col] );

	// Scan the first column and fill the 'Elim' and 'Factor' vectors.
	//
	for( LP.GetColumn( *col, AC, row, clen ); clen; ++AC, ++row, --clen )
	{
		Real_T F;

		if( ExcludeRows[*row] || *row == piv ) continue;

		F = - *AC / *AR;

		if( Q( F ) < fMax ) 
		{
			Factor[ *row ] = F;
			Elim[ *row ]++;
		}
		ElimCnt++;
	}
	Int_T ElimMax = 1;

	//	If the first column was a singleton, we skip the row.
	//
	if( ElimCnt <= 0 ) return 0;

	//
	//	End of analysis of the first column.
	//==========================================================================


	//--------------------------------------------------------------------------
	//	Loop on the remaining columns.
	//
	for( ++AR, ++col, --rlen; rlen && ElimCnt > 0; ++AR, ++col, --rlen )
	{
		if( ExcludeCols[*col] || ColLen[*col] == 0 ) continue;

		//----------------------------------------------------------------------
		//	Abort the whole search if there is a singleton column in the row.
		//
		if( ColLen[*col] == 1 ) return 0;

		for( LP.GetColumn( *col, AC, row, clen ); clen; ++AC, ++row, --clen )
		{
			//	Skip pivot row and removed rows.
			//
			if( ExcludeRows[*row] || *row == piv )
				continue;

			//	If in this iteration we detect, that a row is not a superset of
			//	the non-zero pattern of the pivot row's pattern, we decrement
			//	the 'ElimCnt' counter.
			//
			if( Elim[*row] < ElimMax )
			{
				if( Elim[*row] > 0 )
				{
					Elim[*row] = 0;
					ElimCnt--;
				}
				Factor[*row] = 0.0;
				continue;
			}

			//	Increment row's coinciding non-zero counter.
			//
			Elim[*row]++;

			//	Compute the eliminatin factor for the current eliminated row
			//	and column combination. Also compute its numerical quality
			//	estimate 'fQ'.
			//
			Real_T F	= - *AC / *AR,
				fQ		= Q( F ),
				&fac	= Factor[ *row ];

			if( fQ > fMax ) continue;

			//	Compute numerical quality estimate of the previous value of the
			//	factor for this row 'facQ'.
			//

			//	If there was no acceptable old factor found, store the new
			//	factor. Otherwise remember the better (more stable) factor of
			//	the two).
			//
			if( fQ < Q( fac ) )
				fac = F;
		}

		//	'ElimMax' has to be equal (in the next iteration) the number of
		//	pivot row's non-zeros minus one.
		//
		ElimMax++;
	}
	// End of loop on potential pivot row's columns.
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	//	It is possible, that the search for possible eliminations will be ended
	//	with some (or all) 'Factor' entries having non-zero value, although the
	//	corresponding rows are not supersets of pivot row's non-zero pattern.
	//	In the following loop the 'Factor' vector is thoroughly examined in
	//	order to detect and correct such situations.
	//
	Int_T FactorCnt, i;
	for( ElimCnt = 0, FactorCnt = 0, i = 0; i < m; i++ )
		if( Elim[i] == ElimMax )
		{
			FactorCnt++;
			ElimCnt++;
		}
		else
			Factor[i] = 0.0;

	return FactorCnt;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::NumElimRows( Int_T piv, const Array<Real_T> &Factor )

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

Int_T Presolver::NumElimRows( Int_T piv, const Array<Real_T> &Factor )
{
	Ptr<Real_T> a, apiv;
	Ptr<Int_T> col;
	Int_T i, Trunc, len;
	Int_T ElimNZ = 0;

	//--------------------------------------------------------------------------
	//	Unpack pivot row 'piv' into a work array.
	//
	WorkVector<Int_T> PivRow( n );
	PivRow.Fill( -1, n );

	for( i = 0, LP.GetRow( piv, apiv, col, len ); len; ++col, ++i, --len )
		if( !ExcludeCols[ *col ] ) PivRow[ *col ] = i;

	//--------------------------------------------------------------------------
	//	Loop on entries of 'Factor' array. Perform numerical eliminations.
	//
	for( i = 0; i < m; i++ )
	{
		if( IsZero( Factor[i] ) ) continue;

		assert( i != piv );

		//------------------------------------------------------------------
		//	Update the right and left hand side vectors 'bl' and 'bu' of the
		//	eliminated row.
		//
		if( IsNonZero( bl[piv] ) )
		{
			Real_T shift = Factor[i] * bl[piv];

			if( ( rt[i] & VT_LO ) && IsZero( bl[i] += shift ) )
				bl[i] = 0.0;
			if( ( rt[i] & VT_UP ) && IsZero( bu[i] += shift ) )
				bu[i] = 0.0;
		}

		LP.GetRow( i, a, col, len );

		//------------------------------------------------------------------
		//	Loop on eliminated row entries.
		//
		for( Trunc = 0; len; ++a, ++col, --len )
		{
			//	If column already deleted or no corresponding non-zero in pivot
			//	row, skip current one.
			//
			if( ExcludeCols[ *col ] || PivRow[ *col ] < 0 ) continue;

			assert( *col >= 0 && *col < n );

			if( IsZero( *a += Factor[i] * apiv[ PivRow[ *col ] ] ) )
			{
				*a		= 0.0;
				*col	= n;
				Trunc++;
			}
		}
		//
		//	End of loop on eliminated row entries.
		//------------------------------------------------------------------

		if( Trunc )
		{
			RowList.TruncateLength( i, Trunc );
			ElimNZ += Trunc;
		}
	}
	//
	//	End of loop on 'Factor' entries (elimination factors).
	//----------------------------------------------------------------------

	return ElimNZ;
}


/*------------------------------------------------------------------------------

	Int_T Presolver::NumElimCols( Int_T piv, const Array<Real_T> &Factor,
		Array<Bool_T> &Scan, Int_T &ScanCnt )

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

Int_T Presolver::NumElimCols( Int_T piv, const Array<Real_T> &Factor, // )
	Array<Bool_T> &Scan, Int_T &ScanCnt )
{
	Ptr<Real_T> apiv, a;
	Ptr<Int_T> col, row;
	Int_T Trunc, rlen, clen;
	Int_T ElimNZ = 0;

	//--------------------------------------------------------------------------
	//	Loop on pivot row entries.
	//
	for( LP.GetRow( piv, apiv, col, rlen ); rlen; ++apiv, ++col, --rlen )
	{
		if( ExcludeCols[ *col ] ) continue;

		LP.GetColumn( *col, a, row, clen );

		//----------------------------------------------------------------------
		//	Elimination loop on column entries.
		//
		for( Trunc = 0; clen; ++a, ++row, --clen )
		{
			if( ExcludeRows[ *row ] || IsZero( Factor[ *row ] ) )
				continue;

			assert( *row >= 0 && *row < m );

			if( IsZero( *a += Factor[ *row ] * *apiv ) )
			{
				//	Mark (if necessary) appropriate 'Scan' entry, if the row
				//	needs to be scanned tested as a pivot row candidate.
				//
				if( !Scan[ *row ] &&
					rt[ *row ] == VT_FIXED && RowLen[ *row ] >= 2 )
				{
					Scan[ *row ] = True;
					ScanCnt++;
				}

				//	Remove the eliminated non-zero.
				//
				*a		= 0.0;
				*row	= m;
				Trunc++;
			}
		}
		//
		//	End of loop on column entries.
		//----------------------------------------------------------------------

		if( Trunc )
		{
			ElimNZ += Trunc;
			ColList.TruncateLength( *col, Trunc );
		}
	}
	//
	//	End of loop on pivot row entries.
	//--------------------------------------------------------------------------

	return ElimNZ;
}
