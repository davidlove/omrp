/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines for sparse matrices.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	inverse.cpp
CREATED:			1991.12.29 by prof. A. Ruszczynski
LAST MODIFIED:		1996.02.06

DEPENDENCIES:		smartptr.h, stdtype.h, std_tmpl.h, error.h, inverse.h
					std_math.h, vec_pool.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Inverse::Inverse()
	Inverse::~Inverse()
	Inverse::luclean()
	Inverse::luaddcol()

STATIC FUNCTIONS:
	None.

STATIC DATA:
	None.

------------------------------------------------------------------------------*/


#include <math.h>

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __VEC_POOL_H__
#	include "vec_pool.h"
#endif

#ifndef __INVERSE_H__
#	include "inverse.h"
#endif
#ifndef __INVAUX_H__
#	include "invaux.h"
#endif


/*------------------------------------------------------------------------------

	Inverse::Inverse( Int_T n )

PURPOSE:
	Basis inverse constructor - allocates and initializes necessary storage.

PARAMETERS:
	Int_T n
		Dimension of the square matrix that will be factorized.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Inverse::Inverse( Int_T _n )
	: n( _n ), ia( Int_T(5*n) ), a( ia ), irow( ia ), jcol( ia ),
	U_len( 0 ), M_len( 0 ), rowend( 0 ), colend( 0 ), cmprs( 0 ),
	rptr( n ), cptr( n ), rlen( n ), rlst( n ), clen( n ), clst( n ),
	stabrow( 0.1 ), stabglo( 1.0e-10 ), alert( 0 ),
	A_max( 0.0 ), U_max( 0.0 ), maxties( 20 ), wNz( 0 ),
	Diag( _n, 0.0 ), DiagUsed( _n, -1 )
{
	assert( _n > 0 );

	rlen.Fill( 0, n );
	clen.Fill( 0, n );
	ResetStatisticCounters();
}


/*------------------------------------------------------------------------------

	Inverse::~Inverse( void )

PURPOSE:
	Basis inverse destructor - destroys any work vectors that may be still
attached.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Inverse::~Inverse( void )
{
	WorkVectorPool::Destroy( FTRANL_LABEL );
}


/*------------------------------------------------------------------------------

	void Inverse::Resize( Int_T NewN )

PURPOSE:
	Changes the basis dimension and clears data structures.

PARAMETERS:
	Int_T NewN
		The new basis matrix dimension. Should be positive.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::Resize( Int_T NewN )
{
	assert( NewN > 0 );

	n = NewN;

	rptr.Resize( n );
	rlen.Resize( n );
	rlst.Resize( n );
	cptr.Resize( n );
	clen.Resize( n );
	clst.Resize( n );

	if( ia < 5 * n )
	{
		ia = Int_T( 5 * n );

		a.Resize( ia );
		irow.Resize( ia );
		jcol.Resize( ia );
	}

	Clear();
}


/*------------------------------------------------------------------------------

	void Inverse::Clear( void )

PURPOSE:
	Clears data structures.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::Clear( void )
{
	rlen.Fill( 0, n );
	clen.Fill( 0, n );
	U_len = M_len = 0;
	A_max = 0.0;
}


/*------------------------------------------------------------------------------

	void Inverse::AddCol( Int_T cnb, const Ptr<Real_T> &col,
		const Ptr<Int_T> &rnb )

PURPOSE:
	Add column 'cnb' to matrix

PARAMETERS:
	Int_T cnb
		Number of column of the square matrix.

	const Real_T *col, const Int_T *rnb, Int_T collen
		Column given as a sparse vector (non-zeros in 'col', their row numbers
		in 'rnb', their number in 'collen').

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::AddCol( Int_T cnb, const Ptr<Real_T> &col, // )
	const Ptr<Int_T> &rnb, Int_T len )
{
	const Int_T ul = U_len;

	assert( cnb >= 0 && cnb < n );
	assert( len > 0 && len <= n );

	while( U_len + len > ia )
	{
		ia += ia;
		a.Resize( ia );
		irow.Resize( ia );
		jcol.Resize( ia );
	}

	for( Int_T k = 0, l = len; k < l; k++ )
	{
		Real_T x = fabs( col[k] );

		if( IsZero( x ) ) continue;

		if( x > A_max ) A_max = x;
		a[ U_len ]		= col[k];
		irow[ U_len ]	= rnb[k];
		jcol[ U_len ]	= cnb;
		rlen[ rnb[k] ]++;
		U_len++;
	}

	clen[ cnb ] = Int_T( U_len - ul );

	U_max = A_max;
}


/*------------------------------------------------------------------------------

	void Inverse::ResetStatisticCounters( void )

PURPOSE:
	Resets the counters associated with the inverse representation.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::ResetStatisticCounters( void )
{
	RefactCnt = UpdateCnt = SparseFTRAN_Cnt = DenseFTRAN_Cnt =
		SparseBTRAN_Cnt = DenseBTRAN_Cnt = 0;
}


/*------------------------------------------------------------------------------

	Int_T Inverse::ReadStatisticCounter( int n )
		const

PURPOSE:
	Reads the value of one of the counters.

PARAMETERS:
	int n
		The number of the counter to read.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Inverse::ReadStatisticCounter( int i )
	const
{
	switch( i )
	{
	case Refact:		return RefactCnt;
	case Upd:			return UpdateCnt;
	case SpFTRAN:		return SparseFTRAN_Cnt;
	case DenFTRAN:		return DenseFTRAN_Cnt;
	case SpBTRAN:		return SparseBTRAN_Cnt;
	case DenBTRAN:		return DenseBTRAN_Cnt;
	default:			abort(); return 0;
	}
}


/*------------------------------------------------------------------------------

	void Inverse::RegisterDiagonal( const Array<Real_T> &d, Int_T len )

PURPOSE:
	Registers a diagonal matrix for possible basis recovery.

PARAMETERS:
	const Array<Real_T> &d, Int_T len
		A diagonal matrix and the length of the diagonal (supposed to be equal
		to 'n').

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::RegisterDiagonal( const Array<Real_T> &d, // )
#ifdef NDEBUG
	Int_T )
#else
	Int_T len )
#endif
{
	assert( len == n );
#ifndef NDEBUG
	for( Int_T ii = 0; ii < n; ii++ )
		assert( IsNonZero( d[ii] ) );
#endif

	DiagChanged = False;
	if( DiagPresent )
		for( Int_T i = 0; i < n; i++ )
			if( DiagUsed[i] >= 0 )
			{
				DiagChanged = True;
				break;
			}

	Diag.Copy( d, n, n );
	DiagPresent = True;
}


/*------------------------------------------------------------------------------

	void Inverse::RemoveDiagonal( void )

PURPOSE:
	Disables basis recovery.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::RemoveDiagonal( void )
{
	if( DiagPresent )
		DiagChanged = True;

	DiagPresent = False;
}


/*------------------------------------------------------------------------------

	Int_T Inverse::PositionInFile( Int_T item, Int_T data, Ptr<Int_T> ptr,
		Ptr<Int_T> len, Ptr<Int_T> arr )

PURPOSE:
	A primitive (not meant to be used directly) called by functions
"PositionInRowFile()" and "PositionInColumnFile()". Does their job on the file
defined by the three last arguments (see below).

PARAMETERS:
	Int_T item, Int_T data
		'data' to be found in the item number 'item' in the file.
	
	Ptr<Int_T> ptr, Ptr<Int_T> len, Ptr<Int_T> arr
		The file. Item to be searched for the first occurrence of 'data' starts
		at position 'ptr[item]' in array 'arr' and is 'len[item]' long.

RETURN VALUE:
	The position (in range from 'ptr[item]' to 'ptr[item]' + 'len[item]') at
which 'data' was found. If it is not found -1 is returned.

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Int_T Inverse::PositionInFile( Int_T item, Int_T data, Ptr<Int_T> ptr, // )
	Ptr<Int_T> len, Ptr<Int_T> arr )
{
	assert( item >= 0 && item < n );
	assert( data >= 0 && data < n );

	Int_T start	= ptr[item],
		end		= Int_T( start + len[item] );

	assert( start >= 0 );
	assert( start <= Max( rowend, colend ) );
	assert( end >= start );
	assert( end <= Max( rowend, colend ) );

	arr += start;
	for( Int_T i = start; i < end; ++i, ++arr )
	{
		if( *arr == data ) return i;

		assert( *arr >= 0 && *arr < n );
	}

	return -1;
}
