/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines for sparse matrices.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski - initial version,
					Artur Swietanowski - ammendments & all sparse versions of
					functions

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio

--------------------------------------------------------------------------------

SOURCE FILE NAME:	invsolve.cpp
CREATED:			1993.05.28
LAST MODIFIED:		1996.02.06

DEPENDENCIES:		smartptr.h, stdtype.h, std_tmpl.h, error.h, inverse.h,
					std_math.h, work_vec.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	The above mentioned functions perform solves with factorized basis matrix
and its transpose. The first word of the name indicates whether the input right
hand side of the equation is given as dense or as sparse vector. FTRAN
corresponds to 'forward transformation', that is a solve with basis. BTRAN
stands for 'backward transformation' and solves eqation with basis transpose.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	inverse::SparseFTRAN()
	inverse::DenseFTRAN()
	inverse::DenseBTRAN()
	inverse::SparseBTRAN()

STATIC FUNCTIONS:
	None.

STATIC DATA:
	None.

------------------------------------------------------------------------------*/


#include <string.h>

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
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

#ifndef __INVERSE_H__
#	include "inverse.h"
#endif
#ifndef __INVAUX_H__
#	include "invaux.h"
#endif


/*------------------------------------------------------------------------------

	void Inverse::SparseFTRAN( Array<Real_T> &b, Array<Int_T> &bInd,
		Int_T &bNz )

PURPOSE:
	Solves a matrix equation with basis matrix (stored as LU factors) and a
given right hand side passed to the function as a sparse vector.  Algorithm
used takes advantage of the right hand side's sparsity and skips - if possible
- a large part of computations. The idea was first introduced by Reid ("A
sparsity exploiting variant of LU factorization").

PARAMETERS:
	Array<Real_T> &b, Array<Int_T> &bInd, Int_T &bNz
		These three paramaters correspond to a single logical parameter: packed
		sparse right hand side vector (non-zeros, their row numbers and
		non-zeros' number). The same vector will be used to return the result.
		It is assumed that 'b' and 'bInd' are capable of holding at least 'n'
		non-zeros, where 'n' is the dimension of the basis matrix.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Intermediate result - after 'FTRANL' - is stored in a named sparse work
vector ('w', 'wInd', 'wNz' and 'wMark') and later used by the basis update
procedure.

------------------------------------------------------------------------------*/

void Inverse::SparseFTRAN( Array<Real_T> &b, Array<Int_T> &bInd, Int_T &bNz )
{
	SparseFTRAN_Cnt++;

	Int_T i, k;
	Int_T j;
	WorkVector<Int_T> bMark( n ),
		wMark( n, FTRANL_LABEL ),
		wInd( n );
	WorkVector<Real_T> w( n, FTRANL_LABEL );

	//--------------------------------------------------------------------------
	//	Erase the 'bMark' vector.
	//
	bMark.Fill( -1, n );

	//--------------------------------------------------------------------------
	//	Rewrite argument column into work data 'w', 'wInd', 'wMark' and 'wNz'.
	//
	w.Copy( 	b,		n, n, bNz );
	wInd.Copy(	bInd,	n, n, bNz );
	wMark.Fill( -1, n );
	wNz = bNz;
	for( j = 0; j < wNz; j++ ) wMark[ wInd[j] ] = j;
	bNz = 0;

	//--------------------------------------------------------------------------
	//	First compute 'w := M^(-1)b' (perform sparse FTRANL).
	//
	for( i = Int_T( ia - 1 ), k = M_len; k; k--, i-- )
	{
		Int_T ind1 = wMark[ jcol[i] ],
			ind2 = wMark[ irow[i] ];

		if( ind1 >= 0 && ind2 >= 0 )		// Element update.
		{
			w[ ind1 ] += a[i] * w[ ind2 ];
			if( IsZero( w[ ind1 ] ) )
			{
				assert( wNz > 0 );

				wMark[ wInd[ --wNz ] ]		= ind1;
				wInd[ ind1 ]				= wInd[ wNz ];
				w[ ind1 ]					= w[ wNz ];
				wMark[ jcol[i] ]	= -1;
			}
		}
		else if( ind2 >= 0 )				// Fill in.
		{
			Real_T x = a[i] * w[ ind2 ];
			if( IsNonZero( x ) )
			{
				wMark[ jcol[i] ]	= wNz;
				wInd[ wNz ]					= jcol[i];
				w[ wNz++ ]					= x;

				assert( wNz <= n );
			}
		}
	}
	//
	//	End of sparse FTRANL loop.
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	//	Sparse FTRANU.
	//
	//
	//	wMark codes:
	//		>=  0	- non-zero; no need to scan matrix row.
	//		   -1   - zero; matrix row scan unnecessary.
	//		   -2   - zero; martrix row scan necessary.
	//		<= -3   - non-zero; row scan necessary.
	//
	for( j = Int_T( n - 1 ); j >= 0; j-- )
	{
		Int_T r			= rlst[j],
			&wr			= wMark[r];

		if( wr == -1 ) continue;

		Real_T x = 0.0e0;
		Int_T rowStart	= rptr[r],
			rowEnd		= Int_T( rowStart + rlen[r] );

		if( wr >= 0 )
			x = w[ wr ];
		else if( wr <= -3 )
			x = w[ - wr - 3 ];
		else if( wr != -2 )
			abort();

		if( wr <= -2 )
		{
			Int_T bc;
			for( i = Int_T( rowStart + 1 ); i < rowEnd; i++ )
				if( ( bc = bMark[ jcol[i] ] ) >= 0 )
					x -= a[i] * b[ bc ];
			wr = Int_T( (wr == -2) ? -1 : - wr - 3 );
		}

		x /= a[ rowStart ];

		Int_T c	= jcol[ rowStart ],
			br	= bMark[c];

		//----------------------------------------------------------------------
		//	Write down / remove element.
		//
		if( IsNonZero( x ) )		// Non-zero result.
		{
			if( br >= 0 )			// Element update.
				b[ br ] = x;
			else					// New element creation.
			{
				assert( bNz < n );

				bMark[c]	= bNz;
				bInd[ bNz ]	= c;
				b[ bNz++ ]	= x;
			}
		}
		else if( br >= 0 )			// Non-zero cancellation.
		{
			assert( bNz > 0 );

			bMark[ bInd[ --bNz ] ]	= br;
			bInd[ br ]				= bInd[ bNz ];
			b[ br ]					= b[ bNz ];
			bMark[c]				= -1;
		}

		//----------------------------------------------------------------------
		//  Mark rows that will need to be scanned.
		//
		Int_T colEnd;
		for( i = cptr[c], colEnd = Int_T( i + clen[c] ); i < colEnd; i++ )
		{
			Int_T &wr2 = wMark[ irow[i] ];
				if( wr2 == wr && wr2 != -1 ) continue;

			if( wr2 == -1 )
				wr2 = -2;
			else if( wr2 >= 0 )
				wr2 = Int_T( -3 - wr2 );
			// Otherwise the row is already marked.
		}
	}
	//
	// End of FTRANU loop.
	//--------------------------------------------------------------------------

	for( j = 0; j < n; j++ )
	{
		Int_T &wr = wMark[j];

		if( wr == - 2 )
			wr = -1;
		else if( wr <= -3 )
			wr = Int_T( -wr - 3 );
	}

	w.Detach();
	wMark.Detach();
}


/*------------------------------------------------------------------------------

	void Inverse::DenseFTRAN( Array<Real_T> &b )

PURPOSE:
	Solution of a matrix equation with the basis matrix represented by the
current object of 'Inverse' class. 'b' is the right hand side vector. Solution
will overwrite the RHS vector.

PARAMETERS:
	Array<Real_T> &b
		The dense right hand side vector. Assumed to be allocated and filled
		properly with numerical data.

RETURN VALUE:
	None (data is written into the RHS vector 'b').

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::DenseFTRAN( Array<Real_T> &b )
{
	DenseFTRAN_Cnt++;

	WorkVector<Real_T> w( n );

	w.Copy( b, n, n, n );
	b.Fill( 0.0, n );

	//--------------------------------------------------------------------------
	//	First compute 'w := M^(-1)w' (perform sparse FTRANL).
	//
	for( Int_T i = Int_T( ia - 1 ), k = Int_T( ia - M_len ); i >= k; i-- )
		if( IsNonZero( w[ irow[i] ] ) )
			w[ jcol[i] ] += a[i] * w[ irow[i] ];

	//--------------------------------------------------------------------------
	//	Dense FTRANU.
	//
	for( Int_T ii = Int_T( n - 1 ); ii >= 0; ii-- )
	{
		Int_T l		= rlst[ii];
		Real_T aa	= w[l];

		Int_T j, k;
		for( j = Int_T(rlen[l]-1), k = Int_T(rptr[l]+1); j; j--, k++ )
			aa -= a[k] * b[ jcol[k] ];

		if( IsZero( aa ) ) continue;

		j		= jcol[k = rptr[l]];
		b[j]	= aa / a[k];
	}

	for( Int_T j = 0; j < n; j++ )
		if( IsZero( b[j] ) ) b[j] = 0.0;
}


/*------------------------------------------------------------------------------

	void Inverse::DenseBTRAN( Array<Real_T> &b )

PURPOSE:
	Name stands for "Dense Backward Transformation". Function performs a dense
solve of matrix equation with basis transpose. Basis is stored in form of LU
factors.

PARAMETERS:
	Array<Real_T> &b
		The dense right hand side vector. Assumed to be allocated and filled
		properly with numerical data.

RETURN VALUE:
	None (data is written into the RHS vector 'b').

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::DenseBTRAN( Array<Real_T> &b )
{
	DenseBTRAN_Cnt++;

	Int_T kp, k;
	Int_T i, j, ii;
	WorkVector<Real_T> w( n );

	w.Copy( b, n, n, n );
	b.Fill( 0.0, n );

	//--------------------------------------------------------------------------
	//	BTRANU
	//
	for( ii = 0; ii < n; ii++ )
	{
		i	= clst[ ii ];

		Real_T am = w[i];

		if( IsNonZero( am ) )
		{
			j		= rlst[ ii ];
			kp		= rptr[j];				// Pivotal row has pivot in front.
			am		= am / a[kp];
			b[j]	= am;
			for( kp++, k = Int_T( rlen[j] - 1 ); k; k--, kp++ )
				w[ jcol[kp] ] -= am * a[kp];
		}
	}

	//--------------------------------------------------------------------------
	//	BTRANL
	//
	for( k = Int_T( ia - M_len ); k < ia; k++ )
		b[ irow[k] ] += b[ jcol[k] ] * a[k];

	for( i = 0; i < n; i++ )
		if( IsZero( b[i] ) ) b[i] = 0.0e0;
}


/*------------------------------------------------------------------------------

	void Inverse::SparseBTRAN( Array<Real_T> &b, Array<Short_T> &mark )

PURPOSE:
	Name stands for "Sparse Backward Transformation". Function performs a semi-
sparse solve of a matrix equation with basis transpose. Basis is stored in form
of LU factors.

PARAMETERS:
	Array<Real_T> &b
		The dense right hand side vector. Assumed to be allocated and filled
		properly with numerical data.

	Array<Short_T> &mark
		The vector of non-zero marks. mark[i] = 0 if b[i] = 0.0 and
		mark[i] != 0 if b[i] != 0.0

RETURN VALUE:
	None (data is written into the RHS vector 'b').

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Inverse::SparseBTRAN( Array<Real_T> &b, Array<Int_T> &mark )
{
	SparseBTRAN_Cnt++;

	Int_T kp, k;
	Int_T i, j, ii;
	WorkVector<Real_T> w( n );
	Real_T am;

	w.Copy( b, n, n, n );
	b.Fill( 0.0, n );

	//--------------------------------------------------------------------------
	//	BTRANU
	//
	for( ii = 0; ii < n; ii++ )
		if( mark[ i = clst[ii] ] )
		{
			j		= rlst[ii];
			kp		= rptr[j];				// Pivotal row has pivot in front.
			am		= w[i];
			b[j]	= am /= a[kp];
			for( kp++, k = Int_T( rlen[j] - 1 ); k; k--, kp++ )
			{
				Int_T row = jcol[kp];

				if( IsZero( w[ row ] -= am * a[kp] ) )
				{
					w[ row ]	= 0.0;
					mark[ row ]	= 0;
				}
				else
					mark[ row ]	= 1;
			}
		}

	for( i = 0; i < n; i++ )
		if( IsZero( b[i] ) )
		{
			b[i]	= 0.0e0;
			mark[i]	= 0;
		}
		else
			mark[i]	= 1;

	//--------------------------------------------------------------------------
	//	BTRANL
	//
	for( k = Int_T( ia - M_len ); k < ia; k++ )
		if( mark[ i = jcol[k] ] )
		{
			j = irow[k];
			if( IsZero( b[j] += b[i] * a[k] ) )
			{
				b[j]	= 0.0;
				mark[j]	= 0;
			}
			else
				mark[j]	= 1;
		}
}
