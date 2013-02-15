/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines for sparse matrices.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	invaux.cpp
CREATED:			1991.12.29 by prof. A. Ruszczynski
LAST MODIFIED:		1995.09.12

DEPENDENCIES:		smartptr.h, stdtype.h, error.h, invaux.h
					<math.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	lupack()
	reorder()

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

------------------------------------------------------------------------------*/

#include <math.h>
#include <assert.h>

#ifndef __INVAUX_H__
#	include "invaux.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif


/*------------------------------------------------------------------------------

	void Pack( Array<Real_T> &a, Array<Int_T> &ind, Array<Int_T> &ptr,
		Int_T n, Array<Int_T> &len, Bool_T reals, Int_T &lfile );

PURPOSE:
	Compress file of positive integers. Entry 'j' starts at 'ind[ ptr[j] ]' and
contains 'len[j]' integers, 'j=0..n-1'. Other components of 'ind' are set to
'-1'. If 'reals != 0' array 'a' contains a real file associated with 'ind' and
this is compressed too.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

void Pack( Array<Real_T> &a, Array<Int_T> &ind, Array<Int_T> &ptr, // )
	Int_T n, Array<Int_T> &len, Bool_T reals, Int_T &lfile )
{
#ifndef NDEBUG
	{
		Int_T j, k;
		for( Int_T i = 0; i < n; i++ )
			for( j = ptr[i], k = Int_T( j + len[i] ); j < k; j++ )
				assert( ind[ j ] >= 0 && ind[ j ] < n );
	}
#endif

	//--------------------------------------------------------------------------
	//	Store the last element of entry 'j' in 'ptr[j]' then overwrite it by 'j'
	//
	Int_T knew, kend = lfile, j, k;
	for( j = 0; j < n; j++ )
		if( len[j] )
		{
			k	 				= Int_T( ptr[j] + len[j] - 1 );
			ptr[j]				= ind[ k ];
			ind[ k ]	= j;
		}
		else
			ptr[j]	= -1;

	//--------------------------------------------------------------------------
	//	Skip elements that need not be moved.
	//
	for( knew = 0 ; ind[ knew ] >= 0 && knew < kend; knew++ );

	//--------------------------------------------------------------------------
	//	Now shift entries to the left; 'knew' is the new position.
	//
	if( reals )
		for( k = knew; k < kend; k++ )
		{
			if( ind[ k ] >= 0 )
			{
				a[ knew ]		= a[ k ];
				ind[ knew++ ]	= ind[ k ];
			}
		}
	else
		for( k = knew; k < kend; k++ )
			if( ind[ k ] >= 0 )
				ind[ knew++ ]	= ind[ k ];

	lfile = knew--;

	//--------------------------------------------------------------------------
	//	Set new pointers by jumping over ends of entries. Note that number 'j'
	//	of entry is stored at the last position, so we can calculate where the
	//	end of previous entry is placed.
	//
	while( knew >= 0 )
	{
		j				= ind[ knew ];
		ind[ knew ]		= (Int_T)ptr[j];
		knew			-= len[j];
		ptr[j]			= Int_T( knew + 1 );
	}

	//--------------------------------------------------------------------------
	//	Now we need to set the pointers to zero length entries (rows or columns)
	//	to point to some resonable positions (up until now they were set to -1).
	//
	for( j = Int_T( n - 1 ); j >= 0; j-- )		// Scan the pointers backwards.
		if( ptr[j] < 0 || ptr[j] >= lfile )
		{
 			assert( ptr[j] == -1 && len[j] == 0 );
			ptr[j] = lfile;
		}
}



/*------------------------------------------------------------------------------

	void Reorder( Int_T n, Int_T U_len, Array<Real_T> &a, Array<Int_T> &jcol,
		Array<Int_T> &rptr, Array<Int_T> &irow)

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

void Reorder( Int_T n, Int_T U_len, Array<Real_T> &a, Array<Int_T> &jcol, // )
	Array<Int_T> &rptr, Array<Int_T> &irow )
{
	Int_T	i, ice, icep, j, jce, jcep;
	Real_T	ace,acep;
	Int_T	k, kr, loc;

	for ( j=0 ; j<n ; j++ ) rptr[j]=0;

	/* count the number of elements in each column */
	for( k = 0; k < U_len; k++ )
		rptr[ irow[ k ] ]++;

	/* set the rptr array */
	for ( k=0,j=0 ; j<n ; j++ ) {
		kr			= Int_T( k + rptr[j] );
		rptr[ j ]	= k;
		k			= kr;
	}

	/* reorder the elements into column order by in-place sort		*/

	for ( i=0; i<U_len ; i++ ) {
		/* establish the current entry */
		if ( (jce=irow[i]) < 0 ) continue;
		ace=a[i];
		ice=jcol[i];
		/* clear the location vacated */
		irow[i]=-1;
		/* chain from current entry to store items */
		do {
			/* compute new position for current entry */
			loc = rptr[jce]++;
			/* save contents of that location */
			acep=a[ loc ];
			icep=jcol[ loc ];
			jcep=irow[ loc ];
			/* put current entry */
			a[ loc ]=ace;
			jcol[ loc ]=ice;
			irow[ loc ]=-1;
			/* make saved entry current */
			ace=acep;
			ice=icep;
			jce=jcep;
		} while ( jcep >= 0 );
	}

	/* reset rptr vector */
	Int_T ja, jb;
	for( ja = 0, j = 0; j<n ; j++ )
	{
		jb		= rptr[j];
		rptr[j]	= ja;
		ja		= jb;
	}
}
