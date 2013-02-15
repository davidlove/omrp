/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines for sparse matrices.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	invupd.cpp
CREATED:			1991.12.29 by prof. A. Ruszczynski
LAST MODIFIED:		1996.02.06

DEPENDENCIES:		smartptr.h, stdtype.h, std_tmpl.h, error.h, inverse.h,
					invaux.h, std_math.h, work_vec.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Inverse::luupdate()

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

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

	Short_T Inverse::Update( Int_T mm )

PURPOSE:
	Update factorization of a square basis matrix ('a') of dimension 'n' after
exchange of column 'mm'.

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	Return codes:
 		1	:	success,
		-2	:	basis is singular.

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Short_T Inverse::Update( Int_T mm )
{
	UpdateCnt++;

	Short_T code	= 1;
	Real_T M_entry, U_entry;
	Int_T jm	= mm,
		mcp		= 20,	// Limits the number of compresses before reallocation.
		winbeg, winend, ij, jns, ins, ipiv, jpiv, nz,
		i = 0, j, l, ir, ii, im, m = 0, found, bmpend = 0;
	size_t indent, indent1;
	Int_T kp, kl, kr, kj, kk, k1, k2, krl, knp = 0, kq, kpl, ks, k;
		
	WorkVector<Real_T> w( n, FTRANL_LABEL );
	WorkVector<Int_T> mark( n ), wMark( n, FTRANL_LABEL );
	WorkVector<Int_T> adr( n );

	//--------------------------------------------------------------------------
	//	Remove column 'jm'.
	//
	U_len	-= clen[ jm ];
	im		= irow[ cptr[ jm ] ];
	kp		= cptr[ jm ];
	kl		= Int_T( kp + clen[ jm ] );
	for ( k = kp; k < kl; k++ )
	{
		i					= irow[ k ];
		irow[ k ]	= -1;
		for( kk = rptr[i]; jcol[ kk ] != jm; kk++ )
			assert( kk < rptr[i] + rlen[i] );

		krl					= Int_T( rptr[i] + --rlen[i] );
		a[ kk ]		= a[ krl ];
		jcol[ kk ]	= jcol[ krl ];
		jcol[ krl ]	= -1;
	}
	clen[ jm ] = 0;

	//----------------------------------------------------------------------
	//	A new column will be added at the end of column file. If necessary
	//	pack column file or reallocate it to increase its length. It is the only
	//	place where column file may be reallocated. Therefore, we need to take
	//	into account all possible situations. In particular we presume that a
	//	column of the maximum length ('n') will be moved (in the file) at some
	//	later time. Thus the second condition and the new value of 'ia' (see
	//	below).
	//
	if( colend + M_len + wNz > ia || U_len + M_len + wNz + n > ia )
	{
		while( cmprs > mcp || U_len + U_len/5 + M_len + wNz + n > ia )
		{
			k1 = Int_T( ia - 1 );
			ia = Max( Int_T( 3*ia/2 ),
				Int_T( U_len + U_len/5 + M_len + wNz + n ) );
			k2 = Int_T( ia - 1 );
			cmprs = 0;
			a.Resize(  ia );
			irow.Resize(  ia );
			jcol.Resize(  ia );
			for( kk = 0; kk < M_len; kk++, k1--, k2-- )
			{
				a[ k2 ]		= a[ k1 ];
				irow[ k2 ]	= irow[ k1 ];
				jcol[ k2 ]	= jcol[ k1 ];
			}
		}
		Pack( a, irow, cptr, n, clen, False, colend );
		cmprs++;
	}
	assert( colend + M_len + wNz <= ia && U_len + M_len + wNz + n <= ia );

	//--------------------------------------------------------------------------
	//	Insert new column.
	//
	cptr[ jm ] = colend;				// New column start (length already 0 ).
	for( ii = 0; ii < n; ii++ )
	{
		if( ( i = rlst[ ii ] ) == im ) m = ii;
		if( wMark[i] == -1 || IsZero( w[ wMark[i] ] ) ) continue;
		U_len++;
		bmpend = ii;

		irow[ colend++ ] = i;	// Store non-zero's row number.
		clen[ jm ]++;					// Update column's length.

		//----------------------------------------------------------------------
		//	Add the new element to row file. If possible add it at end of
		//	an existing row. Otherwise pack or reallocate the row file.
		//
		nz	= rlen[i];
		kpl	= Int_T( rptr[i] + nz );
		if( kpl < rowend)
		{
			if( jcol[ kpl ] == -1 ) goto AddAtEnd;
		}
		else if( kpl == rowend && rowend + M_len + 1 < ia )
		{
			rowend++;
			goto AddAtEnd;
		}

		//----------------------------------------------------------------------
		//	Check available space to see if reallocation is required.
		//
		if( rowend + M_len + 2*nz + 1 > ia || U_len + M_len + nz + n + 1  > ia )
		{
			Int_T minLen = Int_T( U_len + U_len/5 + M_len + nz + n + 1 );
			while( cmprs > mcp ||  minLen > ia )
			{
				k1 = Int_T( ia - 1 );
				ia = Max( Int_T( 3*ia/2 ), Int_T( minLen ) );
				k2 = Int_T( ia - 1 );
				cmprs = 0;
				a.Resize(  ia );
				irow.Resize(  ia );
				jcol.Resize(  ia );
				for( kk = 0; kk < M_len; kk++, k1--, k2-- )
				{
					a[ k2 ]		= a[ k1 ];
					irow[ k2 ]	= irow[ k1 ];
					jcol[ k2 ]	= jcol[ k1 ];
				}
			}
			Pack( a, jcol, rptr, n, rlen, True, rowend);
			cmprs++;
		}
		//----------------------------------------------------------------------
		//	Now the row has been moved to the end of row file. Therefore we
		//	don't need to take into account the growing of the space it will
		//	leave. Thus we replaced '2 * nz' from the previous condition with
		//	'nz'.
		//
		assert( rowend + M_len + nz + 1 <= ia );

		//----------------------------------------------------------------------
		//	Move the whole row, or just remember its previous beginning and
		//	store the new one.
		//
		kp		= rptr[i];
		rptr[i]	= rowend;

		if( nz )
			for( k = kp, kpl = Int_T( kp + nz ); k < kpl; k++, rowend++ )
			{
				a[ rowend ]		= a[ k ];
				jcol[ rowend ]	= jcol[ k ];
				jcol[ k ]		= -1;
			}
 		kpl = rowend++;

		//----------------------------------------------------------------------
		//	Finally write down the element in the row file.
AddAtEnd:
 		rlen[i]++;
		a[ kpl ]	= w[ wMark[i] ];
 		jcol[ kpl ]	= jm;
	}
	//
	//	End of column insertion loop.
	//--------------------------------------------------------------------------


	//--------------------------------------------------------------------------
	//	First singularity check (based on non-zero structure only).
	//
	if( rlen[ im ] == 0 || clen[ jm ] == 0 || m > bmpend )
		goto singular;

	//--------------------------------------------------------------------------
	//	Find column singletons. Non-singletons (and the spike) are marked with
	//	'mark[j] = 1'.
	//	Only 'rlst[]' is revised, 'clst[]' is used for workspace.
	//
	for( j = 0; j < n; j++ ) mark[j] = 0;
	ins	= winbeg = m;
	mark[ jm ] = 1;
	for( ii = m; ii <= bmpend; ii++ )
	{
		i = rlst[ ii ];
		if( mark[ clst[ ii ] ] )	//	Non-singleton is stored in 'clst'.
 		{
			indent = rptr[i];
			for( k = rlen[i]; k--; )
				mark[ jcol[ indent++ ] ] = 1;
			clst[ ins++ ] = i;
		} else						//	Singleton can be placed in 'rlst'.
			rlst[ winbeg++ ] = i;
	}

	//--------------------------------------------------------------------------
 	//	Place non-singletons after singletons.
 	//	Place spike at end.
	//
	indent =  m + 1;
	for( ii = winbeg; ii < bmpend; ii++ )
		rlst[ ii ] = clst[ indent++ ];
	rlst[ bmpend ] = im;

	//--------------------------------------------------------------------------
 	//	Find row singletons. Non-singletons and spike row are marked with
	//	'mark[i] = 2'.
	//
	winend = jns = bmpend;
	mark[ im ] = 2;
	j = jm;
	for( ii = bmpend; ii >= winbeg; ii-- )
	{
		i=rlst[ii];
		if( mark[i] == 2 )
		{
 			if( ii != bmpend )
 				j = jcol[ rptr[i] ];
 			indent = cptr[j];
			for( k = clen[j]; k--; )
				mark[ irow[ indent++ ] ] = 2;
			clst[ jns-- ] = i;
		}
		else
			rlst[ winend-- ] = i;
	}

	indent = jns + 1;
	for( ii = winbeg; ii <= winend; ii++, indent++)
	{
		mark[ clst[ indent ] ]	= 3;
		rlst[ ii ]				= clst[ indent ];
	}

	//--------------------------------------------------------------------------
 	//	Deal with singleton spike column. 
	//	NOTE that bump rows are marked by mark[i]=3.
	//
	for( ii = winbeg; ii <= winend; ii++ )
	{
		kp = cptr[ jm ];
		kl = Int_T( kp + clen[ jm ] );
		found = 0;
		for( k = kp; k < kl; k++ )
		{
			l = irow[ k ];
			if( mark[l]==3 )
				if( found )
				{
 					indent =  winbeg;
 					for( ij = ii; ij <= winend; ij++ )
 						clst[ ij ] = rlst[ indent++ ];
					goto test;
 				}
 				else
 				{
					i		= l;
					knp		= k;
					found	= 1;
 				}
		}
		if( !found ) goto singular;
		
		//----------------------------------------------------------------------
 		//	Make (i,jm) a pivot.
 		//
		irow[ knp ]	= irow[ kp ];
		irow[ kp ]	= i;
		kp					= rptr[i];
		for( k = kp; jcol[ k ] != jm; k++ )
			assert( k < kp + rlen[i] );

		M_entry				= a[ kp ];
		a[ kp ]		= a[ k ];
 		a[ k ]		= M_entry;
 		jcol[ k ]	= jcol[ kp ];
		jcol[ kp ]	= jm;
 		jm					= jcol[ k ];
		clst[ ii ]			= i;
		mark[i]				= 2;
	}
	ii = winend;

test:
	//--------------------------------------------------------------------------
	//	'winbeg == winend' means that triangularity has been restored.
	//
	if( winbeg == winend ) goto restore;

	for( i = winbeg; i < winend; i++ ) rlst[i] = clst[i];
	winbeg = ii;
	if( winbeg == winend ) goto restore;

	//--------------------------------------------------------------------------
 	//	'adr' will store incremented indices to pivot row entries. It is cleared
	//	now.
	//
	adr.Fill( 0, n );

	//==========================================================================
	//
	//							E L I M I N A T I O N
	//
	//==========================================================================

 	ir = rlst[ winend  ];
	for ( ii = winbeg; ii <= winend; ii++ )
	{
		ipiv	= rlst[ ii ];
		kp		= rptr[ ipiv ];
		kr		= rptr[ ir ];
		jpiv	= jcol[ kp ];
		if( ii == winend )
			jpiv = jm;

		//----------------------------------------------------------------------
 		//	Search non-pivot row for element to be eliminated (the one in the
 		//	pivot column).
 		//
		krl	= Int_T( kr + rlen[ ir ] );
		knp	= -1;
		for( kk = kr; kk < krl; kk++ )
			if( jcol[ kk ] == jpiv )
			{
				knp = kk;
				break;
			}

		if( knp == -1 )
			if ( ii==winend )
				goto singular;
			else
				continue;

		//----------------------------------------------------------------------
 		//	Bring element to be eliminated to front of its row.
 		//
		M_entry				= a[ knp ];
		a[ knp ]	= a[ kr ];
		a[ kr ]		= M_entry;
		jcol[ knp ]	= jcol[ kr ];
		jcol[ kr ]	= jpiv;
		
		//----------------------------------------------------------------------
		//	Perform pairwise pivoting:
		//		Choose a pivot that will guarantee a multiplier not greater,
		//		than '1/stabrow'. If both the pivot ('ipiv','jpiv') and element
		//		to be eliminated ('ir','jpiv') satisfy this condition, then
		//		choose on sparsity grounds (take pivot with shorter row).
		//		If both rows have the same length, choose larger of the two
		//		pivot candidates.
		//
		if( fabs( M_entry ) > stabrow * fabs( a[ kp ] ) )
			if( ii == winend ||
				fabs( a[ kp ] ) < stabrow * fabs( M_entry ) ||
				rlen[ ir ] < rlen[ ipiv ] ||
				( fabs( M_entry ) > fabs( a[ kp ] ) &&
				rlen[ ir ] == rlen[ ipiv ] ) )
			{
				//--------------------------------------------------------------
				//	Interchange rows 'ir' and 'ipiv'.
				//
				rlst[ winend ]	= ipiv;
				rlst[ ii ]		= ir;
				ir				= ipiv;
				ipiv			= rlst[ ii ];
				k				= kr;
				kr				= kp;
				kp				= k;
 				kj				= cptr[ jpiv ];
				for( k = kj; irow[ k ] != ipiv; k++ )
					assert( k >= 0 && k < kj + clen[ jpiv ] );

				irow[ k ]		= irow[ kj ];
				irow[ kj ]		= ipiv;
 			}

		//----------------------------------------------------------------------
		//	If the resulting pivot is too small, exit and report singularity.
		//
		if( IsZero( a[ kp ] ) ) goto singular;
		
		if( ii == winend ) break;
		M_entry = -a[ kr ] / a[ kp ];

		//----------------------------------------------------------------------
 		//	Compress row file to make room for new row and 'M_entry'.
 		//
		if( rowend + rlen[ ir ] + rlen[ ipiv ] + M_len + 1 > ia )
		{
			while( cmprs > mcp ||
				U_len + U_len/5 + rlen[ ir ] + rlen[ ipiv ] + M_len + 1 > ia )
			{
				cmprs	= 0;
				k1		= Int_T( ia - 1 );
				ia		+= ia/2;              
				k2		= Int_T( ia - 1 );
				a.Resize(  ia );
				irow.Resize(  ia );
				jcol.Resize(  ia );
				for( kk = 0; kk < M_len; kk++, k1--, k2-- )
				{
					a[ k2 ]		= a[ k1 ];
					irow[ k2 ]	= irow[ k1 ];
					jcol[ k2 ]	= jcol[ k1 ];
				}
			}
			Pack( a, jcol, rptr, n, rlen, True, rowend );
			cmprs++;
			kp = rptr[ ipiv ];
			kr = rptr[ ir ];
		}
		assert( rowend + rlen[ ir ] + rlen[ ipiv ] + M_len + 1 <= ia );

		krl	= Int_T( kr + rlen[ ir ] );

		kq	= Int_T( kp + 1 );
		kpl	= Int_T( kp + rlen[ ipiv ] );
		
		//----------------------------------------------------------------------
		//	Place pivot row pattern (excluding pivot) in 'adr'.
		//
		for( k = kq; k < kpl; k++ ) adr[ jcol[ k ] ] = k;

		//----------------------------------------------------------------------
		//	Scan modified row and update non-zeros. Transfer the non-zeros to
		//	the end of row file.
		//
		jcol[ kr++ ]	= -1;
		rptr[ ir ]				= rowend;
	 	for( ks = kr; ks < krl; ks++ )
	 	{
			j			= jcol[ ks ];

			jcol[ ks ]	= -1;
			U_entry				= a[ ks ];

			if( adr[j] )
			{
				U_entry += M_entry * a[ adr[j] ];
				adr[j] = 0;
			}

			if( IsNonZero( U_entry ) )	// Update element.
			{	
				U_max						= Max( U_max, fabs( U_entry ) );
				a[ rowend ]			= U_entry;
				jcol[ rowend++ ]	= j;
			}
			else						// Remove element from column clique.
			{	
				U_len--;
				clen[j]--;
				kl = Int_T( cptr[j] + clen[j] );
				for( kk = cptr[j]; irow[ kk ] != ir; kk++ )
					assert( kk >= 0 && kk < kl );

				irow[ kk ] = irow[ kl ];
				irow[ kl ] = -1;
			}
		}

		//----------------------------------------------------------------------
 		//	Scan pivot row for fills. Also put them at the end of row file.
		//
		for( ks = kq; ks < kpl; ks++ )
		{
			j = jcol[ ks ];
			if( !adr[j] ) continue;		// Skip elements that will not need to
										// be updated.
			
			U_entry					= M_entry * a[ adr[j] ];
			a[ rowend ]		= U_entry;
			jcol[ rowend ]	= j;
			rowend++;
			U_len++;
			adr[j]			= 0;

			//------------------------------------------------------------------
			//	Add new element's row number 'ir' to column clique 'j'.
			//
			nz	= clen[j];
			k	= cptr[j];
			kl	= Int_T( k + nz );

			//------------------------------------------------------------------
			//	If possible place new element at the end of present entry.
			//
			if( kl == colend )
			{
				if( colend + M_len + 1 > ia ) goto compress;
				colend++;
			} 
			else if( irow[ kl ] >= 0 ) goto compress;
			irow[ kl ] = ir;
			goto noneed;

			/* new entry has to be created. */
compress:
			//------------------------------------------------------------------
			//	Compress column file if there is no room for new entry.
			//
			if( colend + M_len + nz + 1 > ia )
			{
  				Pack( a, irow, cptr, n, clen, False, colend );
				cmprs++;
				k	= cptr[j];
				kl	= Int_T( k + nz );
			}
			assert( colend + M_len + nz + 1 <= ia);

			//------------------------------------------------------------------
			//	Transfer old column to its new position - the end of column
			//	file.
			//
			indent	= cptr[j];
			indent1	= colend;
			cptr[j]	= colend;
			kk		= nz;
			while( kk-- )
			{
				irow[ indent1++ ]	= irow[ indent ];
				irow[ indent++ ]	= -1;
			}
			colend	+= nz;
			
			//------------------------------------------------------------------
			//	Add new element.
			//
			irow[ colend++ ]=ir;

noneed:
 			U_max	= Max( fabs( U_entry ), U_max );
			clen[j]	= Int_T( nz + 1 );
		}
		//
 		//	End of fill-in loop.
 		//----------------------------------------------------------------------
 		
 		//----------------------------------------------------------------------
 		//	As the eliminated row is now at the end of the row file, we can
 		//	calculate its length using 'rowend' and 'rptr[ir]'. If the length
 		//	is equal to zero, the matrix is singular.
 		//	(if the length is subzero, we're in trouble)
 		//
		rlen[ ir ] = Int_T( rowend - rptr[ ir ] );
		if( rlen[ ir ] == 0 )
			goto singular;
		else
			assert( rlen[ir] >= 0 );

		//----------------------------------------------------------------------
		//	Store multiplier. Make an entry for multiplier 'M_entry' at the
		//	other end of 'a' and 'ind'.
		//	Compress column file if necessary.
		//
		if( M_len + colend + 1 > ia )
		{
			Pack( a, irow, cptr, n, clen, False, colend );
			cmprs++;
		}
		assert( M_len + colend + 1 <= ia);

		M_len++;
		k					= Int_T( ia - M_len );
		a[ k ]		= M_entry;
		irow[ k ]	= ipiv;
		jcol[ k ]	= ir;
		
		//----------------------------------------------------------------------
		//	Create blank in pivotal column.
		kp			= cptr[jpiv];
		clen[ jpiv ]--;
		kl			= Int_T( kp + clen[jpiv] );
		for( k = kp; irow[ k ] !=ir; k++ )
			assert( k >= 0 && k < kl );

		irow[ k ]	= irow[ kl ];
		irow[ kl ]	= -1;
		U_len--;
	}

	//==========================================================================
	//
	//	Trangularity restored.
	//
	//==========================================================================
restore:
	//--------------------------------------------------------------------------
	//	Construct column permutation and store it in clst[.]
	//
	for( ii=m; ii<=bmpend; ii++ )
	{
		ir			= rlst[ ii ];
		clst[ ii ]	= jcol[ rptr[ ir ] ];
	}

	return code;


	//==========================================================================
	//
	//	Singular matrix detected in the update procedure.
	//
	//==========================================================================
singular:
	code		= -2;
	return code;
}
