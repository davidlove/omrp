/* ======================================================================= */
/*                           FILE Q6                                       */
/*     AUXILIARY SUBROUTINES FOR PERFORMING SIMPLE ARRAY OPERATIONS        */
/* ======================================================================= */
/*																		   */
/*	q6.f -- translated by f2c (version 19940705.1).
	You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
/* ----------------------------------------------------------------------- */
/*	Modified by Artur Swietanowski.										   */
/* ----------------------------------------------------------------------- */



#include <assert.h>

#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif



Real_T dnorm2_( Int_T n, Real_T *a )
{
	assert( n >= 0 );

	Real_T ret_val;
	for( ret_val = 0.0; n; --n, ++a )
		ret_val += *a * *a;
	return ret_val;
}


Real_T ddot_( Int_T n, Real_T *a, Real_T *b )
{
	assert( n >= 0 );

	Real_T ret_val;
	for( ret_val = 0.0; n; --n, ++a, ++b )
		ret_val += *a * *b;
	return ret_val;
}


Real_T ddots_( Int_T n, Real_T *a, Int_T *irn, Real_T *b )
{
	assert( n >= 0 );

	--b;

	Real_T ret_val;
	for( ret_val = 0.0; n; --n, ++a, ++irn )
	{
		assert( *irn > 0 );

		ret_val += *a * b[*irn];
	}
	return ret_val;
}


void dzero_( Int_T n, Real_T *a )
{
	assert( n >= 0 );

	for( ; n; --n, ++a )
		*a = 0.0;
}


void dmult_( Int_T n, Real_T *a, Real_T *t )
{
	assert( n >= 0 );

	for( ; n; --n, ++a )
		*a *= *t;
}


void dcopy_( Int_T n, Real_T *a, Real_T *b )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++b )
		*b = *a;
}


void dpack_( Int_T n, Real_T *a, Real_T *b, Int_T *irn )
{
	assert( n >= 0 );

	for( --a; n; --n, ++b, ++irn )
	{
		assert( *irn > 0 );

		*b = a[*irn];
	}
}


void dunpk_( Int_T n, Real_T *a, Int_T *irn, Real_T *b )
{
	assert( n >= 0 );

	for( --b; n; --n, ++a, ++irn)
	{
		assert( *irn > 0 );

		b[*irn] = *a;
	}
}


void dunne_( Int_T n, Real_T *a, Int_T *irn, Real_T *b )
{
	assert( n >= 0 );

	for( --b; n; --n, ++a, ++irn )
	{
		assert( *irn > 0 );

		b[*irn] = - *a;
	}
}


void dcopne_( Int_T n, Real_T *a, Real_T *b )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++b )
		*b = - *a;
}


void icopy_( Int_T n, Int_T *a, Int_T *b )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++b )
		*b = *a;
}


void dsuma_( Int_T n, Real_T *a, Real_T *b )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++b )
		*a += *b;
}


void ddifr_( Int_T n, Real_T *a, Real_T *b )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++b )
		*a -= *b;
}


void dstep_( Int_T n, Real_T *a, Real_T *d, Real_T tau )
{
	assert( n >= 0 );

	for( ; n; --n, ++a, ++d )
		*a += *d * tau;
}


void dstpck_( Int_T n, Real_T *a, Int_T *irn, Real_T *d, Real_T tau )
{
	assert( n >= 0 );

	for( --d; n; --n, ++a, ++irn )
	{
		assert( *irn > 0 );

		*a += d[*irn] * tau;
	}
}


void dstunp_( Int_T n, Real_T *a, Real_T *d, Int_T *irn, Real_T tau )
{
	assert( n >= 0 );

	for( --a; n; --n, ++d, ++irn )
	{
		assert( *irn > 0 );

		a[*irn] += *d * tau;
	}
}
