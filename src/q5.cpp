/*------------------------------------------------------------------------------
                          FILE Q5
			Routines for advanced matrix operations

     Q5ADRW - PLANE REFLECTIONS TO ADD A ROW,
     Q5DLCO - PLANE REFLECTIONS TO DELETE A COLUMN,
     Q5DLRW - PLANE REFLECTIONS TO DELETE A ROW, 
     Q5GVNS - TO COMPUTE THE GIVENS REFLECTOR, 
     Q5REFL - TO APPLY THE REFLECTOR TO TWO VECTORS, 
     Q5REXP - TO ADD A COLUMN TO A TRIANGULAR MATRIX STORED BY ROWS, 
     Q5RPCK - TO DELETE A COLUMN AND THE FOLLOWING DIAGONAL ELEMENTS 
              FROM A MATRIX STORED BY ROWS, 
     Q5TRIS - TO SOLVE A TRIANGULAR SYSTEM OF EQUATIONS. 
--------------------------------------------------------------------------------
     WRITTEN BY A. RUSZCZYNSKI , INSTITUT FUER OPERATIONS RESEARCH, 
     UNIVERSITAET ZUERICH, FEBRUARY 1985. 
     DATE LAST MODIFIED: OCTOBER 1985. 
------------------------------------------------------------------------------*/
/*   q5.f -- translated by f2c (version 19940705.1).
     You must link the resulting object file with the libraries:
     -lf2c -lm   (in that order)
--------------------------------------------------------------------------------  
     Modified by Artur Swietanowski.											*/


#include <math.h>
#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif


/* Table of constant values */

static Int_T q5_c_1 = 1;




int q5adrw_( Int_T *m, Real_T *r, Real_T *z, Real_T *w, Real_T *pi, // )
	Real_T *zlast, Real_T *wlast )
{
	/* Local variables */
	Real_T c;
	Int_T i;
	Real_T p, s;
	Int_T ir, len;

	/*--------------------------------------------------------------------------
	PURPOSE:
	    APPLIES REFLECTIONS TO ROWS OF R AND TO PI TO ANNULATE PI
	    WHILE PRESERVING UPPER TRIANGULARITY OF R.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--pi;
	--w;
	--z;
	--r;

	/* Function Body */
	if (*m < 1)
		return 0;

	ir = 1;
	len = *m;
	for (i = 1; i <= *m; ++i)
	{
		q5gvns_(&r[ir], &pi[i], &c, &s, &p);
		++ir;
		--len;
		if (s != 0.)
		{
			q5refl_(&len, &r[ir], &pi[i + 1], &c, &s, &p);
			q5refl_(&q5_c_1, &z[i], zlast, &c, &s, &p);
			q5refl_(&q5_c_1, &w[i], wlast, &c, &s, &p);
		}
		ir += len;
	}
	return 0;
}
/* -----END OF Q5ADRW----------------------------------------------------- */


int q5dlco_( Int_T *m, Int_T *idel, Real_T *r, Real_T *z, Real_T *w )
{
	/* Local variables */
	Int_T imax;
	Real_T c;
	Int_T i;
	Real_T p, s;
	Int_T ix, iy, lenrow;

	/*--------------------------------------------------------------------------
	SUBROUTINES CALLED
	    Q5GVNS,Q5RPCK,Q5REFL.

	PURPOSE:
	    DELETES COLUMN IDEL FROM A QR FACTORIZATION BY DELETING
	    COLUMN IDEL FROM R AND RESTORING UPPER TRIANGULAR FORM OF R
	    BY GIVENS REFLECTIONS (ROW OPERATIONS). REFLECTIONS ARE APPLIED
	    ALSO TO VECTORS Z AND W AND TOP PARTS OF Z AND W
	    STILL SATISFY EQUATIONS Z=Q'C AND R'W=B (B,C - FIXED ).
	    NOTE THAT R IS STORED BY ROWS AND Q IS NOT STORED AT ALL.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--w;
	--z;
	--r;

	/* Function Body */
	if (*m < 1)
		return 0;

	if (*idel != *m)
	{
		iy = Int_T( (*idel - 1) * ((*m << 1) - *idel + 2) / 2 + 2 );
		lenrow = Int_T( *m - *idel );
		imax = Int_T( *m - 1 );
		for (i = *idel; i <= imax; ++i)
		{
			ix = iy;
			iy += lenrow;
			--lenrow;
			q5gvns_(&r[ix], &r[iy], &c, &s, &p);
			++iy;
			if (s != 0.)
			{
				q5refl_(&lenrow, &r[ix + 1], &r[iy], &c, &s, &p);
				q5refl_(&q5_c_1, &z[i], &z[i + 1], &c, &s, &p);
				q5refl_(&q5_c_1, &w[i], &w[i + 1], &c, &s, &p);
			}
		}
	}

	q5rpck_(m, &r[1], idel);
	return 0;
}
/* ----END OF Q5DLCO----------------------------------------------------- */


int q5dlrw_( Int_T *m, Real_T *r, Real_T *z, Real_T *w, Real_T *pi, // )
	Real_T *col, Real_T *rho, Real_T *zlast, Real_T *wlast )
{
	/* Local variables */
	Real_T c;
	Int_T i;
	Real_T p, s;
	Int_T ir, len;

	/*--------------------------------------------------------------------------
	PURPOSE:
	    APPLIES REFLECTIONS TO COL AND RHO TO MAKE COL ZERO
	    AND RHO ONE WHILE PRESERVING UPPER TRIANGULARITY OF R.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--col;
	--pi;
	--w;
	--z;
	--r;

	/* Function Body */
	if (*m < 1)
		return 0;

	dzero_(*m, &pi[1]);
	i = Int_T( *m + 1 );
	ir = Int_T( *m * i / 2 + 1 );
	for (len = 1; len <= *m; ++len)
	{
		ir -= len;
		--i;
		q5gvns_(rho, &col[i], &c, &s, &p);
		if (s != 0.)
		{
			q5refl_(&len, &pi[i], &r[ir], &c, &s, &p);
			q5refl_(&q5_c_1, zlast, &z[i], &c, &s, &p);
			q5refl_(&q5_c_1, wlast, &w[i], &c, &s, &p);
		}
	}
	return 0;
}
/* -----END OF Q5DLRW----------------------------------------------------- */


int q5gvns_( Real_T *x, Real_T *y, Real_T *c, Real_T *s, Real_T *p )
{
	/* Local variables */
	Real_T t;

	/*--------------------------------------------------------------------------
	PURPOSE:
	    COMPUTES PARAMETERS FOR THE GIVENS MATRIX G FOR WHICH
	    (X,Y)G = (Z,0). REPLACES (X,Y) BY (Z,0).
	--------------------------------------------------------------------------*/

	if (*y == 0.)
	{
		*c = 1.;
		*s = 0.;
		*p = 0.;
	}
	else if (*x == 0.)
	{
		*c = 0.;
		*s = 1.;
		*p = 1.;
		*x = *y;
		*y = 0.;
	}
	else
	{
		/* Computing 2nd power */
		t = sqrt( *x * *x + *y * *y );
		t = ( *x > 0.0 ) ? t : -t;
		*c = *x / t;
		*s = *y / t;
		*p = *y / (t + *x);
		*x = t;
		*y = 0.;
	}
	return 0;
}
/* -----END OF Q5GVNS----------------------------------------------------- */


int q5refl_( Int_T *n, Real_T *x, Real_T *y, Real_T *c, Real_T *s, Real_T *p )
{
	/* Local variables */
	Int_T i;
	Real_T u;

	/*--------------------------------------------------------------------------
	PURPOSE:
	    REPLACES THE TWO-COLUMN MATRIX (X,Y) BY (X,Y)G ,
	    WHERE G IS THE GIVENS MATRIX FROM SUBROUTINE Q5GVNS.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--y;
	--x;

	/* Function Body */
	if (*n < 1)
		return 0;

	if (*c == 0.)
		for (i = 1; i <= *n; ++i)
		{
			u = x[i];
			x[i] = y[i];
			y[i] = u;
		}
	else
		for (i = 1; i <= *n; ++i)
		{
			u = x[i];
			x[i] = u * *c + y[i] * *s;
			y[i] = (x[i] + u) * *p - y[i];
		}
	return 0;
}
/* -----END OF Q5REFL----------------------------------------------------- */


int q5rexp_( Int_T *m, Real_T *r, Real_T *c )
{
	/* Local variables */
	Int_T iold, inew, j, jdummy;

	/*--------------------------------------------------------------------------
	PURPOSE:
	    AUGMENTS THE UPPER TRIANGULAR MATRIX R STORED BY ROWS,
	    R11,R12,...,R1M,R22,...,R2M,...,RMM,
	    BY ADDING THE COLUMN C, WHICH YIELDS THE NEW R OF THE FORM
	    R11,...,R1M,C(1),R22,...,R2M,C(2),R33,...,RMM,C(M),C(M+1).
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--c;
	--r;

	/* Function Body */
	j = Int_T( *m + 1 );
	iold = Int_T( *m * j / 2 + 1 );
	inew = Int_T( iold + *m );
	r[inew] = c[j];

	while( j != 1 )
	{
		--j;
		--inew;
		r[inew] = c[j];
		for (jdummy = j; jdummy <= *m; ++jdummy)
		{
			--inew;
			--iold;
			r[inew] = r[iold];
		}
	}
	return 0;
}
/* -----END OF Q5REXP----------------------------------------------------- */


int q5rpck_( Int_T *m, Real_T *r, Int_T *idel )
{
	Int_T iold, inew;
	Int_T lbreak, len;

	/*--------------------------------------------------------------------------
	SUBROUTINES CALLED
	    DCOPY

	PURPOSE
	    DELETES COLUMN IDEL FROM AN UPPER TRIANGULAR MATRIX R
	    STORED BY ROWS: R11, ..., R1M, R22, ...,R2M, ..., RMM.
	    DELETES ALSO DIAGONAL ELEMENTS FOR I GREATER THAN IDEL.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--r;

	/* Function Body */
	inew = *idel;
	iold = Int_T( *idel + 1 );
	lbreak = Int_T( *m - *idel );

	for( len = Int_T( *m - 2 ); len >= lbreak; len-- )
	{
		dcopy_(len, &r[iold], &r[inew]);
		iold = Int_T( inew + *m );
		inew += len;
	}

	for( len = lbreak; len >= 1; len-- )
	{
		dcopy_(len, &r[iold], &r[inew]);
		iold = Int_T( iold + len + 1 );
		inew += len;
	}

	return 0;
}
/* -----END OF Q5RPCK----------------------------------------------------- */


int q5tris_( Int_T *m, Real_T *u, Real_T *x, Bool_T *trans )
{
	/* Local variables */
	Int_T lrow, i;
	Real_T s;
	Int_T iu;

	/*--------------------------------------------------------------------------
	SUBROUTINES CALLED
	    DSTEP

	PURPOSE 
	    SUBROUTINE FOR SOLVING A SYSTEM OF M LINEAR EQUATIONS
	    WITH AN UPPER TRIANGULAR MATRIX U, WHICH IS STORED BY ROWS,
	    I.E. U11, U12, ..., U1M, U22, ..., U2M, U33, ..., UMM.
	    (U' LOWER TRIANGULAR AND STORED BY COLUMNS).
	    X IS THE RIGHT HAND SIDE AND IS OVERWRITTEN BY THE SOLUTION.
	    IF TRANS IS .TRUE. THE SYSTEM WITH U(TRANSPOSE) IS SOLVED.
	--------------------------------------------------------------------------*/

	/* Parameter adjustments */
	--x;
	--u;

	/* Function Body */
	if (*m < 1)
		return 0;

	if( ! *trans )
	{
		//     MULTIPLY X BY THE INVERSE OF U.
		//
		s = 0.;
		for( iu = Int_T( *m * (*m + 1) / 2 ), i = *m, lrow = 1; ;
			iu--, i--, lrow++ )
		{
			x[i] = (x[i] - s) / u[iu];
			if( i == 1 )
				break;

			iu -= lrow;
			s = ddot_(lrow, &u[iu], &x[i]);
		}
	}
	else
	{
		//     MULTIPLY X BY THE INVERSE OF THE TRANSPOSE OF U.
		//
		for( iu = 1, i = 1, lrow = Int_T( *m - 1 ); ; lrow-- )
		{
			x[i] /= u[iu];
			if( i == *m )
				break;

			s = x[i];
			++i;
			++iu;
			dstep_(lrow, &x[i], &u[iu], -s);
			iu += lrow;
		}
	}

	return 0;
}
/* -----END OF Q5TRIS----------------------------------------------------- */
