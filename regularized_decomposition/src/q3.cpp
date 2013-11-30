/* ========================================================================	*/
/*                           FILE Q3										*/
/*				Routines for Solving Regularized Master						*/
/*																			*/
/*     Q3CHBD - TO CHOOSE A VIOLATED BOUND,									*/
/*     Q3CHCT - TO CHOOSE A VIOLATED CUT,									*/
/*     Q3CORR - TO MAKE A STEP IN THE DUAL VARIABLES,						*/
/*     Q3GETY - TO RECOVER THE PRIMAL SOLUTION,								*/
/*     Q3LDEP - TO FIND COEFFICIENTS OF LINEAR DEPENDENCE,					*/
/*     Q3PRIC - TO FIND DUAL VARIABLES FOR THE NEW ACTIVE SET,				*/
/*     Q3RESP - TO RESET DUAL VARIABLES SO THAT THEY SUM TO ONE,			*/
/*     Q3RESR - TO RESET THE WHOLE FACTORIZATION,							*/
/*     Q3RESZ - TO RESET FACTORIZATION PARTS DEPENDENT ON					*/
/*              THE REGULARIZING POINT AND ON THE PENALTY.					*/
/*																			*/
/* ------------------------------------------------------------------------	*/
/*     WRITTEN BY A. RUSZCZYNSKI , INSTITUT FUER OPERATIONS RESEARCH,		*/
/*     UNIVERSITAET ZUERICH, MARCH 1985.									*/
/*     DATE LAST MODIFIED: NOVEMBER 1985.									*/
/*		q3.f -- translated by f2c (version 19940705.1).						*/
/*--------------------------------------------------------------------------*/
/*	  Modified by Artur Swietanowski.										*/
/* ========================================================================	*/


#include <assert.h>
#include <math.h>

#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif

static Int_T q3_c_0 = 0;
static Int_T q3_c_1 = 1;



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3CHBD														*/
/*  Purpose   :  VERIFIES MARKED BOUNDS ON VARIABLES.						*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master), Q1CMTE (the upper level algorithm)	*/
/* Subroutines called:  None												*/
/* ------------------------------------------------------------------------	*/


Real_T q3chbd_( Int_T nfree, Real_T *xmin, Real_T *xmax, Real_T *y, // )
	Int_T *irn, Int_T *istat, Int_T *jmax, Real_T tolcut, Int_T istch )
{
		// Parameter adjustments
		--istat; --irn; --y; --xmax; --xmin;


		Real_T gmax = tolcut;

		*jmax = 0;
		for( Int_T j = 1; j <= nfree; ++j )
		{
			Int_T jj = irn[j];

			if( istat[jj] != istch ) continue;

			if( y[jj] < xmin[jj] - gmax )
			{
				*jmax = Int_T( -j );
				gmax = xmin[jj] - y[jj];
			}
			else if( y[jj] > xmax[jj] + gmax )
			{
				*jmax = j;
				gmax = y[jj] - xmax[jj];
			}
		}

		return gmax;
}
/* -----END OF Q3CHBD----------------------------------------------------- */



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3CHCT														*/
/*  Purpose   :  VERIFIES MARKED OBJECTIVE AND FEASIBILITY CUTS.			*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master)										*/
/* Subroutines called:  ddots, fabs											*/
/* ------------------------------------------------------------------------	*/


int q3chct_( Int_T *n, Int_T *mg, Real_T *g, Real_T *a, Real_T *y, // )
	Real_T *v, Int_T *iblock, Int_T *icheck, Int_T *inew, Real_T *gmax,
	Real_T tolcut, Int_T *ieq )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T i, istop;
	Real_T gi;
	Int_T ibi;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--icheck;
	--iblock;
	--a;
	--y;
	--v;
	--ieq;

	/* Function Body */
	*inew = 0;
	istop = 0;

L10:
	i_1 = *mg;
	for (i = 1; i <= i_1; ++i) {
		if (icheck[i] > 0) {
			goto L20;
		}
		if (icheck[i] == -1) {
			icheck[i] = 0;
		}
		goto L30;

L20:
		gi = ddot_(*n, &g[i * g_dim1 + 1], &y[1]) + a[i];
		ibi = iblock[i];
		if (ibi > 0) {
			gi -= v[ibi];
		}
		if (ieq[i] == 0) {
			goto L21;
		}
		gi = fabs(gi);
		if (icheck[i] == -3) {
			gi = 0.;
		}

L21:
		if (gi > tolcut) {
			goto L25;
		}
		icheck[i] = -1;
		goto L30;
L25:
		if (gi <= *gmax) {
			goto L30;
		}
		*gmax = gi;
		*inew = i;
L30:
		;
	}

	if( *gmax <= tolcut)
		goto L35;

	if (ieq[*inew] == 1) {
		icheck[*inew] = -3;
	}
	goto L50;

L35:
	if (istop > 0) {
		goto L50;
	}
	istop = 1;
	i_1 = *mg;
	for (i = 1; i <= i_1; ++i) {
		if (icheck[i] == 0) {
			icheck[i] = 1;
		}
		/* L40: */
	}
	goto L10;

L50:
	return 0;
}
/* -----END OF Q3CHCT----------------------------------------------------- */



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3CORR														*/
/*  Purpose   :CORRECTS DUAL VARIABLES DECREASING NONBASICS IN DIRECTION PI	*/
/*			  UNTIL A NONNEGATIVITY BOUND IS HIT OR (IN CASE OF RELAXATION)	*/
/*			  MAXIMUM STEPSIZE OF ONE IS USED.								*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master)										*/
/* Subroutines called:  DDIFR,DSTEP,DSUMA									*/
/* ------------------------------------------------------------------------	*/

int q3corr_(Int_T *nfix, Int_T *mtot, Int_T *l, Int_T *idel, Int_T *ldel, // )
	Int_T *ifdep, Real_T *pricnb, Real_T *pricba, Real_T *pi, Real_T *d,
	Int_T *inonba, Int_T *iblock, Int_T *ieq )
{
	/* System generated locals */
	Int_T i_1;
	Real_T d_1;

	/* Local variables */
	Int_T i, k, found;
	Real_T theta;
	Int_T ip;
	Real_T pivmax;

	/* Parameter adjustments */
	--ieq;
	--iblock;
	--inonba;
	--d;
	--pi;
	--pricba;
	--pricnb;

	/* Function Body */
	*idel = 0;
	*ldel = 0;
	found = 1;
	theta = 1.;
	if (*ifdep > 0) {
		found = 0;
	}

	/*     FIND THE FIRST BASIC PRICE THAT REACHES ZERO, IF ANY. */

	dzero_(*l, &d[1]);
	ip = Int_T( *nfix + 1 );
	i_1 = *mtot;
	for (i = ip; i <= i_1; ++i) {
		k = inonba[i];
		if (k < 1) {
			goto L5;
		}
		k = iblock[k];
		if (k > 0) {
			d[k] += pi[i];
		}
L5:
		;
	}
	i_1 = *l;
	for (i = 1; i <= i_1; ++i) {
		if (d[i] <= 0.) {
			goto L10;
		}
		if ( found )
			if (pricba[i] >= theta * d[i]) goto L10;
		found = 1;
		theta = pricba[i] / d[i];
		*ldel = i;
L10:
		;
	}

	/*     FIND A NONBASIC PRICE THAT REACHES ZERO NOT LATER. */

	if (*ldel > 0)
		*idel = Int_T( *mtot + 1 );

	pivmax = 0.;
	i_1 = *mtot;
	for (i = 1; i <= i_1; ++i) {
		if (pi[i] >= 0.) continue;
		if ( (i > *nfix) && (ieq[inonba[i]] == 1) ) continue;
		if ( !found ) goto L12;
		if ( (d_1 = pricnb[i] + theta * pi[i]) < 0.) goto L12;
		if ( (d_1 == 0) && (pivmax < -pi[i]) ) goto L12;
		continue;
L12:
		found = 1;
		pivmax = -pi[i];
		theta = pricnb[i] / pivmax;
		*idel = i;
	}
	if (*ifdep > 0) {
		goto L17;
	}
	if (theta < 1.0) {
		goto L25;
	}
	if (*idel == 0) {
		goto L20;
	}
	*idel = 0;
	goto L25;
L17:
	if (*idel > 0) {
		goto L25;
	}
	/*     UNBOUNDEDNESS. */
	*ifdep = -1;
	goto L50;
	/*     STEPSIZE EQUAL ONE IS GOOD. */
L20:
	dsuma_(*mtot, &pricnb[1], &pi[1]);
	ddifr_(*l, &pricba[1], &d[1]);
	goto L50;
	/*     STEPSIZE EQUALS THETA. */
L25:
	dstep_(*mtot, &pricnb[1], &pi[1], theta);
	dstep_(*l, &pricba[1], &d[1], -theta);
	if (*idel > *mtot) {
		goto L28;
	}
	*ldel = 0;
	goto L50;
	/*     FIND THE FIRST NONBASIC PRICE FROM BLOCK LDEL. */
L28:
	i_1 = *mtot;
	for (i = ip; i <= i_1; ++i) {
		k = inonba[i];
		if (iblock[k] == *ldel) {
			goto L35;
		}
		/* L30: */
	}
L35:
	*idel = i;

	assert( *idel <= *mtot );
	assert( *ldel > 0 );

L50:
	return 0;
}
/*-----END OF Q3CORR----------------------------------------------------------*/



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3GETY														*/
/*  Purpose   :  RECOVERS THE PRIMAL SOLUTION FROM THE DUAL.				*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master)										*/
/* Subroutines called:  DMULT, DSTPCK, DZERO								*/
/* ------------------------------------------------------------------------	*/


int q3gety_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *y, Real_T *yb, // )
	Real_T *pi, Real_T *g, Real_T *pricnb, Int_T *inonba, Int_T *irn,
	Real_T *penlty )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;
	Real_T d_1;

	/* Local variables */
	Int_T i, nfree;
	Int_T ii, ip;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--y;
	--yb;
	--pi;
	--pricnb;
	--inonba;
	--irn;

	/* Function Body */
	nfree = Int_T( *n - *nfix );
	if (nfree < 1) {
		return 0;
	}
	dzero_(nfree, &pi[1]);
	if (*m < 1) {
		goto L15;//L16 would work better actually!
	}

	//pi = g*pricnb[of free vars]
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		ip = Int_T( *nfix + i );
		ii = inonba[ip];
		dstpck_(nfree, &pi[1], &irn[1], &g[ii * g_dim1 + 1], pricnb[ip]);
	}
L15:
	//pi = pi*(1/pen)
	if (*penlty != 1.) {
		d_1 = 1. / *penlty;
		dmult_(nfree, &pi[1], &d_1);
	}
L16:    
	//y = yb - pi
	i_1 = nfree;
	for (i = 1; i <= i_1; ++i) {
		ii = irn[i];
		y[ii] = yb[ii] - pi[i];
		/* L20: */
	}
	return 0;
} 
/* -----END OF Q3GETY---------------------------------------------------------*/


/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3LDEP														*/
/*  Purpose   :  CALCULATES THE COEFFICIENTS PI OF LINEAR DEPENDENCE		*/
/*				OF THE CUTS IN THE ACTIVE SET.								*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master)										*/
/* Subroutines called:  DCOPNE, DSTPCK, DZERO, Q5TRIS						*/
/* ------------------------------------------------------------------------	*/


int q3ldep_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *r, Real_T *pi, // )
	Real_T *col, Real_T *g, Int_T *inonba, Int_T *istat, Int_T *inew )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T i;
	Int_T ii;
	Int_T ip;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--r;
	--pi;
	--col;
	--inonba;
	--istat;

	/* Function Body */
	dcopne_(*m, &col[1], &pi[*nfix + 1]);
	q5tris_(m, &r[1], &pi[*nfix + 1], (Bool_T*)&q3_c_0);
	if (*nfix < 1) {
		goto L60;
	}
	if (*inew > 0) {
		goto L10;
	}
	dzero_(*nfix, &pi[1]);
	goto L30;
L10:
	i_1 = *nfix;
	for (i = 1; i <= i_1; ++i) {
		ii = inonba[i];
		pi[i] = g[ii + *inew * g_dim1];
		/* L20: */
	}
L30:
	if (*m < 1) {
		goto L45;
	}
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		ip = Int_T( *nfix + i );
		ii = inonba[ip];
		dstpck_(*nfix, &pi[1], &inonba[1], &g[ii * g_dim1 + 1], pi[ip]);
		/* L40: */
	}
L45:
	i_1 = *nfix;
	for (i = 1; i <= i_1; ++i) {
		ii = inonba[i];
		if (istat[ii] < 0) {
			pi[i] = -pi[i];
		}
		/* L50: */
	}
L60:
	ii = Int_T( *nfix + *m + 1 );
	pi[ii] = 1.;
	return 0;
}
/*-----END OF Q3LDEP----------------------------------------------------------*/


/* ------------------------------------------------------------------------	*/
/*  Function  :  Q3PRIC														*/
/*  Purpose   :  CALCULATES MULTIPLIERS FOR THE NEW ACTIVE SET.				*/
/*																			*/
/* Called by  :  Q2MSTR (Solve Master)										*/
/* Subroutines called:  DCOPY, DMULT, DPACK, DSTPCK, DSUMA, Q5TRIS			*/
/* ------------------------------------------------------------------------	*/

int q3pric_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *y, Real_T *yb, 
	Real_T *g, Real_T *r, Real_T *pi, Real_T *z, Real_T *w, Real_T *penlty,
    Int_T *inonba, Int_T *istat )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T i;
	Int_T ii;
	Int_T ip;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--y;
	--yb;
	--r;
	--pi;
	--z;
	--w;
	--inonba;
	--istat;

	/* Function Body */

	//calculate pi[free]
	dcopy_(*m, &z[1], &pi[*nfix + 1]);      //pi = z
	dsuma_(*m, &pi[*nfix + 1], &w[1]);      //pi = z + w 
	q5tris_(m, &r[1], &pi[*nfix + 1], (Bool_T*)&q3_c_0); //Solve for H by RH=pi, and set pi=H

	//calculate pi[fix]
	dpack_(*nfix, &yb[1], &pi[1], &inonba[1]);   //pi[fixed] = yb[fixed var]
	if (*m < 1) {
		goto L30;
	}
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		ip = Int_T( *nfix + i );
		ii = inonba[ip];		//pi[fix]=pi[fix] - g[fixed var]*pi[free]
		dstpck_(*nfix, &pi[1], &inonba[1], &g[ii * g_dim1 + 1], -pi[ip]);
	}
L30:
	if (*nfix < 1) {
		goto L50;
	}
	i_1 = *nfix;
	for (i = 1; i <= i_1; ++i) {
		ii = inonba[i];
		pi[i] = y[ii] - pi[i];	//pi[fix]= y[fix]-pi[fix]
		if (istat[ii] < 0) {
			pi[i] = -pi[i];
		}
		/* L40: */
	}
	//end of pi[fix]
L50:
	if (*penlty != 1.)		//adjust pi[everything] by the penalty
		dmult_( Int_T( *nfix + *m ), &pi[1], penlty);

	return 0;
}
/*-----END OF Q3PRIC----------------------------------------------------------*/


int q3resp_(Int_T *l, Int_T *m, Real_T *pricba, Real_T *pricnb, // )
	Real_T *dpb, Int_T *inonba, Int_T *iblock, Real_T *weight, Int_T *ieq )
{
	/* System generated locals */
	Int_T i_1;
	Real_T d_1;

	/* Local variables */
	Int_T inbi, i, k;

	/* LOCAL VARIABLES :     DOUBLE PRECISION SUMK	*/
	/* FUNCTIONS :     NONE							*/
	/* SUBROUTINES CALLED :NONE						*/
	/* PURPOSE: */
	/*     NORMALIZES MULTIPLIERS TO SUM TO WEIGHT(K) IN BLOCK K. */

	/* Parameter adjustments */
	--ieq;
	--weight;
	--iblock;
	--inonba;
	--dpb;
	--pricnb;
	--pricba;

	/* Function Body */
	i_1 = *l;
	for (k = 1; k <= i_1; ++k) {
		dpb[k] = pricba[k] - weight[k];
	}
	if (*m < 1) {	//there are no non-basic cuts (that are not simple bounds)
		goto L25;	//there is no need to normalize the multipliers. They should equal to weight
	}
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		inbi = inonba[i];
		k = iblock[inbi];
		if (k <= 0) {	//nb cut does not belong to a subproblem, or a feasibility cut, do nothing
			goto L20;
		}
		if (pricnb[i] <= 0.) {	//??
			goto L20;
		}
		dpb[k] += pricnb[i];
L20:
		;
	}
L25:
	i_1 = *l;
	for (k = 1; k <= i_1; ++k) {
		/* Computing MAX */
		d_1 = pricba[k] - dpb[k];  
		pricba[k] = Max(d_1,0.);   
		/* L30: */
	}
	/*     DO 30 K = 1,L */
	/*       SUMK = PRICBA(K) - WEIGHT(K) */
	/*       NK = 1 */
	/*       IF ( M .LT. 1 ) GOTO 15 */
	/*       DO 10 I = 1,M */
	/*         INBI = INONBA(I) */
	/*         IF ( IBLOCK(INBI) .NE. K ) GOTO 10 */
	/*         IF ( PRICNB(I) .LE. 0.0D0 ) GOTO 10 */
	/*         SUMK = SUMK + PRICNB(I) */
	/*         NK = NK + 1 */
	/*  10   CONTINUE */
	/*  15   IF ( DABS(SUMK) .LT. SMALL ) GOTO 30 */
	/*       SUMK = SUMK / FLOAT(NK) */
	/*       PRICBA(K) = DMAX1( PRICBA(K) - SUMK, 0.0D0 ) */
	/*       IF ( M .LT. 1 ) GOTO 30 */
	/*       DO 20 I = 1,M */
	/*         INBI = INONBA(I) */
	/*         IF ( IBLOCK(INBI) .NE. K ) goto 20 */
	/*         PRICNB(I) = DMAX1( PRICNB(I) - SUMK, 0.0D0 ) */
	/*  20   CONTINUE */
	/*  30 CONTINUE */
	return 0;
}
/*-----END OF Q3RESP----------------------------------------------------------*/


int q3resr_( Int_T *n, Int_T *nfix, Int_T *l, Int_T * m, Real_T *g, // )
	Real_T *a, Real_T *x, Real_T *y, Real_T *yb, Real_T *weight,
	Real_T *penlty, Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *col,
    Real_T *pi, Real_T *pricba, Real_T *pricnb, Int_T *iblock, Int_T *ibasic,
	Int_T *inonba, Int_T *irn )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;
	Real_T d_1;

	/* Local variables */
	Int_T inew, i, ifdep;
	Int_T ii, jj;

	/* PURPOSE */
	/*     RESETS QR FACTORS FOR THE CURRENT ACTIVE SET */

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--a;
	--x;
	--y;
	--yb;
	--weight;
	--q;
	--r;
	--z;
	--w;
	--col;
	--pi;
	--pricba;
	--pricnb;
	--iblock;
	--ibasic;
	--inonba;
	--irn;

	/* Function Body */
	dzero_(*n, &yb[1]);
	i_1 = *l;
	for (i = 1; i <= i_1; ++i) {
		jj = ibasic[i];
		dstep_(*n, &yb[1], &g[jj * g_dim1 + 1], weight[i]);
	}
	d_1 = -1. / *penlty;
	dmult_(*n, &yb[1], &d_1);
	dsuma_(*n, &yb[1], &x[1]);
	i = 0;
L20:
	if (i == *m) {
		return 0;
	}
	/*     ADD THE SUCCESSIVE CUT. */
	ii = Int_T( *nfix + i + 1 );
	inew = inonba[ii];
	ifdep = 1;
	q4adct_(n, nfix, &i, &g[g_offset], &a[1], &q[1], &r[1],
		&z[1], &w[1], &y[1], &yb[1], &col[1], &pi[1], &pricnb[1],
		&inonba[1], &irn[1], & ifdep, (Bool_T*)&q3_c_0, &inew);
	jj = iblock[inew];
	if (ifdep == 0) {
		goto L20;
	}
	/*     AUGMENTATION FAILED. SKIP THIS CUT. */
	if (jj < 1) {
		goto L30;
	}
	pricba[jj] += pricnb[ii];
	jj = ibasic[jj];
	dsuma_(*n, &g[inew * g_dim1 + 1], &g[jj * g_dim1 + 1]);
	a[inew] += a[jj];
L30:
	dcopy_( Int_T( *m - i ), &pricnb[ii + 1], &pricnb[ii]);
	icopy_( Int_T( *m - i ), &inonba[ii + 1], &inonba[ii]);
	--(*m);
	--i;
	goto L20;
}
/*-----END OF Q3RESR----------------------------------------------------------*/


int q3resz_( Int_T *n, Int_T *nfix, Int_T *l, Int_T *m, Real_T *x, // )
	Real_T *yb, Real_T *g, Real_T *r, Real_T *z, Real_T *pi, Int_T *ibasic,
	Int_T *inonba, Int_T *irn, Real_T *weight, Real_T *penlty, Bool_T *newyb )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;
	Real_T d_1;

	/* Local variables */
	Int_T i;
	Int_T jj;

	/* FUNCTIONS */
	/* PURPOSE */
	/*     RESETS THE AUXILIARY VECTORS YB AND Z. */

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--x;
	--yb;
	--r;
	--z;
	--pi;
	--ibasic;
	--inonba;
	--irn;
	--weight;

	/* Function Body */
	if (! *newyb)
		goto L20;

	//recalculate yb
	
	dzero_(*n, &yb[1]);
	
	// multiplies each basic cut by its weight and adds each term by term
	// yb[i] = sum(l, g[l-th cut's ith position]*weight[l])
	i_1 = *l;
	for (i = 1; i <= i_1; ++i) { 
		jj = ibasic[i];
		dstep_(*n, &yb[1], &g[jj * g_dim1 + 1], weight[i]);
	}
	
	//adjust yb by the penalty and add x term by term
	d_1 = -1. / *penlty;
	dmult_(*n, &yb[1], &d_1);
	dsuma_(*n, &yb[1], &x[1]);
L20:
	if (*m < 1) {
		return 0;
	}
	//recalculate z

	//  pi[i] = yb[ith free variable]  
	//  e.g. first free variable is the 3rd variable, then pi[1] = yb[3] 
	dpack_( Int_T( *n - *nfix ), &yb[1], &pi[1], &irn[1]);
	
	//calculates m values of z
	//z[i] = dot_prod (pi, g[free variables])
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];
		z[i] = ddots_( Int_T( *n - *nfix ), &pi[1], &irn[1],
			&g[jj * g_dim1 + 1]);
		/* L30: */
	}
	
	//solves a triangular system to get z
	q5tris_(m, &r[1], &z[1], (Bool_T*)&q3_c_1);
	return 0;
}
/* -----END OF Q3RESZ----------------------------------------------------- */
