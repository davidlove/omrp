/* =======================================================================	*/
/*                          FILE Q4											*/
/*		Routines for Active Set Strategy used in Solving Reg Master			*/
/*																			*/
/*     Q4ADBD - TO ADD A BOUND TO THE ACTIVE SET,							*/
/*     Q4ADCT - TO ADD A CUT TO THE ACTIVE SET,								*/
/*     Q4DTBD - TO DELETE A BOUND FROM THE ACTIVE SET,						*/
/*     Q4DTCT - TO DELETE A CUT FROM THE ACTIVE SET,						*/
/*     Q4ORTG - TO ORTHOGONALIZE A NEW CUT									*/
/*              WITH RESPECT TO THE ACTIVE SET,								*/
/*     Q4SBCL - TO UPDATE FACTORS WHEN THE BASIS CHANGES.					*/
/* ------------------------------------------------------------------------	*/
/*     WRITTEN BY A. RUSZCZYNSKI , INSTITUT FUER OPERATIONS RESEARCH,		*/
/*     UNIVERSITAET ZUERICH, MARCH 1985.									*/
/*     DATE LAST MODIFIED: March 1996										*/
/* =======================================================================	*/
/*       Translated from FORTRAN by f2c (version 19940705.1).				*/
/*       Modified by Artur Swietanowski.                                    */
/* =======================================================================	*/


#include <math.h>
#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif


static Int_T c_1 = 1;
static Int_T c_0 = 0;


/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4ADBD														*/
/*  Purpose   :  ADDS A BOUND ON VARIABLE TO THE ACTIVE SET BY DELETING		*/
/*               THE CORRESPONDING ROW FROM THE QR FACTORIZATION.			*/
/*																			*/
/* Functions  :  DOUBLE PRECISION DNORM2									*/
/* Subroutines:  DCOPY,DCOPNE,DZERO,DMULT,DSTPCK,DSTUNP,Q5TRIS,Q5DLRW		*/
/* ------------------------------------------------------------------------	*/

int q4adbd_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *xmin, Real_T *xmax, // )
	Real_T *g, Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *y,
	Real_T *yb, Real_T *col, Real_T *pi, Real_T *pricnb, Int_T *inonba,
	Int_T *irn, Int_T *istat, Int_T *ifdep, Bool_T *nodel, Int_T *inew )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T irow, mtot, i, inabs, nfree;
	Real_T wlast, prnew = 0.0;
	Real_T zlast;
	Int_T jj;
	Real_T rho;


	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--xmin;
	--xmax;
	--q;
	--r;
	--z;
	--w;
	--y;
	--yb;
	--col;
	--pi;
	--pricnb;
	--inonba;
	--irn;
	--istat;

	/* Function Body */
	inabs = Int_T( (*inew >= 0) ? *inew : - *inew );
	irow = irn[inabs];
	nfree = Int_T( *n - *nfix );
	dzero_(nfree, &q[1]);
	q[inabs] = 1.;
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];
		pi[i] = g[irow + jj * g_dim1];
		/* L10: */
	}
	rho = 1.;

	/*     ORTHOGONALIZE THE UNIT VECTOR WITH RESPECT TO G*R(INV). */

	q4ort1_(n, nfix, &nfree, m, &g[g_offset], &q[1], &r[1],
		&col[1], &rho, & pi[1], &inonba[1], &irn[1]);

	if (rho == 0.) {
		goto L100;
	}

	if (*inew > 0) {
		goto L20;
	}
	istat[irow] = 0;
	goto L25;
L20:
	istat[irow] = -1;
L25:

	/*     END OF ORTHOGONALIZATION. INSERT THE NEW PRICE. */

	/* L35: */
	++(*nfix);
	mtot = Int_T( *nfix + *m );
	if (*ifdep > 0) {
		prnew = pricnb[mtot];
	}
	if (*m < 1) {
		goto L42;
	}
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( mtot - i );
		pricnb[jj + 1] = pricnb[jj];
		inonba[jj + 1] = inonba[jj];
		/* L40: */
	}
L42:
	inonba[*nfix] = irow;
	if (*ifdep == 0) {
		goto L45;
	}
	pricnb[*nfix] = prnew;
	goto L46;
L45:
	pricnb[*nfix] = 0.;

	/*     AUGMENT R, Z AND W, AND CORRECT Y. */

L46:
	if (*inew < 0) {
		goto L50;
	}
	y[irow] = xmax[irow];
	goto L60;
L50:
	y[irow] = xmin[irow];
L60:
	zlast = ddots_(nfree, &q[1], &irn[1], &yb[1]) / rho;
	wlast = -(y[irow] + ddot_(*m, &w[1], &col[1])) / rho;
	dcopy_( Int_T( nfree - inabs ), &q[inabs + 1], &q[inabs]);
	icopy_( Int_T( nfree - inabs ), &irn[inabs + 1], &irn[inabs]);
	--nfree;
	if (*nodel)
		dstunp_(nfree, &y[1], &q[1], &irn[1], -(wlast + zlast) / rho);

	if (*inew > 0) {
		goto L70;
	}
	rho = -rho;
	zlast = -zlast;
	wlast = -wlast;
L70:

	/*     REORDER CUTS SO THAT THE BOUND BECOMES FIRST  */
	/*     AND DELETE THE BOUND FROM EXPLICIT FACTORS.   */

	q5dlrw_(m, &r[1], &z[1], &w[1], &pi[1], &col[1], &rho, &
		zlast, &wlast);
	*ifdep = 0;
	return 0;

	/*     LINEAR DEPENDENCE. SET DATA FOR THE EXCHANGE PROCEDURE. */

L100:
	*ifdep = 1;
	jj = Int_T( *nfix + *m + 1 );
	inonba[jj] = 0;
	pricnb[jj] = 0.;
	if (*inew < 0)
		dcopne_(*m, &col[1], &col[1]);

	return 0;
}
/* -----END OF Q4ADBD----------------------------------------------------- */


	
/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4ADCT														*/
/*  Purpose   :  ADD A NEW CUT TO THE ACTIVE SET AND AUGMENT QR FACTORS.	*/
/*																			*/
/* Called by  :  Q2MSTR_  and Q3RESR_										*/
/* Subroutines called:  DPACK,DSTUNP,Q5REXP,Q4ORTG							*/
/* ------------------------------------------------------------------------	*/


int q4adct_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *a, 
	Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *y, Real_T *yb,
	Real_T *col, Real_T *pi, Real_T *pricnb, Int_T *inonba, Int_T *irn,
	Int_T *ifdep, Bool_T *nodel, Int_T *inew )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T mtot, i;
	Real_T alpha;
	Int_T ii;
	Real_T rho;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--a;
	--q;
	--r;
	--z;
	--w;
	--y;
	--yb;
	--col;
	--pi;
	--pricnb;
	--inonba;
	--irn;

	/* Function Body */
	dpack_( Int_T( *n-*nfix ), &g[*inew * g_dim1 + 1], &q[1], &irn[1]);  //q[free] = g[free] (inew'th cut)

	/*     ORTHOGONALIZE. */

	i_1 = Int_T( *n - *nfix );
	q4ortg_(n, nfix, &i_1, m, &g[g_offset], &q[1], &r[1], 
		&col[1], &rho, &pi[1], &inonba[1], &irn[1]);
	++(*m);
	mtot = Int_T( *nfix + *m );
	if (*ifdep > 0) {
		goto L10;
	}
	inonba[mtot] = *inew;
	pricnb[mtot] = 0.;
L10:
	if (rho > 0.) {
		goto L20;
	}
	*ifdep = 1;
	return 0;
L20:
	*ifdep = 0;

	/*     AUGMENT R, Z AND W, AND CORRECT Y. */

	col[*m] = rho;
	i_1 = Int_T( *m - 1 );
	q5rexp_(&i_1, &r[1], &col[1]);
	alpha = 0.;
	if (*nfix < 1) {
		goto L60;
	}
	i_1 = *nfix;
	for (i = 1; i <= i_1; ++i) {
		ii = inonba[i];
		alpha += g[ii + *inew * g_dim1] * y[ii];
		/* L50: */
	}
L60:
	alpha += a[*inew];
	z[*m] = ddots_( Int_T( *n-*nfix ), &q[1], &irn[1], &yb[1]) / rho;
	w[*m] = (alpha - ddot_( Int_T( *m-1 ), &w[1], &col[1])) / rho;
	if (*nodel)
		dstunp_( Int_T( *n - *nfix ), &y[1], &q[1], &irn[1],
			-(z[*m] + w[*m]) / rho);
	return 0;
}
/* -----END OF Q4ADCT----------------------------------------------------- */



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4DTBD														*/
/*   Purpose  :  DELETE A BOUND ON VARIABLE FROM THE ACTIVE SET				*/
/*               BY ADDING A ROW TO THE QR FACTORIZATION.					*/
/*																			*/
/* Called by  :  Q2MSTR														*/	
/* Subroutines:  DCOPY, ICOPY, Q5ADRW										*/
/* ------------------------------------------------------------------------	*/


int q4dtbd_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *y, // )
	Real_T *yb, Real_T *r, Real_T *z, Real_T *w, Real_T *pi, Real_T *pricnb,
	Int_T *inonba, Int_T *irn, Int_T *istat, Int_T *idel, Int_T *ifdep )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Int_T i;
	Real_T wlast, zlast;
	Int_T ii, jj;


	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--y;
	--yb;
	--r;
	--z;
	--w;
	--pi;
	--pricnb;
	--inonba;
	--irn;
	--istat;

	/* Function Body */
	jj = inonba[*idel];

	/*     PUT THE BOUND IN FRONT OF THE ACTIVE SET. */

	if (*m < 1) {
		goto L20;
	}
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		ii = Int_T( *nfix + i );
		ii = inonba[ii];
		int a = 5;
		pi[i] = g[jj + ii * g_dim1];
		/* L10: */
	}
L20:
	zlast = yb[jj];
	wlast = -y[jj];

	/*     REORDER THE ACTIVE SET SO THAT THE BOUND BECOMES LAST */
	/*     AND DELETE THE BOUND. */

	q5adrw_(m, &r[1], &z[1], &w[1], &pi[1], &zlast, &wlast);
	ii = Int_T( *nfix + *m - *idel + *ifdep );
	dcopy_(ii, &pricnb[*idel + 1], &pricnb[*idel]);
	icopy_(ii, &inonba[*idel + 1], &inonba[*idel]);
	--(*nfix);
	istat[jj] = 1;
	ii = Int_T( *n - *nfix );
	irn[ii] = jj;
	return 0;
}
/* -----END OF Q4DTBD----------------------------------------------------- */


/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4DTCT														*/
/* Purpose  :  DELETES A CUT FROM THE ACTIVE SET AND UPDATES THE QR FACTORS */
/*																			*/
/* Called by  :  Q2MSTR														*/
/* Subroutines:  DSUMA,DDIFR,DSTEP,DCOPY,ICOPY,Q4SBCL						*/
/* ------------------------------------------------------------------------	*/

int q4dtct_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *a, // )
	Real_T *yb, Real_T *weight, Real_T *penlty, Int_T *iblock, Int_T *ibasic,
	Int_T *inonba, Int_T *icheck, Real_T *r, Real_T *z, Real_T *w,
	Real_T *pricnb, Real_T *pricba, Real_T *col, Int_T * idel, Int_T *ldel,
	Int_T *ifdep )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;
	Real_T d_1;

	/* Local variables */
	Int_T i, ibold;
	Int_T ibnew;
	Int_T m1;
	Int_T ii;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--a;
	--yb;
	--weight;
	--iblock;
	--ibasic;
	--inonba;
	--icheck;
	--r;
	--z;
	--w;
	--pricnb;
	--pricba;
	--col;

	/* Function Body */
	if (*ldel > 0) {
		goto L10;
	}
	ii = inonba[*idel];
	icheck[ii] = -1;
	*ldel = iblock[ii];
	if (*ldel < 1) {
		goto L40;
	}
	/*       RESTORE THE CUT. */
	ibold = ibasic[*ldel];
	dsuma_(*n, &g[ii * g_dim1 + 1], &g[ibold * g_dim1 + 1]);
	a[ii] += a[ibold];
	goto L40;         //delete non-basic cut

	/*     CHANGE THE BASIS AND UPDATE YB. */

L10:
	ibold = ibasic[*ldel];
	ibnew = inonba[*idel];
	dstep_(*n, &yb[1], &g[ibnew * g_dim1 + 1], -weight[*ldel] / *penlty);
	m1 = Int_T( *m + *ifdep );
	ii = Int_T( *nfix + m1 );
	if (inonba[ii] == 0) {
		m1 = *m;
	}
	/*     CORRECT REDUCED CUTS FROM BLOCK LDEL. */
	if (m1 < 1) {
		goto L25;
	}
	i_1 = m1;
	for (i = 1; i <= i_1; ++i) {
		ii = Int_T( *nfix + i );
		if (ii == *idel) {
			goto L20;
		}
		ii = inonba[ii];
		if (iblock[ii] != *ldel) {
			goto L20;
		}
		ddifr_(*n, &g[ii * g_dim1 + 1], &g[ibnew * g_dim1 + 1]);
		a[ii] -= a[ibnew];
L20:
		;
	}
L25:
	dsuma_(*n, &g[ibnew * g_dim1 + 1], &g[ibold * g_dim1 + 1]);
	a[ibnew] += a[ibold];
	ibasic[*ldel] = ibnew;
	pricba[*ldel] = pricnb[*idel];
	icheck[ibold] = -1;
	if (*idel - *nfix <= *m) {
		goto L30;
	}

	/*     IF THE NEW CUT ADVANCED, ONLY Z CHANGES. */

	dstep_(*m, &z[1], &col[1], -weight[*ldel] / *penlty);
	*ifdep = 0;
	goto L50;

	/*     IF A NONBASIC CUT ADVANCED, UPDATE R AND Z. */

L30:
	i_1 = Int_T( *idel - *nfix );
	d_1 = weight[*ldel] / *penlty;
	q4sbcl_(&r[1], &z[1], m, &i_1, &iblock[1], &inonba[*nfix + 1],
		&d_1);

	/*     DELETE THE NONBASIC CUT. */

L40:
	i_1 = Int_T( *idel - *nfix );
	q5dlco_(m, &i_1, &r[1], &z[1], &w[1]);
	ii = Int_T( *nfix + *m - *idel + *ifdep );
	//moves everything up in pricba and inonba
	dcopy_(ii, &pricnb[*idel + 1], &pricnb[*idel]);
	icopy_(ii, &inonba[*idel + 1], &inonba[*idel]);
	--(*m);
L50:
	return 0;
}
/* -----END OF Q4DTCT----------------------------------------------------- */

/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4ORTG														*/
/*   Purpose  :  ORTHOGONALIZES Q WITH RESPECT TO COLUMNS					*/
/*				 OF THE ORTHOGONAL MATRIX G*R(INV) WHICH IS NOT STORED.		*/
/*																			*/
/* Called by  :  Q2MSTR														*/
/* Subroutines:  DCOPY,DDIFR,DSUMA,DMULT,DSTPCK,DZERO,Q5TRIS				*/
/* ------------------------------------------------------------------------	*/



int q4ortg_( Int_T *n, Int_T *nfix, Int_T *nfree, Int_T *m, Real_T *g, // )
	Real_T *q, Real_T *r, Real_T * col, Real_T *rho, Real_T *pi, Int_T *inonba,
	Int_T *irn )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Real_T gnor;
	Int_T i;
	Real_T rhold;
	Int_T jj;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--q;
	--r;
	--col;
	--pi;
	--inonba;
	--irn;


	/* Function Body */
	*rho = dnorm2_(*nfree, &q[1]);
	gnor = *rho + 1.;
	if (*m < 1) {
		goto L20;
	}
	dzero_(*m, &col[1]);
L5:
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];
		pi[i] = ddots_(*nfree, &q[1], &irn[1], &g[jj * g_dim1 + 1]);
	}
	q5tris_(m, &r[1], &pi[1], (Bool_T*)&c_1);
	dsuma_(*m, &col[1], &pi[1]);
	if (*m == *nfree) {
		goto L22;
	}
	q5tris_(m, &r[1], &pi[1], (Bool_T*)&c_0);
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];
		dstpck_(*nfree, &q[1], &irn[1], &g[jj * g_dim1 + 1],
			-pi[i]);
		/* L15: */
	}
	rhold = *rho;
	*rho = dnorm2_(*nfree, &q[1]);
	if (*rho < gnor * 1.0e-20) {
		goto L22;
	}
	//if (*rho < gnor * 1.0e-8) {
		//Print("rho/gnor = %G\n", *rho/gnor);
	//}
	if (*rho + *rho < rhold) {
		goto L5;
	}
L20:
	*rho = sqrt(*rho);
	return 0;
L22:
	*rho = 0.;
	return 0;
}
/* -----END OF Q4ORTG----------------------------------------------------- */



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4ORT1														*/
/*   Purpose  :  ORTHOGONALIZES A UNIT VECTOR Q WITH RESPECT TO COLUMNS		*/
/*				 OF THE ORTHOGONAL MATRIX G*R(INV) WHICH IS NOT STORED.		*/
/*																			*/
/* Called by  :  Q4ADBD														*/
/* Subroutines:  DCOPY,DDIFR,DSUMA,DMULT,DSTPCK,DZERO,Q5TRIS				*/
/* ------------------------------------------------------------------------	*/


int q4ort1_( Int_T *n, Int_T *nfix, Int_T *nfree, Int_T *m, Real_T *g, // )
	Real_T *q, Real_T *r, Real_T * col, Real_T *rho, Real_T *pi,
	Int_T *inonba, Int_T *irn )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;

	/* Local variables */
	Real_T gnor;
	Int_T i;
	Real_T rhold;
	Int_T jj;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--q;
	--r;
	--col;
	--pi;
	--inonba;
	--irn;

	/* Function Body */
	gnor = *rho;
	if (*m < 1) {
		return 0;
	}
	dzero_(*m, &col[1]);
	goto L11;

L5:
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];
		pi[i] = ddots_(*nfree, &q[1], &irn[1], &g[jj * g_dim1 + 1]);
	}
L11:
	q5tris_(m, &r[1], &pi[1], (Bool_T*)&c_1);
	dsuma_(*m, &col[1], &pi[1]);
	if (*m == *nfree)
		goto L22;

	q5tris_(m, &r[1], &pi[1], (Bool_T*)&c_0);
	i_1 = *m;
	for (i = 1; i <= i_1; ++i) {
		jj = Int_T( *nfix + i );
		jj = inonba[jj];

		dstpck_(*nfree, &q[1], &irn[1], &g[jj * g_dim1 + 1], -pi[i]);
	}
	rhold = *rho;
	*rho = dnorm2_(*nfree, &q[1]);
	if (*rho < gnor * 1.0e-20) {
		goto L22;
	}
	if (*rho + *rho < rhold) {
		goto L5;
	}
	*rho = sqrt(*rho);
	return 0;
L22:
	*rho = 0.;
	return 0;
} 
/* -----END OF Q4ORT1----------------------------------------------------- */



/* ------------------------------------------------------------------------	*/
/*  Function  :  Q4SBCL														*/
/*   Purpose  :  SUBTRACTS COLUMN IDEL OF AN UPPER TRIANGULAR MARTIX R		*/
/*				 FROM SUBSEQUENT COLUMNS ASSOCIATED WITH THE SAME BLOCK		*/
/*				 AND ( MULTIPLIED BY WEIGHT ) FROM Z. R IS STORED BY ROWS.  */
/*																			*/
/* Called by  :  Q4DTCT														*/
/* Subroutines Called:  None												*/
/* ------------------------------------------------------------------------	*/


int q4sbcl_(Real_T *r, Real_T *z, Int_T *m, Int_T *idel, Int_T *iblock, // )
	Int_T *inonba, Real_T *weight)
{
	/* System generated locals */
	Int_T i_1;

	/* Local variables */
	Int_T inbj, i, j, ib, ic, ir;

	/* Parameter adjustments */
	--inonba;
	--iblock;
	--z;
	--r;

	/* Function Body */
	ib = inonba[*idel];
	ib = iblock[ib];
	ic = *idel;
	i_1 = *idel;
	for (i = 1; i <= i_1; ++i) {
		ir = ic;
		j = *idel;
L5:
		if (j == *m) {
			goto L8;
		}
		++j;
		++ir;
		inbj = inonba[j];
		if (iblock[inbj] == ib) {
			r[ir] -= r[ic];
		}
		goto L5;
L8:
		z[i] -= *weight * r[ic];
		ic = Int_T( ic + *m - i );
		/* L10: */
	}
	return 0;
}
/* -----END OF Q4SBCL----------------------------------------------------- */
