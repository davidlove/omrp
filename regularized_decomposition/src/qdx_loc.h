/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	qdx_loc.h
CREATED:			1994.07.28
LAST MODIFIED:		1996.02.27

DEPENDENCIES:		stdtype.h, std_tmpl.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __QDX_LOC_H__
#define __QDX_LOC_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif


//==============================================================================
//
//	Declarations / definitions for QDX solver (translated from Fortran to C
//	by f2c automatic converter).
//
//==============================================================================

#define FALSE_ 0
#define TRUE_ 1

class RD_SubproblemManager;

/*--- Q1 ---------------------------------------------------------------------*/
int q1slct_( Int_T *n, Int_T *l, Int_T *m, Int_T *mg, Real_T *g, Real_T *a,
	Int_T *iblock, Int_T *ibasic, Int_T *inonba, Int_T *icheck, Int_T *ieq,
	Int_T *status, Int_T *drow );

void CheckConstraintsInSection( Int_T *marks, Int_T j1, Int_T j2,
	Int_T *iptr, Real_T *dmat, Int_T *jcol, Real_T *y, Int_T *status,
	Real_T *bmin, Real_T *bmax, Real_T tolcut, Real_T &gsec, Int_T &imax );

Bool_T CheckAllConstraints( Int_T &j1, Int_T &j2, Int_T mdmat, Int_T n,
	Int_T nsec, Int_T &mg, Int_T *iblock, Int_T *icheck, Int_T g_dim1,
	Real_T *g, Int_T *drow, Real_T *a, Int_T *ieq, Int_T *marks, Int_T *iptr,
	Real_T *dmat, Int_T *jcol, Real_T *y, Int_T *status, Real_T *bmin,
	Real_T *bmax, Real_T tolcut, Int_T nwcmax, Real_T &gmax, Int_T &inew,
	Int_T &jnew, Int_T &ifeas1 );

//@BEGIN------------------------------------------------------------------------
//
//		I am adding the variable istop in the parameter list to 
//		get the "value"s of cuts and "f2"
//		Also, ExpC and ExpCost to store them.. 

void SolveBlocks( Real_T &f1, Real_T &f2, Real_T &fx, Int_T &mg, Int_T l,
	Int_T n, Real_T *y, Real_T &value, Real_T *g, Real_T *weight, Int_T &ifeas2,
	Int_T &iphase, Int_T *iblock, Int_T *icheck, Real_T tolcut, Real_T *a,
	Int_T *ieq, Int_T &inew, Int_T &jnew, Real_T &gmax, Int_T g_dim1, Real_T *v,
	Real_T *x, RD_SubproblemManager &SubMan, Int_T *istop, Real_T *ExpC, Real_T& ExpCost, 
	Real_T *ExpC2, Real_T& ExpCost2);
//
//@END--------------------------------------------------------------------------

void DetermineStepType( Int_T &iphase, Int_T ifeas1, Int_T ifeas2,
	Int_T &ifeas, Real_T f1, Real_T f2, Real_T vsum,
	Int_T mg, Int_T l, Int_T n, Int_T g_dim1, Int_T *icheck, Int_T *iblock,
	Real_T *g, Real_T *a, Int_T iter, Int_T &index, Real_T &funold,
	Real_T &funmin, Int_T &istart, Int_T &iser, Int_T levprt, Real_T *y,
	Real_T *x, Real_T *pi, Real_T *yb, Int_T *istat, Int_T *irn,
	Real_T &penlty, Real_T initpen, Bool_T &newpen,
	Real_T penmax, Int_T &nulinf, Int_T &nserex, Int_T &nulf,
	Int_T &nserap, Real_T tolcut, Real_T gamma, Int_T inew, Int_T jnew, 
	bool lshaped);

void CountCriticalScenarios( Int_T l, Int_T m, Int_T nfix, Int_T *iblock, // )
	Int_T *inonba );

void MakeStep( Int_T n, Real_T *x, Real_T *y, Real_T *yb, Real_T *pi,
	Real_T &funold, Real_T fun, Int_T &index );

void CompressCommittee( Int_T &mg, Int_T n, Int_T l, Int_T m, Int_T nfix, // )
	Int_T *istat, Int_T *inonba, Int_T *icheck, Int_T *ieq, Int_T *status,
	Int_T *drow, Real_T *g, Real_T *a, Int_T *iblock, Int_T *ibasic,
	Int_T g_offset );

/*--- Q2 ---------------------------------------------------------------------*/
int q2mstr_( Int_T n, Int_T *nfix, Int_T *m, Int_T *mg, Int_T *l, Real_T *g,
	Real_T *a, Real_T *xmin, Real_T *xmax, Real_T *x, Real_T *y, Real_T *yb,
	Real_T *v, Real_T *weight, Real_T *penlty, Bool_T *newpen, Int_T *iblock,
	Int_T *ibasic, Int_T *inonba, Int_T *icheck, Int_T *ieq, Int_T *irn,
	Int_T *istat, Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *pricnb,
	Real_T *pricba, Real_T *pi, Real_T *col, Real_T *dpb, Int_T *index,
	Real_T *gmax, Int_T *inew, Int_T *jnew, Int_T *iter,
	Real_T *initpen, Real_T &tolcut );

/*--- Q3 ---------------------------------------------------------------------*/
Real_T q3chbd_( Int_T nfree, Real_T *xmin, Real_T *xmax, Real_T *y,
	Int_T *irn, Int_T *istat, Int_T *jmax, Real_T tolcut, Int_T istch );

int q3chct_( Int_T *n, Int_T *mg, Real_T *g, Real_T *a, Real_T *y,
	Real_T *v, Int_T *iblock, Int_T *icheck, Int_T *inew, Real_T *gmax,
	Real_T tolcut, Int_T *ieq );

int q3corr_( Int_T *nfix, Int_T *mtot, Int_T *l, Int_T *idel, Int_T *ldel,
	Int_T *ifdep, Real_T *pricnb, Real_T *pricba, Real_T *pi, Real_T *d,
	Int_T *inonba, Int_T *iblock, Int_T *ieq );

int q3gety_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *y, Real_T *yb,
	Real_T *pi, Real_T *g, Real_T *pricnb, Int_T *inonba, Int_T *irn,
	Real_T *penlty );

int q3ldep_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *r, Real_T *pi,
	Real_T *col, Real_T *g, Int_T *inonba, Int_T *istat, Int_T *inew );

int q3pric_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *y, Real_T *yb,
	Real_T *g, Real_T *r, Real_T *pi, Real_T *z, Real_T *w, Real_T *penlty,
	Int_T *inonba, Int_T *istat );

int q3resp_( Int_T *l, Int_T *m, Real_T *pricba, Real_T *pricnb,
	Real_T *dpb, Int_T *inonba, Int_T *iblock, Real_T *weight, Int_T *ieq );

int q3resr_( Int_T *n, Int_T *nfix, Int_T *l, Int_T *m, Real_T *g,
	Real_T *a, Real_T *x, Real_T *y, Real_T *yb, Real_T *weight,
	Real_T *penlty, Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *col,
	Real_T *pi, Real_T *pricba, Real_T *pricnb, Int_T *iblock, Int_T *ibasic,
	Int_T *inonba, Int_T *irn );

int q3resz_( Int_T *n, Int_T *nfix, Int_T *l, Int_T *m, Real_T *x,
	Real_T *yb, Real_T *g, Real_T *r, Real_T *z, Real_T *pi, Int_T *ibasic,
	Int_T *inonba, Int_T *irn, Real_T *weight, Real_T *penlty, 
	Bool_T *newyb );

/*--- Q4 ---------------------------------------------------------------------*/
int q4adbd_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *xmin, Real_T *xmax,
	Real_T *g, Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *y,
	Real_T *yb, Real_T *col, Real_T *pi, Real_T *pricnb, Int_T *inonba,
	Int_T *irn, Int_T *istat, Int_T *ifdep, Bool_T *nodel, Int_T *inew );

int q4adct_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *a,
	Real_T *q, Real_T *r, Real_T *z, Real_T *w, Real_T *y, Real_T *yb,
	Real_T *col, Real_T *pi, Real_T *pricnb, Int_T *inonba, Int_T *irn,
	Int_T *ifdep, Bool_T *nodel, Int_T *inew );

int q4dtbd_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *y,
	Real_T *yb, Real_T *r, Real_T *z, Real_T *w, Real_T *pi, Real_T *pricnb,
	Int_T *inonba, Int_T *irn, Int_T *istat, Int_T *idel, Int_T *ifdep );

int q4dtct_( Int_T *n, Int_T *nfix, Int_T *m, Real_T *g, Real_T *a,
	Real_T *yb, Real_T *weight, Real_T *penlty, Int_T *iblock, Int_T *ibasic,
	Int_T *inonba, Int_T *icheck, Real_T *r, Real_T *z, Real_T *w,
	Real_T *pricnb, Real_T *pricba, Real_T *col, Int_T *idel, Int_T *ldel,
	Int_T *ifdep );

int q4ortg_( Int_T *n, Int_T *nfix, Int_T *nfree, Int_T *m, Real_T *g,
	Real_T *q, Real_T *r, Real_T *col, Real_T *rho, Real_T *pi, Int_T *inonba,
	Int_T *irn );

int q4ort1_( Int_T *n, Int_T *nfix, Int_T *nfree, Int_T *m, Real_T *g,
	Real_T *q, Real_T *r, Real_T *col, Real_T *rho, Real_T *pi, Int_T *inonba,
	Int_T *irn );

int q4sbcl_( Real_T *r, Real_T *z, Int_T *m, Int_T *idel, Int_T *iblock,
	Int_T *inonba, Real_T *weight );

/*--- Q5 ---------------------------------------------------------------------*/
int q5adrw_( Int_T *m, Real_T *r, Real_T *z, Real_T *w, Real_T *pi,
	Real_T *zlast, Real_T *wlast );
int q5dlco_( Int_T *m, Int_T *idel, Real_T *r, Real_T *z, Real_T *w );
int q5dlrw_( Int_T *m, Real_T *r, Real_T *z, Real_T *w, Real_T *pi,
	Real_T *col, Real_T *rho, Real_T *zlast, Real_T *wlast );
int q5gvns_( Real_T *x, Real_T *y, Real_T *c, Real_T *s, Real_T *p );
int q5refl_( Int_T *n, Real_T *x, Real_T *y, Real_T *c, Real_T *s, Real_T *p );
int q5rexp_( Int_T *m, Real_T *r, Real_T *c );
int q5rpck_( Int_T *m, Real_T *r, Int_T *idel );
int q5tris_( Int_T *m, Real_T *u, Real_T *x, Bool_T *trans );

/*--- Q6 ---------------------------------------------------------------------*/

Real_T dnorm2_( Int_T n, Real_T *a );
Real_T ddot_( Int_T n, Real_T *a, Real_T *b );
Real_T ddots_( Int_T n, Real_T *a, Int_T *irn, Real_T *b );
void dzero_( Int_T n, Real_T *a );
void dmult_( Int_T n, Real_T *a, Real_T *t );
void dcopy_( Int_T n, Real_T *a, Real_T *b );
void dpack_( Int_T n, Real_T *a, Real_T *b, Int_T *irn );
void dunpk_( Int_T n, Real_T *a, Int_T *irn, Real_T *b );
void dunne_( Int_T n, Real_T *a, Int_T *irn, Real_T *b );
void dcopne_( Int_T n, Real_T *a, Real_T *b );
void icopy_( Int_T n, Int_T *a, Int_T *b );
void dsuma_( Int_T n, Real_T *a, Real_T *b );
void ddifr_( Int_T n, Real_T *a, Real_T *b );
void dstep_( Int_T n, Real_T *a, Real_T *d, Real_T tau );
void dstpck_( Int_T n, Real_T *a, Int_T *irn, Real_T *d, Real_T tau );
void dstunp_( Int_T n, Real_T *a, Real_T *d, Int_T *irn, Real_T tau );

//==============================================================================
//
//	End of Declarations / definitions for QDX.
//
//==============================================================================

#endif
