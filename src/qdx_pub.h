/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	qdx_pub.h
CREATED:			1994.07.28
LAST MODIFIED:		1996.03.01

DEPENDENCIES:		stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __QDX_PUB_H__
#define __QDX_PUB_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

class RD_SubproblemManager;

//==============================================================================
//
//	Definitions for QDX Fortran solver (translated to C by f2c converter).
//


void q1cmte_( const Int_T n, Int_T l, Real_T *x, Real_T *y, Real_T *yb,
	Real_T *xmin, Real_T *xmax, Real_T *v, Real_T *weight, Int_T mdmat,
	Real_T *dmat, Int_T *jcol, Int_T *iptr, Real_T *bmin, Real_T *bmax,
	Int_T *marks, Int_T *status, Real_T *g, Real_T *a, Int_T *iblock,
	Int_T *icheck, Int_T *ieq, Int_T *drow,
	Int_T *ibasic, Real_T *pricba, Int_T *inonba, Int_T *irn, 
	Int_T *istat, Real_T *pricnb, Real_T *pi, Real_T *q, Real_T *r, 
	Real_T *z, Real_T *w, Real_T *col, Real_T *dpb, Int_T itmax, 
	Int_T *istop, RD_SubproblemManager &SubMan, Int_T levprt, Real_T initpen,
	Real_T *ExpC, Real_T& ExpCost, Real_T *ExpC2, Real_T& ExpCost2);


//@BEGIN-----------------------------------------------------------------------
//      Added the last line for testing sol. quality
//@END-------------------------------------------------------------------------


//
//	End of definitions for QDX.
//
//==============================================================================

#endif
