/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines - header file
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method for large
					scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. A. P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	invaux.h
CREATED:			1993.10.31
LAST MODIFIED:		1996.02.06

DEPENDENCIES:		smartdcl.h, stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __INVAUX_H__
#define __INVAUX_H__

#ifndef __SMARTDCL_H__
#	include "smartdcl.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif


#define FTRANL_LABEL ("FTRANL-INTERMEDIATE")

//==============================================================================
//
//	Auxiliary functions' prototypes.
//
//==============================================================================

void Pack( Array<Real_T> &a, Array<Int_T> &ind, Array<Int_T> &ptr, Int_T n,
	Array<Int_T> &len, Bool_T reals, Int_T &lfile );

void Reorder( Int_T n, Int_T U_len, Array<Real_T> &a, Array<Int_T> &jcol, // )
	Array<Int_T> &rptr, Array<Int_T> &irow );

//==============================================================================
//
//	End of auxiliary functions' prototypes.
//
//==============================================================================

#endif
