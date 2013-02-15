/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, simplex type basis construction.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr. Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	crash.h
CREATED:			1994.12.16
LAST MODIFIED:		1994.12.16

DEPENDENCIES:		stdtype.h

--------------------------------------------------------------------------------

HEADER FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __CRASH_H__
#define __CRASH_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif


//==============================================================================
//
//	Local type and symbolic constant definitions.
//
//==============================================================================

enum PREF_SET {
	PS_FREE	= 0,
	PS_NORM	= 1,
	PS_REST	= 2,
	PS_NONE	= 3
};

//------------------------------------------------------------------------------
//--						Structure used in crash
//------------------------------------------------------------------------------
struct PS			// Preference sets
{
	Int_T next,		// Next variable number.
		previous,	// Previous variable number.
		len,		// Constraint matrix column working length.
		origLen;	// Total column length.
	PREF_SET type;	// Variable type.
	Real_T price;	// Variable price (not cost!).
};

//==============================================================================
//
//	End of local type and symbolic constant definitions.
//
//==============================================================================

#endif
