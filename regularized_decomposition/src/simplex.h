/*------------------------------------------------------------------------------
MODULE TYPE:		Simplex algorithm global data declaration.
PROJECT CODE:		Simplex.
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	simplex.h
CREATED:			1993.10.07
LAST MODIFIED:		1995.10.27

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains global declarations and definitions for an
implementation of the simplex method.

------------------------------------------------------------------------------*/

#ifndef __SIMPLEX_H__
#define __SIMPLEX_H__

//==============================================================================
//
//	Beginning of global definitions for simplex.
//
//==============================================================================

//------------------------------------------------------------------------------
//	Enumerated data type and external (global) variable declaration. The
//	variable holds the value corresponding to preferred verbosity level.
//
enum VerbLevel
{
	V_NONE = 0,
	V_LINE,
	V_LINE_HEAD,
	V_LOW,
	V_HIGH
};

//------------------------------------------------------------------------------
//	This definition limits length of variable and constraint labels, as well as
//	lengths of right hand side (RHS) vector, range vector and bounds vectors'
//	names. Length of LP problem name is also assumed to be the same.
//
//	The value of 8 conforms to IBM's MPSX user manual, where MPS file format
//	used in that package is described.
//
#define LAB_LEN	((unsigned)8)

//------------------------------------------------------------------------------
//	Pricing scheme enumeration is presented here. It is used by a simplex
//	optimizer.
//
enum PricingScheme {
	PRS_RC,
	PRS_SE,
	PRS_ASE
};

//==============================================================================
//
//	End of global definitions for simplex.
//
//==============================================================================

#endif
