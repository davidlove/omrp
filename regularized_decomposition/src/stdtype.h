/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	stdtype.h
CREATED:			1993.09.16
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains the most commonly used type definitions. Among
these: definitions of integer (Int_T), long (Long_T), short (Short_T), real
(Real_T) and boolean (Bool_T) types that affect almost all source files.
	There were several machine-dependent sections of type definitions. Now
the types are the same for all machines, but we still check whether the
an appropriate machine type is provided.

	NULL is defined. Also enumerated constants are declared for boolean "True"
and "False" values. Symbolic name INFINITY is defined as a very large double
number.

------------------------------------------------------------------------------*/

#ifndef __STDTYPE_H__
#define __STDTYPE_H__

//==============================================================================
//
//	The most general constants' definitions.
//
//==============================================================================

#ifdef NULL
#	undef NULL
#endif
#define NULL		0

#ifdef INFINITY
#	undef INFINITY
#endif
#define INFINITY	(1.0e+30)

//	End of the general definitions.
//
//==============================================================================

//==============================================================================
//
//	Data type definitions:
//
//	1.	Short_T is used to hold signed data, that is never longer than 8 bits.
//		RECOMMENDED: signed 8 bit or longer integral type.
//
//	2.	Int_T is used for all indexing operations. Thus its range limits the
//		maximum number of rows, columns, non-zeros etc. It HAS TO be a signed
//		type. Has to be a signed type, as sometimes negative numbbers are used
//		to mark something or other. As it is indexing type it should not
//		be longer than 'size_t' type.
//		RECOMMENDED: at least 16 bits signed integral type.
//
//	3.	Long_T is used only occasionally to e.g., hold counters, which may take
//		very large values. Has to be a signed type.
//		RECOMMENDED: at least 32 bits signed floating point type.
//
//	4.	Real_T is a type for all floating point arithmetics.
//		RECOMMENDED: at least 64 bits signed floating point type.
//
//	5.	Bool_T is a type for logical data. It is actually independent of the
//		machine used. It should always be an enumeration - thus an 'int'.
//

enum Bool_T { False = 0, True = 1 };
typedef signed char		Short_T;	// 8 bits
typedef int				Int_T;		// 32 bits
typedef long			Long_T;		// 32 bits
typedef double			Real_T;		// 64 bits

#endif
