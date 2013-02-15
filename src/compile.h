/*------------------------------------------------------------------------------
MODULE TYPE:		Compile-time configuration file.
PROJECT CODE:		-------------------
PROJECT FULL NAME:	-------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	compile.h
CREATED:			1993.09.11
LAST MODIFIED:		modified to reconfigure compilation

DEPENDENCIES:		no local dependencies
					no global dependencies

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains a set of macro definitions that are supposed to
configure at compile time the contents of some (or all) source files of the
project in current directory. Some of the definitions are "standard" (i.e. they
are expected to be defined here by some library modules), others are specific
only to this project.  Most projects require existence of the "compile.h" file
in the current directory (it might be empty though). The macros not covered in
this introductory part are described when they are defined.
	The last part of the header (called "FIXES") is meant to perform macro
definitions' completness and consistency check-up. If, for example, you decide,
that in your project only one of several macros from this file may be defined at
a time, than you should put a piece of macro-code in that last section that
would generate a compile-time error (using "#error" directive) if more than one
happens to be actually defined.

------------------------------------------------------------------------------*/

#ifndef __COMPILE_H__
#define __COMPILE_H__


//==============================================================================
//
//	Module specific macro definitions.
//
//==============================================================================

//------------------------------------------------------------------------------
//	File "error.cc"
//
#define COMP_ERROR_ABORT			// Call 'abort()' on fatal error exit.

//------------------------------------------------------------------------------
//	File "read_lp.cc"
//
#define COMP_READLP_1RHS
#define COMP_READLP_1RANGE
#define COMP_READLP_1BOUND

//------------------------------------------------------------------------------
//	File "solv_lp.cc"
//
#define COMP_TOSTD_LOBND	(1.0e-10)
#define COMP_TOSTD_UPBND	(1.0e+10)


//==============================================================================
//
//	End of module specific macro definitions.
//
//==============================================================================

#endif
