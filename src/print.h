/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose header
PROJECT CODE:		----------------------
PROJECT FULL NAME:	----------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	----------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	print.h
CREATED:			1995.08.12
LAST MODIFIED:		1995.10.31

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A prototype of a function which adds output buffer flushing to the
functionality of a <stdio.h> function "printf()".

------------------------------------------------------------------------------*/

#ifndef __PRINT_H__
#define __PRINT_H__

//==============================================================================
//
//	A prototype.
//
//==============================================================================

#include <stdio.h>

void SetDefaultPrintOutput( FILE *fp );
void Print( const char *fmt, ... );

#endif
