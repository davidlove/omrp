/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose header
PROJECT CODE:		----------------------
PROJECT FULL NAME:	----------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	----------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	print.cpp
CREATED:			1995.08.12
LAST MODIFIED:		1995.10.31

DEPENDENCIES:		print.h
					<stdio.h>, <assert.h>, <stdarg.h>

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void Print( const char *fmt, ... );

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	NDEBUG			- as in <assert.h>

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <assert.h>
#include <stdarg.h>

#ifndef __PRINT_H__
#	include "print.h"
#endif


static FILE *output = stdout;


/*------------------------------------------------------------------------------

	void Print( const char *fmt, ... );

PURPOSE:
	A function which adds output buffer flushing to the functionality of a
<stdio.h> function "printf()". Usage identical to that of printf.

PARAMETERS:
	Same as "printf()".

RETURN VALUE:
	None.

SIDE EFFECTS:
	Standard output buffer is flushed.

------------------------------------------------------------------------------*/

void Print( const char *fmt, ... )
{
	if( output == NULL ) return;

	va_list params;
	va_start( params, fmt );
	vfprintf( output, fmt, params );
	va_end( params );

	fflush( output );
}


void SetDefaultPrintOutput( FILE *fp )
{
	output = fp;
}
