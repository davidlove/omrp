/*------------------------------------------------------------------------------
MODULE TYPE:		Set of functions.
PROJECT CODE:		General purpose module.
PROJECT FULL NAME:	-----------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	error.h
CREATED:			1993.09.11
LAST MODIFIED:		1993.09.11

DEPENDENCIES:		no local dependencies,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains definitions and prototypes for a set of
error message output functions. The module provides some simple error
message output procedures. There are 3 distinct categories of errors:
*	warnings,
*	ordinary errors and
*	fatal errors.

	This is a list of all publicly available functions of the module
(which prototypes are included in this header file):
	void DeclareErrorMessages( int num, const char const * const * messages )
	void SetErrorOutputStream( FILE *fp = stderr )
	int Warning( int error_code, ... )
	int Warning( const char *format, ...)
	int Error( int error_code, ... )
	int Error( const char *format, ...)
	void FatalError( int error_code, ... )
	void FatalError( const char *format, ...)
	int WarningCount( void )
	int ErrorCount( void )
	void ResetWarningCount( void )
	void ResetErrorCount( void )
	void MaxWarnCount( int num )
	void MaxErrCount( int num )
	
------------------------------------------------------------------------------*/

#ifndef __ERROR_H__
#define __ERROR_H__

#include <stdio.h>

//==============================================================================
//
//	Function prototypes
//
//==============================================================================

//
//	Module configuration functions.
//
void DeclareErrorMessages( int num, const char ** messages );
void SetErrorOutputStream( FILE *fp = stderr );

//
//	Actual error message output.
//
int Warning( int error_code, ... );
int Warning( const char *format, ...);
int Error( int error_code, ... );
int Error( const char *format, ...);
void FatalError( int error_code, ... );
void FatalError( const char *format, ...);

//
//	Reading and reseting the warning/error counters; Reacting to exceedingly
//	large numbers of errors.
//
int WarningCount( void );
int ErrorCount( void );

void ResetWarningCount( void );
void ResetErrorCount( void );

void MaxWarnCount( int num = 0 );
void MaxErrCount( int num = 0 );

//==============================================================================
//
//	End of function prototypes
//
//==============================================================================

#endif
