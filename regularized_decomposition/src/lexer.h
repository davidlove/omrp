/*------------------------------------------------------------------------------
MODULE TYPE:		Simple general purpose lexical analyser.
PROJECT CODE:		-------------
PROJECT FULL NAME:	-------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	lexer.h
CREATED:			1993.10.10
LAST MODIFIED:		1995.10.23

DEPENDENCIES:		stdtype.h,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __LEXER_H__
#define __LEXER_H__

#include <stdio.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

//==============================================================================
//
//	Declaration of class "Lexer". 
//
//==============================================================================

class Lexer
{
private:
	static FILE *Input;

public:
	enum { IO_BUF_LEN = 512 };

	static char FileName[ IO_BUF_LEN ];
	static int LineNumber;
	static void *IO_OK;
	static char Buf[ IO_BUF_LEN ];
	static char *BufPtr;
	static double Number;
	static char *LabelPtr[2];
	static const char EmptyString[1];

public:
	static void SetInputStream( const char *name, FILE *fp );

	static Bool_T GetSpace( void );
	static Bool_T GetKeyword( const char *Keyword,
		Bool_T CaseSensitive = False );
	static Bool_T GetLabelFree( int lbl, int = 0, int = 0 );
	static Bool_T GetLabelFixed( int lbl, int StartCol, int EndCol );
	static Bool_T GetNumeric( void );
	static Bool_T GetNewline( Bool_T Force = False );
};

//==============================================================================
//
//	End of declaration of class "Lexer". 
//
//==============================================================================

#endif
