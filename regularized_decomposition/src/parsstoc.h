/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data file parser.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	parsstoc.h
CREATED:			1994.08.17
LAST MODIFIED:		1995.07.23

DEPENDENCIES:		lexer.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER FILE CONTENTS:
	Public interface of the stochastic data file parser. Contains one reader
function and a set of functions which set up semantic actions.

------------------------------------------------------------------------------*/

#ifndef __PARSSTOC_H__
#define __PARSSTOC_H__

#include <stdio.h>


//==============================================================================
//
//	Publicly available functions of the parser - prototypes.
//	Type definitions for pointers to functions.
//
//==============================================================================

class Scenarios;

Bool_T GetStochFile( const char *FileName, FILE *fp, Scenarios *sc );

//==============================================================================
//
//	End of prototypes.
//
//==============================================================================

#endif
