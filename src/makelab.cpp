/*------------------------------------------------------------------------------
MODULE TYPE:		Tool-specific header
PROJECT CODE:		SCE-to-STO converter
PROJECT FULL NAME:	Multistage stochastic problem converter

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. A. Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	makelab.cpp
CREATED:			1996.02.06
LAST MODIFIED:		1996.03.19

DEPENDENCIES:		makelab.h, stdtype.h, error.h
					<assert.h>, <stdio.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	char *MakeLabel( Int_T j, Int_T scen, char *dst )

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __MAKELAB_H__
#	include "makelab.h"
#endif


/*------------------------------------------------------------------------------

	char *MakeLabel( Int_T j, Int_T scen, char *dst )

PURPOSE:
	This function is used for label generation for the deterministic of
two-stage equivalen problem output routine. For each triple of input arguments
a unique label is constructed.
	WARNING: unless an alternative destination (argument "dst") is given, the
label is created in static storage and is overwritten by the next call.

PARAMETERS:
	Int_T j, Int_T scen
		Row/column number and scenario number.

RETURN VALUE:
	Generated label.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

char *MakeLabel( Int_T j, Int_T scen, char *dst )
{
	static char l[9];

	assert( j >= 0 );

	if( dst == NULL ) dst = l;

	if( j > 0xffff || scen > 0xffff )
		FatalError( "Problem too big for the current labeling scheme." );

	sprintf( dst, "%04x%04x", (int)scen, (int)j );
	return dst;
}
