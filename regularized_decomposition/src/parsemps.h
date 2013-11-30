/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	parsemps.h
CREATED:			1993.09.21
LAST MODIFIED:		1993.09.22

DEPENDENCIES:		stdtype.h,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains prototypes of functions which constitute an
interface to a parser for fixed and free MPS file. These functions read an
input file from a stream declared with "SetInputStream" call. The type of file
(one of: free MPS, fixed MPS and our proprietary binary format) is recognized
automatically.
	Then a pointer to an MPS_LP object in which data read from the file is to
be stored must be passed. Finally the file type has to be recognized (by a call
to "GetNameLine"). Text file may be read by "GetMPS_Body".
	Binary file has to be processed elsewhere.

------------------------------------------------------------------------------*/

#ifndef __PARSEMPS_H__
#define __PARSEMPS_H__

#include <stdio.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif

//==============================================================================
//
//	Publicly accessible functions' prototypes.
//
//==============================================================================

void SetInputStream( const char *name, FILE *fp );
void SetLP_TargetObject( MPS_LP *_lp );
FF GetNameLine( void );
Bool_T GetMPS_Body( void );

//==============================================================================
//
//	End of publicly accessible functions' prototypes.
//
//==============================================================================

#endif
