/*------------------------------------------------------------------------------
MODULE TYPE:		File input routine.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	read_tim.h
CREATED:			1994.08.16
LAST MODIFIED:		1994.08.16

DEPENDENCIES:		mps_lp.h, stdtype.h,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	Public interface to time file reader, which does both parsing and performs
the semantic actions (namely reading the numbers of first and second stage
starting rows and columns).

------------------------------------------------------------------------------*/

#ifndef __READ_TIM_H__
#define __READ_TIM_H__

#include <stdio.h>

#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif


//==============================================================================
//
//	Prototypes of functions used for readin the time ('*.tim') file.
//
//==============================================================================

Bool_T ReadTimes( FILE *fp, MPS_LP &LP, Int_T &Stage1Row,
	Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col,
	VerbLevel Verbosity );

//==============================================================================
//
//	End of prototypes.
//
//==============================================================================

#endif
