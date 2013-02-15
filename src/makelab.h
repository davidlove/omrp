/*------------------------------------------------------------------------------
MODULE TYPE:		Tool-specific header
PROJECT CODE:		SCE-to-STO converter
PROJECT FULL NAME:	Multistage stochastic problem converter

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. A. Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	makelab.h
CREATED:			1996.02.06
LAST MODIFIED:		1996.02.06

DEPENDENCIES:		stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __MAKELAB_H__
#define __MAKELAB_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

//==============================================================================
//
//	Function prototype
//

char *MakeLabel( Int_T j, Int_T scen, char *dst = NULL );

//==============================================================================

#endif

