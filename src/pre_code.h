/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	pre_code.h
CREATED:			1995.10.27
LAST MODIFIED:		1995.10.28

DEPENDENCIES:		none
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains a declaration of all enumerations and definitions
used by the presolve procedure.

------------------------------------------------------------------------------*/

#ifndef __PRE_CODE_H__
#define __PRE_CODE_H__

enum LP_STATUS {
	LPS_UNKNOWN,
	LPS_INFEASIBLE,
	LPS_UNBOUNDED,
	LPS_SOLVED
	};


enum LP_REDUCTIONS {
	LPR_NONE			= 0,
	LPR_EMPTY_ROWS		= 0x0001,
	LPR_EMPTY_COLS		= 0x0002,
	LPR_ORIG_FIXED		= 0x0004,
	LPR_NUM_ELIM		= 0x0008,
	LPR_SINGL_ROWS		= 0x0010,
	LPR_SINGL_COLS		= 0x0020,
	LPR_FORC_DOM_CONSTR	= 0x0040,
	LPR_DOM_COLS		= 0x0080,
	LPR_EXPLICIT_SLACKS	= 0x0100,

	LPR_MIN				= LPR_EMPTY_ROWS | LPR_EMPTY_COLS | LPR_ORIG_FIXED,

	LPR_SIMPLE			= LPR_MIN | LPR_SINGL_ROWS,

	LPR_PRIMAL			= LPR_SIMPLE | LPR_SINGL_COLS | LPR_FORC_DOM_CONSTR
						/* | LPR_EXPLICIT_SLACKS */,

	LPR_DUAL			= LPR_MIN | LPR_DOM_COLS | LPR_EMPTY_COLS,

	LPR_ALL				= LPR_PRIMAL | LPR_DUAL | LPR_NUM_ELIM
						/* | LPR_EXPLICIT_SLACKS */
	};


#define DEFAULT_FEASIBILITY_TOL	(1.0e-8)

#endif
