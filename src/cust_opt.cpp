/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	cust_opt.cpp
CREATED:			1995.10.27
LAST MODIFIED:		1995.10.27

DEPENDENCIES:		stdtype.h, smartptr.h, my_defs.h, option.h, cust_opt.h,
					presolve/pre_code.h
					<stdlib.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Static data for classes representing different types of command line options
used in the context of the linear (incl. stochastic) programming.

------------------------------------------------------------------------------*/

#ifndef __CUST_OPT_H__
#	include "cust_opt.h"
#endif


const char *VERBOSITY_Option::labels[5] =
	{ "none", "low", "high", "line", "lineh" };

const int VERBOSITY_Option::states[5] =
	{ V_NONE, V_LOW, V_HIGH, V_LINE, V_LINE_HEAD };

const char *PRESOLVE_MODE_Options::labels[14] =
	{ "on", "all", "off", "none", "min", "simple", "primal", "dual",
	"sr", "fdr", "dc", "fsc", "ne", "es" };

const int PRESOLVE_MODE_Options::states[14] = 
	{ LPR_ALL, LPR_ALL, LPR_NONE, LPR_NONE, LPR_MIN, LPR_SIMPLE, LPR_PRIMAL,
	LPR_DUAL, LPR_SINGL_ROWS, LPR_FORC_DOM_CONSTR, LPR_DOM_COLS,
	LPR_SINGL_COLS, LPR_NUM_ELIM, LPR_EXPLICIT_SLACKS };

const char *PRICING_MODE_Option::labels[3] = { "rc", "se", "ase" };

const int PRICING_MODE_Option::states[3] = { PRS_RC, PRS_SE, PRS_ASE };
