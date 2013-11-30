/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - program configuration file reading.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	parsespc.h
CREATED:			1992.09.29
LAST MODIFIED:		1995.10.27

DEPENDENCIES:		compile.h, simplex.h, stdtype.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header contains class "Spc" declaration. Object of this class holds
all runtime configuration information of the linear programming simplex solver.

------------------------------------------------------------------------------*/

#ifndef __PARSESPC_H__
#define __PARSESPC_H__

#ifndef __COMPILE_H__
#	include "compile.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif


//==============================================================================
//
//	Class "Spc" declaration.
//
//==============================================================================

class Spc
{
public:
	//--------------------------------------------------------------------------
	//	File names.
	//
	//	NOTE: Always allocate memory for these fields. Use 'malloc' type
	//	functions or 'strdup', but never 'new'. The memory will be freed
	//	upon object destruction (with calls to 'free').
	//
	char *MPS_File;
	char *OutputFile;
	char *SolutionFile;

	//--------------------------------------------------------------------------
	//	Verbosity level.
	//
	VerbLevel Verbosity;

	//--------------------------------------------------------------------------
	//	Values for numerical algorithms.
	//
	Real_T LU_PIVOT_TOL,
		FEASIBILITY_TOL,
		OPTIMALITY_TOL,
		MIN_STEP_LENGTH,
		PIVOT_TOL,
		GROWTH_FACTOR,
		LENGTH_FACTOR,
		GOOD_RESID,
		SATISF_RESID,
		POOR_RESID,
		ALARM_RESID;

	Int_T REFACT_FREQ,
		RESID_CHECK_FREQ;

	//--------------------------------------------------------------------------
	//	Some configuration data.
	//
	PricingScheme Pricing;

public:
	Spc( void );			//	Sets reasonable defaults for all data members.
	~Spc( void );			//	Deallocates memory (if necessary).
};

//==============================================================================
//
//	End of class "Spc" declaration.
//
//==============================================================================

#endif
