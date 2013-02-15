/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code header.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	solvtol.h
CREATED:			1993.11.14
LAST MODIFIED:		1996.04.12

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains a number of numerical tolerances for use by a
simplex type linear optimizer.

------------------------------------------------------------------------------*/

#ifndef __SOLVTOL_H__
#define __SOLVTOL_H__

//==============================================================================
//
//	Defaults, lower and upper bounds on settings read in 'spc.cc'.
//
//==============================================================================

#define SMALL_ELEM_LO			(1.0e-12)
#define SMALL_ELEM_HI			(1.0e-8)
#define SMALL_ELEM_DEF			(1.0e-10)

#define LU_PIVOT_TOL_LO			(1.0e-2)
#define LU_PIVOT_TOL_HI			(0.9999)
#define LU_PIVOT_TOL_DEF		(1.0e-1)

#define FEASIBILITY_TOL_LO		(1.0e-6)
#define FEASIBILITY_TOL_HI		(1.0e-7)
#define FEASIBILITY_TOL_DEF		(1.0e-7)

#define OPTIMALITY_TOL_LO		(1.0e-8)
#define OPTIMALITY_TOL_HI		(1.0e-6)
#define OPTIMALITY_TOL_DEF		(1.0e-8)

#define MIN_STEP_LENGTH_DEF		(1.0e-10)

#define PIVOT_TOL_LO			(1.0e-8)
#define PIVOT_TOL_HI			(1.0e-6)
#define PIVOT_TOL_DEF			(1.0e-6)

#define GROWTH_FACTOR_DEF		(1.0e+3)
#define LENGTH_FACTOR_DEF		(2.0)

#define GOOD_RESID_DEF			(1.0e-10)
#define SATISF_RESID_DEF		(1.0e-6)
#define POOR_RESID_DEF			(1.0e-1)
#define ALARM_RESID_DEF			(1.0e+1)

#define REFACT_FREQ_LO			(5)
#define REFACT_FREQ_HI			(200)
#define REFACT_FREQ_DEF			(100)

#define RESID_CHECK_FREQ_LO		(5)
#define RESID_CHECK_FREQ_HI		(500)
#define RESID_CHECK_FREQ_DEF	(100)

//==============================================================================
//
//	End of defaults, lower and upper bounds on settings read in 'spc.cc'.
//
//==============================================================================

#endif
