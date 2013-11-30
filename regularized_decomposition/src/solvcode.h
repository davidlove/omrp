/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	solvcode.h
CREATED:			1993.11.01
LAST MODIFIED:		1995.08.27

DEPENDENCIES:		None

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SOLVCODE_H__
#define __SOLVCODE_H__

//==============================================================================
//
//	Enumerations for simplex-type solver.
//
//==============================================================================

//------------------------------------------------------------------------------
//	These are the codes used in simplex solver hash tables.
//
enum
{
	A2B_LO		= -1,				// Standard codes (used by the solver)
	A2B_UP		= -2,				// which denote the position of non-basic
	A2B_IN		= -3,				// variables: on their lower or upper
	A2B_UNDEF	= -4,				// bounds, or between the bounds.

	A2B_BASIC	= 1,				// Used only for basis construction.

	B2A_UNDEF	= -1				// 
};


//------------------------------------------------------------------------------
//	Possible return codes of the simplex solver are eumerated here. 'Solve'
//	procedure returns one of these values.
//
enum SOLVE_RESULT
{
	SR_UNKNOWN,					// Used for initializing variables
	SR_UNINITIALIZED,
	SR_OK,

	//--------------------------------------------------------------------------
	//	Problem successfuly solved.
	//
	SR_OPTIMUM,					// Optimal solution found.
	SR_INFEASIBLE,				// Primal problem is infeasible.
	SR_UNBOUNDED,				// Dual problem is infeasible (primal is
								// unbounded).

	//--------------------------------------------------------------------------
	//	Problem not solved.
	//
	SR_CANNOT_FIND_INIT_SOL,	// Cannot find initial solution.
	SR_CANNOT_BACKTRACK,		// More backtracking impossible - end of list.
	SR_STALLED,					// Stalling occured.
	SR_RUNNING					// Iteration limit reached.
};


//------------------------------------------------------------------------------
//	These constants are used as bit masks (possibly OR'ed together) denoting
//	the type of residuals' checks to be performed by 'CheckResiduals' function.
//
enum
{
	CHK_INF		= 0x01,
	CHK_PRIM	= 0x02,
	CHK_DUAL	= 0x04,
	CHK_ALL		= CHK_INF | CHK_PRIM | CHK_DUAL
};


//------------------------------------------------------------------------------
//	Enumeration of codes that the pricing function ('FindColumnCandidates') may
//	return.
//
enum SLV_FCC
{
	//--------------------------------------------------------------------------
	//	Pricing procedure result.
	//
	SLV_OPTIMUM				= -1,
	SLV_CANDIDATE_FOUND		= -2,
	SLV_PROBLEM_INFEASIBLE	= -3
};


//------------------------------------------------------------------------------
//	Other symbolic constants used by the solver (in numerous functions) are
//	listed below. The enumerations are divided into subsets corresponding to
//	their meaning.
//
enum
{
	//--------------------------------------------------------------------------
	//	Pivoting procedure results.
	//
	SLV_UNBOUNDED			= -1,
	SLV_TO_BND				= -2,

	//--------------------------------------------------------------------------
	//	Basis update procedure result.
	//
	SLV_BASIS_UPDATED		= 1,
	SLV_BASIS_REFACTORIZED	= 2,
	SLV_UPDATE_FAILED		= 3,

	//--------------------------------------------------------------------------
	//	History handling functions codes.
	//
	SLV_HIST_LIST_EMPTY		= -1
};

enum SOLVER_ERROR
{
	SE_OK,
	SE_DUAL_RESID_POOR,
	SE_DUAL_RESID_ALARM,
	SE_INF_RESID_POOR,
	SE_INF_RESID_ALARM,
	SE_PRIM_RESID_POOR,
	SE_PRIM_RESID_ALARM,
	SE_SINGULAR_BASIS
};


//==============================================================================
//
//	End of enumerations for simplex-type solver.
//
//==============================================================================

#endif
