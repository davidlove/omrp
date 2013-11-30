/*------------------------------------------------------------------------------
MODULE TYPE:		Subproblem solution routines.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	sub_man.h
CREATED:			1994.12.07
LAST MODIFIED:		1996.02.13

DEPENDENCIES:		stdtype.h, smartptr.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SUB_MAN_H__
#define __SUB_MAN_H__


#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __RD_SOLV_H__
#	include "rd_solv.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif


//==============================================================================
//
//	Class "SubproblemManager" declaration.
//
//==============================================================================

class SolvableLP;
class RD_SubproblemLP;
class SolverStateDump;

class RD_SubproblemManager
{
private:
	//--------------------------------------------------------------------------
	//	The linear problem data and subproblem solver object.
	//
	Array<Real_T> FirstStageCost;	// The first stage problem cost Vector.

	RD_SubproblemSolver SubproblemSolver;
									// Solver object responsible for solving all
									// the subproblems.

	RD_SubproblemLP &SubproblemLP;	// Subproblem LP data structure.

	Int_T n1st,						// Number of the first stage variables.
		m2st;						// Number of the second stage constraints.

	//--------------------------------------------------------------------------
	//	The scenarios.
	//
	const Scenarios *Scen;			// Scenario repository. Each scenario is a
									// list of changes to 'T_base', 'd_base'
									// and 'q_base'.
	Int_T NumberOfScenarios;		// Number of scenarios.

	Array<SolverStateDump *> SolverState;
									// Each subproblem solution will be
									// restarted from a solution of some earlier
									// subproblem. The solver states are stored
									// here.

	//--------------------------------------------------------------------------
	//	Other data.
	//
	enum ObjState { EMPTY, INITIALIZED, READY, IN_SOLUTION } ObjectState;
	enum { FIRST_CALL = -2, CALL_AFTER_Y_CHANGED = -1 };

	Int_T PreviousBlockNumber;

	VerbLevel Verbosity;

public:
	enum RestartMode { TREE, SELF, RANDOM };

private:
	RestartMode Restart;

public:
	RD_SubproblemManager( const SolvableLP &A, RD_SubproblemLP &Sub,
		const Spc &spc, RestartMode Restart = SELF );
	~RD_SubproblemManager( void );

	void SetScenarios( const Scenarios &sc );
	void SetVerbosity( VerbLevel v );

	Bool_T SolveSubproblem( Int_T block, Int_T yn, const Real_T *y, Real_T &val,
    	Real_T *grad );

	RD_SubproblemSolver &GetSubproblemSolver( void );

	//@BEGIN----------------------------------------
	void SetNumOfScenarios( int s ); 
	void ChangePreviousBlockNumber ( int pbn );
	void ReSizeSolverState ( int nos, int InitScen );

	void SetObjState ( void ); 
	//@END------------------------------------------

private:
	Bool_T SolveStage1ObjectiveSubproblem( Int_T yn, const Real_T *y,
		Real_T &val, Real_T *grad );
};

//==============================================================================
//
//	End of class "SubproblemManager" declaration.
//
//==============================================================================

//==============================================================================
//
//	Inline function's definitions.
//
//==============================================================================

inline
RD_SubproblemSolver &RD_SubproblemManager::GetSubproblemSolver( void )
{ return SubproblemSolver; }


inline
void RD_SubproblemManager::SetVerbosity( VerbLevel v )
{ Verbosity = v; }

//@BEGIN------------------
inline 
void RD_SubproblemManager::SetNumOfScenarios( int s )
{ NumberOfScenarios = s; }

inline
void RD_SubproblemManager::ChangePreviousBlockNumber ( int pbn )
{ PreviousBlockNumber = pbn; }

inline 
void RD_SubproblemManager::SetObjState ( void )
{ ObjectState = INITIALIZED; }

//@END---------------------


#endif
