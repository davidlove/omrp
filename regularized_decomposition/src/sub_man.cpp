/*------------------------------------------------------------------------------
MODULE TYPE:		Subproblem solution routines.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	sub_man.cpp
CREATED:			1994.12.10
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		sub_man.h, rd_sublp.h, smplx_lp.h, solv_lp.h, mps_lp.h,
					scenario.h, stdtype.h, smartptr.h, print.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __RD_SUBLP_H__
#	include "rd_sublp.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif



//------------------------------------------------------------------------------
//	Explicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )
	template class SmartPointerBase<SolverStateDump *>;
	template class Array<SolverStateDump *>;
	template class Ptr<SolverStateDump *>;

	template SolverStateDump **MALLOC( SolverStateDump **& Table, size_t len );
	template SolverStateDump **REALLOC( SolverStateDump **& Table, size_t len );
	template void FREE( SolverStateDump **& Table );
#endif
//
//------------------------------------------------------------------------------



/*------------------------------------------------------------------------------

	RD_SubproblemManager::RD_SubproblemManager( const SolvableLP &A, 
		RD_SubproblemLP &Sub, const Spc &spc, RestartMode restart )
	RD_SubproblemManager::~RD_SubproblemManager( void )

PURPOSE:
	Constructor: initializes all the linear problem data and the subproblem
solver.
	Destructor: removes the array of solver states which is possibly left over
from the subproblem solution.

PARAMETERS:
	const SolvableLP &A, RD_SubproblemLP &Sub
		The first stage constraints and the second stage constraints resp.

	const Spc &spc
		Subproblem solver configuration structure.

	RestartMode restart
		Subproblem restart mode (one of: RANDOM, SELF, TREE ).

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

RD_SubproblemManager::RD_SubproblemManager( const SolvableLP &A, // )
	RD_SubproblemLP &Sub, const Spc &spc, RestartMode restart )
	: 
	//
	//	The linear problem data and subproblem solver object.
	//
	FirstStageCost( A.GetN(), 0.0 ),
	SubproblemSolver( Sub, spc ), SubproblemLP( Sub ),
	n1st( A.GetN() ), m2st( Sub.GetM() ),

	//
	//	The scenarios.
	//
	Scen( NULL ), NumberOfScenarios( -1 ), SolverState(),

	//
	//	Other data.
	//
	ObjectState( INITIALIZED ), PreviousBlockNumber( FIRST_CALL ),
	Verbosity( V_LOW ), Restart( restart )
{
	for( Int_T n = A.GetN(), j = 0; j < n; j++ )
		FirstStageCost[j] = A.GetC( j );
}


RD_SubproblemManager::~RD_SubproblemManager( void )
{
	for( Int_T i = 0; i < NumberOfScenarios; i++ )
		if( SolverState[i] != NULL )
		{
			delete SolverState[i];
			SolverState[i] = NULL;
		}
}


/*------------------------------------------------------------------------------

	void RD_SubproblemManager::SetScenarios( const Scenarios &sc )

PURPOSE:
	The subproblem manager holds a pointer to a scenario repository object of
class "Scenarios". This function is used to initialize this pointer.

PARAMETERS:
	const Scenarios &sc 
		A reference to an object which represents a scenario repository.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void RD_SubproblemManager::SetScenarios( const Scenarios &sc )
{
	assert( ObjectState == INITIALIZED );

	Scen = &sc;
	NumberOfScenarios = Scen->GetNumberOfScenarios();

	SolverState.Resize( NumberOfScenarios );
	SolverState.Fill( NULL, NumberOfScenarios );

	ObjectState = READY;
}


/*------------------------------------------------------------------------------

	Bool_T RD_SubproblemManager::SolveSubproblem( Int_T block, Int_T yn,
		const Real_T *y, Real_T &val, Real_T *grad )

PURPOSE:
	Solves a linear subproblem by modifying the linear problem according to a
specific scenario and then applying a simplex optimizer.

PARAMETERS:
	Int_T block
		Subproblem number (if equal to the Number of scenarios, then the first
		stage objective is turned into an optimality cut).

	Int_T yn, const Real_T *y
		Dimension and value of the first stage variable vector.

	Real_T &val, Real_T *grad
		On return: feasibility or optimality cut (see below).
		On entry: "grad" should be allocated. Values are irrelevant.

RETURN VALUE:
	"True" if subproblem feasible (and grad/val store an optimality cut) or
"False" otherwise (and grad/val store a feasibility cut).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T RD_SubproblemManager::SolveSubproblem( Int_T block, Int_T yn, // )
	const Real_T *y, Real_T &val, Real_T *grad )
{
	//--------------------------------------------------------------------------
	//	Some static data needed for gathering statistical data from subproblems.
	//
	static Long_T TotalIterCnt	= 0;
	static Int_T NumOptimal		= 0,
		NumInfeasible			= 0,
		SubsSolved				= 0;
	

	//--------------------------------------------------------------------------
	//	Object state check.
	//
	assert( ( PreviousBlockNumber == FIRST_CALL && ObjectState == READY ) ||
		ObjectState == IN_SOLUTION );

	ObjectState = IN_SOLUTION;

	//--------------------------------------------------------------------------
	//	Argument check.
	//
	block--;
	assert( block >= 0 && block <= NumberOfScenarios );
	assert( yn == n1st );
	assert( y != NULL );
	assert( grad != NULL );

	if( block == NumberOfScenarios )
	{
		//----------------------------------------------------------------------
		//	Output statistics from the current iteration.
		//
		if( Verbosity == V_LOW )
			Print( "\t%d subs (%d feas., %d infeas.) in %ld iterations.\n",
				(int)SubsSolved, (int)NumOptimal, (int)NumInfeasible,
				(long)TotalIterCnt );
		TotalIterCnt = 0;
		NumOptimal = NumInfeasible = SubsSolved = 0;

		//----------------------------------------------------------------------
		//	"Solve" first stage objective subproblem.
		//
		PreviousBlockNumber = CALL_AFTER_Y_CHANGED;
		return SolveStage1ObjectiveSubproblem( yn, y, val, grad );
	}

	//--------------------------------------------------------------------------
	//	Modify subproblem with the scenario data and the trial point.
	//
	if( PreviousBlockNumber < 0 )
		//This has "y" changed so recomputes the right hand side T*y
		SubproblemLP.ApplyScenario( (*Scen)[block], True, yn, y );
	else
		//Similarly, "y" has not been changed and and T*y is not computed again 
		SubproblemLP.ApplyScenario( (*Scen)[block], False, yn, y );

	//--------------------------------------------------------------------------
	//
	//	S U B P R O B L E M   S O L U T I O N !
	//
	SOLVE_RESULT sr = SR_UNKNOWN;

	//	If this is the first time a subproblem is solved, then solve the problem
	//	from the very beginning. Otherwise restart the solver from the previous
	//	solution.
	//
	if( PreviousBlockNumber == FIRST_CALL )
	{
		if( Verbosity >= V_LOW )
			Print( "\nSubproblem restart mode: %s.\n",
				( Restart == TREE ) ? "TREE" :
					( Restart == RANDOM ) ? "RANDOM" :
					( Restart == SELF ) ? "SELF" : "????" );

		if( Verbosity >= V_HIGH )
			Print( "    %3s  %10s  %6s  %10s  %10s  %10s  %10s\n",
				"NO", "STATUS", "ITER", "PRIM RESID", "DUAL RESID", "INFEAS",
				"RESULT" );

		if( !SubproblemSolver.InitializeSolver( V_NONE ) )
			FatalError( "Unable to solve the first subproblem." );
		sr = SubproblemSolver.Solve( V_NONE );
	}
	//
	//	If the scenario has been changed, we solve the node at the root of the
	//	tree. We restart from the solution of the same node for the previous
	//	scenario.
	//
	else if( PreviousBlockNumber == CALL_AFTER_Y_CHANGED )
	{
		assert( block == 0 );
		assert( Restart == RANDOM || SolverState[0] != NULL );

		sr = ( Restart == SELF || Restart == TREE ) ?
			SubproblemSolver.RestartAndSolve( V_NONE, SolverState[0] ) :
			SubproblemSolver.RestartAndSolve( V_NONE, -1, True );
	}
	//
	//	Otherwise we're in the middle of the bunch of scenarios.
	//
	else
	{
		switch( Restart )
		{
		case SELF:
			sr = ( SolverState[block] != NULL ) ?
				SubproblemSolver.RestartAndSolve( V_NONE, SolverState[block] ) :
				SubproblemSolver.RestartAndSolve( V_NONE, -1, True );
			break;

		case RANDOM:
			sr = SubproblemSolver.RestartAndSolve( V_NONE, -1, True );
			break;

		case TREE:
			Int_T pred = Scen->PreviousBlockNumber( block );

			assert( pred >= 0 && pred < NumberOfScenarios );
			sr = SubproblemSolver.RestartAndSolve( V_NONE, SolverState[pred] );
			break;
		}
	}

	//--------------------------------------------------------------------------
	//	Store the solver state. If needed, delete the previous state.
	//
	if( Restart == SELF || Restart == TREE )
	{
		if( SolverState[block] != NULL ) delete SolverState[block];
		SolverState[block] = SubproblemSolver.GetSolverStateDump();
		assert( SolverState[block] != NULL );
	}
	PreviousBlockNumber = block;

	if( Verbosity >= V_HIGH )
		Print( "SUB %3d  %10s  %6d  %10.2e  %10.2e  %10.2e  %10.2e\n",
			block,
			(sr == SR_OPTIMUM) ? "OPTIMAL" : (sr == SR_INFEASIBLE) ?
				"INFEASIBLE" : "???",
			(int) SubproblemSolver.GetNumberOfIterations(),
			SubproblemSolver.GetPrimalResiduals(),
			SubproblemSolver.GetDualResiduals(),
			SubproblemSolver.GetBoxConstraintViolation(),
			SubproblemSolver.GetResult()
	);

	SubsSolved++;
	TotalIterCnt += SubproblemSolver.GetNumberOfIterations();

	switch( sr )
	{
	case SR_OPTIMUM:
		NumOptimal++;
		SubproblemSolver.GetOptimalityCut( val, grad, yn, (*Scen)[block] );
		break;

	case SR_INFEASIBLE:
		NumInfeasible++;
		SubproblemSolver.GetFeasibilityCut( val, grad, yn, (*Scen)[block] );
		break;

	case SR_UNBOUNDED:
		FatalError( "Subproblem unbounded!" );
		break;

	default:
#ifndef NDEBUG
		abort();
#endif
		break;
	}

	//--------------------------------------------------------------------------
	//	Return 'True' if the subproblem was feasible, 'False' otherwise.
	//
	return ( sr == SR_OPTIMUM ) ? True : False;
}


/*------------------------------------------------------------------------------

	Bool_T RD_SubproblemManager::SolveStage1ObjectiveSubproblem( Int_T yn,
		const Real_T *y, Real_T &val, Real_T *grad )

PURPOSE:
	Since the first stage objective function is treated just like a subproblem
in the RD method, this function behaves as a solver for this subproblem.
Naturally, the subproblem is always feasible and provides an optimality cut.

PARAMETERS:
	Int_T yn, const Real_T *y
		First stage constraints vector with specified length.
	
	Real_T &val, Real_T *grad
		Optimality cut.

RETURN VALUE:
	Always "True" since the subproblem is always feasible.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T RD_SubproblemManager::SolveStage1ObjectiveSubproblem( Int_T yn, // )
	const Real_T *y, Real_T &val, Real_T *grad )
{
	val = 0.0;

	for( Int_T i = 0; i < yn; i++ )
	{
		grad[i] = FirstStageCost[i];
		val += grad[i] * y[i];
	}

	return True;
}


//@BEGIN------------------------------------------------------------------------
// Purpose:  To Resize the SolverState after new scenarios has been added. 


void RD_SubproblemManager::ReSizeSolverState ( int nos, int InitScen )
{
	SolverState.Resize( nos );

	for (Int_T i = InitScen+1; i <= nos; i++)
		SolverState.FillOne( NULL, i );
}

//@END--------------------------------------------------------------------------
