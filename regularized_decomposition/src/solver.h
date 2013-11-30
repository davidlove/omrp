/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code header.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	solver.h
CREATED:			1992.10.02
LAST MODIFIED:		1996.09.16

DEPENDENCIES:		smartptr.h, stdtype.h, solv_lp.h, inverse.h, solvcode.h,
					parsespc.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	Declaration of class "Solver" - the class representing a penalty-based 
implementation of a revised simplex method for linear programming. The class
can solve a problem, preform a given number of primal simplex iterations,
restart solution process (possibly after modifications to the problem) etc.

------------------------------------------------------------------------------*/

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <stdio.h>
#include <assert.h>

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

//#ifndef __HISTORY_H__
//#	include "history.h"
//#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SOLVCODE_H__
#	include "solvcode.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif


//==============================================================================
//
//	Now some macros.
//

//
//	Macros for setting default starting point and initial basis finding.
//
#define DEFAULT_STARTING_POINT	FindStartingPoint_Basic
#define DEFAULT_INITIAL_BASIS	FindInitialBasis_Bixby
#define EXTERNAL_STARTING_POINT	FindStartingPoint_External
#define EXTERNAL_INITIAL_BASIS	FindInitialBasis_External

#define MAX_FEAS_RESTORE		4

//
//	End of the macros.
//
//==============================================================================

//class History;
class SimplexLP;
class Inverse;

class Solver;

typedef Bool_T (Solver::* FIB)( const Real_T LU_PIVOT_TOL,
	VerbLevel Verbosity );
typedef Bool_T (Solver::* FSP)( VerbLevel Verbosity );


//------------------------------------------------------------------------------
//
//	Class SolverStateDump
//		Used for storing the state of a simplex optimizer class "Solver".
//
class SolverStateDump
{
private:
	Int_T N;				//size of x

	Array<Int_T> A2B;		//Ax = b; but A=[B N] where B is the basic part and 
							//            and N is the nonbasic part
							//A2B gives which columns are basic 

	Array<Real_T> X;		//solution x
	
	SolverStateDump( Int_T n, const Array<Int_T> &a2b, const Array<Real_T> &x,
		Int_T AllocN );		//AllocN = Size of tables which store column-oriented data.

public:
	~SolverStateDump( void );

	friend class Solver;
};

//------------------------------------------------------------------------------
//
//	Class Solver
//		Holds all data needed to solve an LP problem.
//		Method Solve is used to solve the problem.
//
class Solver
{
protected:
	SimplexLP &LP;			// LP class: holds original problem.
	Inverse *B;				// Class inverse handles factorization and
							// matrix equations' solving.
	const Spc &S;			// Configuration database. (specifications)
	Int_T Alen, Ulen;		// Length of factors before and after factorization.
	Real_T Amax, Umax;		// Maximum absolute value of factors before and
							// after factorization.

	Int_T N, M;				// Current problem dimensions.
	const Int_T AllocN;		// Size of tables which store column-oriented data.

	Array<Int_T> A2B, B2A;	// Index tables - tell which variables belong into
							// basis and which positions they occupy.
	Array<Real_T> x,		// Dense vector of primal variables.
		u,					// Upper bounds on variables.
		y,					// Dense vector of dual variables (rarely used).
		z,					// Vector of reduced costs.
		gamma,				// Vector of reduced cost weights (used by the
							// steepest edge algorithm).
		alpha, beta,		// Work vectors for reduced cost updates and SE.
		y_t, y_x;			// Work vectors for split pricing.
	Array<Short_T> VarType;	// Array of variable types description.

	Array<Real_T> w1, w2,	// Two sparse ('w1' and 'w2') and one dense ('w3')
		w3;					// work vectors. Indice for sparse work vectors.
	Array<Int_T> w1Ind, w2Ind;
	Int_T w1Len;

	SOLVE_RESULT Status;
	Real_T Result;			// Objective value (updated on every iteration).
	Real_T PrimalResiduals,
		DualResiduals,
		BoxConstraintViolation;

	//--------------------------------------------------------------------------
	//	Statistics of the 'Solver' object.
	//
	Long_T IterCnt,			// Iteration count. 
		TotalIterCnt,		// Total number of iterations (including restarts).
		RC_FaultCnt,		// Invalid reduced cost count.
		PenaltyAdjustCnt,	// Counter of penalty parameter adjustments.
		ResidCheckCnt,		// Number of times the residuals were checked.
		SE_ResetCnt,		// Number of steepest edge resets.
		PrimVarComputeCnt,	// Number of times the primal var's were computed.
		DualVarComputeCnt,	// -------||---------- dual   -------||----------
		AltPricCnt,			// Number of times the pricas were split.
		InfeasMinCnt;		// Number of invocations of infeasibility
							// minimization routine.

	Int_T ArtifIncrease;	// A counter of artificial variables' increases.

	//--------------------------------------------------------------------------
	//	Iteration history storage and recovery data structures.
	//
//	History H;

	//--------------------------------------------------------------------------
	//	Values for numerical algorithms - copied from Spc instance by SetSpc().
	//
	//	One control variable for assuring, that solver has been initialized.
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
		RESID_CHECK_FREQ,
		COLUMN_REJECT_CNT,
		DEGEN_CNT,
		CYCLE_CNT;

	PricingScheme Pricing;

	Bool_T Initialized;

	Array<Int_T> ExternalA2B;
	Array<Real_T> ExternalX;
	Bool_T UseExternalBasis, UseExternalSolution;

	Int_T FeasibilityRestoreCount, UpdateCnt, ResidCheckClock;

//------------------------------------------------------------------------------
//	Public interface.
//	
public:
	Solver( SimplexLP &l, const Spc &spc );
	~Solver( void );

	SOLVE_RESULT Solve( VerbLevel Verbosity, Long_T IterLimit = -1 );
	SOLVE_RESULT RestartAndSolve( VerbLevel Verbosity, Long_T IterLimit = -1,
		Bool_T ComputeDuals = True );
	SOLVE_RESULT RestartAndSolve( VerbLevel Verbosity,
		const SolverStateDump *dump, Long_T IterLimit = -1 );

	SolverStateDump *GetSolverStateDump( void );

	enum CNT { Iter = 200, RC_Fault, PenaltyAdjust, TotalIter, ResidCheck,
		SE_Reset, PrimVarCompute, DualVarCompute, AltPric, InfeasMin };

	void ResetStatisticCounters( void );
	Long_T ReadStatisticCounter( int n ) const;

	Real_T DualityGap( void ) const;

protected:
	void SetSpc( const Spc &spc );

	SOLVER_ERROR UpdatePrimalVars( Int_T q, Real_T z_q, Real_T theta );
	void UpdateReducedCosts( Int_T p, Int_T pp, Real_T z_q, Real_T gamma_q,
		Real_T pivot );
	Bool_T PeriodicalResidCheck( SOLVER_ERROR &ErrorCode, VerbLevel Verbosity );

	Bool_T FactorizeBasis( void );
	int UpdateBasis( Int_T p = -1 );
	void FixArtificialVariables( void );
	void ComputePrimalVariables( void );
	void ComputeDualVariables( Bool_T CountArtif = True );
	void ComputeDualVariables( Real_T penalty );
	void ComputeSplitDualVariables( void );
	void ComputeResult( void );
	void ComputeReducedCosts( void );
	void ComputeReducedCosts( Int_T Start );
	void ResetWeights( void );
	Real_T MinimumPenaltyToMove( const Bool_T DependsOnArtif );
	Real_T MinimumPenaltyToStayBounded( Int_T j, const Bool_T DependsOnArtif );
	Bool_T VerifyUnboundedness( Int_T j );
	Real_T CheckResiduals( int &Mode );
	void MakeSolutionFeasible( Real_T ResidLevel, VerbLevel Verbosity,
		Bool_T MinimizeInfeasibility = True );
	Bool_T MinimizeLambda( Array<Real_T> &lambda, VerbLevel Verbosity );
	Bool_T MinimizeInfOneVariable( Int_T col, Array<Real_T> &lambda,
		Array<Bool_T> &NonZero, Int_T MaxFillIn, Int_T &InfCnt,
		Real_T &InfNorm );
	Bool_T HandleNumericalDifficulties( SOLVER_ERROR Err, VerbLevel Verbosity );
	Int_T BackTrack( Int_T Steps );

	void InGoingColumnNumber( Int_T &q, VerbLevel Verbosity );
	SLV_FCC FindColumnCandidates( Int_T &qMin );
	void UnMark( void );
	Real_T OutGoingColumnNumber( Int_T q, Real_T z_q, Int_T &p, Real_T &theta,
		Int_T &bound );
	Bool_T NonZeroArtificials( void );
	Bool_T BasicArtificials( void );
	void ToBounds( void );

	void ResetTolerances( VerbLevel Verbosity );
	void HardLPTolerances( VerbLevel Verbosity );
	void EasyLPTolerances( VerbLevel Verbosity );
	void FinalTolerances( VerbLevel Verbosity );

private:
#ifndef NDEBUG
	void CheckA2B_Consistency( void );
#endif

	//--------------------------------------------------------------------------
	//	Function "InitializeSolver()" is called by the user of the class before
	//	the solution may begin. "FSP" and "FIB" arguments are pointers to member
	//	functions responsible for finding an initial basis and the initial
	//	solution.
	//	The other functions are possible choices for initial basis and initial
	//	starting point computation.
	//
public:
	Bool_T InitializeSolver( VerbLevel Verbosity,
		FSP FindStartingPoint = (FSP)NULL,
		FIB FindInitialBasis = (FIB)NULL,
		Bool_T BasisFirst = True );

//	Bool_T SetExternalBasis( const Array<Int_T> &ExternA2B, Int_T Len );
//	void ClearExternalBasis( void );

	Bool_T SetExternalSolution( const Array<Real_T> &ExternX, Int_T Len );
	void ClearExternalSolution( void );

	//
	//	Initial basis (functions corresp. to "FIB" argument in
	//	"InitializeSolver()" above).
	//
	Bool_T FindInitialBasis_Bixby( const Real_T LU_PIVOT_TOL,
		VerbLevel Verbosity );
//	Bool_T FindInitialBasis_External( const Real_T LU_PIVOT_TOL,
//		VerbLevel Verbosity );

	//
	//	Initial basis (functions corresp. to "FSP" argument in
	//	"InitializeSolver()" above).
	//
	Bool_T FindStartingPoint_Basic( VerbLevel Verbosity );
	Bool_T FindStartingPoint_AllOnes( VerbLevel Verbosity );
	Bool_T FindStartingPoint_AllNonBasicOnes( VerbLevel Verbosity );
	Bool_T FindStartingPoint_External( VerbLevel Verbosity );

	//--------------------------------------------------------------------------
public:
	const Array<Real_T>	&GetPrimalVariables( void );
	const Array<Real_T>	&GetDualVariables( void );

	Real_T	GetResult( void );
	Real_T	GetPrimalResiduals( void );
	Real_T	GetDualResiduals( void );
	Real_T	GetBoxConstraintViolation( void );
	Long_T	GetNumberOfIterations( void );
	Long_T	GetNumberOfInvalidRC( void );
	Long_T	GetNumberOfPenaltyAdjust( void );

	Solution * GetSolution( int mask, Bool_T LP_Valid = True );

	void OutputBasisToFile( FILE *fp );
};


//==============================================================================
//
//	Inline definitions of some functions.
//
//==============================================================================


inline
const Array<Real_T> & Solver::GetPrimalVariables( void )
	{ return x; }


inline
const Array<Real_T> & Solver::GetDualVariables( void )
	{ return y; }

inline
Real_T Solver::GetResult( void )
	{ return Result; }


inline
Real_T Solver::GetPrimalResiduals( void )
	{ return PrimalResiduals; }


inline
Real_T Solver::GetDualResiduals( void )
	{ return DualResiduals; }


inline
Real_T Solver::GetBoxConstraintViolation( void )
	{ return BoxConstraintViolation; }


inline
Long_T Solver::GetNumberOfIterations( void )
	{ return IterCnt; }


inline
Long_T Solver::GetNumberOfInvalidRC( void )
	{ return RC_FaultCnt; }

inline
Long_T Solver::GetNumberOfPenaltyAdjust( void )
	{ return PenaltyAdjustCnt; }


inline
void Solver::ClearExternalSolution( void )
	{ UseExternalSolution = False; }



inline
SolverStateDump::~SolverStateDump( void )
{}

#endif
