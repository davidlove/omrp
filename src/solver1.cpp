/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solver1.cpp
CREATED:            1992.09.29
LAST MODIFIED:		1996.09.16

DEPENDENCIES:       smartptr.h, stdtype.h, error.h, solver.h, smplx_lp.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This module contains the whole simplex LP solver driver. Almost all the
routines are called from the 'Solve' function, which constitutes the driver.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/


#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __SOLVCODE_H__
#	include "solvcode.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __INVERSE_H__
#	include "inverse.h"
#endif


/*------------------------------------------------------------------------------

	Solver::Solver( SimplexLP &l, Spc &spc )
	Solver::~Solver( void )

PURPOSE:
	Solver constructor. Allocates memory for all arrays (using array
constructors) and initializes other data. Also copies specification.

	Solver destructor. Deallocates the basis inverse representation.

PARAMETERS:
	SimplexLP &l
		Reference to linear problem.

	const Spc &spc
		Reference to problem specification class object.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Solver::Solver( SimplexLP &l, const Spc &spc )
	: LP( l ), B( NULL ), S( spc ),
	Alen( 0 ), Ulen( 0 ), Amax( 0.0 ), Umax( 0.0 ),
	N( Int_T( LP.GetStructN() + LP.GetSlackN() ) ),
	M( LP.GetM() ), AllocN( Int_T( N + M ) ),
	A2B( AllocN, A2B_UNDEF ), B2A( M, B2A_UNDEF ),

	x( AllocN ), u( AllocN ), y( M ), z( AllocN ), gamma( AllocN ),
	alpha( AllocN ), beta( AllocN ), y_t( M ), y_x( M ), VarType( AllocN ),
	Status( SR_UNINITIALIZED ), Result( 0.0 ),

	PrimalResiduals( 0.0 ), DualResiduals( 0.0 ), BoxConstraintViolation( 0.0 ),

	ArtifIncrease( 0 ),

//	H( LP.GetM() ),

	DEGEN_CNT( 100 ), CYCLE_CNT( 5000 ), Pricing( PRS_RC ),

	Initialized( False ),

	ExternalA2B( LP.GetStructN() ), ExternalX( LP.GetStructN() ),
	UseExternalBasis( False ), UseExternalSolution( False ),

	FeasibilityRestoreCount( 0 ), UpdateCnt( 0 ), ResidCheckClock( 0 )
{
	B = new Inverse( LP.GetM() );
	if( !B )
		FatalError( "Not enough memory to allocate inverse representation." );
	SetSpc( spc );
	ResetStatisticCounters();

	w1.Resize( M );
	w1Ind.Resize( M );
	w2.Resize( M );
	w2Ind.Resize( M );
	w3.Resize( M );

	w1Len = 0;
}


Solver::~Solver( void )
{
	if( B ) delete B;
}


/*------------------------------------------------------------------------------

	void Solver::SetSpc( const Spc &spc )

PURPOSE:
	Copies the initial configuration information from the class Spc.  Creates
data objects dependent on information read from the .spc file.  These are:
columns' lists and history list.

PARAMETERS:
	const Spc *spc
		Reference to an object holding all essential numerical tolerances.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solver::SetSpc( const Spc &spc )
{
	LU_PIVOT_TOL		= spc.LU_PIVOT_TOL;

	FEASIBILITY_TOL		= spc.FEASIBILITY_TOL;

	OPTIMALITY_TOL		= spc.OPTIMALITY_TOL;
	MIN_STEP_LENGTH		= spc.MIN_STEP_LENGTH;
	PIVOT_TOL			= spc.PIVOT_TOL;
	GROWTH_FACTOR		= spc.GROWTH_FACTOR;
	GOOD_RESID			= spc.GOOD_RESID;
	SATISF_RESID		= spc.SATISF_RESID;
	POOR_RESID			= spc.POOR_RESID;
	ALARM_RESID			= spc.ALARM_RESID;

	REFACT_FREQ			= spc.REFACT_FREQ;
	RESID_CHECK_FREQ	= spc.RESID_CHECK_FREQ;
	LENGTH_FACTOR		= spc.LENGTH_FACTOR;

	Pricing				= spc.Pricing;
}


/*------------------------------------------------------------------------------

	Int_T Solver::ReadStatisticCounter( CNT n )
		const

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

Long_T Solver::ReadStatisticCounter( int n )
	const
{
	switch( n )
	{
	case Iter:				return Iter;
	case RC_Fault:			return RC_Fault;
	case PenaltyAdjust:		return PenaltyAdjustCnt;
	case TotalIter:			return TotalIterCnt;
	case ResidCheck:		return ResidCheckCnt;
	case SE_Reset:			return SE_ResetCnt;
	case PrimVarCompute:	return PrimVarComputeCnt;
	case DualVarCompute:	return DualVarComputeCnt;
	case AltPric:			return AltPricCnt;
	case InfeasMin:			return InfeasMinCnt;

	default:
		if( n == Inverse::Refact || n == Inverse::Upd ||
			n == Inverse::SpFTRAN || n == Inverse::DenFTRAN ||
			n == Inverse::SpBTRAN || n == Inverse::DenBTRAN )
			return B->ReadStatisticCounter( n );
#ifndef NDEBUG
		else
			abort();
#endif
	}
	return 0;
}


/*------------------------------------------------------------------------------

	void Solver::ResetStatisticCounters( void )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

void Solver::ResetStatisticCounters( void )
{
	IterCnt = RC_FaultCnt = PenaltyAdjustCnt = TotalIterCnt = ResidCheckCnt =
		SE_ResetCnt = PrimVarComputeCnt = DualVarComputeCnt = AltPricCnt =
		InfeasMinCnt = 0;

	B->ResetStatisticCounters();
}
