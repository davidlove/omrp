/*------------------------------------------------------------------------------
MODULE TYPE:		Subproblem solution routines.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	rd_solv.h
CREATED:			1994.12.07
LAST MODIFIED:		1996.02.13

DEPENDENCIES:		stdtype.h, smartptr.h, rd_sublp.h, smplx_lp.h, solv_lp.h,
					mps_lp.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __RD_SOLV_H__
#define __RD_SOLV_H__

#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __RD_SUBLP_H__
#	include "rd_sublp.h"
#endif


class RD_SubproblemLP;
class Scenario;

//------------------------------------------------------------------------------
//
//	Declaration of the "RD_SubproblemSolver" class.
//
//------------------------------------------------------------------------------

class RD_SubproblemSolver : public Solver
{
private:
	RD_SubproblemLP &LP;

public:
	RD_SubproblemSolver( RD_SubproblemLP &LP, const Spc &spc );
	~RD_SubproblemSolver( void );

	void GetFeasibilityCut( Real_T &value, Real_T *grad, Int_T n,
		const Scenario &Scen );
	void GetOptimalityCut( Real_T &value, Real_T *grad, Int_T n,
		const Scenario &Scen );

private:
	void CalculateGradient( Real_T *grad, Int_T n, Array<Real_T> &pi,
		const Scenario &Scen ) const;
};


//------------------------------------------------------------------------------
//
//	Inline definitions of some functions.
//
//------------------------------------------------------------------------------

inline
RD_SubproblemSolver::RD_SubproblemSolver( RD_SubproblemLP &lp, const Spc &spc )
	: LP( lp ), Solver( lp, spc )
{}

inline
RD_SubproblemSolver::~RD_SubproblemSolver( void )
{}


#endif
