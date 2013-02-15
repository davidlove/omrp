/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski
--------------------------------------------------------------------------------

HEADER FILE NAME:	rd_sublp.h
CREATED:			1995.07.27
LAST MODIFIED:		1996.02.13

DEPENDENCIES:		smplx_lp.h, scenario.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x
------------------------------------------------------------------------------*/

#ifndef __RD_SUBLP_H__
#define __RD_SUBLP_H__

#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif

//==============================================================================
//
//	Class RD_SubproblemLP declaration.
//
//==============================================================================

class RD_SubproblemLP : public SimplexLP
{
public:
	//--------------------------------------------------------------------------
	//	Some data needed in order to 
	//
	const MPS_LP &T_base;		// The common (base) value of the technology
								// matrix.

private:
	enum { UNINIT, INIT } State;// Object state (set to "INIT" after the
								// following data are initially allocated and
								// filled.

	Array<Real_T> d_base,		// The common (base) value of the right hand
		q_base;					// side/cost vectors for all scenarios.

	Int_T n1st,					// Number of first stage variables.
		m2st;					// Number of second stage constraints.

	//--------------------------------------------------------------------------
	//	Work vectors needed to solve the subproblems.
	//
	Array<Real_T> T_base_y,		// Work vector; calculated once per each scan
								// of scenarios; equals T_base * y.

		h, q;					// Work vectors; calculated anew each time
								// subproblem solution is requested; holds the
								// current right hand side and cost vector.

public:
	RD_SubproblemLP( const MPS_LP &T );
	virtual ~RD_SubproblemLP( void );

	void InitializeRD_Subproblem( void );
	void ApplyScenario( const Scenario &Sc, Bool_T NewTrialPoint, Int_T nn,
		const Real_T *TrialPoint );
};


//==============================================================================
//
//	End of class RD_SubproblemLP declaration.
//
//==============================================================================


inline
RD_SubproblemLP::RD_SubproblemLP( const MPS_LP &T )
	: SimplexLP(), T_base( T ), State( UNINIT ), d_base( b ),
	n1st( T.GetN() ), m2st( T.GetM() ),
	T_base_y( T.GetM(), 0.0 ), h( T.GetM(), 0.0 )
{}

inline
RD_SubproblemLP::~RD_SubproblemLP( void )
{}

#endif
