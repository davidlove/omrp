/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	stochsol.h
CREATED:			1996.06.21
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		solution.h,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __STOCHSOL_H__
#define __STOCHSOL_H__

#ifndef __SOLUTION_H__
#	include "solution.h"
#endif

//==============================================================================
//
//	Class "StochSolution" declaration.
//
//==============================================================================

class StochSolution : public Solution
{
public:
	enum StochSolutionMask {
		Weights	= 0x40,

		PrimalWeights = Primal | Res | Weights
	};

public:
	Array<Real_T> weights, f;
	Int_T Scenarios;

public:
	StochSolution( Int_T scen, Int_T nn );
	virtual ~StochSolution( void );

	virtual void WriteText( FILE *fp, int cntnts ) const;
};


inline
StochSolution::StochSolution( Int_T scen, Int_T nn )
	: Solution( 0, nn, PrimalWeights ), Scenarios( scen )
{
	assert( nn > 0 );
	assert( scen > 0 );
	weights.Resize( scen );
	f.Resize( scen );
}


inline
StochSolution::~StochSolution( void )
{}


//==============================================================================
//
//	End of class "StochSolution" declaration.
//
//==============================================================================

#endif
