/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	determlp.h
CREATED:			1995.07.31
LAST MODIFIED:		1996.03.28

DEPENDENCIES:		stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __DETERMLP_H__
#define __DETERMLP_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif


class Scenario;
class Scenarios;

//==============================================================================
//
//	Class "DeterministicLP" declaration.
//
class DeterministicLP : public MPS_LP
{
private:
	Int_T CostRow;
	Int_T ActualStage2Rows;

	Array<Bool_T> IncludeRows, IncludeCols;

	const char *ObjName;

	Real_T FixedAdjustment;
	Array<Real_T> Shift;

public:
	DeterministicLP( void );
	virtual ~DeterministicLP( void );

	Int_T GetCostRow( void ) const;
	Int_T ScaleObjective( Scenarios &scen );
	Bool_T CheckStructure( Int_T Stage2Row, Int_T Stage2Col );
	const char *RevealObjectiveName( void ) const;

	void ZeroRandomData( const Scenario &Scen );
	void RenumberIndiceInScenarios( Scenarios &Scen, Int_T Stage2Col );

	Int_T DivideIntoStages( MPS_LP &A, MPS_LP &T, MPS_LP &W, Int_T Stage2Col,
		const Scenario &Scen );
	Int_T GetActualStage2Rows( void );

};
//
//	End of class "DeterministicLP" declaration.
//
//==============================================================================


//==============================================================================
//
//	Inline functions
//
inline
DeterministicLP::DeterministicLP( void )
	: CostRow( -1 ), ActualStage2Rows( -1 ), IncludeRows(), IncludeCols(),
	ObjName( NULL ), FixedAdjustment( 0.0 ), Shift()
{}


inline
DeterministicLP::~DeterministicLP( void )
{}


inline
Int_T DeterministicLP::GetActualStage2Rows( void )
{ return ActualStage2Rows;  }

//
//	End of inline functions.
//
//==============================================================================

#endif
