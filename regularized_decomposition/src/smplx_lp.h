/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	smplx_lp.h
CREATED:			1994.03.11
LAST MODIFIED:		1995.10.02

DEPENDENCIES:		smartptr.h, solv_lp.h,

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/


#ifndef __SMPLX_LP_H__
#define __SMPLX_LP_H__

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif


//==============================================================================
//
//	Class "SimplexLP" declaration.
//
//==============================================================================

class SimplexLP : public SolvableLP
{
private:
	Bool_T Standard;		// Boolean value stating wether the LP was converted
							// to standard form.
	Array<Real_T> SlackU;	// Upper bounds on slacks (meaningful when RANGES
							// occur).
	Array<Real_T> Slack,	// Slacks and artificial variables.
		Lambda;				// Super artificial variables (always 'm' of those,
							// but some may be fixed at 0.0).
	Real_T M;				// Penalty value for all artificial variables.
	Int_T SlackLen;			// Number of (non-artificial) slacks added.
	Array<Int_T> SlackRow;	// Rows of consecutive slack variables.
	Array<Short_T> SlackVT,	// Slack variable type (bounded in range rows,
							// normal otherwise).
		LambdaVT;			// Variable type for the lambda vector.
	Array<Int_T> LambdaRow;	// Needed only to pass information through the
							// 'GetColumn' function. Almost constant.

	Array<Int_T> SlackPosByRows, SlackColByRows;
							// Helps in fast search for slack columns when
							// performing row access.

public:
	SimplexLP( void );
	virtual ~SimplexLP( void );

private:
	void ComputePenalty( void );
	void ShiftLowerBoundsToZero( void );

public:
	void GetPenaltyEstimates( Real_T &MinM, Real_T &MaxM );

	void ToStandard( VerbLevel Verbosity );
	void UndoStandard( void );
	virtual void ProcessSolution( Solution &sol );

	//--------------------------------------------------------------------------
	//	Some functios for infeasibility handling.
	//
	//	The first function accepts a lambda vector (a sparse vector with
	//	implicit lower bound of 0.0 and infinite upper bound) and stores it in
	//	data structures. The old lambda vector is not preserved.
	//
	//	The second one returns a reference to the current infeasibility vector.
	//
	//	The last two functions handle the penalty value.
	//
	void CreateLambda( const Array<Real_T> &v );
	const Array<Real_T> &GetLambda( void ) const;

	Real_T GetPenalty( void ) const;
	void SetPenalty( const Real_T NewPenalty );
	void FixLambda( Int_T j );

	//--------------------------------------------------------------------------
	//	This function constructs a simplex-type initial basis. It is a crash
	//	algorithm (for further details - see "crash.cpp").
	//
	void InitialBasis( Real_T PivotTol, Array<Int_T> &A2B,
		VerbLevel Verbosity );

	//--------------------------------------------------------------------------
	//	Inline implementations of the below listed functions - see the end of
	//	the header file. These are all virtual functions (the virtuality is
	//	inherited from "MPS_LP" ans "SolvableLP" classes).
	//
	//	Only the two last functionsare not inherited (and not virtual).
	//
	Int_T	GetN( void )				const;
	Int_T	GetNZ( void )				const;
	Real_T	GetL( Int_T j )				const;
	Real_T	GetU( Int_T j )				const;
	Real_T	GetC( Int_T j )				const;
	Short_T	GetVarType( Int_T j )		const;

	Int_T GetStructN( void )			const;
	Int_T GetSlackN( void )				const;

	void GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row, Int_T &len )
		const;
	void GetRow( Int_T Row, Ptr<Real_T> &A, Ptr<Int_T> &Col, Int_T &len,
		Int_T &SlackCol, Real_T &Slack, Int_T &LambdaCol, Real_T &Lambda )
		const;
};


//==============================================================================
//
//	End of class "SimplexLP" declaration.
//
//==============================================================================


//==============================================================================
//
//	Inline definitions of some functions.
//


//------------------------------------------------------------------------------
inline
SimplexLP::~SimplexLP( void )
{}


//------------------------------------------------------------------------------
inline
Int_T SimplexLP::GetN( void )
const
	{ return Int_T( MPS_LP::GetN() + SlackLen + m ); }


//------------------------------------------------------------------------------
inline
Int_T SimplexLP::GetStructN( void )
const
	{ return MPS_LP::GetN(); }


//------------------------------------------------------------------------------
inline
Int_T SimplexLP::GetSlackN( void )
const
	{ return SlackLen; }


//------------------------------------------------------------------------------
inline
Int_T SimplexLP::GetNZ( void )
const
	{ return Int_T( MPS_LP::GetNZ() + SlackLen + m ); }


//------------------------------------------------------------------------------
inline
const Array<Real_T> & SimplexLP::GetLambda( void )
const
{
	assert( Standard );
	return Lambda;
}

//------------------------------------------------------------------------------
inline
Real_T SimplexLP::GetPenalty( void )
const
	{ return M; }


//------------------------------------------------------------------------------
inline
void SimplexLP::SetPenalty( const Real_T NewPenalty )
{
	assert( NewPenalty >= M );
	M = NewPenalty;
}


//------------------------------------------------------------------------------
inline
void SimplexLP::FixLambda( Int_T j )
{
	j -= n + SlackLen;

	assert( j >= 0 && j < m );

	LambdaVT[j]	= VT_FIXED | VT_ARTIF;
}

//
//	End of inline definitions of functions.
//
//==============================================================================

#endif
