/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	presolve.h
CREATED:			1993.10.07
LAST MODIFIED:		1995.10.26

DEPENDENCIES:		simplex.h, error.h, smartptr.h, memblock.h, smartdcl.h,
					sptr_deb.h, stdtype.h, lp_codes.h, solv_lp.h, std_tmpl.h,
					compile.h, mps_lp.h, sort_lab.h, myalloc.h, std_math.h,
					cl_list.h, pre_code.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains a declaration of "PresolvedLP" class, which is
derived from "SolvableLP". It holds data necessary for:
a)	presolving LP previously stored in an underlying "SolvableLP" object and
b)	undoing the presolve operations once the solution to the presolved (and
	thus simplified) problem is known, so that a solution to standard form
	problem is recovered.

	All the presolving operations are performed by "PresolvedLP::Presolve"
function.

------------------------------------------------------------------------------*/

#ifndef __PRESOLVE_H__
#define __PRESOLVE_H__

#include <stdio.h>

#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __PRE_CODE_H__
#	include "pre_code.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __CL_LIST_H__
#	include "cl_list.h"
#endif

class Postsolver;

class Presolver
{
private:
	SolvableLP &LP;

	//--------------------------------------------------------------------------
	//	Linear problem data.
	//
	Int_T m, n;				// Current dimensions of the constraint matrix.

	Array<Int_T> RowLen,	// Current row and column lengths.
		ColLen;				//

	LP_STATUS Status;		// Problem solution status (unknown, solved,
							// infeasible, unbounded).
	Real_T PrimalInf,		// Primal infeasibility or dual infeasibility
		DualInf;			// detected during presolve phase.
	Real_T ftol;			// Feasibility tolerance (primal and dual).

	//--------------------------------------------------------------------------
	// Presolver data.
	//
	Array<Bool_T> ExcludeRows, ExcludeCols, ExplSlackRemoved;
	Array<Real_T> bl, bu;	// Lower/upper bounds on row activity.
	Array<Short_T> rt;		// Row type coded as variable types are.

	Real_T f;				// Fixed adjustment to the objective function.
		
	//
	// List of columns/rows ordered by their lengths.
	//
	ListByLength ColList, RowList;

	//--------------------------------------------------------------------------
	//	Presolver statistics - how many basic operations of each kind were
	//	performed.
	//
	Int_T EliminatedNonZeros,
		EliminatedRows,
		EliminatedCols,
		OrigFixedVars,
		//
		//	Free singleton columns and other singleton column reductions.
		//
		OrigFreeSingletonCols,
		ImpliedFreeSingletonCols,
		RelaxedConstraints,
		//
		//	Row analysis.
		//
		ForcingRows,
		DominatedRows,
		VariableBoundsTightened;

	//--------------------------------------------------------------------------
	//	Postsolver and communication with it.
	//
	Postsolver *PostSolve;

public:
	Presolver( SolvableLP &lp, Postsolver *postsolve = NULL );

	void ConnectToPostsolve( Postsolver *postsolve = NULL );

	void Presolve( const int Mode = LPR_ALL, VerbLevel Verbosity = V_HIGH );
	void UpdateLP_AfterReductions( void );
	void ReleaseWorkMemory( void );

	LP_STATUS ProblemStatus( void ) const;
	Real_T GetPrimalInfeasibility( void ) const;
	Real_T GetDualInfeasibility( void ) const;

	void SetFTOL( Real_T v );
	Real_T GetFTOL( void ) const;

private:
	//--------------------------------------------------------------------------
	//	Auxiliary functions for internal data management.
	//
	void InitializePresolverData( void );

	//--------------------------------------------------------------------------
	//	Auxiliary error-checking functions.
	//
#ifndef NDEBUG
	void CheckMatrixIntegrity( void );
	void CheckColumnStructure( void );
#endif

	//--------------------------------------------------------------------------
	//	Categories of presolve techniques (called by 'Presolve' member
	//	function).
	//
	Int_T EliminateEmptyRows( void );
	Int_T EliminateEmptyColumns( void );
	Int_T EliminateFixedVariables( void );
	Int_T EliminateSingletonRows( void );
	Int_T DealWithSigletonColumns( const int Mode );

	//--------------------------------------------------------------------------
	//	Functions for numerical eliminations.
	//
	Int_T NumericalEliminations( Real_T PivMax );
	Int_T IsNonZeroSubset( Int_T piv, Array<Real_T> &Factor, const Real_T tol );
	Int_T NumElimRows( Int_T piv, const Array<Real_T> &Factor );
	Int_T NumElimCols( Int_T piv, const Array<Real_T> &Factor,
		Array<Bool_T> &Scan, Int_T &ScanCnt );

	//--------------------------------------------------------------------------
	//	Functions for dealing with singleton columns.
	//
	Bool_T EliminateExplicitSlack( Int_T i, Int_T j, Real_T a_ij );
	Bool_T IsImpliedFreeSingleton( Int_T i, Int_T j, Real_T a_ij );

	//--------------------------------------------------------------------------
	//	Functions for detecting forcing and dominated rows.
	//
	Int_T ForcingAndDominatedRows( void );
	Bool_T IsDominatedRow( Int_T row, Real_T AL, Real_T AU, Short_T RT );
	Bool_T IsForcingRow( Int_T row, Real_T AL, Real_T AU, Short_T RT );
	Int_T TightenBounds( Int_T row, Real_T AL, Real_T AU, Short_T RT );

	void ComputeImpliedBounds( Int_T i, Int_T j, Real_T a_ij,
		Real_T &il, Real_T &iu, Short_T &vt );
	void ComputeRowActivityLimits( Int_T row, Real_T &al, Real_T &au,
		Short_T &rt, Int_T xcol = -1 );

	//--------------------------------------------------------------------------
	//	Functions that remove one col/row.
	//
	void FixVariable( Int_T j, int vt );
	void FixVariable( Int_T j, Real_T val );
	void RemoveColumn( Int_T j );
	void RemoveRow( Int_T i );
	void EliminateFreeSingletonColumn( Int_T i, Int_T j, Real_T a_ij );

	//--------------------------------------------------------------------------
	//	Functions for tightening primal variables' bounds and freeing variables.
	//
	Bool_T SetU( Int_T j, Real_T val );
	Bool_T SetL( Int_T j, Real_T val );
	void FreeVar( Int_T j );

	//--------------------------------------------------------------------------
	//	Functions for converting RHS and range vector (e.g. read from MPS file)
	//	to lower and upper bounds on row activity. And the other way around,
	//	too.
	//
	void BLU2RHS( void );
	void RHS2BLU( void );
};


//==============================================================================
//
//	Inline implementations of some of the functions.
//
//==============================================================================

inline
LP_STATUS Presolver::ProblemStatus( void )
const
{ return Status; }


inline
Real_T Presolver::GetPrimalInfeasibility( void )
const
{ return PrimalInf; }


inline
Real_T Presolver::GetDualInfeasibility( void )
const
{ return DualInf; }


inline
void Presolver::ConnectToPostsolve( Postsolver *postsolve )
{ PostSolve = postsolve; }


inline
void Presolver::SetFTOL( Real_T v )
{
	assert( v >= 0.0 && v < 1.0 );
	ftol = v;
}


inline
Real_T Presolver::GetFTOL( void )
const
{ return ftol; }

#endif
