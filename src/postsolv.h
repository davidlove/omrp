/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem postsolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	postsolv.h
CREATED:			1994.05.08
LAST MODIFIED:		1996.04.16

DEPENDENCIES:		my_defs.h, smartptr.h, error.h, memblock.h, smartdcl.h,
					sptr_deb.h, stdtype.h, simplex.h, myalloc.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __POSTSOLV_H__
#define __POSTSOLV_H__


#include <stdio.h>

#ifndef __MY_DEFS_H__
#	include "my_defs.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif


#define POST_SOLV_MIN_LEN	100
#define FEAS_TOL (1.0e-8)

//==============================================================================
//
//	Class "Postsolver" declaration.
//

class PresolverAction;
class Solution;

class Postsolver
{
private:
	Int_T n;					// Linear problem dimension before presolve.	
	Array<Bool_T> ExcludeCols;	// Boolean arrays of markers: "True" if the
								// corresp. column was reduced; "False"
								// otherwise.

	Int_T Len, MaxLen;			// Length and total allocated length of
								// the array of actions.
	Array<PresolverAction *>	Actions;

	Real_T FTOL;				// Feasibility tolerance.

public:
	Postsolver( void );
	Postsolver( Int_T nn );
	~Postsolver( void );

	void SetFTOL( Real_T tol );
	Real_T GetFTOL( void );
	void SetSize( Int_T nn );

	void FreeMemory( void );

	Bool_T Undo( Solution &solution );

	void FixedAdjustment( Real_T val );
	void VariableFixing( Int_T ind, Real_T val );
	void ExplicitSlackRemoval( Int_T ind,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
		Short_T vartype, Real_T ls, Real_T us,
		Short_T rowtype, Real_T bl, Real_T bu );
	void FreeSingletonColumnRemoval( Int_T ind,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
		Real_T af, Real_T b, Real_T slOpt );

private:
	void ExtendActionArray( void );

public:
	Bool_T ReadActions( const char *fname );
	Bool_T ReadActions( FILE *fp );
	Bool_T WriteActions( const char *fname );
	Bool_T WriteActions( FILE *fp );

private:
	//--------------------------------------------------------------------------
	//	Functions used for action file reading.
	//
//	Bool_T GetRowsSection( void );
	Bool_T GetColumnsSection( void );
	Bool_T GetExclusionSection( Array<Bool_T> &Exclude, Int_T k,
		const char *Section );
	Bool_T GetAdjustment( void );
	Bool_T GetAction( Int_T &nn );
	Bool_T GetFixAction( Int_T &nn, Int_T var );
	Bool_T GetRowData( Real_T &val, Int_T &Len, Array<Real_T> &a,
		Array<Int_T> &ind );
	Bool_T GetExplicitSlackAction( Int_T &nn, Int_T var );
	Bool_T GetFreeSingletonAction( Int_T &nn, Int_T var );
	void SkipToNextSection( void );
};


inline
Postsolver::Postsolver( void )
	: n( 0 ), ExcludeCols(), Len( 0 ), MaxLen( 0 ), Actions(),
	FTOL( FEAS_TOL )
	{}


inline
Postsolver::Postsolver( Int_T nn )
	: n( nn ), ExcludeCols( nn+1, False ), Len( 0 ), MaxLen( 0 ), Actions(),
	FTOL( FEAS_TOL )
	{ assert( nn >= 0 ); ExcludeCols[nn] = True; }


inline
void Postsolver::SetFTOL( Real_T tol )
	{ assert( tol > 0.0 && tol < 1.0 ); FTOL = tol; }


inline
Real_T Postsolver::GetFTOL( void )
	{ return FTOL; }

//
//	End of class "Postsolver" declaration.
//
//==============================================================================


//==============================================================================
//
//	Class hierarchy for representing presolver actions.
//

abstract class PresolverAction
{
public:
	virtual ~PresolverAction( void );

	virtual Bool_T Write( FILE *fp ) pure;
	virtual Bool_T Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols ) pure;
};


class FixedAdjustment : public PresolverAction
{
private:
	Real_T Val;

public:
	FixedAdjustment( Real_T val );
	~FixedAdjustment( void );

	virtual Bool_T Write( FILE *fp );
	virtual Bool_T Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols );
};


class VariableFixing : public PresolverAction
{
private:
	Int_T j;
	Real_T Val;

public:
	VariableFixing( Int_T jj, Real_T val );
	~VariableFixing( void );

	virtual Bool_T Write( FILE *fp );
	virtual Bool_T Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols );
};


class ExplicitSlackRemoval : public PresolverAction
{
private:
	Int_T j;
	Array<Real_T> A;
	Array<Int_T> Col;
	Int_T Len;
	Real_T As, Ls, Us, bL, bU;
	Short_T VarType, RowType;

public:
	ExplicitSlackRemoval( Int_T jj, const Array<Bool_T> *ExcludeCols,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
		Short_T vartype, Real_T ls, Real_T us,
		Short_T rowtype, Real_T bl, Real_T bu );
	~ExplicitSlackRemoval( void );

	virtual Bool_T Write( FILE *fp );
	virtual Bool_T Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols );
};


class FreeSingletonColumnRemoval : public PresolverAction
{
private:
	Int_T j;
	Array<Real_T> A;
	Array<Int_T> Col;
	Int_T Len;
	Real_T Af, B, SlOpt;

public:
	FreeSingletonColumnRemoval( Int_T jj, const Array<Bool_T> *ExcludeCols,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
		Real_T af, Real_T b, Real_T slOpt );
	~FreeSingletonColumnRemoval( void );

	virtual Bool_T Write( FILE *fp );
	virtual Bool_T Undo( Solution &solution, Real_T FTOL,
		Array<Bool_T> &ExcludeCols );
};

//
//	End of declarations of classes belonging to the hierarchy of
//	"PresolverActions".
//
//==============================================================================



//==============================================================================
//
//	INLINE definitions of some of the functions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	inline:

	Postsolver::~Postsolver( void )
	PresolverAction::~PresolverAction( void )
	FixedAdjustment::FixedAdjustment( Real_T val )
	VariableFixing::~VariableFixing( void )
	ExplicitSlackRemoval::~ExplicitSlackRemoval( void )
	FreeSingletonColumnRemoval::~FreeSingletonColumnRemoval( void )

PURPOSE:
	Destructors for classes holding the information about the presolve steps.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Postsolver::~Postsolver( void )
	{ FreeMemory(); }

inline
PresolverAction::~PresolverAction( void )
	{}

inline
FixedAdjustment::FixedAdjustment( Real_T val )
	: Val( val )
	{}

inline
FixedAdjustment::~FixedAdjustment( void )
	{}

inline
VariableFixing::VariableFixing( Int_T jj, Real_T val )            
	: j( jj ), Val( val )          
	{ assert( jj >= 0 ); }

inline
Bool_T VariableFixing::Undo( Solution &solution, Real_T, // )
	Array<Bool_T> &ExcludeCols )
{
	assert( ExcludeCols[j] );
	solution.x[j] = Val;
	ExcludeCols[j] = False;
	return True;
}

inline
VariableFixing::~VariableFixing( void )
	{}

inline
ExplicitSlackRemoval::~ExplicitSlackRemoval( void )
	{}

inline
FreeSingletonColumnRemoval::~FreeSingletonColumnRemoval( void )
	{}


#endif
