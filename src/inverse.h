/*------------------------------------------------------------------------------
MODULE TYPE:		Factorization routines - header file
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		prof. Andrzej Ruszczynski (original version),
					Artur Swietanowski (revisions).

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	inverse.h
CREATED:			1991.12.11
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		smartptr.h, stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header contains declaration of class "inverse", which holds a
representation of basis inverse (a LU factorization) in simplex solver.

------------------------------------------------------------------------------*/

#ifndef __INVERSE_H__
#define __INVERSE_H__

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

//
//	Only for VERY thorough self-debugging of factorization routine.
//
//#define FACTOR_DEBUG

//==============================================================================
//
//	Data structures' definition.
//
//==============================================================================

class Inverse
{
private:
	Int_T n,				// Dimension of the matrix to factorize.
		ia;					// Length of arrays a, irow, jcol

	Array<Real_T> a;		// Array of entries of U and M.
	Array<Int_T> irow,		// Row #s for U and pivot rows for M.
		jcol;				// Col #s for U and eliminated rows for M .

	Int_T U_len, M_len,		// Number of nonzeros of U, M.
		rowend,				// Length of row collection for U.
		colend,				// Length of file of column cliques for U.
		cmprs;				// Number of compresses.
	Array<Int_T> rptr,		// Pointers to rows of U in a and jcol.
		cptr,				// Pointers to columns of U in irow.
		rlen,				// Row lengths in U.
		rlst,				// Row list in pivotal order.
		clen,				// Column lengths in U.
		clst;				// Column list in pivotal order.

	Real_T stabrow;			// Safe ratio of pivot to maximum in row
							// (LU update tolerance).
	Real_T stabglo;			// Safe ratio to max in A (when alert) 
	int	alert;				// 0 use stabrow, 1 use stabglo and stabrow.

	Real_T A_max, U_max;	// Maximum entry of A, U.
	Int_T maxties;			// Maximum number of ties allowed in factor.

	Int_T wNz;				// Keeps the length of a sparse work vector
							// (local to some member functions).

	Array<Real_T> Diag;		// An optional diagonal matrix for basis recovery.
	Bool_T DiagPresent;		// Flag allowing/disallowing basis recovery.
	Array<Int_T> DiagUsed;	// Array marking all diagonal entries that are
							// currently being used (stores the numbers of the
							// columns that have been replaced).
	Bool_T DiagChanged;		// Flag marking when the diagonal has been changed
							// and factorization will be needed instead of
							// update.

	//--------------------------------------------------------------------------
	//	Statistic counters.
	//
	Int_T RefactCnt,		// Number of refactorizations.
		UpdateCnt,			// Number of updates.
		SparseFTRAN_Cnt,	// Number of respective types of solve's.
		DenseFTRAN_Cnt,		// 
		SparseBTRAN_Cnt,	// 
		DenseBTRAN_Cnt;		// 

public:
	//--------------------------------------------------------------------------
	//	Basis creation, destruction and management.
	//
	Inverse( Int_T n );						// constructor
	~Inverse( void );						// destructor

	void RegisterDiagonal( const Array<Real_T> &d, Int_T len );
											// Passes a diagonal matrix for
											// basis recovery.
	void RemoveDiagonal( void );			// Disables basis recovery.

	void Clear( void );						// Clean data structures.
	void Resize( Int_T NewN );				// Change basis dimension.

	//--------------------------------------------------------------------------
	//	Basis creation, factorization and update.
	//
	void AddCol( Int_T cnb, const Ptr<Real_T> &col, const Ptr<Int_T> &rnb,
		Int_T len );						// Add column 'cnb' to basis matrix.

	Short_T Update( Int_T m );				// Update factors (column m changed)
	Short_T Factor( void );					// Factorize matrix.

	//--------------------------------------------------------------------------
	//	Inverse representation stability data access.
	//
	Int_T GetFactorLen( void );				// Get length of factors.
	Real_T GetUMax( void );					// Get current stabilty parameters.

	//--------------------------------------------------------------------------
	//	Sparse and dense FTRAN / BTRAN solves.
	//
	void SparseFTRAN( Array<Real_T> &b, Array<Int_T> &bInd, Int_T &bNz );
	void DenseFTRAN( Array<Real_T> &b );
	void SparseBTRAN( Array<Real_T> &b, Array<Int_T> &mark );
	void DenseBTRAN( Array<Real_T> &b );

	//--------------------------------------------------------------------------
	//	Numerical stability parameters management.
	//
	void SetLU_Tol( Real_T tol );
	void SetStab( Real_T stabglo );			// set global stabilty parameter
	void SetAlert( int alert );				// set alert

	//--------------------------------------------------------------------------
	//	Statistic counters management.
	//
	enum CNT { Refact = 100, Upd, SpFTRAN, DenFTRAN, SpBTRAN, DenBTRAN };

	void ResetStatisticCounters( void );
	Int_T ReadStatisticCounter( int n ) const;

private:
	//--------------------------------------------------------------------------
	//	Auxiliary functions for basis factorization.
	//
	Int_T PositionInRowFile( Int_T row, Int_T col );
	Int_T PositionInColumnFile( Int_T row, Int_T col );

	Int_T PositionInFile( Int_T item, Int_T data, Ptr<Int_T> ptr,
		Ptr<Int_T> len, Ptr<Int_T> arr );

	void ReplaceEmptyRowsAndColumns( void );
	void CountRowsAndCols( void );
	void LargestInAllRowsToFront( void );
	void ConstructColumnCliques( void );
	void ConstructEqualLengthLists( Array<Int_T> &rpre, 
		Array<Int_T> &rsuc, Array<Int_T> &cpre, Array<Int_T> &csuc );

	Bool_T FindPivot( Int_T &ipiv, Int_T &jpiv, const Ptr<Int_T> rsuc,
		const Ptr<Int_T> csuc );
	Bool_T FindPivotInCol( const Int_T col, const Int_T nz, Int_T &ipiv,
		Int_T &jpiv, Real_T &bestCost, Int_T &ties, Real_T &bestRatio );
	Bool_T FindPivotInRow( const Int_T row, const Int_T nz, Int_T &ipiv,
		Int_T &jpiv, Real_T &bestCost, Int_T &ties, Real_T &bestRatio );

	void RemoveRowsFromList( const Int_T col, Ptr<Int_T> rsuc,
		Ptr<Int_T> rpre );
	void RemoveColsFromList( const Int_T row, Ptr<Int_T> csuc,
		Ptr<Int_T> cpre );

	void MakeRoomInTheRowFile( Int_T maxLen, Int_T increase );
	Int_T ElimPivotRowFromColumnFile( Int_T ipiv, Int_T jpiv );
	void StoreMultiplier( Int_T row, Int_T col, Real_T val );
	void CompressColumnFile( Int_T NewLength );
	Real_T LargestInRowToFront( Int_T row );
	void MoveRowToEndOfFile( Int_T row );
	void MoveCliqueToEndOfFile( Int_T col );
	void UnpackPivotRowIndice( Int_T row, Array<Int_T> &adr );

	Real_T EliminateElementInPivotColumn( Int_T ipiv, Int_T jpiv, Int_T ielim );
	void EliminateNoFillIn( Int_T ielim, Real_T factor, Array<Int_T> &adr );
	void EliminateWithFillIn( Int_T ipiv, Int_T ielim, Real_T factor,
		Array<Int_T> &adr );

#ifdef FACTOR_DEBUG
	void CheckIntegrity( const Ptr<Int_T> rpre, const Ptr<Int_T> cpre );
#endif
};


//==============================================================================
//
//	End of data structures' definition.
//
//==============================================================================


//==============================================================================
//
//	Inline implementations of some functions.
//
//==============================================================================

inline
Int_T Inverse::GetFactorLen( void )
	{ return Int_T( U_len + M_len ); }

//------------------------------------------------------------------------------
inline
Real_T Inverse::GetUMax( void )
	{ return U_max; }

//------------------------------------------------------------------------------
inline
void Inverse::SetStab( Real_T _stabglo )
	{ stabglo = _stabglo; }

//------------------------------------------------------------------------------
inline
void Inverse::SetAlert( int _alert )
	{ alert = _alert; }

//------------------------------------------------------------------------------
inline
void Inverse::SetLU_Tol( Real_T _tol )
	{ stabrow = _tol; }

//------------------------------------------------------------------------------
inline
Int_T Inverse::PositionInRowFile( Int_T row, Int_T col )
{
	return PositionInFile( row, col, rptr, rlen, jcol );
}

//------------------------------------------------------------------------------
inline
Int_T Inverse::PositionInColumnFile( Int_T row, Int_T col )
{
	return PositionInFile( col, row, cptr, clen, irow );
}


//==============================================================================
//
//	End of inline implementations.
//
//==============================================================================

#endif
