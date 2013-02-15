/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	solv_lp.h
CREATED:			1993.09.27
LAST MODIFIED:		1996.03.14

DEPENDENCIES:		compile.h, stdtype.h, stdtmpl.h, mps_lp.h, lp_codes.h,
					smartptr.h
					<assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains declaration of class "SolvableLP" and some inline
functions' definitions. This class is derived from "MPS_LP" class (see: file
"mps_lp.h").
	An object of this class defines a solvable linear programming (LP)
minimization problem. It differs from the MPS_LP in the following:
*	fixed adjustment is added to objective function,
*	first free row found in MPS file is interpreted as an objective function; it
	is removed from the column file; number of rows "m" is reduced by 1,
*	all remaining free rows are also removed, "m" is reduced accordingly,
*	row file representing the LP is added (to the column file already available
	in the "MPS_LP" class),
*	"GetRow" method is provided, which gives access to the row file.

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	SCAL_ACT_COND		-	Lowest value of condition estimator that would
							trigger on scaling.
	SCAL_ACT_MIN		-	Lowest non-zero absolute value allowed when scaling
							is to be omitted.
	SCAL_ACT_MAX		-	Greatest non-zero absolute value allowed when
							scaling is to be omitted.

------------------------------------------------------------------------------*/

#ifndef __SOLV_LP_H__
#define __SOLV_LP_H__

#include <assert.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

#ifndef __COMPILE_H__
#	include "compile.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif
#ifdef SUPPORT_LP_DIT
#	include "lp_dit.h"
#endif


#define SCAL_LP_PASSES		(2)
#define SCAL_ACT_COND		(2.0e0)
#define SCAL_ACT_MIN		(1.0e-1)
#define SCAL_ACT_MAX		(1.0e+1)
#define SCAL_LP_COLS
#define SCAL_LP_ROWS
#define SCAL_LP_EQUILIBRATE


class Solution;

//==============================================================================
//
//	Declaration of class SolvableLP.
//
//==============================================================================

class SolvableLP : public MPS_LP
{
	//==========================================================================
	//	INTERNAL DATA ADDED TO MPS_LP OBJECT: ROW FILE, OBJECTIVE VECTOR, FIXED
	//	ADJUSTMENT, OBJECTIVE NAME
	//
protected:
	Lbl Obj;				// Objective function name.
	Array<Real_T> ar;		// Constraint matrix held by rows (excl. free).
	Real_T f;				// Fixed adjustment to objective function.
	Array<Real_T> c;		// Dense vector representing the objective function.
	Array<Int_T> RowStart;	// Indices of matrix row beginings.
	Array<Int_T> Col;		// Column numbers (in matrix held by rows).

	Bool_T RowsPresent;		// Flag which enables/disables row operations.
	Bool_T CreateRows;		// Flag which encourages row structure creation.

	//==========================================================================
	//	SCALING FACTORS AND SCALING STATE FLAG
	//
protected:
	Array<int> ScaleCol,	// Tables of scaling factors (all scaling factors
		ScaleRow;			// are integer powers of 2, that are stored).
	int CostScale,			// Objective function and RHS/range scaling factors
		RHS_Scale;			// stored as integer powers of two.
	Bool_T Scaled;			// 'True' if scaling was performed.

	//==========================================================================
	//	CONSTRUCTION / DESTRUCTION
	//
	//--------------------------------------------------------------------------
	//	Constructor (initilizes an empty object, does not allocate memory) and
	//	destructor (deallocates memory).
	//
public:
	SolvableLP( void );
	virtual ~SolvableLP( void );

	//==========================================================================
	//	INTERNAL DATA MANIPULATION ROUTINES
	//
	//--------------------------------------------------------------------------
	//	This procedure is called after some rows and/or columns have been
	//	removed (or rather marked as redundant) from the problem.
	//
	void UpdateAfterReduction( const Array<Bool_T> *ExcludeRows,
		const Array<Bool_T> *ExcludeCols = NULL );

	//--------------------------------------------------------------------------
	//	Function for creating row file as a copy of the column file.
	//
	void CreateRowStructure( Bool_T flag );
	Bool_T UpdateRowStructure( void );
	void RemoveRowStructure( void );

	//==========================================================================
	//	LINEAR PROBLEM DATA FILE INPUT AND OUTPUT.
	//
	//--------------------------------------------------------------------------
	//	Method "ReadLP" first calls its inherited counterpart and then processes
	//	the LP problem into standard form. Processing may or may not include
	//	scaling the LP (depending on the value of 'Scale' argument).
	//
	virtual Bool_T ReadAndScaleLP( const char *FileName,
		VerbLevel Verbosity, const Bool_T DIT = False,
		const Bool_T DoScale = False );

	virtual Bool_T ReadAndScaleLP( const char *FileName, FILE *LP_File,
		VerbLevel Verbosity, const Bool_T DIT = False,
		const Bool_T DoScale = False );

	virtual void ProcessSolution( Solution &sol );

	virtual Bool_T CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin,
		Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow );
	virtual Bool_T CreateAsSubmatrix( const MPS_LP &Src, 
		const Array<Bool_T> &IncludeRows, const Array<Bool_T> &IncludeCols );

	virtual Bool_T WriteMPS( FILE *MPS_File, VerbLevel Verbosity ) const;

	virtual Bool_T ScaleLP( VerbLevel Verbosity );
	virtual void UnScaleLP( Bool_T UnScaleMatrix = False );

	virtual Array<Int_T> GetRowLen( void );

	//==========================================================================
	//	ACCESS functions
	//
	//--------------------------------------------------------------------------
	//	This function gives access to rows stored in the row file. It returns
	//	pointers INTO the tables which store the LP problem. Thos is not a
	//	virtual function!
	//
	void GetRow( Int_T Row, Ptr<Real_T> &A, Ptr<Int_T> &Col, Int_T &Len ) const;

	//--------------------------------------------------------------------------
	//  These functions are to be used with extreme caution! They give access
	//  to the whole constraint matrix in raw form. They breach the security
	//  system of 'Array'/'Ptr' template classes. they are only needed when
	//  interfacing to some obscure languages, like Fortran.
	//
public:
	const Real_T *GetNonZerosByRows( void )	const;
	const Int_T *GetColumnNumbers( void )	const;
	const Int_T *GetRowStarts( void )		const;

	//--------------------------------------------------------------------------
	//	Inline implementations of the below listed functions - see the end of
	//	the header file. Both functions are not virtual!
	//
public:
	Real_T	GetC( Int_T j )	const;
	Real_T	GetF( void )	const;

	const char *RevealObjectiveName( void ) const;

	//--------------------------------------------------------------------------
	//	The functions listed below are used to change the cost vector data and
	//	some lower bounds on variables.
	//
	void    SetUnscaledL( Int_T j, Real_T l_j, Bool_T Finite = True );
	void    SetUnscaledU( Int_T j, Real_T u_j, Bool_T Finite = True );

	Bool_T	SetC( Int_T j, Real_T val );
	void	SetF( Real_T fx );
	void	IncrementF( Real_T fx );

	//==========================================================================
	//	LINEAR PROBLEM SCALING / UNSCALING.
	//
	//--------------------------------------------------------------------------
	//	Functions for scaling / unscaling the constraint matrix. Scales only the
	//	column file.
	//
private:
	Bool_T	Scale( VerbLevel Verbosity );
	Bool_T	ScaleCostAndBounds( void );
	void	UnScaleCostAndBounds( void );
	Bool_T	ScaleRangeAndRHS( void );
	void	UnScaleRangeAndRHS( void );

#ifdef SCAL_LP_COLS
	void	ScaleColumns( void );
#endif

#ifdef SCAL_LP_ROWS
	void	ScaleRows( void );
#endif

#ifdef SCAL_LP_EQUILIBRATE
	void	EquilibrateColumns( void );
#endif

#if defined( SCAL_LP_COLS ) || defined( SCAL_LP_EQUILIBRATE ) || \
	defined( SCAL_LP_ROWS )
	void UnScaleMatrix( void );
#endif

	//--------------------------------------------------------------------------
	//	A section of public functions to be used for LP-DIT input/output
	//	(have to be public to be accessed through C code).
	//
#ifdef SUPPORT_LP_DIT

public:

	//
	//	Linear problem output functions.
	//
	virtual void lpi_par( LP_HEAD &h );
	virtual void lpi_cols( LP_HEAD &h, LP_VAR *v );
	virtual void lpi_rows( LP_HEAD &h, LP_VAR *v );
	virtual void lpi_vect( LP_HEAD &h, LP_MAT *v );
#endif

protected:
	void FindFreeRows( Array<Bool_T> &ExcludeRows, Bool_T FindObj = False );
	void PutObjectiveInC( void );
	void DetectFixedVariables( void );

	Bool_T CreateFromMPS_LP( Bool_T DoScale, VerbLevel Verbosity );
};

//==============================================================================
//
//	End of declaration of class SolvableLP.
//
//==============================================================================


//==============================================================================
//
//	Inline definitions of some functions.
//
//==============================================================================


//------------------------------------------------------------------------------
inline
SolvableLP::~SolvableLP( void )
{}


//------------------------------------------------------------------------------
inline
Real_T SolvableLP::GetF( void )
const
	{ return f; }


//------------------------------------------------------------------------------
inline
Real_T SolvableLP::GetC( Int_T j )
const
{
	assert( j >= 0 && j < n );
	return c[j];
}


inline
Bool_T SolvableLP::SetC( Int_T j, Real_T val )
{
	assert( j >= 0 && j < n );
	assert( fabs( val ) < +INFINITY );

	if( Scaled && IsNonZero( val ) )
		val = ldexp( val, -ScaleCol[j] );
	c[j] = val;

	return True;
}


inline
void SolvableLP::SetF( Real_T fx )
{
	assert( fabs( fx ) < +INFINITY );
	f = fx;
}


inline
void SolvableLP::IncrementF( Real_T fx )
{
	f += fx;
	assert( fabs( f ) < +INFINITY );
}

/*------------------------------------------------------------------------------

	inline
	void GetRow( Int_T Row, Ptr<Real_T> &A, Ptr<Int_T> &Col, Int_T &Len ) const

PURPOSE:
	These functions give access to rows of the LP problem. They return pointers
INTO the tables which store the LP problem. That's why the "A", "Row" and "Col"
are pointers to "const" data in the public version of the function. They are,
however, non-const in the protected version.

PARAMETERS:
	Int_T i
		Row number (in range 0 : "m" ) of the row to be accessed.

	Ptr<Real_T> &A
		Array of row "i" non-zeros. It's value on entry is irrelevant.

	Ptr<Int_T> &Col
		Array of row "i" non-zero pattern. It's value on entry is irrelevant.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
void SolvableLP::GetRow( Int_T i, Ptr<Real_T> &a, Ptr<Int_T> &col, Int_T &Len )
const
{
	assert( i >= 0 && i < m );
	assert( RowsPresent );

	Len = RowStart[i+1];
	size_t Start = RowStart[i];

	a.ExtractFragment(   ar,  Start, (size_t) Len );
	col.ExtractFragment( Col, Start, (size_t) Len );

	Len -= Start;
}


/*------------------------------------------------------------------------------

	(inline)
	const Real_T *SolvableLP::GetNonZerosByRows( void ) const
	const Int_T *SolvableLP::GetColumnNumbers( void ) const
	const Int_T *SolvableLP::GetRowStarts( void ) const

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

inline
const Real_T *SolvableLP::GetNonZerosByRows( void )
const
{
	assert( RowsPresent );
	return ( nz > 0 ) ? & ar[0] : (const Real_T *)NULL;
}


inline
const Int_T *SolvableLP::GetColumnNumbers( void )
const
{
	assert( RowsPresent );
	return ( nz > 0 ) ? & Col[0] : (const Int_T *)NULL;
}


inline
const Int_T *SolvableLP::GetRowStarts( void )
const
{
	assert( RowsPresent );
	return & RowStart[0];
}


inline
const char *SolvableLP::RevealObjectiveName( void )
const
{ return Obj; }


inline
void SolvableLP::CreateRowStructure( Bool_T flag )
{ CreateRows = flag; }


//==============================================================================
//
//	End of inline definitions of functions.
//
//==============================================================================

#endif
