/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - simplex-type solvers.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solv_lp.cpp
CREATED:			1993.09.27
LAST MODIFIED:		1996.04.16

DEPENDENCIES:		stdtype.h, stdtmpl.h, mps_lp.h, solv_lp.h,
					myalloc.h, lp_codes.h, simplex.h, std_math.h,
					smartptr.h, work_vec.h, solution.h, print.h
					<math.h>, <stdlib.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This file contains definitions of non-inline member functions of
"SolvableLP" class.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	SolvableLP::SolvableLP()
	SolvableLP::~SolvableLP()
	void SolvableLP::ReadAndScaleLP()
	Bool_T SolvableLP::UpdateRowStructure( void )
	void SolvableLP::UpdateAfterReduction( const Array<Bool_T> *ExcludeRows,
		const Array<Bool_T> *ExcludeCols )
	Bool_T SolvableLP::Scale( VerbLevel Verbosity )
	void SolvableLP::EquilibrateColumns( void )
	void SolvableLP::ScaleColumns( void )
	void SolvableLP::ScaleRows( void )

STATIC FUNCTIONS:
	static Real_T ConditionEstimator( MPS_LP *LP, Real_T &aMin,
		Real_T &aMax, Int_T n )
	static void FindColumnRanges( MPS_LP *LP, Array<Real_T> &cMin,
		Array<Real_T> &cMax, Int_T n )
	static void FindRowRanges( MPS_LP *LP, Array<Real_T> &rMin,
		Array<Real_T> &rMax, Int_T n )

STATIC DATA:
	static const Real_T M_SQRT

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	TOSTD_LOBND,
	TOSTD_UPBND			-	HAVE TO BE DEFINED! Those macros denote
							respectively lower and upper bound on admissible
							non-zero absolute value in linear problem.

------------------------------------------------------------------------------*/

#include <stdlib.h>

#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif


#define LN_10 (2.302585092994)


/*------------------------------------------------------------------------------

	SolvableLP::SolvableLP( void )

PURPOSE:
	"SolvableLP" class constructor initializes all pointers to NULL (destructor
depends on it) and sets all local problem dimensions to 0. 

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SolvableLP::SolvableLP( void )
	: f( 0.0 ), RowsPresent( False ), CreateRows( True ),
	CostScale( 0 ), RHS_Scale( 0 ), Scaled( False )
{ *Obj = '\0'; }


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::ReadAndScaleLP( const char *FileName, // )
		VerbLevel Verbosity, const Bool_T DIT, const Bool_T DoScale )

	Bool_T SolvableLP::ReadAndScaleLP( const char *FileName, FILE *LP_File, //
		VerbLevel Verbosity, const Bool_T DIT, const Bool_T DoScale )

PURPOSE:
	This function reads LP problem from the indicated file (see: parameters)
using inherited method "MPS_LP::ReadLP" and then converts it to standard
(solvable) form. By standard form we denote a non-equality linear problem with a
single objective function, a right hand side vector and two vectors of simple
bounds on decision variables:

			min c^T x

			Ax >=< b
			0 <= x <= u

	The problem in MPS form is converted to this form by a sequence of simple
transformations:
*	free rows are removed from the constraint matrix (problem dimension is
	reduced accoringly),
*	first free row is interpreted as objective function,
*	fixed variables which have not been explicitly declared fixed have their
	type changed to 'VT_FIXED'.

PARAMETERS:
	const char *FileName
		Input file name (any valid string will do, since it will only be used
		for printing reports). It is passed to "MPS_LP::ReadLP" and not used
		directly.
		If NULL, the input file is not read.

	FILE *LP_File
		Stream descriptor of the input file. It has to represent a valid input
		stream from which LP problem data will be read. It is passed to
		"MPS_LP::ReadLP" and not used directly.
		If NULL, the input file is not read.

	VerbLevel Verbosity
		Specifies the verbosity level.

	const Bool_T DIT
		If 'True', then the linear problem is read in LP-DIT format.

	const Bool_T DoScale
		If 'True' - the problem is scaled.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	File pointer in the input stream will move.

------------------------------------------------------------------------------*/

Bool_T SolvableLP::ReadAndScaleLP( const char *FileName, // )
	VerbLevel Verbosity, const Bool_T DIT, const Bool_T DoScale )
{
	assert( FileName != NULL && *FileName != '\0' );

	if( DIT )
		return ReadAndScaleLP( FileName, (FILE*)NULL, Verbosity, DIT, DoScale );
	else
	{
		FILE *fp = fopen( FileName, "rt" );
		if( fp == NULL ) return False;

		Bool_T result = ReadAndScaleLP( FileName, fp, Verbosity, DIT, DoScale );
		fclose( fp );
		return result;
	}
}


Bool_T SolvableLP::ReadAndScaleLP( const char *FileName, FILE *LP_File, // )
	VerbLevel Verbosity, const Bool_T DIT, const Bool_T DoScale )
{
	//--------------------------------------------------------------------------
	//	Read the input file unless arguments 'FileName' and 'LP_File' are
	//	non-null. Their zero value is considered an error.
	//
	assert( (void *)FileName != NULL );
	assert( DIT || LP_File != NULL );
		
#ifdef SUPPORT_LP_DIT
	if( DIT )
	{
		if( !MPS_LP::ReadLP_DIT( FileName, Verbosity ) )
			return False;
	}
#else
	if( DIT )
	{
		Error( "LP-DIT not supported in the current configuration." );
		return False;
	}
#endif
	else if( !MPS_LP::ReadLP( FileName, LP_File, Verbosity ) )
		return False;

	return CreateFromMPS_LP( DoScale, Verbosity );
}


Bool_T SolvableLP::ScaleLP( VerbLevel Verbosity )
{ return ( Scale( Verbosity ) && UpdateRowStructure() ) ? True : False; }


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin,
		Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow )

	Bool_T SolvableLP::CreateAsSubmatrix( const MPS_LP &Src,
		const Array<Bool_T> &IncludeRows, const Array<Bool_T> &IncludeCols )

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

Bool_T SolvableLP::CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin, // )
	Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow )
{
	//--------------------------------------------------------------------------
	//	Create the lower level object's data by calling appropriate function.
	//
	if( !MPS_LP::CreateAsSubmatrix( Src, RowMin, RowMax, ColMin, ColMax,
		CostRow ) )
		return False;

	return CreateFromMPS_LP( False, V_NONE );
}


Bool_T SolvableLP::CreateAsSubmatrix( const MPS_LP &Src, // )
	const Array<Bool_T> &IncludeRows, const Array<Bool_T> &IncludeCols )
{
	//--------------------------------------------------------------------------
	//	Create the lower level object's data by calling appropriate function.
	//
	if( !MPS_LP::CreateAsSubmatrix( Src, IncludeRows, IncludeCols ) )
		return False;

	return CreateFromMPS_LP( False, V_NONE );
}


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::UpdateRowStructure( void )

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

Bool_T SolvableLP::UpdateRowStructure( void )
{
	if( !CreateRows ) return True;

	//--------------------------------------------------------------------------
	//	Allocate the arrays of non-zeros and their column numbers for the row
	//	file.
	//
	ar.Resize( (size_t) nz );
	Col.Resize( (size_t) nz );

	RowStart.Resize( m + 1 );
	RowStart.Fill( 0, m + 1 );

	Int_T i, j, len;
	Ptr<Int_T> row;
	Ptr<Real_T> a;

	//--------------------------------------------------------------------------
	//	Scan WHOLE column file (i.e. including free rows) to compute row
	//	lengths in 'RowStart[]' array.
	//
	for( j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++row )
		{
			assert( *row >= 0 && *row < m );

			RowStart[ *row ]++;
		}

	//--------------------------------------------------------------------------
	//	Make 'RowStart' entries point just beyond the last positions of row
	//	non-zeros.
	//
	for( i = 1; i <= m; i++ )
		RowStart[i] += RowStart[i-1];

	assert( RowStart[m] == nz );

	//--------------------------------------------------------------------------
	//	Copy matrix representation by columns into matrix by rows.
	//
	for( j = Int_T( n - 1 ); j >= 0; j-- )
		for( GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			Int_T rs = --RowStart[ *row ];

			ar[ (size_t) rs ]	= *a;
			Col[ (size_t) rs ]	= j;
		}

	RowsPresent = True;

	//--------------------------------------------------------------------------
	//	Compute control sums by rows and by columns and compare them.
	//
#ifndef NDEBUG
	Real_T RowCI	= 0.0,
		RowCR		= 0.0,
		ColCI		= 0.0,
		ColCR		= 0.0;
	Ptr<Int_T> col;

	for( i = 0; i < m; i++ )
		for( GetRow( i, a, col, len ); len; --len, ++a, ++col )
		{
			RowCI	+= *col + i;
			RowCR	+= *a;
		}

	for( j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			ColCI	+= *row + j;
			ColCR	+= *a;
		}

	assert( IsEqual( ColCI, RowCI ) && IsEqual( ColCR, RowCR ) );
#endif

	return True;
}


/*------------------------------------------------------------------------------

	void SolvableLP::UpdateAfterReduction( const Array<Bool_T> &ExcludeRows,
		const Array<Bool_T> &ExcludeCols );

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

void SolvableLP::UpdateAfterReduction( const Array<Bool_T> *ExcludeRows, // )
	const Array<Bool_T> *ExcludeCols )
{
	Int_T i, j;
	Bool_T Cols = ( ExcludeCols ) ? True : False;

	//--------------------------------------------------------------------------
	//	Renumber the rows (use a work array to compute and hold the row
	//	permutation). Let (-1) denote rows, which have been removed and thus
	//	do not need a new number.
	//
	Array<Int_T> Num( m, -1 );

	{
		Int_T CurNum = 0;

		for( i = 0; i < m; i++ )
			if( !(*ExcludeRows)[i] )
				Num[i] = CurNum++;
	}

	//--------------------------------------------------------------------------
	//	Pack the column file (remove all non-zeros corresponding to erased
	//	columns or rows). At the same time update the 'ColStart' table.
	//	Permute row numbers during this scan (use 'RowNum' array).
	//
	//	Access column file directly.
	//
	Int_T OldColStart	= 0,
		NewColStart		= 0;
	Int_T OldLen,
		NewN			= 0;

	for( j = 0; j < n; j++ )
	{
		if( Cols && (*ExcludeCols)[j] ) continue;

		OldColStart		= ColStart[j];
		OldLen			= Int_T( ColStart[j+1] - OldColStart );

		ColStart[NewN]	= NewColStart;

		for( ; OldLen; OldLen--, OldColStart++ )
		{
			Int_T row = Row[ (size_t) OldColStart ];

			if( (*ExcludeRows)[row] ) continue;

			ac[ (size_t) NewColStart ]	= ac[ (size_t) OldColStart ];
			Row[ (size_t) NewColStart ]	= Num[ row ];
			NewColStart++;
		}

		NewN++;
	}
	Int_T New_nz = ColStart[NewN] = NewColStart;

	//--------------------------------------------------------------------------
	//	Optionally update the column-related data (pack objective vector, simple
	//	bounds and variable type vectors). Also the column ordering data is
	//	packed.
	//
	Int_T New_n = n;

	if( Cols )
	{
		//----------------------------------------------------------------------
		//	'i' will point to the new position in the vectors, 'j' - to the old
		//	one.
		//
		for( i = j = 0; j < n && !(*ExcludeCols)[j]; j++, i++ );

		for( ; j < n; j++ )
			if( !(*ExcludeCols)[j] )
			{
				c[i]			= c[j];
				l[i]			= l[j];
				u[i]			= u[j];
				VarType[i]		= VarType[j];

				i++;
			}

		New_n = i;

		ColLabels.RemoveLabels( *ExcludeCols );
	}

	//--------------------------------------------------------------------------
	//	Update the row-related data (row activity bounds, row type and range
	//	vectors).
	//
	//	'j' will point to the new position in the vectors, 'i' - to the old one.
	//
	for( i = 0; !(*ExcludeRows)[i]; i++ )
		;

	for( j = i; i < m; i++ )
		if( !(*ExcludeRows)[i] )
		{
			b[j]			= b[i];
			r[j]			= r[i];
			RowType[j]		= RowType[i];

			j++;
		}
	Int_T New_m = j;

	RowLabels.RemoveLabels( *ExcludeRows );

	//--------------------------------------------------------------------------
	//	Resize 'Num' and use it to renumber columns.
	//
	if( Cols )
	{
		Num.Resize( n );

		Int_T CurCol = 0;

		for( j = 0; j < n; j++ )
			Num[j] = Int_T( ExcludeCols[0][j] ? -1 : CurCol++ );
	}

	Num.Resize( 0 );

	//--------------------------------------------------------------------------
	//	Update the problem dimensions. Reallocate the column file.
	//
	if( nz > 1.3 * New_nz )
	{
		if( New_nz )
		{
			ar.Resize( (size_t) New_nz );
			Col.Resize( (size_t) New_nz );
		}
		else
		{
			ar.Resize( 0 );
			Col.Resize( 0 );
		}
	}
	nz = New_nz;

	if( n > New_n )
	{
		ColStart.Resize( New_n + 1 );
		c.Resize( New_n );
		u.Resize( New_n );
		l.Resize( New_n );
		VarType.Resize( New_n );
	}
	n = New_n;

	if( m > New_m )
	{
		b.Resize( New_m );
		r.Resize( New_m );
		RowType.Resize( New_m );
	}
	m = New_m;

	//--------------------------------------------------------------------------
	//	Mark the row structure as invalid.
	//
	RowsPresent = False;
}


/*------------------------------------------------------------------------------

	void MPS_LP::SetUnscaledL( Int_T j, Real_T l_j, Bool_T Finite )
	void MPS_LP::SetUnscaledU( Int_T j, Real_T u_j, Bool_T Finite )

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SolvableLP::SetUnscaledL( Int_T j, Real_T l_j, Bool_T Finite )
{
	if( Finite && Scaled && IsNonZero( l_j ) )
		l_j = ldexp( l_j, ScaleCol[j] );
	MPS_LP::SetL( j, l_j, Finite );
}

 
void SolvableLP::SetUnscaledU( Int_T j, Real_T u_j, Bool_T Finite )
{
	if( Finite && Scaled && IsNonZero( u_j ) )
		u_j = ldexp( u_j, ScaleCol[j] );
	MPS_LP::SetU( j, u_j, Finite );
}

 
/*------------------------------------------------------------------------------

	Bool_T SolvableLP::ScaleCostAndBounds( void )

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

static const Real_T M_SQRT = sqrt( 2.0 );

Bool_T SolvableLP::ScaleCostAndBounds( void )
{
	CostScale = 0;						// Reset the scaling factor.

	//--------------------------------------------------------------------------
	//	Scale the objective vector.
	//

	//
	//	Pass 1: Apply scaling factors found during column scaling to the
	//	objective function. At the same time find the largest cost coefficient.
	//
	Real_T cMax = 0.0;
	int pwr;
	Int_T j;

	for( j = 0; j < n; j++ )
		if( ( pwr = ScaleCol[j] ) == 0 )
			cMax = Max( fabs( c[j] ), cMax );
		else if( IsNonZero( c[j] ) )
			cMax = Max( fabs( c[j] = ldexp( c[j], -pwr ) ), cMax );

	if( IsZero( cMax ) )				// Nothing to scale.
		return False;

	if( cMax > sqrt( +INFINITY ) )
		FatalError( "The cost vector holds huge numbers (%10.2E).\n",
			cMax );

	//
	//	Pass 2: Determine the objective scaling factor as an approximation of
	//	the largest cost by an integer power of two. Divide the whole cost
	//	vector by this factor.
	//
	frexp( cMax / M_SQRT, &CostScale );

	if( CostScale != 0 )
		for( j = 0; j < n; j++ )
			if( IsNonZero( c[j] ) )
				c[j] = ldexp( c[j], -CostScale );

	//--------------------------------------------------------------------------
	//	Recalculate bound vectors. Use scaling factors found in column scaling.
	//
	for( j = 0; j < n; j++ )
	{
		if( ( pwr = ScaleCol[j] ) == 0 ) continue;

		if( VarType[j] & VT_HAS_UP_BND && IsNonZero( u[j] ) )
			u[j] = ldexp( u[j], pwr );

		if( ( VarType[j] & VT_HAS_LO_BND ) && IsNonZero( l[j] ) )
			l[j] = ldexp( l[j], pwr );
	}

	return True;
}


/*------------------------------------------------------------------------------

	void SolvableLP::UnScaleCostAndBounds( void )

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

void SolvableLP::UnScaleCostAndBounds( void )
{
	if( !Scaled && CostScale == 0 ) return;

	//--------------------------------------------------------------------------
	//	Un-scale the objective vector and variable simple bounds.
	//
	//	Apply scaling factors found during column scaling and the objective
	//	scaling factor to the objective function.
	//
	for( Int_T j = 0; j < n; j++ )
	{
		int pwr = ScaleCol[j];

		if( IsNonZero( c[j] ) )
			c[j] = ldexp( c[j], pwr + CostScale );

		//----------------------------------------------------------------------
		//	Recalculate bound vectors. 
		//
		if( pwr )
		{
			if( VarType[j] & VT_HAS_UP_BND && IsNonZero( u[j] ) )
				u[j] = ldexp( u[j], -pwr );

			if( ( VarType[j] & VT_HAS_LO_BND ) && IsNonZero( l[j] ) )
				l[j] = ldexp( l[j], -pwr );
		}
	}
}


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::ScaleRangeAndRHS( void )

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

Bool_T SolvableLP::ScaleRangeAndRHS( void )
{
	RHS_Scale = 0;						// Reset the scaling factor.

	//--------------------------------------------------------------------------
	//	Scale the right hand side vector and range vector.
	//

	//
	//	Pass 1: Apply scaling factors found in row scaling. Find the largest
	//	RHS element.
	//
	Real_T bMax = 0.0;
	int pwr;
	Int_T i;

	for( i = 0; i < m; i++ )
	{
		if( ( pwr = ScaleRow[i] ) != 0 )
		{
			if( IsNonZero( b[i] ) )
				b[i] = ldexp( b[i], -pwr );
			if( GetRowType( i ) & RT_RNG && IsNonZero( r[i] ) )
				r[i] = ldexp( r[i], -pwr );
		}
		bMax = Max( fabs( b[i] ), bMax );
	}

	if( IsZero( bMax ) )				// Empty RHS -> nothing to scale.
		return False;

	if( bMax > sqrt( +INFINITY ) )
		FatalError( "The right hand side vector holds huge numbers (%10.2E).\n",
			bMax );

	//
	//	Pass 2: Compute the RHS/range vector scaling bound. Scale both vectors.
	//	Also scale the lower/upper bounds on variables.
	//
	frexp( bMax * M_SQRT, &RHS_Scale );

	if( RHS_Scale )
	{
		for( i = 0; i < m; i++ )
		{
			if( IsNonZero( b[i] ) )
				b[i] = ldexp( b[i], -RHS_Scale );
			if( ( GetRowType( i ) & RT_RNG ) && IsNonZero( r[i] ) )
				r[i] = ldexp( r[i], -RHS_Scale );
		}

		for( i = 0; i < n; i++ )
		{
			if( ( VarType[i] & VT_HAS_UP_BND ) && IsNonZero( u[i] ) )
				u[i] = ldexp( u[i], -RHS_Scale );

			if( ( VarType[i] & VT_HAS_LO_BND ) && IsNonZero( l[i] ) )
				l[i] = ldexp( l[i], -RHS_Scale );
		}
	}

	return True;
}


/*------------------------------------------------------------------------------

	void SolvableLP::UnScaleRangeAndRHS( void )

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

void SolvableLP::UnScaleRangeAndRHS( void )
{
	if( !Scaled && RHS_Scale == 0 ) return;

	//--------------------------------------------------------------------------
	//	Un-scale the right hand side vector and range vector.
	//
	//	Apply scaling factors found in row scaling.
	//
	Int_T i;
	for( i = 0; i < m; i++ )
	{
		int pwr = ScaleRow[i];

		if( IsNonZero( b[i] ) )
			b[i] = ldexp( b[i], pwr + RHS_Scale );
		if( GetRowType( i ) & RT_RNG && IsNonZero( r[i] ) )
			r[i] = ldexp( r[i], pwr + RHS_Scale );
	}

	//
	//	Un-scale the lower/upper bounds on variables.
	//
	if( RHS_Scale )
	{
		for( i = 0; i < n; i++ )
		{
			if( ( VarType[i] & VT_HAS_UP_BND ) && IsNonZero( u[i] ) )
				u[i] = ldexp( u[i], RHS_Scale );

			if( ( VarType[i] & VT_HAS_LO_BND ) && IsNonZero( l[i] ) )
				l[i] = ldexp( l[i], RHS_Scale );
		}
	}
}


//==============================================================================
//
//	Static function prototypes - matrix scaling.
//
//==============================================================================

static Real_T ConditionEstimator( MPS_LP *LP, Real_T &aMin, Real_T &aMax,
	Int_T n );

#ifdef SCAL_LP_COLS
static void FindColumnRanges( MPS_LP *LP, Array<Real_T> &cMin,
	Array<Real_T> &cMax, Int_T n );
#endif

#ifdef SCAL_LP_ROWS
static void FindRowRanges( MPS_LP *LP, Array<Real_T> &rMin,
	Array<Real_T> &rMax, Int_T n );
#endif

//==============================================================================
//
//	End of static function prototypes.
//
//==============================================================================


//==============================================================================
//
//	Implementation of all scaling functions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::Scale( VerbLevel Verbosity )

PURPOSE:
	This function scales a linear problem with integer powers of two. A typical
scaling routine would first do two sequential passes in which it would scale
columns and rows of the constraint matrix. The purpose of the scaling is to make
all row / column non-zeros as close to one as possible (on logarithmic scale).
To this end we use a geometrical mean of the row's / column's maximum and
minimum non-zero absolute values as a scaling factor.
	Additionally we don't want to introduce any roundoff error. Therefore we
actually approximate the scaling factors with integer powers of two. This allows
us to avoid almost all arithmetic multiplications and divisions - instead we
only adjust scaled numbers' exponents. We also store integer scaling factors -
the powers of two. To obtain the best possible approximation with available
fnuction "frexp", we divide the scaling factor by a square root of two.
	After the passes one additional operation - a so called equilibration - may
be performed. Column equilibration attempts to bring the largest absolute value
of a column non-zero to one. This makes sense only for simplex-type algorithm.
Once again we use approximation with integer powers of two.
	If the initial condition estimator is less then "SCAL_ACT_COND" and the
smallest non-zero absolute value is larger than "SCAL_ACT_MIN" and the largest
non-zero absolute value is smaller than "SCAL_ACT_MAX" the scaling is omitted.
If only the first of the above mentioned conditions is met, there is only one
pass of scaling. Otherwise there are "SCAL_LP_PASSES" passes.

PARAMETERS:
	VerbLevel Verbosity
		Represents report verbosity level required by the calling prcedure.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	The column file will be scaled, the row file will not.

------------------------------------------------------------------------------*/

Bool_T SolvableLP::Scale( VerbLevel Verbosity )
{
	Bool_T ScaleMatrix = True;

	if( Scaled ) return True;

	if( Verbosity >= V_LOW )
		Print( "\nScaling" );

	//--------------------------------------------------------------------------
	//	Declare local variables that will be used during the processing.
	//
	Real_T Cond;			// Condition estimator - will be computed before and
							// after scaling.

	//--------------------------------------------------------------------------
	//	First check if the scaling is necessary.
	//	We consider it unneeded if the largest element in the matrix is no more
	//	than a hundred times larger than the smallest one and if the elements
	//	are in range "SCAL_ACT_MIN" - "SCAL_ACT_MAX".
	//
	Real_T aMin, aMax;		// Minimum / maximum matrix element absolute value.
	Cond = ConditionEstimator( this, aMin, aMax, n );

	if( Cond <= SCAL_ACT_COND && aMin >= SCAL_ACT_MIN && aMax <= SCAL_ACT_MAX )
		ScaleMatrix = False;

	//--------------------------------------------------------------------------
	//	Allocate space for data members "ScaleCol" and "ScaleRow" and initialize
	//	them.
	//
	ScaleCol.Resize( n );		ScaleCol.Fill( 0, n );
	ScaleRow.Resize( m );		ScaleRow.Fill( 0, m );

	int MaxPass = ( Cond > SCAL_ACT_COND ) ? SCAL_LP_PASSES : 1;

	if( Verbosity >= V_LOW )
	{
		if( ScaleMatrix )
		{
			Print(
#ifdef SCAL_LP_EQUILIBRATE
				" (%d pass%s + equilibration):\n"
#else
				" (%d pass%s):\n"
#endif
				"\tInitial condition = %9.3E\n",
				(int) MaxPass, (char *)( MaxPass==1 ? "" : "es" ),
				(double)Cond
			);
		}
		else
		{
			Print(
#ifdef SCAL_LP_EQUILIBRATE
				" (equilibration only):\n"
#else
				" (not needed):\n"
#endif
				"\tInitial condition = %9.3E\n", (double) Cond
			);
		}
	}
	else if( Verbosity >= V_LINE )
	{
		if( ScaleMatrix )
			Print( " %10.3E |", (double) Cond );
		else
			Print( " %10.3E | %10.3E |", (double) Cond, (double) Cond );
	}

	//--------------------------------------------------------------------------
	// Scaling loop ("SCAL_LP_PASSES" passes are performed).
	//
	if( ScaleMatrix )
	{
		for( int pass = 0; pass < MaxPass; pass++ )
		{
#ifdef SCAL_LP_COLS
			ScaleColumns();
#endif

#ifdef SCAL_LP_ROWS
			ScaleRows();
#endif
		}
		Scaled = True;
	}

	//--------------------------------------------------------------------------
	//	Approximate column equilibration (dividing columns by the largest
	//	absolute value of column's non-zero) with integer powers of two.
	//
#ifdef SCAL_LP_EQUILIBRATE
	EquilibrateColumns();
	Scaled = True;
#endif

	if( ScaleCostAndBounds() )
		Scaled = True;

	if( ScaleRangeAndRHS() )
		Scaled = True;

	//--------------------------------------------------------------------------
	//	Compute the condition estimator again.
	//
	if( ScaleMatrix )
	{
		Cond = ConditionEstimator( this, aMin, aMax, n );
		if( Verbosity >= V_LOW )
			Print( "\tFinal condition = %9.3E\n", (double) Cond );
		else if( Verbosity >= V_LINE )
			Print( " %10.3E |", (double) Cond );
	}

	return True;
}


/*------------------------------------------------------------------------------

	void SolvableLP::UnScaleLP( Bool_T UnScaleMatrix )

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

void SolvableLP::UnScaleLP( Bool_T unscaleMatrix )
{
	if( !Scaled ) return;

	if( unscaleMatrix )
	{
#if defined( SCAL_LP_COLS ) || defined( SCAL_LP_EQUILIBRATE ) || \
	defined( SCAL_LP_ROWS )
		UnScaleMatrix();
#endif
	}

	UnScaleCostAndBounds();
	UnScaleRangeAndRHS();
}


#ifdef SCAL_LP_EQUILIBRATE
void SolvableLP::EquilibrateColumns( void )
{
	Ptr<Real_T> a;	// Those three data objects will be used when
	Ptr<Int_T> row;	// accesing the matrix using "GetColumn" function.
	Int_T len;

	for( Int_T j = 0; j < n; j++ )
	{
		//----------------------------------------------------------------------
		//	1.	Find largest absolute value of column non-zero element.
		//
		Real_T aMax = 0.0e0;

		for( GetColumn( j, a, row, len ); len; --len, ++a )
			aMax = Max( fabs( *a ), aMax );

		//----------------------------------------------------------------------
		//	2.	Determine its approximation by an integer power of two.
		//
		int pwr;
		
		frexp( aMax / M_SQRT, &pwr );

		//----------------------------------------------------------------------
		//	3.	If power is non-zero:
		//		1.	store the scaling factor (add it to previously calculated
		//			one),
		//		2.	divide all column non-zeros by that factor.
		//
		if( pwr )
		{
			ScaleCol[j] += pwr;
			for( GetColumn( j, a, row, len ); len; --len, ++a )
				*a = ldexp( *a, -pwr );
		}
	}
}
#endif


#ifdef SCAL_LP_COLS
void SolvableLP::ScaleColumns( void )
{
	//--------------------------------------------------------------------------
	//	Allocate workspace.
	//
	WorkVector<Real_T> cMin( n ),	// These are tables of minimum and maximum
		cMax( n );					// column and non-zero absolute values.

	//--------------------------------------------------------------------------
	//	Find minimum and maximum absolute values of column elements.
	//
	FindColumnRanges( this, cMin, cMax, n );
		
	Ptr<Real_T> a;	// Those three data objects will be used when
	Ptr<Int_T> row;	// accesing the matrix using "GetColumn".
	Int_T len;

	//----------------------------------------------------------------------
	//	Scale:
	//	1.	Compute scaling factors that will not generate a roundoff
	//		error. To assure the most accurate approximation of the "x"
	//		by "2**pwr", pass "x*sqrt(2)" to "frexp".
	//	2.	Scale the constraint matrix.
	//
	for( Int_T j = 0; j < n; j++ )
	{
		Real_T x = sqrt( cMin[j] * cMax[j] ) / M_SQRT;
		int pwr;

		if( IsNonZero( x ) )
		{
			frexp( x, &pwr );
			if( pwr )
			{
				ScaleCol[j] += pwr;

				for( GetColumn( j, a, row, len ); len; --len, ++a )
					*a = ldexp( *a, -pwr );
			}
		}
	}
}
#endif


#ifdef SCAL_LP_ROWS
void SolvableLP::ScaleRows( void )
{
	//--------------------------------------------------------------------------
	//	Allocate workspace.
	//
	WorkVector<Real_T> rMin( m ),	// These are tables of minimum and
		rMax( m );					// maximum row non-zero absolute values.
	WorkVector<int> rStore( m );

	rStore.Fill( 0, m );

	//--------------------------------------------------------------------------
	//	Find minimum and maximum absolute values of row elements.
	//
	FindRowRanges( this, rMin, rMax, n );

	//--------------------------------------------------------------------------
	//	Compute scaling factors that will not generate a roundoff error.
	//	Store factors temporarily in "rStore" vector.
	//
	int pwr;

	for( Int_T i = 0; i < m; i++ )
	{
		Real_T x = sqrt( rMin[i] * rMax[i] ) / M_SQRT;

		if( IsNonZero( x ) )
		{
			frexp( x , &pwr );
			ScaleRow[i] += rStore[i] = (Int_T) pwr;
		}
	}

	//--------------------------------------------------------------------------
	//	Scale the constraint matrix (by columns).
	//
	Ptr<Real_T> a;		// Those three data objects will be used when
	Ptr<Int_T> row;		// accesing the matrix using "GetColumn" function.
	Int_T len;

	for( Int_T j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			pwr = rStore[ *row ];

			if( pwr )
				*a = ldexp( *a, -pwr );
		}
}
#endif


#if defined( SCAL_LP_COLS ) || defined( SCAL_LP_EQUILIBRATE ) || \
	defined( SCAL_LP_ROWS )
void SolvableLP::UnScaleMatrix( void )
{
	Ptr<Real_T> a;	// Those three data objects will be used when
	Ptr<Int_T> row;	// accesing the matrix using "GetColumn".
	Int_T len;

	//----------------------------------------------------------------------
	//	Un-scale:
	//
	for( Int_T j = 0; j < n; j++ )
	{
		int pwr = ScaleCol[j];

		for( GetColumn( j, a, row, len ); len; --len, ++a )
			*a = ldexp( *a, pwr + ScaleRow[ *row ] );
	}
}
#endif


/*------------------------------------------------------------------------------

	f()

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

static Real_T ConditionEstimator( MPS_LP *LP, Real_T &aMin, // )
	Real_T &aMax, Int_T n )
{
	aMin = INFINITY;
	aMax = 0.0e0;

	Ptr<Real_T> a;	// Those three data objects will be used when
	Ptr<Int_T> row;	// accesing the matrix using "GetColumn" function.
	Int_T len;

	//--------------------------------------------------------------------------
	//	The first loop is there only to ensure that both 'aMin' nd 'aMax' will
	//	be initialized with propper values. It finds the first non-empty column
	//	and it scans it.
	//
	Int_T j;
	for( j = 0; j < n; j++ )
	{
		LP->GetColumn( j, a, row, len );

		if( len == 0 ) continue;

		for( ; len; --len, ++a )
		{
			Real_T x = fabs( *a );

			if( x < aMin )	aMin = x;
			if( x > aMax )	aMax = x;
			j++;
		}
		break;
	}

	//--------------------------------------------------------------------------
	//	The second loop searches the whole matrix, starting from the place where
	//	the first one left off.
	//
	for( ; j < n; j++ )
	{
		for( LP->GetColumn( j, a, row, len ); len; --len, ++a )
		{
			Real_T x = fabs( *a );

			if( x < aMin )
				aMin = x;
			else if( x > aMax )
				aMax = x;
		}
	}

	assert( aMax > 0.0 && aMin > 0.0 );

	return log( aMax/aMin )/ LN_10;
}


/*------------------------------------------------------------------------------

	f()

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

#ifdef SCAL_LP_COLS
static void FindColumnRanges( MPS_LP *LP, Array<Real_T> &cMin, // )
	Array<Real_T> &cMax, Int_T n )
{
	Ptr<Real_T> a;	// Those three data objects will be used when
	Ptr<Int_T> row;	// accesing the matrix using "GetColumn" function.
	Int_T len;

	//----------------------------------------------------------------------
	//	Find largest and smallest column elements.
	//
	for( Int_T j = 0; j < n; j++ )
	{
		LP->GetColumn( j, a, row, len );

		if( len )
		{
			Real_T &min = cMin[j],
				&max = cMax[j];
				
			min = max = fabs( *a );

			for( --len, ++a; len; --len, ++a )
			{
				Real_T x = fabs( *a );

				if( x < min )
					min = x;
				else if( x > max )
					max = x;
			}
		}
		else
			cMin[j] = cMax[j] = 0.0;
	}
}
#endif


/*------------------------------------------------------------------------------

	f()

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

#ifdef SCAL_LP_ROWS
static void FindRowRanges( MPS_LP *LP, Array<Real_T> &rMin, // )
	Array<Real_T> &rMax, Int_T n )
{
	//----------------------------------------------------------------------
	//	Initialize the work vectors.
	//
	rMin.Fill( INFINITY,	LP->GetM() );
	rMax.Fill( 0.0e0,		LP->GetM() );

	//----------------------------------------------------------------------
	//	Find largest and smallest row elements (search by columns).
	//
	for( Int_T j = 0; j < n; j++ )
	{
		Ptr<Real_T> a;	// Those three data objects will be used when
		Ptr<Int_T> row;	// accesing the matrix using "GetColumn" function.
		Int_T len;

		for( LP->GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			Real_T x	= fabs( *a ),
				&min	= rMin[ *row ],
				&max	= rMax[ *row ];

			if( x < min )	min = x;
			if( x > max )	max = x;
		}
	}
}
#endif


/*------------------------------------------------------------------------------

	void SolvableLP::ProcessSolution( Solution &sol )

PURPOSE:
	This function takes the solution to a linear problem and (if the problem
was scaled before solving) retrieves the solution to the original problem.

PARAMETERS:
	Solution &sol
		the solution structure, which holds:
		-	solution vector's dimensions,
		-	primal and dual variables,
		-	the objective function value.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SolvableLP::ProcessSolution( Solution &sol )
{
	Array<Real_T> &x	= sol.x,
		&y				= sol.y,
		&z				= sol.z,
		&ra				= sol.ra;

	int contents = sol.GetContents();

	assert( !( contents & Solution::Primal ) || sol.GetN() == n );
	assert( !( contents & Solution::Dual ) || sol.GetM() == m );

	if( Scaled )
	{
		Int_T i, j;

		if( contents & Solution::Primal )
			for( j = 0; j < n; j++ )
				x[j] = ldexp( x[j], RHS_Scale - ScaleCol[j] );

		if( contents & Solution::RC )
			for( j = 0; j < n; j++ )
				z[j] = ldexp( z[j], CostScale + ScaleCol[j] );

		if( contents & Solution::RowAct )
			for( i = 0; i < m; i++ )
				ra[i] = ldexp( ra[i], RHS_Scale - ScaleRow[i] );

		if( contents & Solution::Dual )
			for( i = 0; i < m; i++ )
				y[i] = ldexp( y[i], CostScale - ScaleRow[i] );

		if( contents & Solution::Res )
			sol.result = ldexp( sol.result, CostScale + RHS_Scale );
	}
}


/*------------------------------------------------------------------------------

	Array<Int_T> SolvableLP::GetRowLen( void )

PURPOSE:
	Returns an array of row lengths. Needed during writing an LP-DIT output
file.

PARAMETERS:
	None.

RETURN VALUE:
	The array.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Array<Int_T> SolvableLP::GetRowLen( void )
{
	assert( RowsPresent );

	Array<Int_T> RowLen( m );

	for( Int_T i = 0; i < m; i++ )
		RowLen[i] = Int_T( RowStart[i+1] - RowStart[i] );

	return RowLen;
}


void SolvableLP::FindFreeRows( Array<Bool_T> &ExcludeRows, Bool_T FindObj )
{
	//--------------------------------------------------------------------------
	//	Find free rows. NOTE: The first free row will be considered to be the
	//	objective function. Remember its number and store its row label.
	//
	//	If the objective row was given before (e.g., by LP-DIT input) it will
	//	not be changed.
	//
	for( Int_T i = 0; i < m; i++ )
		if( RowType[i] & RT_FR )
		{
			ExcludeRows[i] = True;
			if( FindObj && ObjRow == -1 )
			{
				ObjRow = i;
				if( OS == OS_FULL )
					strcpy( Obj, RowLabels.FindLabel( i ) );
				else
					*Obj = '\0';
			}
		}
}


void SolvableLP::PutObjectiveInC( void )
{
	assert( ObjRow >= 0 && ObjRow < m );

	c.Resize( n );
	c.Fill( 0.0 , n);

	//--------------------------------------------------------------------------
	//	Fill the objective function vector 'c'.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	for( Int_T j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++a, ++row )
			if( *row == ObjRow )
			{
				c[j] = *a;
				break;
			}
}


void SolvableLP::DetectFixedVariables( void )
{
	//--------------------------------------------------------------------------
	//	Detect FIXED variables (that have not been declared to be fixed).
	//
	for( Int_T j = 0; j < n; j++ )
		if( ( VarType[j] & VT_FIXED ) != VT_FIXED && IsEqual( u[j], l[j] ) )
			VarType[j] = VT_FIXED;
}


Bool_T SolvableLP::CreateFromMPS_LP( Bool_T DoScale, VerbLevel Verbosity )
{
	//--------------------------------------------------------------------------
	//	Find out if inherited MPS_LP object contains LP data. Check wether the
	//	current object has not already been converted to standard form.
	//
	assert( OS == OS_FULL || OS == OS_NO_LABELS );

	//--------------------------------------------------------------------------
	//	Construct a work vector for marking the rows for removal from the
	//	matrix.
	//
	Array<Bool_T> ExcludeRows( m, False );

	//--------------------------------------------------------------------------
	//	Find free rows. The first free row will be considered to be the
	//	objective function. Remember its number and store its row label.
	//
	//	If the objective row was given before (e.g., by LP-DIT input) it will
	//	not be changed.
	//
	//	Fill the objective function vector 'c'. Zero the fixed adjustment.
	//
	FindFreeRows( ExcludeRows, ObjRow == -1 ? True : False );
	PutObjectiveInC();
	f = 0.0;

	//--------------------------------------------------------------------------
	//	Remove free rows from the constraint matrix. If required - scale the
	//	column file. Afterwards create the row file.
	//
	UpdateAfterReduction( &ExcludeRows );
	if( DoScale && !Scale( Verbosity ) )
		return False;
	UpdateRowStructure();

	//--------------------------------------------------------------------------
	//	See what's left of the original problem and protest if there's too
	//	little or too much of it.
	//
	if( ObjRow == -1 )
	{
		Error( "Objective function not found" );
		return False;
	}
#ifndef NDEBUG
	else if( m <= 0 )
	{
		Error( "No constraint rows in linear problem" );
		return False;
	}
#endif

	DetectFixedVariables();

	return True;
}
