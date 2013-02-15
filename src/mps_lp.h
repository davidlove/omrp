/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	mps_lp.h
CREATED:			1993.09.16
LAST MODIFIED:		1996.04.16

DEPENDENCIES:		stdtype.h, smartptr.h, simplex.h, sort_lab.h, std_math.h,
					compile.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains class MPS_LP declaration. This class objects holds
exactly the contents of an MPS file (read by one of the functions). The class
also offers procedures to write (and read) binary format MPS files. Finally text
output of an MPS file is made possible (mainly for interfacing with solvers not
supporting the binary format).
	Reading the file is done by "ReadLP" function. Classes derived from "MPS_LP"
may define their versions of "ReadLP", that will first call "MPS_LP::ReadLP",
and then perform some processing (e.g conversion to standard form, presolving,
scaling). The most important data structures of the class will be preserved (in
form but not with respect to contents), some more data - corresponding to the
operations performed - may be added.

	WARNING: Only the first right hand side (RHS) vector, the first range
vector, the first pair of bounds vectors are read. All subsequent ones are
ignored (see: "COMP_READLP_1xxx" macros defined in "compile.h" and desccribed
in "read_lp.cc").

------------------------------------------------------------------------------*/

#ifndef __MPS_LP_H__
#define __MPS_LP_H__


#include <stdio.h>
#include <assert.h>

#ifndef __COMPILE_H__
#	include "compile.h"
#endif

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif


#ifdef SUPPORT_LP_DIT
#	include "lp_dit.h"
#endif

//------------------------------------------------------------------------------
//	Format of input file type. There are 3 file formats recognised
//	automatically. "MPS_LP::ReadLP_Name" returns this information to "ReadLP".
//
enum FF
{
	FF_FIXED_MPS,	// Fixed format MPS text file.
	FF_FREE_MPS,	// Free format MPS text file.
	FF_UNKNOWN		// Unknown file format.
};


//==============================================================================
//
//	Declaration of class MPS_LP.
//
//==============================================================================

class MPS_LP
{
protected:
	//==========================================================================
	//	Nested type definitions.
	//

	//--------------------------------------------------------------------------
	//	Label / name type.
	//
	typedef char Lbl[ LAB_LEN + 1 ];

	//--------------------------------------------------------------------------
	//	Object status information.
	//
	enum ObjectStatus
	{
		OS_EMPTY,			// Empty object - data may be read (either from MPS
							// or binary file).
		OS_NO_LABELS,		// Object contains LP problem data without variable
							// and constraint labels.
		OS_FULL				// Object contains a complete LP problem.
	};
	
	//
	//	End of nested type definitions.
	//==========================================================================

	//--------------------------------------------------------------------------
	//	General problem information.
	//
	ObjectStatus OS;		// Object status (empty / filled).

	Lbl Name;				// Problem name.
	Int_T n, m;				// Dimensions of the problem (numbers of variables,
	Int_T nz;				// constraints and non-zeros respectively).

	//--------------------------------------------------------------------------
	//	ROWS section.
	//
	SortedArrayOfLabels RowLabels;
	Array<Short_T> RowType;	// Row types.

	//--------------------------------------------------------------------------
	//	COLUMNS section (constraint matrix with free rows stored as a file of
	//	packed colmuns).
	//
	SortedArrayOfLabels ColLabels;
	Array<Real_T> ac;		// Matrix non-zeros.
	Array<Int_T> Row;		// Row numbers of non-zeroes.
	Array<Int_T> ColStart;	// Column starts.

	//--------------------------------------------------------------------------
	//	RHS section.
	//
	Lbl RHS_Name;			// RHS vector name (only the first vector in the RHS
							// section is read and stored).
	Array<Real_T> b;		// RHS vector (dense).

	//--------------------------------------------------------------------------
	//	BOUNDS section.
	//
	Lbl BoundsName;			// Bound vector name (only the first vector in the
							// BOUNDS section is read and stored).
	Array<Short_T> VarType;	// Variable types (bit masks; see "lp_codes.h").
	Array<Real_T> l, u;		// Values of lower and upper bounds on variables.

	//--------------------------------------------------------------------------
	//	RANGES section.
	//
	Lbl RangesName;			// RANGE vector name (only the first vector in the
							// RANGES section is read and stored).
	Array<Real_T> r;		// Range vector.

	//--------------------------------------------------------------------------
	//	Objective row number. Used when the LP-DIT input takes place.
	//
	Int_T ObjRow;

	//--------------------------------------------------------------------------
	//	Some statistics: numbers of rows and variables of some types.
	//
	Int_T mE, mG, mL, mR,	// Rows: equality, greater then, less than, range.
		mF, nPL, nFX, nFR,	// Variables: normal, fixed, free, non-positive,
		nMI, nUP;			// bounded.

public:
	//==========================================================================
	//	CONSTRUCTION / DESTRUCTION
	//
	//--------------------------------------------------------------------------
	//	Constructor fills the object with zeros. Destructor deletes all
	//	previously allocated memory using FREE inline function (see:
	//	"myalloc.h"). It is assumed, that all allocations for this class are
	//	done using "C" functions: "malloc", "calloc", "realloc".
	//
	MPS_LP( void );
	virtual ~MPS_LP( void );

	//--------------------------------------------------------------------------
	//	Create matrix *this as a copy of a fragment of the matrix 'Src'
	//	using rows 'RowMin' to 'RowMax' (inclusive) and columns 'ColMin' tp
	//	'ColMax'. Also copy free row 'CostRow'. Skip all other free rows.
	//
	virtual Bool_T CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin,
		Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow );

	virtual Bool_T CreateAsSubmatrix( const MPS_LP &Src, 
		const Array<Bool_T> &IncludeRows, const Array<Bool_T> &IncludeCols );

	//--------------------------------------------------------------------------
	//	This procedure frees all allocated storage of the MPS_LP object. It is
	//	called by destructor and ReadMPS/ReadBin functions (when they fail).
	//
private:
	void FreeStorage( void );

	//==========================================================================
	//	ACCESS functions
	//
	//--------------------------------------------------------------------------
	//	These functions give (read-only) access to problem dimensions, variable
	//	types and simple bounds and the right hand side vector. All functions
	//	are non-virtual!
	//
public:
	const char *GetName( void )										const;
	Int_T	GetN( void )											const;
	Int_T	GetM( void )											const;
	Int_T	GetNZ( void )											const;
	Real_T	GetL( Int_T j )											const;
	Real_T	GetU( Int_T j )											const;
	Real_T	GetB( Int_T i )											const;
	Real_T	GetR( Int_T i )											const;
	Short_T	GetVarType( Int_T j )									const;
	Short_T	GetRowType( Int_T i )									const;
	void	GetRHS( Array<Real_T> &b )								const;
	void	GetConstraintRange( Int_T i, Real_T &min, Real_T &max )	const;

	//--------------------------------------------------------------------------
	//	These functions allow to change the right hand side vector and
	//	the bounds vectors.
	//
	void SetRHS( const Array<Real_T> &rhs );
	void SetRHS_Elem( Int_T row, Real_T val );
	void SetConstraintRange( Int_T i, Real_T min, Bool_T minFinite, Real_T max,
		Bool_T maxFinite );

	void SetL( Int_T j, Real_T l_j, Bool_T Finite = True );
	void SetU( Int_T j, Real_T u_j, Bool_T Finite = True );

	virtual Bool_T SetC( Int_T j, Real_T val );

	Bool_T SetMatrixElement( Int_T row, Int_T col, Real_T val );

	//--------------------------------------------------------------------------
	//	This function gives access to columns stored in the column file. It
	//	returns arrays, which are not a copy of the tables which store the LP
	//	problem. This is not a virtual function!
	//
public:
	void GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row, Int_T &Len )
		const;

	//--------------------------------------------------------------------------
	//	These functions are to be used with extreme caution! They give access
	//	to the whole constyraint matrix in raw form. They breach the security
	//	system of 'Array'/'Ptr' template classes. they are only needed when
	//	interfacing to some obscure languages, like Fortran.
	//
public:
	const Real_T *GetNonZerosByColumns( void )	const;
	const Int_T *GetRowNumbers( void )			const;
	const Int_T *GetColumnStarts( void )		const;


	//==========================================================================
	//	LINEAR PROBLEM DATA FILE INPUT, OUTPUT AND SCALING
	//
	//--------------------------------------------------------------------------
	//	These four procedures read and write LP's from and to files. MPS and
	//	binary file formats are available. For reading "ReadLP" function is
	//	publicly available. It recognizes FIXED or FREE MPS text file formats
	//	and our own binary file format. Then it either invokes a parser (see:
	//	"parsemps.cc") for an MPS file, or calls "ReadBin" if the file is
	//	binary.
	//
	//	FREE and FIXED MPS file formats conform (in as much as possible) to
	//	input file formats of the IBM's MPSX linear programming package.
	//
public:
	Bool_T ReadLP( const char *FileName, VerbLevel Verbosity );
	Bool_T ReadLP( const char *FileName, FILE *LP_File, VerbLevel Verbosity );
	Bool_T WriteMPS( const char *MPS_File, VerbLevel Verbosity ) const;
	virtual Bool_T WriteMPS( FILE *MPS_File, VerbLevel Verbosity ) const;
	
/*#ifdef SUPPORT_LP_DIT 
	Bool_T ReadLP_DIT( const char *FileName, VerbLevel Verbosity );
	virtual Bool_T WriteLP_DIT( const char *FileName, VerbLevel Verbosity )
		const;
#endif */

	SortedArrayOfLabels &RevealRowLabels( void );
	SortedArrayOfLabels &RevealColumnLabels( void );

	const char *RevealProblemName( void ) const;
	const char *RevealRHS_Name( void ) const;
	const char *RevealBoundsName( void ) const;
	const char *RevealRangesName( void ) const;

protected:
	virtual Array<Int_T> GetRowLen( void );

private:
	void CheckLP( void );

public:
	//==========================================================================
	//	SEMANTIC ACTIONS FOR PROBLEM INPUT FROM A TEXT FILE (MPS FORMAT)
	//
	//--------------------------------------------------------------------------
	//	The functions named below are semantic actions called by the parser
	//	during MPS file processing. The first two actions are called by a
	//	function which reads the first line of the file in order to recognise
	//	the file's format and problem name.
	//	The remaining functions are called by a parser, which processes the MPS
	//	file body. For more detail - see "read_lp.cc".
	//
	//	Attempt to execute them on an object that is not empty will cause the
	//	program to be aborted. Otherwise (i.e. if they are called on an object
	//	already holding data) the behavior of the program is undefined.
	//
	void SetLP_Name( const char *lab0 );

	void BeginRows( void  );
	void NewRow( Short_T RowType, const char *lab0 );
	void BeginColumns( void );
	void NewNonZero( const char *lab0, const char *lab1, double val );
	void BeginRHS( void );
	void NewRHS( const char *lab0, const char *lab1, double val );
	void BeginRanges( void );
	void NewRange( const char *lab0, const char *lab1, double val );
	void BeginBounds( void );
	void NewBound( Short_T BoundType, const char *lab0, const char *lab1,
		double val );
	void Endata( void );

	//--------------------------------------------------------------------------
	//	A section of public functions to be used for LP-DIT input/output
	//	(have to be public to be accessed through C code).
	//
#ifdef SUPPORT_LP_DIT

public:
	//
	//	Linear problem input functions.
	//
	virtual void lpo_par( LP_HEAD &h );
	virtual void lpo_specs( LP_HEAD &h, char **str );
	virtual void lpo_cols( LP_HEAD &h, LP_VAR *v );
	virtual void lpo_rows( LP_HEAD &h, LP_VAR *v );
	virtual void lpo_vect( LP_HEAD &h, LP_MAT *v );

	//
	//	Linear problem output functions.
	//
	virtual void lpi_par( LP_HEAD &h );
	virtual void lpi_specs( LP_HEAD &h, char **str );
	virtual void lpi_cols( LP_HEAD &h, LP_VAR *v );
	virtual void lpi_rows( LP_HEAD &h, LP_VAR *v );
	virtual void lpi_vect( LP_HEAD &h, LP_MAT *v );
#endif
};

//==============================================================================
//
//	End of declaration of class MPS_LP.
//
//==============================================================================

//==============================================================================
//
//	INLINE DEFINITIONS OF SOME OF THE MEMBER FUNCTIONS.
//
//==============================================================================


/*------------------------------------------------------------------------------

	MPS_LP::~MPS_LP( void ) :

PURPOSE:
	This is class MPS_LP destructor. It frees all previously allocated storage
using private FreeStorage function (defined below).

------------------------------------------------------------------------------*/

inline
MPS_LP::~MPS_LP( void )
	{ FreeStorage(); }


//------------------------------------------------------------------------------
inline
const char *MPS_LP::GetName( void )
const
	{ return Name; }


//------------------------------------------------------------------------------
inline
Int_T MPS_LP::GetM( void )
const
	{ return m; }


//------------------------------------------------------------------------------
inline
Int_T MPS_LP::GetN( void )
const
	{ return n; }


//------------------------------------------------------------------------------
inline
Int_T MPS_LP::GetNZ( void )
const
	{ return nz; }


//------------------------------------------------------------------------------
inline
Real_T MPS_LP::GetL( Int_T j )
const
{
	assert( j >= 0 && j < n );

	return l[j];
}


//------------------------------------------------------------------------------
inline
Real_T MPS_LP::GetU( Int_T j )
const
{
	assert( j >= 0 && j < n );

	return u[j];
}


//------------------------------------------------------------------------------
inline
Real_T MPS_LP::GetB( Int_T i )
const
{
	assert( i >= 0 && i < m );

	return b[i];
}


//------------------------------------------------------------------------------
inline
Real_T MPS_LP::GetR( Int_T i )
const
{
	assert( i >= 0 && i < m );

	return r[i];
}


//------------------------------------------------------------------------------
inline
Short_T MPS_LP::GetVarType( Int_T j )
const
{
	assert( j >= 0 && j < n );

	return VarType[j];
}


//------------------------------------------------------------------------------
inline
Short_T MPS_LP::GetRowType( Int_T i )
const
{
	assert( i >= 0 && i < m );

	return RowType[i];
}


/*------------------------------------------------------------------------------

	inline
	void MPS_LP::GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row,
		Int_T &Len )
		const

PURPOSE:
	This function gives access to columns of the LP problem. It returns pointers
INTO the tables which store the LP problem. That's why the "A", "Row" and "Col"
are pointers to "const" data.

PARAMETERS:
	Int_T j
		Column number (in range 0 : "n" ) of column to be accessed.

	Ptr<Real_T> &a
		Array in which column "j" non-zeros start address will be stored. It's
		value on entry is irrelevant.

	Ptr<Int_T> &row
		Array in which column "j" non-zero pattern start address will be stored.
		It's value on entry is irrelevant.


RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
void MPS_LP::GetColumn( Int_T j, Ptr<Real_T> &a, Ptr<Int_T> &row, Int_T &Len )
	const
{
	assert( j >= 0 && j < n );

	Len = ColStart[j+1];
	size_t Start = ColStart[j];

	a.ExtractFragment(   ac,  Start, (size_t) Len );
	row.ExtractFragment( Row, Start, (size_t) Len );

	Len -= Start;
}


/*------------------------------------------------------------------------------

	void MPS_LP::GetRHS( Array<Real_T> &rhs ) const

PURPOSE:
	x

PARAMETERS:
	x

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
void MPS_LP::GetRHS( Array<Real_T> &rhs )
	const
	{ rhs.Copy( b, m, m ); }


/*------------------------------------------------------------------------------

	(inline)
	void MPS_LP::SetRHS( const Array<Real_T> &rhs )
	void MPS_LP::SetRHS_Elem( Int_T row, Real_T val );

PURPOSE:
	x

PARAMETERS:
	x

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
void MPS_LP::SetRHS( const Array<Real_T> &rhs )
	{ b.Copy( rhs, m, m ); }


inline
void MPS_LP::SetRHS_Elem( Int_T row, Real_T val )
{
	assert( row >= 0 && row < m );
	assert( val > -INFINITY && val < INFINITY );

	b[row] = IsNonZero( val ) ? val : 0.0;
}


inline
Bool_T MPS_LP::SetC( Int_T, Real_T )
{ return False; }

/*------------------------------------------------------------------------------

	(inline)
	SortedArrayOfLabels &MPS_LP::RevealRowLabels( void ) 
	SortedArrayOfLabels &MPS_LP::RevealColumnLabels( void )
	const char *MPS_LP::RevealProblemName( void ) const
	const char *MPS_LP::RevealRHS_Name( void ) const
	const char *MPS_LP::RevealBoundsName( void ) const
	const char *MPS_LP::RevealRangesName( void ) const

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
SortedArrayOfLabels &MPS_LP::RevealRowLabels( void )
	{ return RowLabels; }


inline
SortedArrayOfLabels &MPS_LP::RevealColumnLabels( void )
	{ return ColLabels; }


inline
const char *MPS_LP::RevealProblemName( void )
	const
	{ return Name; }


inline
const char *MPS_LP::RevealRHS_Name( void )
	const
	{ return RHS_Name; }


inline
const char *MPS_LP::RevealBoundsName( void )
	const
	{ return BoundsName; }


inline
const char *MPS_LP::RevealRangesName( void )
	const
	{ return RangesName; }


/*------------------------------------------------------------------------------

	(inline)
	Real_T *MPS_LP::GetNonZerosByColumns( void )
	Int_T *MPS_LP::GetRowNumbers( void )
	Int_T *MPS_LP::GetColumnStarts( void )

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
const Real_T *MPS_LP::GetNonZerosByColumns( void )
	const
{
	return ( nz > 0 ) ? & ac[0] : (const Real_T *) NULL;
}


inline
const Int_T *MPS_LP::GetRowNumbers( void )
	const
{
	return ( nz > 0 ) ? & Row[0] : (const Int_T *) NULL;
}


inline
const Int_T *MPS_LP::GetColumnStarts( void )
	const
{
	return & ColStart[0];
}



//==============================================================================
//
//	END OF THE INLINE DEFINITIONS.
//
//==============================================================================

#endif
