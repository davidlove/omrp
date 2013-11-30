/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	read_lp.cpp
CREATED:			1993.09.16
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		stdtype.h, std_tmpl.h, error.h, myalloc.h, mps_lp.h,
					lexer.h, parsemps.h, lp_codes.h, simplex.h, print.h,
					std_math.h
					<stdio.h>, <string.h>, <math.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This source file contains functions used for reading a linear problem (LP)
from an input file. MPS and binary formats are supported. In addition to one
publicly available function ("ReadLP") a large number of other functions (most
of which are semantic actions called by MPS file parser) is also included.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Bool_T MPS_LP::ReadLP( FILE *LP_File, RL ReadLabels = RB_ON )

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS FROM COMPILE.H AND THEIR MEANING:
	COMP_READLP_1RHS		- single RHS vector expected,
	COMP_READLP_1RANGE		- single range vector expected,
	COMP_READLP_1BOUND 		- single bound vector expected,

------------------------------------------------------------------------------*/

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif
#ifndef __PARSEMPS_H__
#	include "parsemps.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif

//==============================================================================
//
//	Static function prototype.
//
//==============================================================================

//------------------------------------------------------------------------------
//	Function for perforimng label comparison during sorting. It takes const
//	void* instead of const Label * (or &), because it is passed to "qsort" and
//	"bsearch" ANSI C functions, which perform sorting and binary search.
//
//static int CompareLabels( const void *l1, const void *l2 );
//static int ComparePositions( const void *l1, const void *l2 );

//==============================================================================
//
//	End of static function prototype.
//
//==============================================================================

static VerbLevel StatVerbosity;

/*------------------------------------------------------------------------------

	Bool_T MPS_LP::ReadLP( const char *FileName, FILE *LP_File,
		VerbLevel Verbosity )

PURPOSE:
	This function recognizes one of three possible file types of the file read 
from an open stream "LP_File". It then calls the "ReadMPS" file reading
procedure which reads the rest of the data and stores it in MPS_LP object. If
the file type is not recognized, procedure returns failure status (False). If
more than 20 errors or more than 20 warnings are encountered, file processing
is aborted with a call to "FatalError".
	This procedure sets object status (MPS_LP::OS) to OS_FULL when it is
successful. Otherwise the object is emptied and object status is set to
OS_EMPTY. You have to take account of the return code.

PARAMETERS:
	FILE *MPS_FIle
	Open stream from which the data will be read. The stream has to be opened in
"read" mode. If it is not, a fatal error is reported.

RETURN VALUE:
	Boolean success status (True on success, False otherwise). When there are
too many errors and warnings, program is aborted and (obviously) no value is
returned.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T MPS_LP::ReadLP( const char *FileName, VerbLevel Verbosity )
{
	FILE *fp = fopen( FileName, "rt" );
	if( fp == NULL )
	{
		if( Verbosity > V_NONE )
			Print( "Unable to open input file '%s'.\n", FileName );
		return False;
	}

	Bool_T result = ReadLP( FileName, fp, Verbosity );
	fclose( fp );
	return result;
}


Bool_T MPS_LP::ReadLP( const char *FileName, FILE *LP_File, // )
	VerbLevel Verbosity )
{
	mE = mG = mL = mF = mR = nPL = nFX = nFR = nMI = nUP = 0;
	StatVerbosity = Verbosity;

	assert( OS == OS_EMPTY );

	ResetWarningCount();
	ResetErrorCount();

	MaxWarnCount( 20 );
	MaxErrCount( 20 );

	Lexer::SetInputStream( FileName, LP_File );
	SetLP_TargetObject( this );

	if( Verbosity >= V_HIGH )
		Print( "Reading input file: %s\n", FileName );

	FF FileType = GetNameLine();

	switch( FileType )
	{
	case FF_FIXED_MPS:			// Fixed or free format MPS file.
	case FF_FREE_MPS:			// Read file in MPS format (use parser).
		if( !GetMPS_Body() || ErrorCount() > 0 )
			goto error;
		if( Verbosity == V_HIGH )
			Print( "Finished reading.\n" );
		break;

	case FF_UNKNOWN:			// Unknown file format.
		goto error;
	}

	CheckLP();
	if( ErrorCount() > 0 ) goto error;

	OS = OS_FULL;
	MaxWarnCount();
	MaxErrCount();

//success:
	if( Verbosity >= V_LOW )
	{
		Print(
			/* FORMAT STRING */
			"\nLinear problem '%s' statistics (as read from '%s').\n"
			"\t%-18s%10d (%d E, %d L, %d G, %d N)\n"
			"\t%-18s%10d\n"
			"\t%-18s%10d (%d PL, %d UP, %d FX, %d FR, %d MI)\n"
			"\t%-18s%10ld\n"
			"\t%-18s%10.4f %%\n",

			/* DATA */
			Name, FileName,
			"No. of rows:",			(int) m, (int) mE, (int) mL, (int) mG,
									(int) mF,
			"No. of ranges:",		(int) mR,
			"No. of variables:",	(int) n, (int) nPL, (int) nUP, (int) nFX,
									(int) nFR, (int) nMI,
			"No. of non-zeros:",	(long int) nz,
			"Density:",				100.0 * double( nz ) / double( m ) /
									double( n )
			);
	}
	else if( Verbosity >= V_LINE )
	{
		Print(
			/* FORMAT STRING */
			"| %-8s | %10d | %10d | %10ld | %10.4f |", 

			/* DATA */
			Name, (int) m, (int) n, (long int) nz,
			100.0 * double( nz ) / double( m ) / double( n )
			);
	}
	return True;

error:
	Error( "Errors encountered while reading input file. Cannot proceed." );
	FreeStorage();
	return False;
}


/*------------------------------------------------------------------------------
	Semantic actions
------------------------------------------------------------------------------*/

#define INIT_ROWS		500
#define INIT_VARS		1000
#define INIT_NON_ZEROS	3 * INIT_VARS

static Int_T mMax, nMax;		// Allocated dimensions of the tables.
static Int_T nzMax;

void MPS_LP::BeginRows( void  )
{
	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tReading ROWS\n" );
	mMax = INIT_ROWS;
	m = 0;
	RowType.Resize( mMax );
}


void MPS_LP::NewRow( Short_T _RowType, const char *lab0 )
{
	assert( OS == OS_EMPTY );

	if( m >= mMax )
	{
		mMax += mMax / 2 + 1;
		assert( mMax > 0 );
		RowType.Resize( mMax );
	}

	if( RowLabels.AddLabel( lab0 ) )
		Error( "Duplicate label: %s", lab0 );

	RowType[ m ]			= _RowType;
	m++;

	assert( m > 0 );

	switch( _RowType )
	{
	case RT_EQ: mE++; break;
	case RT_LE: mL++; break;
	case RT_GE: mG++; break;
	case RT_FR: mF++; break;
	}
}


void MPS_LP::BeginColumns( void )
{
#ifdef NDEBUG
	RowLabels.SortLabels();
#else
	assert( !RowLabels.SortLabels() );
#endif

	if( RowLabels.Duplicates() > 0 )
		FatalError( "Duplicate row labels detected." );

	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tReading COLUMNS\n" );

	//--------------------------------------------------------------------------
	//	If necessary shrink overgrown "RowLab" and "RowType" tables.
	//	Then allocate space for a sortable array of labels.
	//
	if( m < mMax )
	{
		mMax = m;
		RowType.Resize( mMax );
	}

	//--------------------------------------------------------------------------
	//	Now is time to allocate the structures for matrix non-zeros.
	//
	nMax = Max( Int_T( INIT_VARS ), Int_T( 2 * m ) );
	nzMax = Max( Int_T( INIT_NON_ZEROS ), Int_T( 4 * m ) );
	nz = n = 0;

	ac.Resize( (size_t) nzMax );
	Row.Resize( (size_t) nzMax );
	ColStart.Resize( nMax + 1 );
}


void MPS_LP::NewNonZero( const char *ColumnLab, const char *RowLab, double Val )
{
	assert( OS == OS_EMPTY );

	//
	// See if a new column starts here.
	//
	if( !n || ( ColumnLab &&
		strncmp( ColumnLab, ColLabels.FindLabel( n-1 ), LAB_LEN ) ) )
	{
		if( n >= nMax )
		{
			nMax += nMax / 2 + 1;
			assert( nMax > 0 );

			ColStart.Resize( nMax + 1 );
		}

		if( ColLabels.AddLabel( ColumnLab ) )
			Error( "Duplicate column label: %s", ColumnLab );

		ColStart[ n ]			= nz;
		n++;

		assert( n > 0 );
	}

	if( IsZero( Val ) )
		return;

	if( nz >= nzMax )
	{	
		nzMax += nzMax / 2 + 1;
		assert( nzMax > 0 );

		ac.Resize( (size_t) nzMax );
		Row.Resize( (size_t) nzMax );
	}

	Int_T Num = (Int_T) RowLabels.FindLabel( RowLab );

	if( Num < 0 )
	{
		Error( "Row label %s found in COLUMNS missing in ROWS", RowLab );
		return;
	}

	ac[ (size_t) nz ]	= Val;
	Row[ (size_t) nz ]	= Num;
	nz++;

	assert( nz > 0 );
}


void MPS_LP::BeginRHS( void )
{
#ifdef NDEBUG
	ColLabels.SortLabels();
#else
	assert( !ColLabels.SortLabels() );
#endif

	if( RowLabels.Duplicates() > 0 )
		FatalError( "Duplicate column labels detected." );

	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tReading RHS\n" );

	//--------------------------------------------------------------------------
	//	If necessary shrink overgrown "ColLab", "ColStart", "a" and "Row"
	//	tables. Set the last element of "ColStart". Allocate and fill "ColLen".
	//
	nMax = n;
	ColStart.Resize( nMax + 1 );

	nzMax = nz;
	ac.Resize( (size_t) nzMax );
	Row.Resize( (size_t) nzMax );

	ColStart[n] = nz;

	//--------------------------------------------------------------------------
	//	Now is time to allocate and initialize the structures for the right hand
	//	side vector, range and bound vector.
	//
	*RHS_Name = *RangesName = *BoundsName = '\0';
	b.Resize( mMax );			b.Fill( 0.0,					mMax );
	r.Resize( mMax );           r.Fill( INFINITY,				mMax );
	l.Resize( nMax );           l.Fill( 0.0,					nMax );
	u.Resize( nMax );           u.Fill( INFINITY,				nMax );
	VarType.Resize( nMax );		VarType.Fill( VTM_UNDEFINED,	nMax );
}


void MPS_LP::NewRHS( const char *lab0, const char *RowLab, double Val )
{
	assert( OS == OS_EMPTY );

	if( lab0 && !*RHS_Name )				// First RHS non-zero.
	{
		strncpy( RHS_Name, lab0, LAB_LEN );
		RHS_Name[ LAB_LEN ] = '\0';
	}
#ifdef COMP_READLP_1RHS
	else if( lab0 && strncmp( RHS_Name, lab0, LAB_LEN ) )
		//
		// Next (alternative) RHS vector is ignored.
		//
		return;
#endif

	if( IsZero( Val ) )						// Ignore small elements.
		return;

	//--------------------------------------------------------------------------
	//	If we can't ignore it, we'll attempt to store it.
	//
	Int_T Num = (Int_T) RowLabels.FindLabel( RowLab );

	if( Num < 0 )
	{
		Error( "Row label %s found in RHS missing in ROWS", RowLab );
		return;
	}

	b[ Num ] = Val;
}


void MPS_LP::BeginRanges( void )
{
	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tReading RANGES\n" );
}


void MPS_LP::NewRange( const char *lab0, const char *RowLab, double Val )
{
	assert( OS == OS_EMPTY );

	if( lab0 && !*RangesName )				// First range vector.
	{
		strncpy( RangesName, lab0, LAB_LEN );
		RangesName[ LAB_LEN ] = '\0';
	}
#ifdef COMP_READLP_1RANGE
	else if( lab0 && strncmp( RangesName, lab0, LAB_LEN ) )
		//
		// Next (alternative) range vector is ignored.
		//
		return;
#endif

	if( IsZero( Val ) )						// Ignore small elements.
		return;

	//--------------------------------------------------------------------------
	//	If we can't ignore it, we'll attempt to store it.
	//
	Int_T Num = (Int_T) RowLabels.FindLabel( RowLab );

	if( Num < 0 )
	{
		Error( "Row label %s found in RANGES missing in ROWS", RowLab );
		return;
	}

	if( Val < 0.0e0 )
		Warning( "Negative range for row %s; absolute value taken", RowLab );
	r[ Num ]		= fabs( Val );
	RowType[ Num ]	|= RT_RNG;
	mR++;
}


void MPS_LP::BeginBounds( void )
{
	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tReading BOUNDS\n" );
}


void MPS_LP::NewBound( Short_T BoundType, const char *lab0, // )
	const char *lab1, double Val )
{
	assert( OS == OS_EMPTY );

	if( !*BoundsName )						// First bound.
	{
		strncpy( BoundsName, lab0, LAB_LEN );
		BoundsName[ LAB_LEN ] = '\0';
	}
#ifdef COMP_READLP_1BOUND
	else if( strncmp( BoundsName, lab0, LAB_LEN ) )
		//
		// Next (alternative) bound vector is ignored.
		//
		return;
#endif

	//--------------------------------------------------------------------------
	//	If we can't ignore it, we'll attempt to store it.
	//
	Int_T Num = (Int_T) ColLabels.FindLabel( lab1 );

	if( Num < 0 )
	{
		Error( "Column label %s found in BOUNDS missing in COLUMNS", lab1 );
		return;
	}

	Int_T pos = Num;
	Short_T &Bnd = VarType[ pos ];

	//--------------------------------------------------------------------------
	//	In the first place we have to see if the bound just read is
	//	a)	correct and
	//	b)	not in conflict with previously declared bounds.
	//
	//	If this is the first time this bound is mentioned - any bound will do.
	//	Otherwise the new bound is checked against those previously declared.
	//	Bound redefinition is not allowed (because the author does not know
	//	the redefinition rules). So far only redefinition from posession of a
	//	single bound to fixed bound is allowed.
	//
	if( Bnd & VTM_UNDEFINED )
		Bnd = BoundType;
	else
	{
		switch( Bnd )
		{
		case VTM_MI:
			if( BoundType & ( VTM_FR | VTM_MI | VTM_PL | VTM_LO | VTM_FX ) )
				goto error;
			else
				break;

		case VTM_PL:
			if( BoundType & ( VTM_FR | VTM_MI | VTM_PL | VTM_UP | VTM_FX ) )
				goto error;
			else
				break;

		case VTM_LO:
			if( BoundType & VTM_FX )		// Redefinition LO -> FX
				Bnd = 0;
			else if( BoundType & ( VTM_FR | VTM_MI | VTM_LO ) )
				goto error;
			break;

		case VTM_UP:
			if( BoundType & VTM_FX )		// Redefinition UP -> FX
				Bnd = 0;
			else if( BoundType & ( VTM_FR | VTM_PL | VTM_UP ) )
				goto error;
			break;

		case VTM_FX:
		case VTM_FR:	
			goto error;

		default:
			abort();
		}

		Bnd |= BoundType;
	}

	//--------------------------------------------------------------------------
	//	If we got this far, we are sure the bound is OK. We may even store the
	//	numeric values now. We DON'T check infeasibility (i.e. l[j] > u[j] ).
	//
	switch( BoundType )
	{
	case VTM_MI:	if( !( Bnd & VTM_UP ) ) u[ pos ] = 0.0e0;
	case VTM_FR:	l[ pos ] = -INFINITY;
	case VTM_PL:	break;

	case VTM_LO:	l[ pos ] = Val;				break;
	case VTM_UP:	u[ pos ] = Val;				break;
	case VTM_FX:	l[ pos ] = u[ pos ] = Val;	break;
	}

	return;

error:
	Error( "Conflicting or duplicate bound in input stream" );
}


void MPS_LP::Endata( void )
{
	assert( OS == OS_EMPTY );

	if( StatVerbosity >= V_HIGH )
		Print( "\tFound ENDATA\n" );

	//--------------------------------------------------------------------------
	//	Translate "VTM_XX" variable types into "VT_XXXX" variable types, that
	//	will later be used by solvers. A the same time:
	//	*	set bounds that have not been set in BOUNDS section to PL.
	//
	for( Int_T i = 0; i < n; i++ )
		if( VarType[i] & VTM_UNDEFINED )
		{
			nPL++;
			VarType[i] = VT_NORM;
		}
		else if( VarType[i] & VTM_MI )
		{
			nMI++;
			VarType[i] = VT_MI;
		}
		else if( VarType[i] & VTM_FR )			// Declared free.
		{
			nFR++;
			VarType[i] = VT_FREE;
		}
		else if( VarType[i] & VTM_FX )			// Declared fixed.
		{
			nFX++;
			VarType[i] = VT_FIXED;
		}
		else if( VarType[i] & VTM_UP )		// Not fixed, not non-positive
		{										// and has upper bound,
			nUP++;								// then it must be bounded.
			VarType[i] = VT_BOUNDED;
		}
		else									// In all other cases it must
		{										// be a normal variable.
			nPL++;
			VarType[i] = VT_NORM;
		}
}


void MPS_LP::SetLP_Name( const char *LP_Name )
{
	assert( OS == OS_EMPTY );

	strncpy( Name, LP_Name, LAB_LEN );
	Name[ LAB_LEN ] = '\0';
}
