/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem postsolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	read_act.cpp
CREATED:			1995.10.06
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		std_tmpl.h, solution.h, compile.h, smartptr.h,
					error.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h,
					stdtype.h, myalloc.h, simplex.h, postsolv.h, my_defs.h
					<assert.h>, <stdio.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif
#ifndef __POSTSOLV_H__
#	include "postsolv.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T Postsolver::ReadActions( const char *fname )
	Bool_T Postsolver::ReadActions( FILE *fp )

PURPOSE:
	Reads the presolver actions from a file.

PARAMETERS:
	FILE *fp
		Input file (assumed to be opened for reading).

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::ReadActions( const char *fname )
{
	assert( fname != NULL && *fname != '\0' );

	FILE *fp = fopen( fname, "rt" );
	if( fp == NULL ) return False;

	Bool_T result = ReadActions( fp );
	fclose( fp );
	return result;
}


Bool_T Postsolver::ReadActions( FILE *fp )
{
	//--------------------------------------------------------------------------
	//	Initialize the object as an empty object. Set up the read data buffer,
	//	work variables etc. Set up the lexical analysis module.
	//
	FreeMemory();

	Lexer::SetInputStream( "action file", fp );
	ResetErrorCount();
	MaxErrCount( 20 );

	Int_T ReducedCols = 0;

	//--------------------------------------------------------------------------
	//	Read the first line.
	//
	if( !Lexer::GetKeyword( "PRESOLVER" ) || !Lexer::GetSpace() ||
		!Lexer::GetKeyword( "REPORT" ) )
	{
		Error( "Title line incorrect or missing." );
		goto error;
	}

	if( !Lexer::GetNewline() )
	{
		Error( "Extra input after the title line ignored." );
		Lexer::GetNewline( True );
	}

	//--------------------------------------------------------------------------
	//	Read the ROWS section and create the row exclusion array. Then read the
	//	COLUMNS section and do the same.
	//
	if( /* !GetRowsSection() || */ !GetColumnsSection() )
		goto error;

	assert( n >= 0 );

	{ for( Int_T j = 0; j < n; j++ ) if( ExcludeCols[j] ) ReducedCols++; }

	//--------------------------------------------------------------------------
	//	Read the actions in loop. Then read the "ENDATA" and finish processing.
	//
	while( GetAction( ReducedCols ) || GetAdjustment() )
		;

	if( ReducedCols > 0 )
	{
		Error( "Some postsolve actions missing." );
		goto error;
	}

	if( !Lexer::GetKeyword( "ENDATA" ) )
	{
		Error( "ENDATA expected." );
		goto error;
	}

	if( ErrorCount() > 0 )
		goto error;

	return True;

error:
	FreeMemory();
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetRowsSection( void )
	Bool_T Postsolver::GetColumnsSection( void )

PURPOSE:
	Two very similar functions plus a general engine (the third function) for
reading the two initial sections of the presolve action file:
a)	rows,
b)	columns.
Those two sections declare the solution vector dimension as well as which
variables were eliminated and which are supposed to be in present in the
solution vector of the reduced problem.

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

/*
Bool_T Postsolver::GetRowsSection( void )
{
	//--------------------------------------------------------------------------
	//	Read the indicator line for the section and store the number of row
	//	labels.
	//
	if( !Lexer::GetKeyword( "ROWS" ) )
	{
		Error( "ROWS section header missing." );
		return False;
	}

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Number of row labels missing." );
		return False;
	}

	m = (Int_T) Lexer::Number;

	if( !Lexer::GetNewline() )
	{
		Error( "Extra input after the number of row labels ignored." );
		Lexer::GetNewline( True );
	}

	if( m <= 0 )
	{
		Error( "Invalid number of rows (%d).", (int)m );
		return False;
	}

	//--------------------------------------------------------------------------
	//	Construct the row exclusion array.
	//

	return GetExclusionection( RowLabels, m, "ROWS" );
}
*/


Bool_T Postsolver::GetColumnsSection( void )
{
	//--------------------------------------------------------------------------
	//	Read the indicator line for the section and store the number of column
	//	labels.
	//
	if( !Lexer::GetKeyword( "COLUMNS" ) )
	{
		Error( "COLUMNS section header missing." );
		return False;
	}

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Number of column labels missing." );
		return False;
	}

	n = (Int_T) Lexer::Number;

	if( !Lexer::GetNewline() )
	{
		Error( "Extra input after the number of column labels ignored." );
		Lexer::GetNewline( True );
	}

	if( n <= 0 )
	{
		Error( "Invalid number of columns (%d).", (int)n );
		return False;
	}

	//--------------------------------------------------------------------------
	//	Resize the column exclusion array.
	//
	ExcludeCols.Resize( n+1 );
	ExcludeCols.Fill( False, n+1 );
	ExcludeCols[n] = True;

	return GetExclusionSection( ExcludeCols, n, "COLUMNS" );
}


Bool_T Postsolver::GetExclusionSection( Array<Bool_T> &Exclude, Int_T k, // )
	const char *Section )
{
	for( Int_T cnt = 0; Lexer::GetSpace(); cnt++ )
	{
		Lexer::GetNumeric();

		Int_T num = (Int_T)Lexer::Number;

		if( num < 0 || num >= k )
		{
			Error( "Index out of range in %s section: %d", Section, (int)num );
			return False;
		}

		if( Exclude[num] )
		{
			Error( "Number %d repeated in section %s", (int)num, Section );
			return False;
		}

		Exclude[num] = True;
		Lexer::GetNewline( True );
	}

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetAdjustment( void )

PURPOSE:
 	Reads in the section defining the fixed adjustment.

PARAMETERS:
 	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetAdjustment( void )
{
	if( !Lexer::GetKeyword( "ADJUST" ) ) return False;
	Lexer::GetNewline( True );

	if( !Lexer::GetSpace() || !Lexer::GetKeyword( "VALUE" ) ||
		!Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Invalid fixed adjustment note." );
		SkipToNextSection();
		return False;
	}

	FixedAdjustment( Lexer::Number );
	Lexer::GetNewline( True );

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetAction( Int_T &nn )

PURPOSE:
	Reads in a section corresponding to a snigle presolve action.

PARAMETERS:
	Int_T &n, Int_T &nn
		The number of variables read so far an the number remaining.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetAction( Int_T &nn )
{
 	//--------------------------------------------------------------------------
 	//	See if this section actually defines an action.
 	//
	if( !Lexer::GetKeyword( "ACTION" ) ) return False;

 	//--------------------------------------------------------------------------
 	//	Read in the variable label and see if it corresponds to a known label.
 	//
	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Variable index missing in the action description line." );
		SkipToNextSection();
		return False;
	}

	Int_T var = (Int_T)Lexer::Number;

	if( var >= n || var < 0 )
	{
		Error( "Invalid variable index: %d.", (int)var );
		SkipToNextSection();
		return False;
	}
	else if( !ExcludeCols[var] )
	{
		Error( "A postsolve action declared for non-presolved variable." );
		SkipToNextSection();
		return False;
	}

	Lexer::GetSpace();

 	//--------------------------------------------------------------------------
 	//	Read in the section body.
 	//
	if( !GetFixAction( nn, var ) &&
		!GetExplicitSlackAction( nn, var ) &&
		!GetFreeSingletonAction( nn, var ) )
		return False;

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetFixAction( Int_T &nn, Int_T var )

PURPOSE:
 	Reads in a section corresponding to variable fixing.

PARAMETERS:
 	Int_T var
 		Variable number (read in from the section's first line).
 
 	Int_T &n, Int_T &nn
 		The number of variables read so far an the number remaining.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetFixAction( Int_T &nn, Int_T var )
{
	if( !Lexer::GetKeyword( "FIX" ) ) return False;
	Lexer::GetNewline( True );

	nn--;

	if( !Lexer::GetSpace() || !Lexer::GetKeyword( "VALUE" ) ||
		!Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Invalid line after ACTION FIX section header." );
		return False;
	}

	Lexer::GetNewline( True );

	VariableFixing( var, Lexer::Number );
	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetRowData( Real_T &val, Int_T &len,
		Array<Real_T> &a, Array<Int_T> &ind )

PURPOSE:
	Reads a part of the explicit slack or free singleton column section. Reads
the keyword "VALUE" together with the numeric value which follows, and then the
entire row.

PARAMETERS:
	Real_T &val
		The value read (matrix coefficient of the variable in the given row).

	Int_T len, Array<Real_T> &a, Array<Int_T> &ind
		The row: its length, non-zeros and their positions.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetRowData( Real_T &val, Int_T &len, // )
	Array<Real_T> &a, Array<Int_T> &ind )
{
	//--------------------------------------------------------------------------
	//	Read the removed variable's coefficient in the row in question.
	//
	if( !Lexer::GetSpace() )
	{
		Error( "Premature end of section." );
		return False;
	}

	if( !Lexer::GetKeyword( "VALUE" ) || !Lexer::GetSpace() ||
		!Lexer::GetNumeric() )
	{
		Error( "Variable's coefficient expected." );
		return False;
	}

	val = Lexer::Number;

	Lexer::GetNewline( True );

	//--------------------------------------------------------------------------
	//	Read the row's coefficients.
	//
	{
		Int_T max_len = 10;
		a.Resize( max_len);		a.Fill( 0.0, max_len );
		ind.Resize( max_len );	ind.Fill( -1, max_len );

		for( len = 0; Lexer::GetSpace(); len++ )
		{
			if( !Lexer::GetKeyword( "COEFF" ) || !Lexer::GetSpace() )
			{
				Error( "Matrix coefficient expected." );
				return False;
			}

			if( len >= max_len )
			{
				Int_T new_max_len = Int_T( Max( max_len + 10, 3*max_len/2 ) );

				a.Resize( new_max_len );
				ind.Resize( new_max_len );

				a.Fill( 0.0, new_max_len, max_len );
				ind.Fill( -1, new_max_len, max_len );

				max_len = new_max_len;
			}

			if( !Lexer::GetNumeric() || !Lexer::GetSpace() )
			{
				Error( "Variable index expected." );
				return False;
			}
			ind[len]	= (Int_T) Lexer::Number;

			if( !Lexer::GetNumeric() )
			{
				Error( "Matrix coefficient expected." );
				return False;
			}

			a[len]		= Lexer::Number;

			if( ind[len] < 0 || ind[len] >= n )
			{
				Error( "Variable index out of range: %d", ind[len] );
				return False;
			}
			Lexer::GetNewline( True );
		}
	}
	
	return True;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetExplicitSlackAction( Int_T &nn, Int_T var )

PURPOSE:
 	Reads in a section corresponding to explicit slack removal.

PARAMETERS:
 	Int_T var
 		Variable number (read in from the section's first line).
 
RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetExplicitSlackAction( Int_T &nn, Int_T var )
{
	if( !Lexer::GetKeyword( "EXPL_SLACK" ) ) return False;
	Lexer::GetNewline( True );

	nn--;

	Real_T bl	= -INFINITY,
		bu		= +INFINITY;
	Short_T rt	= RT_UNDEFINED;
	Real_T l	= -INFINITY,
		u		= +INFINITY;
	Short_T vt	= VTM_UNDEFINED;
	Real_T sl	= 0.0;

	//--------------------------------------------------------------------------
	//	Read the row type and the row activity bounds.
	//
	if( !Lexer::GetSpace() )
	{
		Error( "Premature end of the ACTION EXPL_SLACK section." );
		return False;
	}

	if( Lexer::GetKeyword( "GE" ) )
		rt = VT_NORM;
	else if( Lexer::GetKeyword( "LE" ) )
		rt = VT_MI;
	else if( Lexer::GetKeyword( "EQ" ) )
		rt = VT_FIXED;
	else if( Lexer::GetKeyword( "RG" ) )
		rt = VT_BOUNDED;
	else
	{
		Error( "Unrecognized row type in presolver action file." );
		goto error;
	}

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Numeric value(s) expected." );
		goto error;
	}

	switch( rt )
	{
	case VT_NORM:		bl = Lexer::Number;			break;
	case VT_MI:			bu = Lexer::Number;			break;
	case VT_FIXED: 		bl = bu = Lexer::Number;	break;
	case VT_BOUNDED:	bl = Lexer::Number;			break;
#ifndef NDEBUG
	default:			abort();
#endif
	}

	if( rt == RT_RNG )
	{
		if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
		{
			Error( "Numeric value(s) expected." );
			goto error;
		}
		bu = Lexer::Number;
	}

	Lexer::GetNewline( True );
	
	//--------------------------------------------------------------------------
	//	Read the variable type and its simple bounds.
	//
	if( !Lexer::GetSpace() )
	{
		Error( "Premature end of section." );
		return False;
	}

	if( Lexer::GetKeyword( "FX" ) )
		vt = VT_FIXED;
	else if( Lexer::GetKeyword( "UP" ) )
		vt = VT_BOUNDED;
	else if( Lexer::GetKeyword( "PL" ) )
		vt = VT_NORM;
	else if( Lexer::GetKeyword( "MI" ) )
		vt = VT_MI;
	else
	{
		Error( "Unrecognized variable type in presolver action file." );
		goto error;
	}

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "Numeric value(s) expected." );
		goto error;
	}

	switch( vt )
	{
	case VT_FIXED:		l = u = Lexer::Number;	break;
	case VT_BOUNDED:	l = Lexer::Number;		break;
	case VT_NORM:		l = Lexer::Number;		break;
	case VT_MI:			u = Lexer::Number;		break;
#ifndef NDEBUG
	default:			abort();
#endif
	}

	if( vt == VT_BOUNDED )
	{
		if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
		{
			Error( "Numeric value(s) expected." );
			goto error;
		}
		u = Lexer::Number;
	}

	Lexer::GetNewline( True );

	//--------------------------------------------------------------------------
	//	Read the matrix row.
	//
	{
		Int_T len;
		Array<Real_T> a;
		Array<Int_T> col;

		if( !GetRowData( sl, len, a, col ) )
			return False;

		ExplicitSlackRemoval( var, a, col, len, sl, vt, l, u, rt, bl, bu );
	}

	return True;

error:
	SkipToNextSection();
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::GetFreeSingletonAction( Int_T &nn, Int_T var )

PURPOSE:
 	Reads in a section corresponding to free singleton variable removal.

PARAMETERS:
 	Int_T var
 		Variable number (read in from the section's first line).
 
RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::GetFreeSingletonAction( Int_T &nn, Int_T var )
{
	if( !Lexer::GetKeyword( "FREE_SINGL" ) ) return False;
	Lexer::GetNewline( True );

	nn--;

	Real_T b	= 0.0,
		slOpt	= 0.0,
		val		= 0.0;

	//--------------------------------------------------------------------------
	//	Read the optional optimal slack value.
	//
	if( !Lexer::GetSpace() )
	{
		Error( "Premature end of the ACTION FREE_SINGL section." );
		return False;
	}

	if( Lexer::GetKeyword( "SLACK" ) )
	{
		if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
		{
			Error( "SLACK value expected." );
			goto error;
		}

		slOpt = Lexer::Number;
	
		Lexer::GetNewline( True );

		if( !Lexer::GetSpace() )
		{
			Error( "Premature end of the ACTION FREE_SINGL section." );
			return False;
		}
	}

	//--------------------------------------------------------------------------
	//	Read the RHS value.
	//
	if( !Lexer::GetKeyword( "RHS" ) || !Lexer::GetSpace() ||
		!Lexer::GetNumeric() )
	{
		Error( "RHS value expected." );
		goto error;
	}

	b = Lexer::Number;

	Lexer::GetNewline( True );

	//--------------------------------------------------------------------------
	//	Read the matrix row.
	//
	{
		Int_T len;
		Array<Real_T> a;
		Array<Int_T> col;

		if( !GetRowData( val, len, a, col ) )
			return False;

		FreeSingletonColumnRemoval( var, a, col, len, val, b, slOpt );
	}

	return True;

error:
	SkipToNextSection();
	return False;
}


/*------------------------------------------------------------------------------

	void Postsolver::SkipToNextSection( void )

PURPOSE:
	In case an error is detected during reading of a section of the action file,
it may be possible to continue reading the file after omittin the rest of the
current section. Of course, the file contents shall be ignored, but perhaps more
than just one error will be detected in one run.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Postsolver::SkipToNextSection( void )
{
	do
	{
		Lexer::GetNewline( True );
	} while( !Lexer::GetSpace() );
}
