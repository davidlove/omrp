/*------------------------------------------------------------------------------
MODULE TYPE:		File input routine.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	periods.cpp
CREATED:			1995.12.27
LAST MODIFIED:		1996.01.30

DEPENDENCIES:		read_tim.h, mps_lp.h, stdtype.h, smartptr.h, print.h,
					<stdio.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Time file reader function (parses and performs the semantic actions).

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

------------------------------------------------------------------------------*/


#ifndef __PERIODS_H__
#	include "periods.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif



/*------------------------------------------------------------------------------

	Periods::Periods( void )

PURPOSE:
	Class "Periods" constructor. Initializes an object.

PARAMETERS:
	None.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Periods::Periods( void )
	: maxLen( 10 ), len( 0 ), PeriodLab( 10 ),
	Row( maxLen, -1 ), Col ( maxLen, -1 ),
	TotalRows( -1 ), TotalCols( -1 )
{
	*ProblemName = '\0';
}


/*------------------------------------------------------------------------------

	Bool_T Periods::ReadPeriodFile( FILE *fp, MPS_LP &LP, VerbLevel Verbosity )

PURPOSE:
	Reads the time file from a stream. 

PARAMETERS:
	FILE *fp
		Input stream.

	MPS_LP &LP
		Linear problem read previously from the core file. It contains (among
		all other data) the labels, that are also used in the time file.

	VerbLevel Verbosity
		Level of verbosity. Directs the amount of printouts during the reading.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	File pointer in the output stream is moved past the "ENDATA" marker or past
	the point in which the first error is detected.

------------------------------------------------------------------------------*/

Bool_T Periods::ReadPeriodFile( FILE *fp, MPS_LP &LP, VerbLevel Verbosity )
{
	Bool_T err = False;

	if( Verbosity >= V_HIGH )
		Print( "\nReading the time file.\n" );

	Lexer::SetInputStream( "Time file", fp );

	//--------------------------------------------------------------------------
	//	Read in the first line.
	//
	if( !Lexer::GetKeyword( "TIME" ) || !Lexer::GetSpace() ||
		!Lexer::GetLabelFree( 0 ) )
	{
		Error( "Error reading the first line of the time file." );
		goto error;
	}

	strncpy( ProblemName, Lexer::LabelPtr[0], LAB_LEN );
	ProblemName[LAB_LEN] = '\0';

	if( Verbosity >= V_HIGH )
		Print( "\t%-20s%-20s\n", "Problem name:", ProblemName );

	if( !Lexer::GetNewline() )
	{
		Error( "Extra input in the first line." );
		goto error;
	}

	//--------------------------------------------------------------------------
	//	Read in the 'PERIODS' keyword and the lines defining actual periods.
	//
	if( !Lexer::GetKeyword( "PERIODS" ) )
	{
		Error( "'PERIODS' expected." );
		goto error;
	}

	Lexer::GetSpace();
	Lexer::GetKeyword( "LP" );

	if( !Lexer::GetNewline( True ) )
		Error( "Unexpected end of text." );

	while( Lexer::GetSpace() )
	{
		if( !Lexer::GetLabelFree( 0 ) || !Lexer::GetSpace() ||
			!Lexer::GetLabelFree( 1 ) )
		{
			Error( "Column and row labels expected." );
			goto error;
		}

		const char *col	= Lexer::LabelPtr[0],
			*row		= Lexer::LabelPtr[1];

		if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) )
		{
			Error( "Period label expected." );
			goto error;
		}

		if( !AddPeriod( LP, Lexer::LabelPtr[0], row, col ) )
		{
			Error( "Labels for %d'th period \'%s' not found.", (int)len,
				Lexer::GetLabelFree( 0 ) );
			err = True;
		}
			
		if( !Lexer::GetNewline( True ) )
			goto error;
	}

	//--------------------------------------------------------------------------
	//	Read in the keyword 'ENDATA' and end reading file.
	//
	if( !Lexer::GetKeyword( "ENDATA" ) )
	{
		Error( "'ENDATA' expected. (Only two period problems understood.)" );
		goto error;
	}

	PeriodLab.SortLabels();

	if( Verbosity >= V_HIGH )
		Print( "Finished reading.\n" );

	//--------------------------------------------------------------------------
	//	Check the structure: see if the periods are specified orderly.
	//
	if( err ) goto error;

	if( Check( LP ) )
	{
		if( Verbosity >= V_HIGH )
			PrintPeriods( LP );
	}
	else
	{
		Error( "Inconsistent data read from the time file." );
		goto error;
	}

	//--------------------------------------------------------------------------
	//	Finally, store the dimensions of the LP.
	//
	TotalRows = LP.GetM();
	TotalCols = LP.GetN();

	assert( Row[len-1] < TotalRows );
	assert( Col[len-1] < TotalCols );

	return True;

error:
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Periods::AddPeriod( MPS_LP &LP, const char *label,
		const char *rowLab, const char *colLab )

PURPOSE:
	Class "Periods" constructor. Initializes an object.

PARAMETERS:
	MPS_LP &LP
		The deterministic linear problem (from which the labels are read).

	const char *label
		Period label.

	const char *rowLab, const char *colLab
		The labels of the first row and column corresponding to the named period.

RETURN VALUE:
	Boolean success status ("False" only if the appropriate labels are not found
in the linear problem).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Periods::AddPeriod( MPS_LP &LP, const char *label, // )
	const char *rowLab, const char *colLab )
{
	assert( label != NULL );
	assert( rowLab != NULL );
	assert( colLab != NULL );

	//--------------------------------------------------------------------------
	//	Adjust the size of the "Row" and "Col" arrays (if necessary).
	//
	if( len >= maxLen )
	{
		Int_T NewMaxLen = Int_T(maxLen + 10);

		Row.Resize( NewMaxLen );
		Row.Fill( -1, NewMaxLen, maxLen );

		Col.Resize( NewMaxLen );
		Col.Fill( -1, NewMaxLen, maxLen );

		maxLen = NewMaxLen;
	}

	//--------------------------------------------------------------------------
	//	Convert the labels into indice. Stor the period label.
	//
	Row[len] = (Int_T) LP.RevealRowLabels().FindLabel( rowLab );
	Col[len] = (Int_T) LP.RevealColumnLabels().FindLabel( colLab );

	PeriodLab.AddLabel( label );

	len++;

	assert( PeriodLab.NumberOfLabels() == len );

	//--------------------------------------------------------------------------
	//	Return success status.
	//
	return ( ( Row[len-1] >= 0 ) && ( Col[len-1] >= 0 ) ) ? True: False;
}


/*------------------------------------------------------------------------------

	Bool_T Periods::Check( MPS_LP &LP )

PURPOSE:
	Checks the object data consistency. Subsequent periods should be specified
in such order that the row and column indice constitute a non-decreasing
sequence bounded from below by zero and from above by appropriate dimension of
the matrix.

PARAMETERS:
	MPS_LP &LP
		The linear problem which the period file refers to. The problem
		dimensions are taken from it.

RETURN VALUE:
	"True" when the object data appears to be consistent, "False" otherwise.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Periods::Check( MPS_LP &LP )
{
	Int_T RowLo	= 0,
		RowUp	= LP.GetM(),
		ColLo	= 0,
		ColUp	= LP.GetN();

	for( Int_T i = 0; i < len; i++ )
	{
		if( Row[i] < RowLo || Row[i] >= RowUp ||
			Col[i] < ColLo || Col[i] >= ColUp )
			return False;
		RowLo = Row[i];
		ColLo = Col[i];
	}

	return True;
}


/*------------------------------------------------------------------------------

	void Periods::PrintPeriods( VerbLevel Verbosity, MPS_LP &LP )

PURPOSE:
	Prints out the periods with their labels.

PARAMETERS:
	MPS_LP &LP
		We take the row and column labels from this LP.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Periods::PrintPeriods( MPS_LP &LP )
{
	SortedArrayOfLabels &RowLab	= LP.RevealRowLabels(),
		&ColLab					= LP.RevealColumnLabels();

	Print(
		"\t%-3s  %-10s  %-10s (%-6s)  %-10s (%-6s)\n"
		"\t%3s  %-10s  %-10s %8s  %10s %8s\n",
		"LP", "PERIOD", "ROW", "ind", "COL", "ind",
		"---", "----------", "----------", "--------", "----------", "--------"
	);

	for( Int_T i = 0; i < len; i++ )
	{
		Print(
			"\t%3d. %-10s  %-10s (%6d)  %-10s (%6d)\n",
			i,								PeriodLab.FindLabel( i ),
			RowLab.FindLabel( Row[i] ),		Row[i],
			ColLab.FindLabel( Col[i] ),		Col[i]
		);
	}

	Print( "\n" );
}


/*------------------------------------------------------------------------------

	Int_T Periods::ColumnInPeriod( Int_T col ) const
	Int_T Periods::RowInPeriod( Int_T row ) const

PURPOSE:
	Calculates the period to which a given row/column belongs.

PARAMETERS:
	Int_T row
	Int_T col
		Row/column number.

RETURN VALUE:
	Period number.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Periods::ColumnInPeriod( Int_T col )
const
{
	assert( col >= 0 && col < TotalCols );

	Int_T p = -1;

	for( Int_T i = 0; i < len; i++ )
	{
		if( col >= Col[i] )
			p = i;
		else if( col < Col[i] )
			break;
	}

	return p;
}


Int_T Periods::RowInPeriod( Int_T row )
const
{
	assert( row >= 0 && row < TotalRows );

	Int_T p = -1;

	for( Int_T i = 0; i < len; i++ )
	{
		if( row >= Row[i] )
			p = i;
		else if( row < Row[i] )
			break;
	}

	return p;
}


/*------------------------------------------------------------------------------

	void Periods::RowRange( Int_T period, Int_T &start, Int_T &end ) const
	void Periods::ColumnRange( Int_T period, Int_T &start, Int_T &end ) const

PURPOSE:
	Calculates the row/column range of a given period.

PARAMETERS:
	Int_T period
		Period number.

	Int_T &start, Int_T &end
		Row/column range

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Periods::RowRange( Int_T period, Int_T &start, Int_T &end )
const
{
	assert( period < len );

	start	= Row[period];
	end		= ( period == len - 1 ) ? TotalRows : Row[period+1];
}


void Periods::ColumnRange( Int_T period, Int_T &start, Int_T &end )
const
{
	assert( period < len );

	start	= Col[period];
	end		= ( period == len - 1 ) ? TotalCols : Col[period+1];
}
