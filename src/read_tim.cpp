/*------------------------------------------------------------------------------
MODULE TYPE:		File input routine.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	read_tim.cpp
CREATED:			1994.08.16
LAST MODIFIED:		1995.08.12

DEPENDENCIES:		read_tim.h, mps_lp.h, stdtype.h, smartptr.h, print.h,
					<stdio.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Time file reader function (parses and performs the semantic actions).

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Bool_T ReadTimes( FILE *fp, MPS_LP &LP, Int_T &Stage1Row,
		Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col,
		VerbLevel Verbosity )

------------------------------------------------------------------------------*/


#ifndef __READ_TIM_H__
#	include "read_tim.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T ReadTimes( FILE *fp, MPS_LP &LP, Int_T &Stage1Row,
		Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col,
		VerbLevel Verbosity )

PURPOSE:
	Reads the time file from a stream. 

PARAMETERS:
	FILE *fp
		Input stream.

	MPS_LP &LP
		Linear problem read previously from the core file. It contains (among
		all other data) the labels, that are also used in the time file.

	Int_T &Stage1Row, Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col
		These are output arguments. The numbers of the rows and columns starting
		the first and second stage data.

	VerbLevel Verbosity
		Level of verbosity. Directs the amount of printouts during the reading.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	File pointer in the output stream is moved past the "ENDATA" marker or past
	the point in which the first error is detected.

------------------------------------------------------------------------------*/

Bool_T ReadTimes( FILE *fp, MPS_LP &LP, Int_T &Stage1Row, // )
	Int_T &Stage1Col, Int_T &Stage2Row, Int_T &Stage2Col, VerbLevel Verbosity )
{
	if( Verbosity >= V_HIGH )
		Print( "\nReading the time file.\n" );

	Lexer::SetInputStream( "Time file", fp );

	const unsigned NAME_LEN = 8;
	typedef char Label[NAME_LEN+1];
	Label ProblemName,
		StageRowLabel[2],
		StageColLabel[2],
		StagePeriodLabel[2];
	SortedArrayOfLabels &Rows	= LP.RevealRowLabels(),
		&Cols					= LP.RevealColumnLabels();
	int i;

	//--------------------------------------------------------------------------
	//	Read in the first line.
	//
	if( !Lexer::GetKeyword( "TIME" ) || !Lexer::GetSpace() ||
		!Lexer::GetLabelFree( 0 ) )
	{
		Error( "Error reading the first line." );
		goto error;
	}
	strncpy( ProblemName, Lexer::LabelPtr[0], NAME_LEN );
	ProblemName[NAME_LEN] = '\0';
	if( Verbosity >= V_HIGH )
		Print( "\t%-20s%-20s\n", "Problem name:", ProblemName );

	if( !Lexer::GetNewline() )
	{
		Error( "Extra input in the first line." );
		goto error;
	}

	//--------------------------------------------------------------------------
	//	Read in the 'PERIODS' keyword and two lines defining actual periods.
	//
	if( !Lexer::GetKeyword( "PERIODS" ) )
	{
		Error( "'PERIODS' expected." );
		goto error;
	}

	Lexer::GetSpace();
	Lexer::GetKeyword( "LP" );

	if( !Lexer::GetNewline( True ) )
		Error( "Unexpected end of text?" );

	for( i = 0; i < 2; i++ )
	{
		if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) ||
			!Lexer::GetSpace() || !Lexer::GetLabelFree( 1 ) )
		{
			Error( "Column and row labels expected." );
			goto error;
		}
		strncpy( StageColLabel[i], Lexer::LabelPtr[0], NAME_LEN );
		strncpy( StageRowLabel[i], Lexer::LabelPtr[1], NAME_LEN );
		StageRowLabel[i][NAME_LEN] = StageColLabel[i][NAME_LEN] = '\0';

		if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) )
		{
			Error( "Period label expected." );
			goto error;
		}
		strncpy( StagePeriodLabel[i], Lexer::LabelPtr[0], NAME_LEN );
		StagePeriodLabel[i][NAME_LEN] = '\0';

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
	if( Verbosity >= V_HIGH )
		Print( "Finished reading.\n" );

	//--------------------------------------------------------------------------
	//	Now is the time to convert labels into numbers and check if they make
	//	sense.
	//
	Stage1Row = Stage1Col = Stage2Col = Stage2Row = -1;

	Stage1Col = (Int_T) Cols.FindLabel( StageColLabel[0] );
	Stage2Col = (Int_T) Cols.FindLabel( StageColLabel[1] );
	Stage1Row = (Int_T) Rows.FindLabel( StageRowLabel[0] );
	Stage2Row = (Int_T) Rows.FindLabel( StageRowLabel[1] );

	if( Verbosity >= V_HIGH )
		Print(
			"\t%10s  %10s  %10s\n"
			"\t%10s  %10s  %10s\n"
			"\t%10s  %10d  %10d\n"
			"\t%10s  %10s  %10s\n"
			"\t%10s  %10d  %10d\n",

			"PERIOD",				"ROW",				"COL",
			StagePeriodLabel[0],	StageRowLabel[0],	StageColLabel[0], 
			"",						(int)Stage1Row,		(int)Stage1Col,
			StagePeriodLabel[1],	StageRowLabel[1],	StageColLabel[1], 
			"",						(int)Stage2Row,		(int)Stage2Col
		);

	assert( Stage1Col >= 0 && Stage1Col < LP.GetN() );
	assert( Stage2Col >= 0 && Stage2Col < LP.GetN() );
	assert( Stage1Row >= 0 && Stage1Row < LP.GetM() );
	assert( Stage2Row >= 0 && Stage2Row < LP.GetM() );

	return True;

error:
	Stage2Row = Stage2Col = -1;
	return False;
}
