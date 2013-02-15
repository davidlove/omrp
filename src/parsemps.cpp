/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	parsemps.cpp
CREATED:			1993.09.17
LAST MODIFIED:		1996.10.02

DEPENDENCIES:		compile.h, error.h, stdtype.h, mps_lp.h, lp_codes.h,
					parsemps.h, lexer.h
					<stdio.h>, <ctype.h>, <stdlib.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This file contains functions that parse a text file in MPS format
containing LP problem description. IBM's MPSX manual has been used as a source
of information on both fixed and free MPS formats.

	WARNING: While syntax checking is almost 100% accurate, some features of
this parser differ from the MPS standard. This is a list of (hopefully) all of
them:
*	the only differences between fixed and free formats is that in fixed format:
	*	labels are allowed to contain spaces and
	*	BOUNDS vector label is optional (may be skipped) - as it happens in one
		of the NETLIB linear problems: GFRD-PNC.
*	comments starting with "$" (dollar) character are not recognized,
*	empty lines are allowed and skipped (see: lexer.cpp),
*	in fixed format column label in COLUMNS section is always required (even
	when it is the same column as in the previous line of input),
*	the range of mantissa or exponent in representation of numeric values is
	not checked (see: lexer.cpp),
*	in fixed format there may be no non-space characters between the fields
	(note that some MPS-reading procedures ignore such characters); such
	occurence will (in most cases) cause an error,
*	all input after the last field meaningful for the MPS parser is totally
	ignored; no warnings are issued.

*	Additionally if macros "COMP_READLP_1RHS", "COMP_READLP_1RANGE" and
	"COMP_READLP_1BOUND" are defined in "compile.h" only single RHS, range and
	bound vectors are expected, and thus their names are not checked to see
	when next RHS / range / bound vectors begin.

	Since the MPS file format is rather simple, the parser may take advantage
of some simplyfying assumptions. This will make analysis and processing much
faster. This is the list of those assumptions and their consequences:

*	In case of an error we ignore the rest of current line and proceed. This
	almost substitutes for intelligent error analysis and recovery  (see:
	lexer.cpp).
*	As this is not a commercial implementation, we will only use logn used LP
	problems, which are syntactically correct. Therefore we don't check for
	most errors. We also read FIXED and FREE MPS formats in almost exactly the
	same manner, although they are different and normally should be
	distinguished.
*	If "COMP_PARSE_UPPERCASE" is defined, we assume all keywords and indicators
	are uppercase (as they indeed ought to be - according to the MPS standard).
*	Column labels in COLUMN section are always present (see also the list
	above).

	The syntax of the free format MPS file is defined as follows:

MPS_File			::= <NameLine> <MPS_Body>
NameLine			::= "NAME" <Space> <ProblemName> <Space> "FREE" <Newline>
MPS_Body			::=
	"ROWS"		<Newline> <RowsSection>
	"COLUMNS"	<Newline> <ColumnsSection>
	"RHS"		<Newline> <RHS_Section>
	[ "RANGES"	<Newline> <RangesSection> ]
	[ "BOUNDS"	<Newline> <BoundsSection> ]
	"ENDATA"
RowsSection			::=
	{ <Space> <RowTypeID> <Space> <RowName> <Newline> }
ColumnsSection		::=
	{ <Space> <ColumName> <Space>
		<RowName> <Space> <NumericConst>
		[ <Space> <RowName> <Space> <NumericConst> ]
	<Newline> }
RHS_Section			::=
	{ <Space> <RHS_VectorName> <Space>
		<RowName> <Space> <NumericConst>
		[ <Space> <RowName> <Space> <NumericConst> ]
	<Newline> }
RangesSection		::=
	{ <Space> <RangeVectorName> <Space>
		<RowName> <Space> <NumericConst>
		[ <Space> <RowName> <Space> <NumericConst> ]
	<Newline> }
BoundsSection		::=
	{
		( <Space> <BoundWithValueID> <Space> <BoundVectorName> <Space>
			<ColumName> <Space> <NumericConst> <Newline>
		) |
		( <Space> <BoundW/O_ValueID> <Space> <BoundVectorName> <Space>
			<ColumName> <Newline>
		)
	}
RowTypeID			::=	"L" | "G" | "E" | "N"
BoundWithValueID	::=	"LO" | "UP" | "FX"
BoundW/O_ValueID	::=	"PL" | "MI" | "FR"
ProblemName			::= <Label>
RowName				::= <Label>
ColumName			::= <Label>
RHS_VectorName		::= <Label>
RangeVectorName		::= <Label>
BoundVectorName		::= <Label>

<Space>, <Newline>, <Label> and <NumericConst> are lexical symbols  (see:
lexer.cpp).

<Space>			denotes any number of space and tab characters.
<Newline>		denotes a newline character, possibly preceded by some space and
				tab characters.
<Label>			denotes a string of 1 to 8 non-space characters.
<NumericConst>	denotes a numeric constant of the format:
				[+|-] {digit}[.{digit}][e[+|-]digit{digit}]

Notation:
	"x"			terminal symbol "x",
	<x>			non-terminal symbol "x",
	{x}			zero or more of "x",
	[x]			zero or one of "x",
	x|y			"x" or "y",
	()			used for grouping.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void SetLP_TargetObject( MPS_LP *_lp )
	Bool_T GetNameLine( void )
	Bool_T GetMPS_Body( void )

LOCAL TYPE DEFINITIONS:
	typedef void (MPS_LP::*NewDataPtr)( const char *lab0, const char *lab1,
		double val );

STATIC FUNCTIONS:
	static void GetRowsSection( void )
	static Bool_T GetColumnTypeLine( const char *Lab0, const char *Lab1,
		NewDataPtr NewData )
	static void GetColumnsSection( void )
	static void GetRHS_Section( void )
	static void GetRangesSection( void )
	static void GetBoundsSection( void )
	static Bool_T GetRowType( void )
	static Bool_T GetBoundType( void )

STATIC DATA:
	static FF FileType
	static int RowType, BoundType
	static MPS_LP *lp
	static Bool_T (*GetLabel)( int lbl, int StartCol, int EndCol );

--------------------------------------------------------------------------------

USED MACROS FROM COMPILE.H AND THEIR MEANING:
	COMP_PARSE_UPPERCASE	-	if defined - we assume all keywords and
								indicators are uppercase (as they indeed ought
								to be - according to the MPS standard).

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
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
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

//==============================================================================
//
//	Static data
//
//==============================================================================

//------------------------------------------------------------------------------
//	File type (binary or MPS) is stored here.
//
static FF FileType = FF_UNKNOWN;

//------------------------------------------------------------------------------
//	Type of the current constraint matrix row and type of the current variable
//	are stored using symbolic values defined in header "lp_codes.h".
//
static int RowType,
	BoundType;

//------------------------------------------------------------------------------
//	Linear problem data structure.
//	Pointer to MPS_LP member function data type definition. This definition
//	will later bee needed when declaring static functions' prototypes.
//
static MPS_LP *lp = NULL;
typedef void (MPS_LP::*NewDataPtr)( const char *lab0, const char *lab1,
	double val );

//------------------------------------------------------------------------------
//	Static pointers to functions that are different in FIXED and FREE MPS
//	formats (now there's only one such function - GetLabel).
//
static Bool_T (*GetLabel)( int lbl, int StartCol, int EndCol ) = NULL;


//==============================================================================
//
//	End of static data.
//
//==============================================================================

//==============================================================================
//
//	Static function prototypes
//
//==============================================================================

static void GetRowsSection( void );
static Bool_T GetColumnTypeLine( const char *Lab0, const char *Lab1,
	NewDataPtr NewData, Bool_T AllowEmptyLabel = False );
static void GetColumnsSection( void );
static void GetRHS_Section( void );
static void GetRangesSection( void );
static void GetBoundsSection( void );
static Bool_T GetRowType( void );
static Bool_T GetBoundType( void );

//==============================================================================
//
//	End of static function prototypes
//
//==============================================================================

//==============================================================================
//
//	Public functions definitions.
//
//==============================================================================


/*------------------------------------------------------------------------------

	void SetLP_TargetObject( MPS_LP *_lp )

PURPOSE:
	Since the parser stores retrieved information in an MPS_LP object, a
pointer to such object is passed. Parser will (when necessary) call public
member functions of MPS_LP class to pass retrieved information. The functions
it calls are refered to as "semantic actions".

PARAMETERS:
	MPS_LP *_lp
		Pointer to target object.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Storing the pointer for further use.

------------------------------------------------------------------------------*/

void SetLP_TargetObject( MPS_LP *_lp )
{
	assert( _lp != NULL );

	lp = _lp;
}


/*------------------------------------------------------------------------------

	FF GetNameLine( void )

PURPOSE:
	This is the first public function that actually attempts to parse the
input stream. It reads the LP problem name and then tries to recognize file
format. The possible formats are: free MPS, fixed MPS and proprietary binary
(see MPS_LP::ReadBin() function description).
	This function invokes two semantic actions - MPS_LP::SetFileType() and
MPS_LP::SetLP_Name() - in order to inform the MPS_LP object of the file format
and the LP problem name.
	When I/O error occurs or the file format is not recognised, function
reports failure.

PARAMETERS:
	None.

RETURN VALUE:
	File type (or FF_UNKNOWN on error).

SIDE EFFECTS:
	Setting

------------------------------------------------------------------------------*/

FF GetNameLine( void )
{
	assert( lp != NULL );

	static const char *Messages[] = {
		"NAME expected",					// 0
		"Extra characters after NAME",		// 1
		"LP problem name expected",			// 2
		"File I/O error encountered"		// 3
		};
	int MsgCode = -1;

	//--------------------------------------------------------------------------
	if( !Lexer::GetKeyword( "NAME" ) )		{ MsgCode = 0; goto error; }

	if( !Lexer::GetNewline() )
	{
		if( !Lexer::GetSpace() )				{ MsgCode = 1; goto error; }
		if( !Lexer::GetLabelFree( 0, 0, 0 ) )	{ MsgCode = 2; goto error; }

		// We have to store the name before a new line is fetched.
		//
		lp->SetLP_Name( Lexer::LabelPtr[0] );

		Lexer::GetSpace();
		if( Lexer::GetKeyword( "FREE" ) )
			FileType = FF_FREE_MPS;
		else
			FileType = FF_FIXED_MPS;
		if( !Lexer::GetNewline( True ) )	{ MsgCode = 3; goto error; }
	}
	else if( ! Lexer::IO_OK )				{ MsgCode = 3; goto error; }
	else
	{
		lp->SetLP_Name( "" );				// Empty problem name string.

		FileType = FF_FIXED_MPS;
	}

	if( FileType == FF_FREE_MPS )
		GetLabel = Lexer::GetLabelFree;
	else if( FileType == FF_FIXED_MPS )
		GetLabel = Lexer::GetLabelFixed;
	else
		GetLabel = NULL;
		
	//--------------------------------------------------------------------------
	//	Success.
	//
	return FileType;

	//--------------------------------------------------------------------------
	//	Error.
	//
error:
	if( MsgCode >= 0 )
		Error( "%s: %s in line %d", Lexer::FileName, Messages[ MsgCode ],
			Lexer::LineNumber );
	lp->SetLP_Name( "" );
	GetLabel = NULL;
	return FF_UNKNOWN;
}


/*------------------------------------------------------------------------------

	Bool_T GetMPS_Body( void )

PURPOSE:
	This is the second publicly available function that reads the input file.
It assumes the first line was read successfully. This function parses the whole
input file (up to ENDATA indicator line) using the remaining functions from the
module. The data read in will be passed into the MPS_LP object using the
semantic actions (see: "mps_lp.h").
	Syntactic and lexical errors will be reported using "Warning", "Error" and
"FatalError" functions. Therefore message output may be externally (from
outside the module) supressed or redirected (see: "error.h", "error.cpp").

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success / failure status.
	WARNING: You have to check the number of errors after this procedure was
called, because this value will be "True" even though there may have been some
minor errors during file readnig. This value will only be "False" in case of
I/O error, or a major error, like missing section, missing ENDATA etc.

SIDE EFFECTS:
	Reads in the file from the input stream (moves the file pointer) and passes
information it reads to MPS_LP object.

------------------------------------------------------------------------------*/

Bool_T GetMPS_Body( void )
{
	assert( lp != NULL );

	static const char *Messages[] = {
		"ROWS expected",					// 0
		"COLUMNS expected",					// 1
		"RHS expected",						// 2
		"ENDATA expected",					// 3
		"Unexpected EOF or file I/O error"	// 4
		};
	int MsgCode = -1;

	//--------------------------------------------------------------------------
	//	ROWS section.
	//
	if( !Lexer::GetKeyword( "ROWS" ) )		{ MsgCode = 0; goto error; }
	if( !Lexer::GetNewline( True ) )		{ MsgCode = 4; goto error; }
	lp->BeginRows();
	GetRowsSection();

	//--------------------------------------------------------------------------
	//	COLUMNS section.
	//
	if( !Lexer::GetKeyword( "COLUMNS" ) )	{ MsgCode = 1; goto error; }
	if( !Lexer::GetNewline( True ) )		{ MsgCode = 4; goto error; }
	lp->BeginColumns();
	GetColumnsSection();

	//--------------------------------------------------------------------------
	//	RHS section.
	//
	if( !Lexer::GetKeyword( "RHS" ) )		{ MsgCode = 2; goto error; }
	if( !Lexer::GetNewline( True ) )		{ MsgCode = 4; goto error; }
	lp->BeginRHS();
	GetRHS_Section();

	//--------------------------------------------------------------------------
	//	RANGES section.
	//
	if( Lexer::GetKeyword( "RANGES" ) )
	{
		if( !Lexer::GetNewline( True ) )	{ MsgCode = 4; goto error; }
		lp->BeginRanges();
		GetRangesSection();
	}

	//--------------------------------------------------------------------------
	//	BOUNDS section.
	//
	if( Lexer::GetKeyword( "BOUNDS" ) )
	{
		if( !Lexer::GetNewline( True ) )	{ MsgCode = 4; goto error; }
		lp->BeginBounds();
		GetBoundsSection();
	}

	//--------------------------------------------------------------------------
	//	ENDATA.
	//
	if( !Lexer::GetKeyword( "ENDATA" ) )	{ MsgCode = 3; goto error; }
	lp->Endata();

	//--------------------------------------------------------------------------
//success:
	return True;

error:
	if( MsgCode >= 0 )
		Error( "%s: %s in line %d", Lexer::FileName, Messages[ MsgCode ],
			Lexer::LineNumber );
	return False;
}

//==============================================================================
//
//	End of public functions' definitions.
//
//==============================================================================


//==============================================================================
//
//	Static functions definitions.
//
//==============================================================================

/*------------------------------------------------------------------------------

	static void GetRowsSection( void )

PURPOSE:
	This function reads ROWS section of an MPS file.

PARAMETERS:
	None.

RETURN VALUE:
	None (this procedure always succeeds).

SIDE EFFECTS:
	Moves forward the file pointer.

------------------------------------------------------------------------------*/

static void GetRowsSection( void )
{
	static const char *Messages[] = {
		"Row type indicator expected",					// 0
		"Extra characters after row type indicator",	// 1
		"Row label expected"							// 2
		};
	static int MsgCode = -1;

	//--------------------------------------------------------------------------
	while( Lexer::GetSpace() )
	{
		if( !GetRowType() )						{ MsgCode = 0; goto error; }
		if( !Lexer::GetSpace() )				{ MsgCode = 1; goto error; }
		if( !GetLabel( 0, 4, 11 ) )				{ MsgCode = 2; goto error; }

		lp->NewRow( (Short_T)RowType, Lexer::LabelPtr[0] );

	error:
		if( MsgCode >= 0 )
		{
			Error( "%s: %s in line %d", Lexer::FileName, Messages[ MsgCode ],
				Lexer::LineNumber );
			MsgCode = -1;
		}
		if( !Lexer::GetNewline( True ) )
		//	Any possible 'SOSROW" markers after row name are ignored.
			return;
	}
	//--------------------------------------------------------------------------
}


/*------------------------------------------------------------------------------

	static Bool_T GetColumnTypeLine( const char *msg0, const char *msg1,
		NewDataPtr NewData )

PURPOSE:
	This function reads one line of COLUMN, RHS or RANGES section. According to
MPS standard specification lines in one of those section have identical format
and intepretation. Therefore a single function is used for reading these
sections' lines.

PARAMETERS:
	const char *msg0
	const char *msg1
		These two strings are used when constructing error messages (e.g. when
		a row label is missing in COLUMNS section, msg0 will be "COLUMN" and
		msg1 - "Row", and together with "label missing" they will create an
		error message).

	NewDataPtr NewData
		This is a pointer to MPS_LP function member, which is a semantic action
		corresponding to the MPS file section being processed:
		*	MPS_LP::NewNonZero	in COLUMNS section,
		*	MPS_LP::NewRHS		in RHS section,
		*	MPS_LP::NewRANGE	in RANGES section.

RETURN VALUE:
	Boolean success status. Failure is reported ONLY on I/O error. Otherwise
success is returned even if errors were encountered and reported.

SIDE EFFECTS:
	Moves forward the file pointer.

------------------------------------------------------------------------------*/

static Bool_T GetColumnTypeLine( const char *msg0, const char *msg1, // )
	NewDataPtr NewData, Bool_T AllowEmptyLabel )
{
	static const char *Messages[] = {
		"Extra characters after",				// 0
		"Numeric expected",						// 1
		"Extra characters after numeric value"	// 2
		};
	int MsgCode = -1;

	//--------------------------------------------------------------------------
	if( !AllowEmptyLabel )
	{
		if( !GetLabel( 0, 4, 11 ) )				{ MsgCode = 0; goto error; }
	}
	else
		GetLabel( 0, 4, 11 );
	Lexer::GetSpace();
	if( !GetLabel( 1, 14, 21 ) )				{ MsgCode = 2; goto error; }

	//
	//	Ignore marker lines (those that containt a word 'MARKER' in the second
	//	field.
	//
	if( strcmp( Lexer::LabelPtr[1], "'MARKER'" ) == 0 )
	{
		if( !Lexer::GetNewline( True ) )		goto IO_error;
	}
	else
	{
		Lexer::GetSpace();
		if( !Lexer::GetNumeric() )				{ MsgCode = 4; goto error; }

		(lp->*NewData)( Lexer::LabelPtr[0], Lexer::LabelPtr[1], Lexer::Number );

		if( !Lexer::GetNewline() )
		{
			if( !Lexer::IO_OK )					goto IO_error;
			if( !Lexer::GetSpace() )			{ MsgCode = 4; goto error; }
			if( !GetLabel( 1, 39, 46 ) )		{ MsgCode = 2; goto error; }
			Lexer::GetSpace();
			if( !Lexer::GetNumeric() )			{ MsgCode = 3; goto error; }

			(lp->*NewData)( NULL, Lexer::LabelPtr[1], Lexer::Number );

			if( !Lexer::GetNewline( True ) )	goto IO_error;
		}
	}
	//--------------------------------------------------------------------------

//success:
	return True;

error:
	switch( MsgCode )
	{
	case 0:
		Error( "%s: %s label expected in line %d", Lexer::FileName, msg0,
			Lexer::LineNumber );
		break;
	case 1:
		Error( "%s: %s %s label in line %d", Lexer::FileName, Messages[0], msg0,
			Lexer::LineNumber );
		break;
	case 2:
		Error( "%s: %s label expected in line %d", Lexer::FileName, msg1,
			Lexer::LineNumber );
		break;
	case 3:
		Error( "%s: %s in line %d", Lexer::FileName, Messages[1],
			Lexer::LineNumber );
		break;
	case 4:
		Error( "%s: %s in line %d", Lexer::FileName, Messages[2],
			Lexer::LineNumber );
		break;
	}
	if( !Lexer::GetNewline( True ) )	return False;
	return True;

IO_error:
	return False;
}


/*------------------------------------------------------------------------------

	static void GetColumnsSection( void )
	static void GetRangesSection( void )
	static void GetRHS_Section( void )

PURPOSE:
	These functions read COLUMNS, RANGES and RHS sections respectively. Both
their implementations and the syntax of MPS file sections they read are very
similar. They all use "GetColumnTypeLine()" function to read a single line of
their sections of input file.

PARAMETERS:
	None.

RETURN VALUE:
	None (they are assumed to always succeed).

SIDE EFFECTS:
	Moves forward the file pointer.

------------------------------------------------------------------------------*/

static void GetColumnsSection( void )
{
	while( Lexer::GetSpace() )
		GetColumnTypeLine( "COLUMN", "Row", &MPS_LP::NewNonZero );
}


static void GetRangesSection( void )
{
	while( Lexer::GetSpace() )
		GetColumnTypeLine( "RANGE vector", "Row", &MPS_LP::NewRange,
			( FileType == FF_FIXED_MPS ) ? True : False );
}


static void GetRHS_Section( void )
{
	while( Lexer::GetSpace() )
		GetColumnTypeLine( "RHS vector", "Row", &MPS_LP::NewRHS,
			( FileType == FF_FIXED_MPS ) ? True : False );
}


/*------------------------------------------------------------------------------

	static void GetBoundsSection( void )

PURPOSE:
	This function reads BOUNDS section of MPS file. It recognizes type of bound
introduced by the input stream line being parsed, and depending on its type may
read a numerical value. Calls semantic action "MPS_LP::NewBound()" to pass read
information.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Moves forward the file pointer.

------------------------------------------------------------------------------*/

static void GetBoundsSection( void )
{
	static const char *Messages[] =
		{
		"Bound type indicator expected",				// 0
		"Extra characters after bound type indicator",	// 1
		"Column label expected",						// 2
		"Numeric value expected",						// 3
		"Variable label expected",						// 4
		};
	static int MsgCode = -1;

	const Bool_T LabelRequired = ( FileType == FF_FREE_MPS ) ? True : False;

	while( Lexer::GetSpace() )
	{
		if( !GetBoundType() )					{ MsgCode = 0; goto error; }
		if( FileType == FF_FREE_MPS &&
			!Lexer::GetSpace() )				{ MsgCode = 1; goto error; }
		if( LabelRequired )
		{
			if( !GetLabel( 0, 4, 11 ) )			{ MsgCode = 2; goto error; }
		}
		else
			GetLabel( 0, 4, 11 );
		Lexer::GetSpace();
		if( !GetLabel( 1, 14, 21 ) )			{ MsgCode = 4; goto error; }
		if( BoundType == VTM_MI || BoundType == VTM_PL || BoundType == VTM_FR )
		{
			lp->NewBound( (Short_T)BoundType, Lexer::LabelPtr[0],
				Lexer::LabelPtr[1], 0.0e0 );
		}
		else if( BoundType == VTM_BV )
		{
			BoundType = VTM_UP;
			lp->NewBound( (Short_T)BoundType, Lexer::LabelPtr[0],
				Lexer::LabelPtr[1], 1.0e0 );
		}
		else
		{
			Lexer::GetSpace();
			if( !Lexer::GetNumeric() )			{ MsgCode = 3; goto error; }

			lp->NewBound( (Short_T)BoundType, Lexer::LabelPtr[0],
				Lexer::LabelPtr[1], Lexer::Number );
		}

	error:
		if( MsgCode >= 0 )
		{
			Error( "%s: %s in line %d", Lexer::FileName, Messages[ MsgCode ],
				Lexer::LineNumber );
			MsgCode = -1;
		}

		if( !Lexer::GetNewline( True ) )	return;
	}
}


/*------------------------------------------------------------------------------

	static Bool_T GetBoundType( void )

PURPOSE:
	Reads one of six possible simple variable bounds. If "COMP_PARSE_UPPERCASE"
is defined in "compile.h" symbols read in from the stream are assumed to be
uppercase.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and static buffer pointer behind the bound type
indicator.

------------------------------------------------------------------------------*/

static Bool_T GetBoundType( void )
{
	register char *BufPtr = Lexer::BufPtr;

#ifndef COMP_PARSE_UPPERCASE
	switch( *BufPtr++ )
	{
	//	UP and UI (for integer variables)
	case 'u':
	case 'U':
		if( *BufPtr == 'P' || *BufPtr == 'p' ||
			*BufPtr == 'I' || *BufPtr == 'i' )
			{ BoundType = VTM_UP; break; }
		else goto error;

	//	LO and LI (for integer variables)
	case 'l':
	case 'L':
		if( *BufPtr == 'O' || *BufPtr == 'o' ||
			*BufPtr == 'I' || *BufPtr == 'i' )
			{ BoundType = VTM_LO; break; }
		else goto error;

	//	FX and FR
	case 'f':
	case 'F':
		if( *BufPtr == 'X' || *BufPtr == 'x' )
			{ BoundType = VTM_FX; break; }
		else if( *BufPtr == 'R' || *BufPtr == 'r' )
			{ BoundType = VTM_FR; break; }
		else goto error;

	//	MI
	case 'm':
	case 'M':
		if( *BufPtr == 'I' || *BufPtr == 'i' )
			{ BoundType = VTM_MI; break; }
		else goto error;

	//	PL
	case 'p':
	case 'P':
		if( *BufPtr == 'L' || *BufPtr == 'l' )
			{ BoundType = VTM_PL; break; }
		else goto error;

	//	BV
	case 'b':
	case 'B':
		if( *BufPtr == 'V' || *BufPtr == 'v' )
			{ BoundType = VTM_BV; break; }
		else goto error;

	default:	
		BoundType = VTM_UNDEFINED;
		goto error;
	}
#else
	switch( *BufPtr++ )
	{
	//	UP and UI (for integer variables)
	case 'U':
		if( *BufPtr == 'P' || *BufPtr == 'I' )
			{ BoundType = VTM_UP; break; }
		else goto error;

	//	LO and LI (for integer variables)
	case 'L':
		if( *BufPtr == 'O' || *BufPtr == 'I' )
			{ BoundType = VTM_LO; break; }
		else goto error;

	//	FX and FR
	case 'F':
		if( *BufPtr == 'X' ) { BoundType = VTM_FX; break; }
		else if( *BufPtr == 'R' ) { BoundType = VTM_FR; break; }
		else goto error;

	//	MI
	case 'M':
		if( *BufPtr == 'I' ) { BoundType = VTM_MI; break; } else goto error;

	//	PL
	case 'P':
		if( *BufPtr == 'L' ) { BoundType = VTM_PL; break; } else goto error;

	//	BV
	case 'b':
	case 'B':
		if( *BufPtr == 'V' ) { BoundType = VTM_BV; break; } else goto error;

	default:	
		BoundType = VTM_UNDEFINED;
		goto error;
	}
#endif

	Lexer::BufPtr = ++BufPtr;
	return True;

error:
	Error( "%s: Unknown bound type in line %d", Lexer::FileName,
		Lexer::LineNumber );
	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetRowType( void )

PURPOSE:
	Reads one of four possible row types.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and static buffer pointer behind the row type
indicator.

------------------------------------------------------------------------------*/

static Bool_T GetRowType( void )
{
	char c = *Lexer::BufPtr;

#ifdef COMP_PARSE_UPPERCASE
	if( islower( c ) ) c += 'A'-'a';
#endif

	switch( c )
	{
	case 'e': case 'E':	RowType = RT_EQ;			break;
	case 'l': case 'L':	RowType = RT_LE;			break;
	case 'g': case 'G':	RowType = RT_GE;			break;
	case 'n': case 'N':	RowType = RT_FR;			break;
	default:			RowType = RT_UNDEFINED;		goto error;
	}

	Lexer::BufPtr++;
	return True;

error:
	Error( "%s: Unknown row type indicator in line %d", Lexer::FileName,
		Lexer::LineNumber );
	return False;
}
