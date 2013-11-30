/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	lexer.cpp
CREATED:			1993.10.10
LAST MODIFIED:		1995.10.23

DEPENDENCIES:		compile.h, error.h, stdtype.h, lexer.h, mps_lp.h
					<stdio.h>, <ctype.h>, <stdlib.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This file contains functions that perform lexical analysis of a text file.
Available lexical symbols are <Space>, <Newline>, <Label>, <Keyword> and
<NumericConst>:
<Space>			denotes any number of space and tab characters.
<Newline>		denotes a newline character, possibly preceded by some space and
				tab characters.
<Label>			denotes a string of 1 to 8 non-space characters (in FREE MPS
				format) or 8 any characters (in FIXED MPS format).
<Keyword>		denotes a single word with no embedded spaces,
<NumericConst>	denotes a floating point numeric constant of the format:
				[+|-] {digit}[.{digit}][e[+|-]digit{digit}]

Notation:
	"x"			terminal symbol "x",
	<x>			non-terminal symbol "x",
	{x}			zero or more of "x",
	[x]			zero or one of "x",
	x|y			"x" or "y",
	()			used for grouping.

	During lexical analysis some more actions are performed, including
stripping off comment lines. Comment lines start with a "*" (dollar) character
(as in MPS standard for IBM's MPSX linear algebra package). They are ignored.
Empty lines are allowed and skipped.
	The range of mantissa or exponent in representation of numeric values is
not checked. Input after the last maningful field of the input file line is
ignored; no warnings are issued.

	Since input file format supported is rather simple, the lexical analyser
(lexer) may take advantage of some simplyfying assumptions. This will make
analysis and processing much faster. This is the list of those assumptions and
their consequences:

*	Input is divided into lines of text. It is read one text line at a time.
	We only need an input buffer capable to hold one line of text. We never have
	to re-read a line of input. In worst case we may sometimes have to go back
	a couple of characters in the current line - a single buffer pointer will
	be capable of holding all information for this limited backtracking.
*	In case of an error we ignore the rest of current line and proceed. This
	almost substitutes for intelligent error analysis and recovery.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void Lexer::SetInputStream( const char *name, FILE *fp )
	Bool_T Lexer::GetSpace( void )
	Bool_T Lexer::GetKeyword( const char *Keyword,
		Bool_T CaseSensitive = False );
	Bool_T Lexer::GetLabelFree( int lbl, int, int )
	Bool_T Lexer::GetLabelFixed( int lbl, int StartCol, int EndCol )
	Bool_T Lexer::GetNumeric( void )
	Bool_T Lexer::GetNewline( Bool_T Force = False )

------------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif

#ifndef __COMPILE_H__
#	include "compile.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif

//==============================================================================
//
//	Static "Lexer" class data (both private and public).
//
//==============================================================================

//------------------------------------------------------------------------------
//	Input stream descriptor pointer (as passed through "SetInputStream" call).
//
FILE *Lexer::Input;

//------------------------------------------------------------------------------
//	These four data store information concerning the input stream: stream
//	descriptor pointer, filename (as passed through "SetInputStream" call)
//	which may be used when issuing error messages, current line number (also
//	for use by error message output) and finally file type (binary or MPS).
//
char Lexer::FileName[ Lexer::IO_BUF_LEN ];
int Lexer::LineNumber;

//------------------------------------------------------------------------------
//	This value will be NULL whenever IO error will occur. Function "GetNewline"
//	which is responsible for reading in new input lines uses ANSI C function
//	"fgets", which returns a NULL pointer on error. The value returned by
//	"fgets" is stored in this variable.
//
void *Lexer::IO_OK;

//------------------------------------------------------------------------------
//	This is file input buffer together with a pointer to current buffer
//	position (the next character to read). The buffer holds one line of input.
//
char Lexer::Buf[ Lexer::IO_BUF_LEN ];
char *Lexer::BufPtr;

//------------------------------------------------------------------------------
//	Values of lexical symbols: numeric constant (stored as double), type of
//	constraint matrix row, type of variable and pointers to labels. The
//	constraint and variable types are stored using symbolic values defined in
//	header "lp_codes.h".
//
double Lexer::Number;
char *Lexer::LabelPtr[2];
const char Lexer::EmptyString[1] = { '\0' };

//==============================================================================
//
//	End of class "Lexer" static data.
//
//==============================================================================


//==============================================================================
//
//	Functions' definitions.
//
//==============================================================================

/*------------------------------------------------------------------------------

	void Lexer::SetInputStream( const char *name, FILE *fp )

PURPOSE:
	By calling this function you tell the lexical analyser what is the file
name of the input file and pass a descriptor of an open input stream from which
the data will be read.

PARAMETERS:
	const char *name
		Input file name (any valid string). This file name is only used for
		error message output.

	FILE *fp
		Stream descriptor for data input. The stream has to be opened before
		calling this function..

RETURN VALUE:
	None.

SIDE EFFECTS:
	Clears all static data and then stores in it both input file name and input
stream descriptor. Also reads the first line and places it in the input buffer.

------------------------------------------------------------------------------*/

void Lexer::SetInputStream( const char *name, FILE *fp )
{
	assert( (void *)name != NULL );
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	Erase all static data.
	//
	Input			= NULL;
	*FileName		= '\0';
	LineNumber		= 0;
	IO_OK			= (void *)1;
	*Buf			= '\0';
	BufPtr			= Buf;

	Number			= 0.0e0;
	LabelPtr[ 0 ]	=
	LabelPtr[ 1 ]	= NULL;

	//--------------------------------------------------------------------------
	//	Copy the name of the input file.
	//
	char *c = FileName;
	while( ( *c++ = *name++ ) != 0 )
		;

	//--------------------------------------------------------------------------
	//	Check and initialize the stream (read the first line).
	//
	Input = fp;

	GetNewline( True );
}


/*------------------------------------------------------------------------------

	Bool_T Lexer::GetSpace( void )

PURPOSE:
	Read one or more space and tab characters. At least one space character is
required.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and buffer pointer behind the last space
character.

------------------------------------------------------------------------------*/

Bool_T Lexer::GetSpace( void )
{
	char *_BufPtr = BufPtr;

	if( *_BufPtr != ' ' && *_BufPtr != '\t' )
		goto error;

	_BufPtr++;
	while( *_BufPtr == ' ' || *_BufPtr == '\t' )
		_BufPtr++;

//success:
	BufPtr = _BufPtr;
	return True;

error:
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Lexer::GetKeyword( const char *Keyword, Bool_T CaseSensitive )

PURPOSE:
	Reads a keyword passed as a "Keyword" parameter.

PARAMETERS:
	const char *Keyword
		Keyword to be read (has to be uppercase).

	Bool_T CaseSensitive
		If "True" comparison is case sensitive (i.e. uppercase is required in
		the input stream).

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and buffer pointer behind the keyword.

------------------------------------------------------------------------------*/

Bool_T Lexer::GetKeyword( const char *Keyword, Bool_T CaseSensitive )
{
	char *_BufPtr = BufPtr;

	if( CaseSensitive )
	{
		while( *Keyword && *_BufPtr && *Keyword++ == *_BufPtr++ )
			;
	}
	else
		for( ; *Keyword && *_BufPtr; _BufPtr++, Keyword++ )
		{
			char c = *_BufPtr;

			if( islower( c ) )
				c += 'A'-'a';
			if( *Keyword != c )
				goto error;
		}

	if( *Keyword )
		goto error;

//success:
	BufPtr = _BufPtr;
	return True;

error:
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Lexer::GetLabelFree( int lbl, int , int )
	Bool_T Lexer::GetLabelFixed( int lbl, int StartCol, int EndCol )

PURPOSE:
	This functions read in a label when FREE / FIXED format MPS file is read. It
is case sensitive. It sets one of the static LabetPtr table fields to point to
the label. The first character behind the label is required to be either a
space, a tab or a newline. It is overwritten with NULL (to terminate the label).
	In FIXED format we allow the label to contain spaces.

PARAMETERS:
	int lbl
		This is the number of LabelPtr table entry, which will store the
		pointer to label string. The pointer will point into static input data
		buffer.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and static buffer pointer behind the appended
NULL character.

------------------------------------------------------------------------------*/

Bool_T Lexer::GetLabelFree( int lbl, int, int )
{
	assert( lbl >= 0 && lbl <= 1 );

	char *_BufPtr = BufPtr;

	LabelPtr[ lbl ] = _BufPtr;
	for( int i = 0; !isspace( *_BufPtr ) && i < (int)LAB_LEN; i++, _BufPtr++ )
		;

	if( *_BufPtr != ' ' &&  *_BufPtr != '\t' && *_BufPtr != '\n' )
		goto error;

	*_BufPtr++ = '\0';

//success:
	BufPtr = _BufPtr;
	return True;

error:
	return False;
}

Bool_T Lexer::GetLabelFixed( int lbl, int StartCol, int EndCol )
{
	assert( lbl >= 0 && lbl <= 1 && StartCol >= 0 && StartCol <= EndCol &&
		EndCol <= 71 );

	char *_BufPtr = BufPtr = Buf + StartCol;
	int i, len = EndCol - StartCol + 1;

	//--------------------------------------------------------------------------
	//	Strip off the initial spaces.
	//
	for( ; *_BufPtr == ' ' && len; _BufPtr++, len-- )
		;
	LabelPtr[ lbl ] = BufPtr = _BufPtr;	// Store the label string start.

	//--------------------------------------------------------------------------
	//	Scan the label (max. 'len' characters, or until the end of the line,
	//	whichever comes first).
	//
	for( i = 0; *_BufPtr && *_BufPtr != '\n' && i < len; i++, _BufPtr++ )
		;

	//--------------------------------------------------------------------------
	//	Scan backwards to remove any trailing spaces.
	//
	while( *( _BufPtr - 1 ) == ' ' && _BufPtr > BufPtr )
		_BufPtr--;

	if( _BufPtr - BufPtr == 0 )				// Only spaces - no label found.
	{
		LabelPtr[ lbl ] = (char *) EmptyString;
		return False;						// Error!
	}
	else
	{
		*_BufPtr = '\0';					// Null-terminate the label string.
		_BufPtr = BufPtr += i + 1;			// Move the buffer pointer forwards.

		return True;						// Label is OK!
	}
}


/*------------------------------------------------------------------------------

	Bool_T Lexer::GetNumeric( void )

PURPOSE:
	This function accepts a floating point numeric value of the format

	[+|-] {digit}[.{digit}][e[+|-]digit{digit}]

where at least one digit (or decimal point) in whole is required in mantissa and
in exponent.  After the numeric constant is successfully parsed, it is
converted into a double (using ANSI C "atof()" function) and stored in a static
variable "Number". Length and precision of the constant is not checked.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	Moves input stream pointer and buffer pointer behind the last
character thet may possibly belong to the nuumertic constant.

------------------------------------------------------------------------------*/

Bool_T Lexer::GetNumeric( void )
{
	char *_BufPtr = BufPtr;
	char SaveChar;
	Bool_T DigitInMantissa = False;

	const char *Messages[] = {
		"Mantissa empty in numeric constant",
		"Exponent empty in numeric constant"
		};
	int MsgCode;

	if( *_BufPtr == '+' || *_BufPtr == '-' )// Sign (unary + or -).
		do
		{
			_BufPtr++;						// Spaces accepted between sign and
											// number.
		} while( *_BufPtr == ' ' || *_BufPtr == '\t' );

	if( isdigit( *_BufPtr ) )
		DigitInMantissa = True;
	while( isdigit( *_BufPtr ) )			// Accept digits.
		_BufPtr++;

	if( *_BufPtr == '.' )					// Accept decimal point ...
	{
		_BufPtr++;
		DigitInMantissa = True;

		while( isdigit( *_BufPtr ) )		// ... and decimal places.
			_BufPtr++;
	}

	if( !DigitInMantissa )					// At least one digit in mantissa
		{ MsgCode = -1; goto error; }		// required.

	if( *_BufPtr == 'E' || *_BufPtr == 'e' )// Accept exponent in number.
	{
		_BufPtr++;							// Scan exponent character ('E').
		
		if( *_BufPtr == '+' || *_BufPtr == '-' )
			_BufPtr++;						// Scan sign in exponent
		
		if( !isdigit( *_BufPtr ) )			// Exponent must have some digits.
			{ MsgCode = 1; goto error; }

		do
		{
			_BufPtr++;						// Scan exponent digits.
		} while( isdigit( *_BufPtr ) );
	}

	SaveChar = *_BufPtr;
	*_BufPtr = '\0';
	Number = atof( BufPtr );
	*_BufPtr = SaveChar;

//success:
	BufPtr = _BufPtr;
	return True;

error:
	if( MsgCode >= 0 )
		Error( "%s: %s in line %d", FileName, Messages[ MsgCode ], LineNumber );
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T Lexer::GetNewline( Bool_T Force )

PURPOSE:
	This function scans the input stream in search of a newline character. It
accepts spaces and tabs on the way. When a newline character is found, new
input line is fetched from the stream and placed in the buffer. Comment lines
(beginning with '*' character in the first column are skipped.
	Static buffer is filled with new line data (if one is read) and the buffer
pointer is set to the first position in buffer.

PARAMETERS:
	Bool_T Force
		If "True", the rest of the current line is not scanned and a new line
		of input is fetched into the buffer immediately.

RETURN VALUE:
	Boolean success status. If "Force" is "True", failure means an I/O error.
Otherwise it may either be that, or some non-spae characters were encountered
before a newline was found. Also when the line is not terminated with a newline
(i.e. input line is too long to hold it in the buffer) "False" is returned.

SIDE EFFECTS:
	If a new line is fetched, the buffer is filled with new line data and
the buffer pointer is set to the first position in buffer.

------------------------------------------------------------------------------*/

Bool_T Lexer::GetNewline( Bool_T Force )
{
	char *_BufPtr = BufPtr;

	//--------------------------------------------------------------------------
	//	If new line search is not forced, we skip only whitespace, and then
	//	attempt to read a new line from the input stream. Otherwise we read a
	//	new line ignoring the remainder of the current one.
	//
	if( !Force )
	{	
		while( *_BufPtr == ' ' || *_BufPtr == '\t' )
			_BufPtr++;
		if( *_BufPtr && *_BufPtr != '\n' )
			goto error;
	}

	//--------------------------------------------------------------------------
	//	We search for first non-empty line that is not commented. We ignore
	//	comment lines (starting with '*' in the first column). 
	//	Finally we skip lines that contain only whitespace.
	//
	do
	{
		do							// Skip comment lines or empty lines
		{
			_BufPtr = Buf;			// Reset global buffer pointer.
			IO_OK = fgets( Buf, IO_BUF_LEN, Input );
			if( ! IO_OK )
			{
				BufPtr = Buf;
				Buf[0] = '\n';
				Buf[1] = '\0';
				goto error;
			}
			else
				LineNumber++;
		} while( *_BufPtr == '*' || *_BufPtr == '\n' );
									// Also skip lines that contain only
									// whitespace characters.
		while( *_BufPtr == ' ' || *_BufPtr == '\t' )
			_BufPtr++;
	} while( *_BufPtr == '\n' );

//success:
	BufPtr = Buf;
	return True;

error:
	return False;
}
