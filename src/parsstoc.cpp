/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data file parser.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	parsstoc.cpp
CREATED:			1994.08.17
LAST MODIFIED:		1995.08.24

DEPENDENCIES:		lexer.h, error.h, parsstoc.h, scenario.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Stochastic data file parser. Contains one reader function and a set of
functions which set up semantic actions.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Bool_T GetStochFile( const char *FileName, FILE *fp )
	void SetPassRowColLabelsProcedure( void_fptr_const_char_ptr_const_char_ptr )
	void SetPassNumericValueProcedure( void_fptr_double )
	void SetPassPeriodLabelProcedure( void_fptr_const_char_ptr )
	void SetPassProbabilityProcedure( void_fptr_double )
	void SetPassEndData( void_fptr_void )

STATIC FUNCTIONS:
	static Bool_T GetSections( void );
	static Bool_T GetIndepSection( void );
	static Bool_T GetBlockSection( void );
	static Bool_T GetScenarioSection( void );
	static Bool_T GetIgnoredSection( void );
	static Bool_T GetIndepDistcreteSection( void );
	static Bool_T GetIndepUniformSection( void );
	static Bool_T GetIndepNormalSection( void );
	static Bool_T GetIndepSubroutineSection( void );
	static Bool_T GetBlocksDiscreteEntries( void );
	static Bool_T GetIndepDiscreteEntries( void );

STATIC DATA:
	static Scenarios *Scen;

------------------------------------------------------------------------------*/


#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif
#ifndef __PARSSTOC_H__
#	include "parsstoc.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif


//==============================================================================
//
//	Function prototypes.
//
//==============================================================================

static Bool_T GetSections( void );
static Bool_T GetIndepSection( void );
static Bool_T GetBlockSection( void );
static Bool_T GetScenarioSection( void );
static Bool_T GetIgnoredSection( void );
static Bool_T GetIndepDistcreteSection( void );
static Bool_T GetIndepUniformSection( void );
static Bool_T GetIndepNormalSection( void );
static Bool_T GetIndepSubroutineSection( void );
static Bool_T GetBlocksDiscreteEntries( void );
static Bool_T GetIndepDiscreteEntries( void );

//==============================================================================
//
//	End of function prototypes.
//
//==============================================================================


//------------------------------------------------------------------------------
//	Static pointer to a scenario repository class.
//
static Scenarios *Scen;
// 
//------------------------------------------------------------------------------


/*------------------------------------------------------------------------------

	Bool_T GetStochFile( const char *FileName, FILE *fp, Scenarios *Scen )

PURPOSE:
	Parses a scenario file, which contains 

PARAMETERS:
	const char *FileName, FILE *fp
		Scenario file name and pointer to the corresponding stream (assumed to
		be opened for reading).
	
	Scenarios *Scen
		Pointer to an object of scenario repository class.

RETURN VALUE:
	Success status ("True" on success, "False" otherwise).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T GetStochFile( const char *FileName, FILE *fp, Scenarios *sc )
{
	assert( sc != NULL );
	Scen = sc;

	Lexer::SetInputStream( FileName, fp );

	if( !Lexer::GetKeyword( "STOCH" ) )
	{
		Error( "File %s: Invalid first line.", Lexer::FileName );
		goto Error;
	}

	Lexer::GetNewline( True );

	while( GetSections() )
		;

	if( !Lexer::GetKeyword( "ENDATA" ) )
	{
		Error( "File %s, line %d: ENDATA expected.",
			Lexer::FileName, Lexer::LineNumber );
		goto Error;
	}
	Scen->EndData();

	return True;

Error:
	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetSections( void )

PURPOSE:
	The stochastic data file consists of a number of sections of different
random variables, groups of variables, scenarios etc.

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Bool_T GetSections( void )
{
	if( !GetIndepSection() && !GetBlockSection() && !GetScenarioSection() )
		return False;
	else
		return True;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepSection( void )

PURPOSE:
	There is a number of possible independent random variable sections. They all
start with the keyword "INDEP" in the first column of the input line. This
function reads an "INDEP" section.

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Bool_T GetIndepSection( void )
{
	if( !Lexer::GetKeyword( "INDEP" ) )
		return False;

	if( !Lexer::GetSpace() )
	{
		Error( "File %s, line %d: Separator missing after INDEP.",
			Lexer::FileName, Lexer::LineNumber );
		return False;
	}

	if( GetIndepUniformSection() || GetIndepNormalSection() ||
		GetIndepSubroutineSection() || GetIndepDistcreteSection() )
		return True;
	else
	{
		Error( "File %s, line %d: Unknown type of INDEP section encoutered.",
			Lexer::FileName, Lexer::LineNumber );
		return False;
	}
}


/*------------------------------------------------------------------------------

	static Bool_T GetBlockSection( void )

PURPOSE:
	There is a number of possible sections of blocks of random variables. They
all start with the "BLOCKS" keyword in the first column of the input line. This
procedure reads in a block section (after recognizing the apprppriate type).

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Bool_T GetBlockSection( void )
{
	if( !Lexer::GetKeyword( "BLOCKS" ) )
		goto Error;

	if( !Lexer::GetSpace() || !Lexer::GetKeyword( "DISCRETE" ) )
	{
		Error( "File %s, line %d: Only BLOCKS DISCRETE section recognized.",
			Lexer::FileName, Lexer::LineNumber );
		goto Error;
	}

	if( !Lexer::GetNewline() )
	{
		Warning( "File %s, line %d: Extra characters after BLOCKS "
			"DISCRETE ignored.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
	}

	if( !GetBlocksDiscreteEntries() )
		goto Error;

	return True;

Error:
	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetScenarioSection( void )

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

static Bool_T GetScenarioSection( void )
{
	if( Lexer::GetKeyword( "SCENARIOS" ) )
	{
		Error( "File %s, line %d: SCENARIOS sec. recognition not implemented.",
			Lexer::FileName, Lexer::LineNumber );

		while( GetIgnoredSection() )
			;
	}

	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIgnoredSection( void )

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

static Bool_T GetIgnoredSection( void )
{
	while( Lexer::GetSpace() && Lexer::GetNewline( True ) )
		;

	return True;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepUniformSection( void )

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

static Bool_T GetIndepUniformSection( void )
{
	if( Lexer::GetKeyword( "UNIFORM" ) )
		Error( "File %s, line %d: INDEP UNIFORM sections not recognized yet.",
			Lexer::FileName, Lexer::LineNumber );

	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepNormalSection( void )

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

static Bool_T GetIndepNormalSection( void )
{
	if( Lexer::GetKeyword( "NORMAL" ) )
		Error( "File %s, line %d: INDEP NORMAL sections not recognized yet.",
			Lexer::FileName, Lexer::LineNumber );

	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepSubroutineSection( void )

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

static Bool_T GetIndepSubroutineSection( void )
{
	if( Lexer::GetKeyword( "SUB" ) )
		Error( "File %s, line %d: INDEP SUB sections not recognized yet.",
			Lexer::FileName, Lexer::LineNumber );

	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepDistcreteSection( void )

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

static Bool_T GetIndepDistcreteSection( void )
{
	if( !Lexer::GetKeyword( "DISCRETE" ) )
		goto Error;

	if( !Lexer::GetNewline() )
	{
		Warning( "File %s, line %d: Extra input ignored after INDEP DISCRETE.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
	}

	while( GetIndepDiscreteEntries() )
		;

	return True;

Error:
	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetIndepDiscreteEntries( void )

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

static Bool_T GetIndepDiscreteEntries( void )
{
	const char *PeriodLabel		= NULL,
		*ColLabel				= NULL,
		*RowLabel				= NULL;
	Real_T Probability		 	= 0.0;
	Real_T Value				= 0.0;

	if( !Lexer::GetSpace() )
		goto Error;

	if( !Lexer::GetLabelFree( 0 ) || !Lexer::GetSpace() ||
		!Lexer::GetLabelFree( 1 ) )
	{
		Error( "File %s, line %d: Invalid line in INDEP DISCRETE section.",
			Lexer::FileName, Lexer::LineNumber );
		goto Error;
	}

	ColLabel = Lexer::LabelPtr[0];
	RowLabel = Lexer::LabelPtr[1];

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "File %s, line %d: Distribution value missing.",
			Lexer::FileName, Lexer::LineNumber );
		goto Error;
	}

	Value = Lexer::Number;

	if( !Lexer::GetSpace() )
	{
		Error( "File %s, line %d: Missing delimiter",
			Lexer::FileName, Lexer::LineNumber );
		goto Error;
	}

	PeriodLabel = "";

	if( !Lexer::GetNumeric() )
	{
		if( !Lexer::GetLabelFree( 0 ) )
		{
			Error( "File %s, line %d: Period label expected",
				Lexer::FileName, Lexer::LineNumber );
			goto Error;
		}

		PeriodLabel =  Lexer::LabelPtr[0];

		if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
		{
			Error( "File %s, line %d: Probability missing.",
				Lexer::FileName, Lexer::LineNumber );
			goto Error;
		}
	}

	Probability = Lexer::Number;

	/* Semantic action */
	Scen->NewIndepDiscreteEntry( RowLabel, ColLabel, PeriodLabel, Value,
		Probability );

	if( !Lexer::GetNewline() )
	{
		Warning( "File %s, line %d: Extra characters after line.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
	}

	return True;

Error:
	return False;
}


/*------------------------------------------------------------------------------

	static Bool_T GetBlocksDiscreteEntries( void )

PURPOSE:
	x

PARAMETERS:
	None.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

static Bool_T GetBlocksDiscreteEntries( void )
{
	Bool_T BL_Found = True;

	//--------------------------------------------------------------------------
	//	The first line of the BLOCK DISCRETE section entry should have the
	//	following form:
	//		<spc> BL <spc> <block name> <spc> <period name> <spc> <probability>
	//
	if( !Lexer::GetSpace() || !Lexer::GetKeyword( "BL" ) ) goto Error;

	while( BL_Found )
	{
		BL_Found = False;

		if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) ||
			!Lexer::GetSpace() || !Lexer::GetLabelFree( 1 ) )
		{
			Error( "File %s, line %d: Block label and period label expected.",
				Lexer::FileName, Lexer::LineNumber );
			goto Error;
		}

		if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
		{
			Error( "File %s, line %d: Block's probability expected.",
				Lexer::FileName, Lexer::LineNumber );
			goto Error;
		}

		/* Semantic action */
		Scen->NewBlocksDiscrete( Lexer::LabelPtr[0], Lexer::LabelPtr[1],
			Lexer::Number );

		if( !Lexer::GetNewline() )
		{
			Warning( "File %s, line %d: Extra characters after block "
				"identifier line ignored", Lexer::FileName, Lexer::LineNumber );
			Lexer::GetNewline( True );
		}

		//----------------------------------------------------------------------
		//	Now, in loop we will input the lines of one block's description.
		//	They should look like this:
		//		<spc> <column label> <spc> <row label> <spc> <value>
		//
		for(;;)
		{
			if( !Lexer::GetSpace() ) break;
			if( Lexer::GetKeyword( "BL" ) )
			{
				BL_Found = True;
				break;
			}

			if( !Lexer::GetLabelFree( 0 ) || !Lexer::GetSpace() ||
				!Lexer::GetLabelFree( 1 ) || !Lexer::GetSpace() ||
				!Lexer::GetNumeric() )
			{
				Error( "File %s, line %d: Random variable expected in block.",
					Lexer::FileName, Lexer::LineNumber );
				goto Error;
			}

			/* Semantic action */
			Scen->NewBlocksDiscreteEntry( Lexer::LabelPtr[1],
				Lexer::LabelPtr[0], Lexer::Number );
	
			if( !Lexer::GetNewline() )
			{
				Warning( "File %s, line %d: Extra characters after a data item "
					"in block ignored", Lexer::FileName, Lexer::LineNumber );
				Lexer::GetNewline( True );
			}
		}
	}

	return True;

Error:
	return False;
}
