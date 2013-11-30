/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	sc_tree.cpp
CREATED:			1996.01.20
LAST MODIFIED:		1996.04.09

DEPENDENCIES:		sc_tree.h, extendar.h, error.h, smartptr.h, smartdcl.h,
					sptr_ndb.h, myalloc.h, std_tmpl.h, sort_lab.h, stdtype.h,
					simplex.h, solv_lp.h, compile.h, lp_codes.h, mps_lp.h,
					std_math.h, periods.h, changelp.h, my_defs.h, lexer.h,
					print.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SC_TREE_H__
#	include "sc_tree.h"
#endif
#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __PERIODS_H__
#	include "periods.h"
#endif
#ifndef __CHANGELP_H__
#	include "changelp.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __MAKELAB_H__
#	include "makelab.h"
#endif


#define INIT_SC		10
#define INIT_BL		100


//------------------------------------------------------------------------------
//	Static data of class "ScenarioTree".
//
Bool_T ScenarioTree::ParseError		= False;
const SolvableLP *ScenarioTree::LP	= NULL;
const Periods *ScenarioTree::PER	= NULL;
//
//	End of the static data of class "ScenarioTree".
//------------------------------------------------------------------------------

/*------------------------------------------------------------------------------

	ScenarioTree::ScenarioTree( void )

PURPOSE:
	Object constructor. Initializes an empty object.

PARAMETERS:
	None.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

ScenarioTree::ScenarioTree( void )
	: ScenNum( 0 ), BlockNum( 0 ), PeriodNum( 0 ),
	ScenLab( INIT_SC ), Structure( 0, INIT_SC ),
	StochasticData( 0, INIT_BL )
{}


/*------------------------------------------------------------------------------

	void ScenarioTree::Clean( void )

PURPOSE:
	Performs all the action of a destructor: frees all allocated memory.
Also restores the object to its initial state.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::Clean( void )
{
	Int_T i;

	for( i = 0; i < ScenNum; i++ )
		if( Structure[i] )
			delete Structure[i];

	for( i = 0; i < BlockNum; i++ )
	{
		ExtendArray<ChangeLP *> *clp = StochasticData[i];

		if( clp != NULL )
		{
			for( Int_T j = 0, l = clp->Len(); j < l; j++ )
				if( (*clp)[j] )
				{
					delete (*clp)[j];
					(*clp)[j] = NULL;
				}

			delete clp;
		}
	}

	ScenNum = BlockNum = PeriodNum = 0;
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::ReadScenarioFile( FILE *fp, const char *fname,
		const SolvableLP &lp, const Periods &per, VerbLevel Verbosity )

PURPOSE:
	Reads the file containing the scenario data.

PARAMETERS:
	FILE *fp, const char *fname
		Input stream and associated file name (the file name is only used for
		messages, so it is irrelevant whether it actually corresponds to a file
		name).

	const SolvableLP &lp, const Periods &per
		CORE and TIME file data (needed for scenario file parsing).

	VerbLevel Verbosity
		Verbosity level - affects the on-screen report.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::ReadScenarioFile( FILE *fp, const char *fname, // )
	const SolvableLP &lp, const Periods &per, VerbLevel Verbosity )
{
	assert( fp != NULL );
	assert( fname != NULL );

	if( Verbosity >= V_LOW )
		Print( "Reading the scenarios from '%s'.\n", fname );

	if( PeriodNum != 0 || BlockNum != 0 || ScenNum != 0 )
		Clean();

	if( ( PeriodNum = per.NumberOfPeriods() ) == 0 )
		return False;

	Lexer::SetInputStream( fname, fp );
	Bool_T success = ReadFile( lp, per );

	if( success )
	{
		ScenLab.SortLabels();
		if( ScenLab.Duplicates() )
		{
			Error( "Duplicate scenario labels found in the scenario tree." );
			success = False;
		}
	}

	if( success ) success = CheckTree();

	if( success )
	{
		CalculateProbabilities();
		PrintInfo();
	}
	else
		Print( "Scenario reading failed!\n" );

	return success;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::PrintInfo( void )

PURPOSE:
	Prints some statistics of the current scenario tree.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::PrintInfo( void )
{
	const Real_T AvgBranch = pow( double(ScenNum), 1.0/double(PeriodNum-1) );

	Print(
		"\t%-30s%10d\n"
		"\t%-30s%10d\n"
		"\t%-30s%10G\n",

		"Number of periods:",			PeriodNum,
		"Number of scenarios:",			ScenNum,
		"Average branching factor:",	(double) AvgBranch
	);

/*
	{
		Print( "\nScenario structure:\n" );

		for( Int_T i = 0; i < ScenNum; i++ )
		{
			for( Int_T p = 0; p < PeriodNum; p++ )
			{
				ScStruct &s = Structure[i][0][p];
				Print( "  %5d", (int)s.scen );
			}
			Print( "\n" );
		}
		Print( "\n" );
	}

	{
		Print( "\nConditional probabilities:\n" );

		for( Int_T i = 0; i < ScenNum; i++ )
		{
			for( Int_T p = 0; p < PeriodNum; p++ )
			{
				ScStruct &s = Structure[i][0][p];

				Print( " %5d(%.5f)", (int)s.scen, (double)s.prob );
			}
			Print( "\n" );
		}
		Print( "\n" );
	}
*/
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::CheckTree( void )

PURPOSE:
	This function checks consistency of the tree structure.
	Additionally, for the converter to work properly, the tree must have a
special structure, which is checked by this procedure. Among the requirements we
need to mention perfect symmetry of the tree as the most important one.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status ("True" if tree is self-consistent).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::CheckTree( void )
{
	Int_T p, i;

	//--------------------------------------------------------------------------
	//	See if there is only one scenario coming out of the root.
	//
	for( i = 1; i < ScenNum; i++ )
	{
		ScStruct &s = Structure[i][0][0];

		if( s.scen != 0 )
		{
			Error( "More than one scenario coming from the ROOT node." );
			return False;
		}
	}

	//--------------------------------------------------------------------------
	//	See if the tree has been generated in such order, that the sequence in
	//	the columns of the "Structure" array is nondecreasing.
	//
	for( p = 0; p < PeriodNum; p++ )
	{
		Int_T LastSc = Structure[0][0][p].scen;

		assert( LastSc >= 0 && LastSc < ScenNum );

		for( i = 1; i < ScenNum; i++ )
		{
			ScStruct &s = Structure[i][0][p];

			if( s.scen < LastSc )
			{
				Error( "Sequence in columns of the scenario tree should be "
					"non-decreasing." );
				return False;
			}
			else if( s.scen > LastSc )
			{
				LastSc = s.scen;
				assert( LastSc >= 0 && LastSc < ScenNum );
			}
		}
	}

	//--------------------------------------------------------------------------
	//	Check whether the tree is symmetric (at each stage the same number of
	//	branches comes out of the node).
	//
	for( p = 1; p < PeriodNum; p++ )
	{
		Int_T LastSc	= Structure[0][0][p].scen,
			same		= 0;

		while( Structure[same][0][p].scen == LastSc && same < ScenNum )
			same++;

		for( i = same; i < ScenNum; )
		{
			if( LastSc == Structure[i][0][p].scen )
			{
				Error( "Tree non-symmetric at period %d.", p );
				return False;
			}

			LastSc = Structure[i][0][p].scen;

			for( Int_T j = 0; j < same; j++, i++ )
				if( Structure[i][0][p].scen != LastSc )
				{
					Error( "Tree non-symmetric at period %d.", p );
					return False;
				}
		}
	}

	return True;
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::ReadFile( const SolvableLP &lp, const Periods &per )

PURPOSE:
	This function parses the whole scenario file and stores the data in the data
structures. 

PARAMETERS:
	const SolvableLP &lp
		The CORE file data associated with the scenarios.

	const Periods &per
		The TIME file data.

	VerbLevel Verbosity
		Verbosity level that affects the output.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::ReadFile( const SolvableLP &lp, const Periods &per )
{
	ParseError	= False;
	LP			= &lp;
	PER			= &per;

	while( GetSections() )
		;

	return ( ParseError ? False : True );
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::GetSections( void )

PURPOSE:
	The stochastic data file consists of a number of sections of different
random variables, groups of variables, scenarios etc.

	So far we only allow oen type of sections: "SCENARIOS DISCRETE"!

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::GetSections( void )
{
	return GetScenarioSection();
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::GetScenarioSection( void )

PURPOSE:
	Gets a single SCENARIOS section.

PARAMETERS:
	None.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::GetScenarioSection( void )
{
	if( Lexer::GetKeyword( "SCENARIOS" ) )
	{
		if( !Lexer::GetSpace() || !Lexer::GetKeyword( "DISCRETE" ) )
			Error( "File %s, line %d: 'DISCRETE' expected after 'SCENARIOS'",
				Lexer::FileName, Lexer::LineNumber );

		Lexer::GetNewline( True );

		//----------------------------------------------------------------------
		//	Read in loop the consecutive scenarios. The first line has to start
		//	with 'SC'.
		//
		Int_T  AttachAtPeriod = -1;

		if( !Lexer::GetSpace() || !Lexer::GetKeyword( "SC" ) )
		{
			Error( "File %s, line %d: section 'SCENARIOS' should begin with "
				"'SC'.", Lexer::FileName, Lexer::LineNumber );
		}
		else if( !GetScenarioIndicator( AttachAtPeriod ) )
			return False;

		while( Lexer::GetSpace() )
			if( Lexer::GetKeyword( "SC" ) )
			{
				if( !GetScenarioIndicator( AttachAtPeriod ) )
					return False;
			}
			else
			{
				if( !GetScenarioEntry( AttachAtPeriod ) )
					return False;
			}

		return True;
	}

	return False;
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::GetScenarioIndicator( Int_T &AttachAtPeriod )

PURPOSE:
	Reads a single scenario (the starting "SC" keyword has been read before
entering the procedure).

PARAMETERS:
	Int_T &AttachAtPeriod
		At which period is this scenario attached to its predecessor(this value
		is returned to the calling function).

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::GetScenarioIndicator( Int_T &AttachAtPeriod )
{
	//--------------------------------------------------------------------------
	//	Since we're quite sure it's a scenario we're processing, we now assume
	//	strict compliance with the syntax.
	//
	//	Still reading the first line.
	//
	if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) ||
		!Lexer::GetSpace() || !Lexer::GetLabelFree( 1 ) )
	{
		Error( "File %s, line %d: Scenario and relative root labels expected.",
			Lexer::FileName, Lexer::LineNumber );
		return False;
	}

	const char *ScLab	= Lexer::LabelPtr[0],
		*RootLab		= Lexer::LabelPtr[1];

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "File %s, line %d: Numeric value expected.",
			Lexer::FileName, Lexer::LineNumber );
		return False;
	}

	const Real_T prob = Lexer::Number;

	if( prob < 0.0 || prob > 1.0 )
	{
		Error( "File %s, line %d: Probability outside <0,1> range.",
			Lexer::FileName, Lexer::LineNumber );
		ParseError = True;
	}

	if( !Lexer::GetSpace() || !Lexer::GetLabelFree( 0 ) )
	{
		Error( "File %s, line %d: Period label expected.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
		return False;
	}

	const char *PeriodLab	= Lexer::LabelPtr[0];

	//-------------------------------------------------------------------------
	//	Now store the data that has just been read.
	//
	AddNewScenario( ScLab, RootLab, PeriodLab, prob, AttachAtPeriod );

	//-------------------------------------------------------------------------
	//	Read in the new line and exit.
	//
	Lexer::GetNewline( True );
	return True;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::IgnoreTillEndOfSection( void )

PURPOSE:
	Ignores the data till the end of the current section.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::IgnoreTillEndOfSection( void )
{
	Lexer::GetNewline( True );
	while( Lexer::GetSpace() )
		Lexer::GetNewline( True );
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::GetScenarioEntry( Int_T AttachAtPeriod )

PURPOSE:
	Reads a single line of the scenario.

PARAMETERS:
	Int_T AttachAtPeriod
		At which period is the random data allowed to be found.

RETURN VALUE:
	Boolean success status. Also returns "False" if a new scenario start is
detected.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::GetScenarioEntry( Int_T AttachAtPeriod )
{
	//-------------------------------------------------------------------------
	//	Read the entry (one line).
	//
	if( !Lexer::GetLabelFree( 0 ) ||
		!Lexer::GetSpace() || !Lexer::GetLabelFree( 1 ) )
	{
		Error( "File %s, line %d: Two labels expected in scenario data line.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
		return False;
	}

	if( !Lexer::GetSpace() || !Lexer::GetNumeric() )
	{
		Error( "File %s, line %d: Value expected in scenario data line.",
			Lexer::FileName, Lexer::LineNumber );
		Lexer::GetNewline( True );
		return False;
	}

	//-------------------------------------------------------------------------
	//	Add the entry to the data structures.
	//
	const char *Column	= Lexer::LabelPtr[0],
		*Row			= Lexer::LabelPtr[1];
	const Real_T Val	= Lexer::Number;

	AddEntry( Column, Row, Val, AttachAtPeriod );

	Lexer::GetNewline( True );

	return True;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::AddNewScenario( const char *ScLab, const char *RootLab,
		const char *PeriodLab, Real_T prob )

PURPOSE:
	Stores the data read from a single scenario indicator line in the
"Structure" array. The data is checked for consistency.

PARAMETERS:
	const char *ScLab, const char *RootLab, const char *PeriodLab
		Current scenario, its predecessor and period labels, respectively.

	Real_T prob

	Int_T &AttachAtPeriod
		At which period is this scenario attached to its predecessor (this value
		is returned to the calling function).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::AddNewScenario( const char *ScLab, const char *RootLab, // )
	const char *PeriodLab, Real_T prob, Int_T &AttachAtPeriod )
{
	//--------------------------------------------------------------------------
	//	Look for the predecessor scenario.
	//
	Int_T RootNum	= -1;

	AttachAtPeriod	= -1;

	if( ScenNum == 0 )
	{
		if( strcmp( RootLab, "ROOT" ) != 0 )
		{
			Error( "File %s, line %d: The first scenario must have ROOT as its "
				"predecessor.", Lexer::FileName, Lexer::LineNumber );
			ParseError = True;
		}
		RootNum			= 0;
		AttachAtPeriod	= 0;
	}
	else
		RootNum = (Int_T) ScenLab.FindLabel( RootLab );

	assert( ScenNum == 0 || RootNum < ScenNum );

	if( RootNum < 0 )
	{
		Error( "File %s, line %d: Predecessor scenario not found.",
			Lexer::FileName, Lexer::LineNumber );
		ParseError = True;
		RootNum = 0;
	}

	//--------------------------------------------------------------------------
	//	Now add the current scenario label to the list.
	//
	ScenLab.AddLabel( ScLab );

	//--------------------------------------------------------------------------
	//	See if the attachment period specified is correct.
	//
	if( AttachAtPeriod < 0 )
		AttachAtPeriod = PER->PeriodNumber( PeriodLab );

	assert( AttachAtPeriod < PeriodNum );

	if( AttachAtPeriod < 0 )
	{
		Error( "File %s, line %d: Invalid period label.",
			Lexer::FileName, Lexer::LineNumber );
		ParseError = True;
		AttachAtPeriod = 0;
	}

	//--------------------------------------------------------------------------
	//	Expand the "Structure" so that it may contain the new scenario data.
	//	Perform the attachment and fill the new entry.
	//
	Structure.Resize( ScenNum+1 );
	Structure[ScenNum] = new ExtendArray<ScStruct>( PeriodNum );
	if( Structure[ScenNum] == NULL ) FatalError( "Out of memory." );

	for( Int_T i = 0; i < PeriodNum; i++ )
	{
		ScStruct &s = Structure[ScenNum][0][i];

		if( i < AttachAtPeriod )
		{
			s = Structure[RootNum][0][i];
		}
		else
		{
			s.scen	= ScenNum;
			s.block	= -1;
			s.prob	= 0.0;
		}
	}

	Structure[ScenNum][0][PeriodNum-1].prob = prob;

	ScenNum++;

	//--------------------------------------------------------------------------
	//	Now expand the "StochasticData" array to accomodate data blocks
	//	corresponding to the new scenario.
	//
	Int_T NewBlockNum = Int_T( BlockNum + PeriodNum - AttachAtPeriod );

	assert( (Int_T)StochasticData.Len() == BlockNum );

	StochasticData.Resize( NewBlockNum );
	StochasticData.FillFragment( (ExtendArray<ChangeLP *> *)NULL, BlockNum,
		NewBlockNum );

	BlockNum = NewBlockNum;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::AddEntry( const char *Column, const char *Row,
		Real_T Val )

PURPOSE:
	Stores the data read from a single scenario data line in the "Structure"
array as well as in the "StochasticData" array.

PARAMETERS:
	const char *Column, const char *Row, Real_T Val
		The radom data (row and column labels and the value).

	Int_T AttachAtPeriod
		At which period is the random data allowed to be found.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::AddEntry( const char *Column, const char *Row, // )
	Real_T Val, Int_T AttachAtPeriod )
{
	//--------------------------------------------------------------------------
	//	Check the arguments. Convert the labels to numbers and create a
	//	"ChangeLP" object which stores the data. See what period does the
	//	"ChangeLP" structure belong to.
	//
	assert( ScenNum >= 1 );
	assert( AttachAtPeriod >= 0 );

	ChangeLP *NewCLP = NULL;

	SortedArrayOfLabels &ColLabs	= ((SolvableLP *)LP)->RevealColumnLabels(),
		&RowLabs					= ((SolvableLP *)LP)->RevealRowLabels();

	Int_T r	= (Int_T) RowLabs.FindLabel( Row ),
		c	= (Int_T) ColLabs.FindLabel( Column ),
		p	= -1;

	assert( r < LP->GetM() );
	assert( c < LP->GetN() );

	if( strcmp( LP->RevealRHS_Name(), Column ) == 0 )
	{
		if( r < 0 ) goto RowLabelNotFound;
		NewCLP = new ChangeRHS( r, Val );

		p = PER->RowInPeriod( r );
	}
	else if( strcmp( LP->RevealObjectiveName(), Row ) == 0 )
	{
		if( c < 0 ) goto ColLabelNotFound;
		NewCLP = new ChangeCost( c, Val );

		p = PER->ColumnInPeriod( c );
	}
	else
	{
		if( r < 0 ) goto RowLabelNotFound;
		if( c < 0 ) goto ColLabelNotFound;
		NewCLP = new ChangeMatrix( r, c, Val );

		p = Max( PER->RowInPeriod( r ), PER->ColumnInPeriod( c ) );
	}

	assert( p >= 0 && p < PeriodNum );
	if( p < AttachAtPeriod  )
		Error( "File %s, line %d: Random data belongs to a period earlier"
			"than declared.", Lexer::FileName, Lexer::LineNumber );

	if( NewCLP == NULL ) FatalError( "Out of memory. ");

	//--------------------------------------------------------------------------
	//	Store the new "ChangeLP" object in the "StochasticData" array.
	//
	//	Notice careful calculation of the new stochastic data block to be
	//	assigned for the current scenario and period!
	//
	//	If necessary, a new stochastic data block is allocated.
	//
	{
		Int_T block = Int_T( BlockNum + p - PeriodNum );
		ExtendArray<ChangeLP *> *&b = StochasticData[block];

		if( b == NULL )
		{
			assert( Structure[ScenNum-1][0][p].block == -1 );

			Structure[ScenNum-1][0][p].block = block;

			b = new ExtendArray<ChangeLP *>( 0, 20 );
			if( b == NULL ) FatalError( "Out of memory." );
		}

		Int_T pos = (Int_T)b->Len();

		b->Resize( pos + 1 );
		(*b)[pos] = NewCLP;
	}

	return;

	//--------------------------------------------------------------------------
	//	Labels for some simple error handling inside the function.
	//
RowLabelNotFound:
	Error( "File %s, line %d: Row label not found.",
		Lexer::FileName, Lexer::LineNumber );
	goto End;

ColLabelNotFound:
	Error( "File %s, line %d: Column label not found.",
		Lexer::FileName, Lexer::LineNumber );
	goto End;

End:
	ParseError = True;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::CalculateProbabilities( void )

PURPOSE:
	When the scenario file is read, only the probabilities of the whole
scecnarios are known. For some purposes (including conversion to a different
number of stages) we need to know the conditional probabilities at all nodes of
the scenario tree.
	This procedure will calculate those probabilities. They will be stored in
the "Structure" array, in "prob" field of the "ScStruct" objects. The root node
will have a conditional probability "one". Then each node "x" which can be
reached by branching from node "y" will hold the conditional probability of
reaching node "x", subject to the fact that the node "y" was reached before.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::CalculateProbabilities( void )
{
	Int_T p;

	//--------------------------------------------------------------------------
	//	Pass one: pass back the probabilities.
	//
	//	The probabilities "flow" from the leaves to the root of the tree.
	//
	for( p = Int_T( PeriodNum - 1 ); p > 0; p-- )
		for( Int_T i = 0; i < ScenNum; i++ )
			Structure[i][0][p-1].prob = Structure[i][0][p].prob;

	//--------------------------------------------------------------------------
	//	Pass two: sum up the probabilities.
	//
	//	The probabilities corresponding to a single node in the tree are
	//	consolidated (added) and stored.
	//
	//	After this phase we have unconditional probability at each node.
	//
	for( p = 0; p < PeriodNum-1; p++ )
		for( Int_T start = 0; start < ScenNum; )
		{
			Real_T sum		= 0.0;
			Int_T LastSc	= Structure[start][0][p].scen;

			Int_T len;
			for( len = 0; start + len < ScenNum &&
				Structure[start + len][0][p].scen == LastSc; len++ )
				sum += Structure[start + len][0][p].prob;

			for( ; len > 0; len-- )
				Structure[start++][0][p].prob = sum;
		}

	//--------------------------------------------------------------------------
	//	Probability errors are checked (but not corrected) here.
	//	Probabilities in all columns of the "Structure" array should add up to
	//	one.
	//
	for( p = 0; p < PeriodNum; p++ )
	{
		Real_T TotalProb	= 0.0;
		Int_T LastSc		= -1;

		for( Int_T i = 0; i < ScenNum; i++ )
			if( LastSc != Structure[i][0][p].scen )
			{
				LastSc		= Structure[i][0][p].scen;
				TotalProb	+= Structure[i][0][p].prob;
			}

		if( fabs( TotalProb - 1.0 ) > 1e-3 )
			Warning( "Probabilities do not sum up to one: %.6f at period %d.",
				TotalProb, p);
	}

	//--------------------------------------------------------------------------
	//	Pass three: scale.
	//
	//	Conditional probabilities are calculated.
	//
	for( p = Int_T( PeriodNum - 1 ); p >= 1; p-- )
		for( Int_T i = 0; i < ScenNum; i++ )
			Structure[i][0][p].prob /= Structure[i][0][p-1].prob;
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::ConvertToTwoStage( SolvableLP &lp,
		const Periods &periods, Int_T CutAt, FILE *fpCOR, FILE *fpSTO,
		FILE *fpTIM, VerbLevel Verbosity )

	Bool_T ScenarioTree::OutputCORE( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt, VerbLevel Verbosity ) const
	Bool_T ScenarioTree::OutputTIME( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt, VerbLevel Verbosity ) const
	Bool_T ScenarioTree::OutputSTOCH( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt, VerbLevel Verbosity ) const


PURPOSE:
	The first procedure expands the stochastic problem to a two-stage equivalent
problem (or, if the division point is set behind the last period, to a
deterministic equivalent). It uses the remaining three procedures to output
consecutive files of the two stage equivalent problem.

	WARNING: the linear problem that is written to the file may become very
large. Tens of megabytes are nothing unusual.

PARAMETERS:
	SolvableLP &lp
		Core file data.

	const Periods &periods
		the structure defining the time periods.

	Int_T CutAt
		At which period the problem is to be divided. The argument may range
		from 1 to "PeriodNum" (inclusive).

	FILE *fpCOR, FILE *fpSTO, FILE *fpTIM
		Files to write COR, STO and TIM to, respectively. If a deterministic
		equivalent is to be produced, the last two streams may be omitted (NULL
		pointers may be given instead of the real names).

	FILE *fp
		One of the above (different in each one of the last three functions).

	VerbLevel Verbosity
		Directs the level of verbosity (detail of possible reportig on the
		progress of the linear problem output routine).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/


Bool_T ScenarioTree::ConvertToTwoStage( SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt, FILE *fpCOR, FILE *fpSTO, FILE *fpTIM,
	VerbLevel Verbosity )
const
{
	assert( CutAt >= 1 && CutAt <= PeriodNum );
	assert( fpCOR != NULL );
	assert( fpSTO != NULL || CutAt == PeriodNum );
	assert( fpTIM != NULL || CutAt == PeriodNum );

	//--------------------------------------------------------------------------
	//	Output the CORE file.
	//
	if( !OutputCORE( fpCOR, lp, periods, CutAt, Verbosity ) )
		return False;

	//--------------------------------------------------------------------------
	//	TIME and STOCH files are to be prodeuced only when a two stage
	//	equivalent is requested. When "CutAt" == "PeriodNum" (i.e., when a
	//	deterministic equivalent is to be produced) the rest of this procedure
	//	is skipped.
	//
	if( CutAt == PeriodNum )
		return True;

	//--------------------------------------------------------------------------
	//	Output the TIME file.
	//
	if( !OutputTIME( fpTIM, lp, periods, CutAt, Verbosity ) )
		return False;

	//--------------------------------------------------------------------------
	//	Output the STOCH file.
	//
	if( !OutputSTOCH( fpSTO, lp, periods, CutAt, Verbosity ) )
		return False;

	return True;
}


Bool_T ScenarioTree::OutputCORE( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt, VerbLevel Verbosity )
const
{
	if( Verbosity >= V_HIGH )
		Print( "\nWriting output CORE file.\n" );

	//--------------------------------------------------------------------------
	//	Output the header line.
	//
	if( fprintf( fp, "%-14s%-8s\n", "NAME", lp.GetName() ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Output the ROWS section.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ROWS\n" );

	if( !OutputROWS( fp, lp, periods, CutAt ) )
		return False;

	//--------------------------------------------------------------------------
	//	Output the COLUMNS section.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting COLUMNS\n" );

	if( !OutputCOLUMNS( fp, lp, periods, CutAt ) )
		return False;

	//--------------------------------------------------------------------------
	//	Section RHS output.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting RHS\n" );

	if( !OutputRHS( fp, lp, periods, CutAt ) )
		return False;

	//--------------------------------------------------------------------------
	//	BOUNDS section.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting BOUNDS\n" );

	if( !OutputBOUNDS( fp, lp, periods, CutAt ) )
		return False;

	//--------------------------------------------------------------------------
	//	Write "ENDATA" and finish processing.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ENDATA\n" );

	if( fprintf( fp, "ENDATA\n" ) == EOF )
		return False;

	return True;
}


Bool_T ScenarioTree::OutputTIME( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt, VerbLevel Verbosity )
const
{
	//--------------------------------------------------------------------------
	//	Output the whole file at once.
	//
	if( Verbosity >= V_HIGH )
		Print( "\nWriting output TIME file.\n" );

	if( fprintf( fp, "TIME          %s\nPERIODS       LP\n", lp.GetName() )
			== EOF ||
		fprintf( fp, "    %s", MakeLabel( 0, 0 ) ) == EOF ||
		fprintf( fp, "  %s  %15s%-8s\n", MakeLabel( 0, 0 ), "", "TIME0" )
			== EOF ||
		fprintf( fp, "    %s", MakeLabel( (int)periods.FirstColumn(CutAt), 0 ) )
			== EOF ||
		fprintf( fp, "  %s  %15s%-8s\n", MakeLabel(
			(int)periods.FirstRow(CutAt), 0 ), "", "TIME1" ) == EOF ||
		fprintf( fp, "ENDATA\n" ) == EOF )
		return False;

	return True;
}


Bool_T ScenarioTree::OutputSTOCH( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt, VerbLevel Verbosity )
const
{
	//--------------------------------------------------------------------------
	//	STOCH file: output the header line
	//
	if( Verbosity >= V_HIGH )
		Print( "\nWriting output STOCH file.\n" );

	if( fprintf( fp, "STOCH         %s\nBLOCKS        DISCRETE\n",
		lp.GetName() ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Loop on scenarios belonging to the first period which are expanded in
	//	the second stage problem.
	//
	for( Int_T i = 0, scen = -1; i < ScenNum; )
	{
		scen = Structure[i][0][CutAt].scen;

		const Real_T prob = TotalProbability( scen, CutAt, PeriodNum );
		const Int_T s0 = i;
		for( ; i < ScenNum && Structure[i][0][CutAt].scen == scen; i++ )
			;

		if( fprintf( fp, " BL BLOCK     TIME1           %f\n", prob ) == EOF )
			return False;

		OutputStochVectors( fp, CutAt, s0, i );
		OutputStochMatrix( fp, lp, periods, CutAt, s0, i );
	}

	//--------------------------------------------------------------------------
	//	STOCH file: output ENDATA
	//
	if( fprintf( fp, "ENDATA\n" ) == EOF )
		return False;

	return True;
}


/*------------------------------------------------------------------------------

	void ScenarioTree::ApplyScenario( SolvableLP &lp, Int_T scen, Int_T period )
	const

PURPOSE:
	As the deterministic or two-stage equivalent problems are written, parts of
the input core file need to be modified according to the data stored in the
scenario file. This routine applies the necessary transformations.

PARAMETERS:
	SolvableLP &lp
		Linear problem (read from the core file).

	Int_T scen, Int_T period
		Scenario and period numbers (coordinates of the selected stochastic data
		block).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ScenarioTree::ApplyScenario( SolvableLP &lp, Int_T scen, Int_T period )
const
{
	Int_T block = Structure[scen][0][period].block;

	assert( block >= -1 && block < BlockNum );

	if( block == -1 ) return;	// No stochastic data block associated.

	ExtendArray<ChangeLP *> *CLP_Array = StochasticData[block];

	assert( CLP_Array != NULL );

	for( Int_T i = 0, len = (Int_T)CLP_Array->Len(); i < len; i++ )
	{
		assert( (*CLP_Array)[i] != NULL );

		CLP_Array[0][i]->Apply( lp );
	}
}


/*------------------------------------------------------------------------------

	Bool_T ScenarioTree::OutputROWS( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt ) const
	Bool_T ScenarioTree::OutputCOLUMNS( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt ) const
	Bool_T ScenarioTree::OutputRHS( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt ) const
	Bool_T ScenarioTree::OutputBOUNDS( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T CutAt ) const

PURPOSE:
	These functions are called by the "ConvertToDeterministic()" function in
order to output the data to a file. "ConvertToDeterministic()" designs the
problem structure, while the above functions output whole blocks of the linear
problem.

PARAMETERS:
	FILE *fp
		Output file.

	SolvableLP &lp
		Core file data.

	const Periods &periods
		Time periods.

	Int_T CutAt
		At which period the problem is to be divided. The argument may range
		from 1 to "PeriodNum" (inclusive).

RETURN VALUE:
	Boolean success status ("False" only on I/O errors).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ScenarioTree::OutputROWS( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt  )
const
{
	//--------------------------------------------------------------------------
	//	The section indicator line and the objective row goes first.
	//
	if( ( fprintf( fp, "ROWS\n" ) == EOF ) ||
		( fprintf( fp, " %c  %-8s\n", 'N', "COST" ) == EOF ) )
		return False;

	//--------------------------------------------------------------------------
	//	Then the remaining rows.
	//
	Int_T ScenMax = ScenNum;

	for( Int_T p = 0; p < PeriodNum; p++ )
	{
		if( p == CutAt ) ScenMax = HowManyNodeDuplicates( CutAt );

		for( Int_T i = 0, scen = -1; i < ScenMax; i++ )
		{
			if( Structure[i][0][p].scen == scen ) continue;
			scen = Structure[i][0][p].scen;

			Int_T start, end;

			periods.RowRange( p, start, end );

			for( Int_T j = start; j < end; j++ )
			{
				char c = ' ';

				switch( lp.GetRowType( j ) & RT_TYPE )
				{
				case RT_FR:	c = 'N'; break;
				case RT_EQ:	c = 'E'; break;
				case RT_GE:	c = 'G'; break;
				case RT_LE:	c = 'L'; break;
				}
				if( fprintf( fp, " %c  %-8s\n", c, MakeLabel( j, scen ) )
					== EOF )
					return False;
			}
		}
	}

	return True;
}


Bool_T ScenarioTree::OutputCOLUMNS( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt  )
const
{
	//--------------------------------------------------------------------------
	//	The section indicator line goes first.
	//
	if( fprintf( fp, "%s\n", "COLUMNS" ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	First handle the 0'th (deterministic) period columns. Those columns are
	//	to be output all at once.
	//
	if( !OutputCOLUMNSFromSubTree( fp, lp, periods, 0, PeriodNum, 0, ScenNum,
		CutAt ) )
		return False;

	//--------------------------------------------------------------------------
	//	Now the columns which will belong to the first stage problem.
	//
	Int_T p;
	for( p = 1; p < CutAt; p++ )
		for( Int_T i = 0, scen = -1; i < ScenNum; )
		{
			scen = Structure[i][0][p].scen;

			const Int_T sStart = i;
			for( ; i < ScenNum && Structure[i][0][p].scen == scen; i++ )
				;

			if( !OutputCOLUMNSFromSubTree( fp, lp, periods, p, CutAt, sStart, i,
				CutAt ) )
				return False;
		}

	//--------------------------------------------------------------------------
	//	Now we output the subtree corresponding to the second stage problem.
	//
	const Int_T ScenMax = HowManyNodeDuplicates( CutAt );

	for( ; p < PeriodNum; p++ )
		for( Int_T i = 0, scen = -1; i < ScenMax; )
		{
			scen = Structure[i][0][p].scen;

			const Int_T sStart = i;
			for( ; i < ScenMax && Structure[i][0][p].scen == scen; i++ )
				;

			if( !OutputCOLUMNSFromSubTree( fp, lp, periods, p, PeriodNum,
				sStart, i, CutAt ) )
				return False;
		}

	return True;
}


Bool_T ScenarioTree::OutputCOLUMNSFromSubTree( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T p0, Int_T p1, Int_T s0, Int_T s1,
	Int_T CutAt )
const
{
	Int_T colStart, colEnd;

	periods.ColumnRange( p0, colStart, colEnd );

	//--------------------------------------------------------------------------
	//	For cost output and evaluation.
	//
	Int_T sc			= Structure[s0][0][p0].scen;
	Real_T probability	= TotalProbability( sc, p0, CutAt );

#ifndef NDEBUG
	{
		for( Int_T s = s0 + 1; s < s1; s++ )
			assert( Structure[s][0][p0].scen == sc );
	}
#endif

	//--------------------------------------------------------------------------
	//	Loop on columns of the period, which is to be output now.
	//
	for( Int_T j = colStart; j < colEnd; j++ )
	{
		char lab[9];
		strcpy( lab, MakeLabel( j, sc ) );

		//----------------------------------------------------------------------
		//	Output cost. It is always the cost of the variables of period 'p0'
		//	at the scenario 'sc' (it is assumed that there is only one scenario
		//	in range 's0-'s1' at period 'p0').
		//
		ApplyScenario( lp, sc, p0 );
		if( IsNonZero( lp.GetC( j ) ) )
			if( fprintf( fp, "    %-8s  %-8s  %12G\n", lab, "COST",
				probability * lp.GetC( j ) ) == EOF )
				return False;

		//----------------------------------------------------------------------
		//	Now output the remaining column data (period by period).
		//
		Int_T sMax = s1;

		for( Int_T p = p0; p < p1; p++ )
		{
			Int_T rowStart, rowEnd;
			periods.RowRange( p, rowStart, rowEnd );

			if( p == CutAt ) sMax = HowManyNodeDuplicates( CutAt );

			for( Int_T i = s0, scen = -1; i < sMax; i++ )
			{
				if( Structure[i][0][p].scen == scen ) continue;
				scen = Structure[i][0][p].scen;

				ApplyScenario( lp, scen, p );

				Ptr<Real_T> a;
				Ptr<Int_T> row;
				Int_T len;

				lp.MPS_LP::GetColumn( j, a, row, len );
				for( ; len; --len, ++a, ++row )
					if( *row >= rowStart && *row < rowEnd )
						if( fprintf( fp, "    %-8s  %-8s  %12G\n", lab,
							MakeLabel( *row, scen ), *a ) == EOF )
							return False;
			}
		}
	}

	return True;
}


Bool_T ScenarioTree::OutputRHS( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt  )
const
{
	//--------------------------------------------------------------------------
	//	The section indicator line goes first.
	//
	if( fprintf( fp, "%s\n", "RHS" ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Then the right hand sides of the whole problem.
	//
	Int_T ScenMax = ScenNum;

	for( Int_T p = 0; p < PeriodNum; p++ )
	{
		if( p == CutAt ) ScenMax = HowManyNodeDuplicates( CutAt );

		for( Int_T i = 0, scen = -1; i < ScenMax; i++ )
		{
			if( Structure[i][0][p].scen == scen ) continue;
			scen = Structure[i][0][p].scen;
			ApplyScenario( lp, scen, p );

			Int_T start, end;

			periods.RowRange( p, start, end );

			for( Int_T j = start; j < end; j++ )
			{
				Real_T b = lp.GetB( j );

				if( IsNonZero( b ) )
					if( fprintf( fp, "    %-8s  %-8s  %12G\n",
						"RHS", MakeLabel( j, scen ), b ) == EOF )
						return False;
			}
		}
	}

	return True;
}


Bool_T ScenarioTree::OutputBOUNDS( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T CutAt  )
const
{
	//--------------------------------------------------------------------------
	//	The section indicator line goes first.
	//
	if( fprintf( fp, "BOUNDS\n" ) == EOF )
		return False;

	//--------------------------------------------------------------------------
	//	Then the bounds for the consecutive columns.
	//	Entries created only when necessary, i.e. for fixed, bounded and free
	//	variables or for variables with non-zero lower bound.
	//

	//
	//	Outer loop on periods.
	//
	Int_T ScenMax = ScenNum;

	for( Int_T p = 0; p < PeriodNum; p++ )
	{
		Int_T colStart, colEnd;
		periods.ColumnRange( p, colStart, colEnd );

		if( p == CutAt ) ScenMax = HowManyNodeDuplicates( CutAt );

		//----------------------------------------------------------------------
		//	A loop on scenarios in the given period.
		//
		for( Int_T i = 0, scen = -1; i < ScenMax; i++ )
		{
			if( Structure[i][0][p].scen == scen ) continue;
			scen = Structure[i][0][p].scen;
			ApplyScenario( lp, scen, p );

			//------------------------------------------------------------------
			//	A loop on columns in the given scenario/period pair.
			//
			for( Int_T j = colStart; j < colEnd; j++ )
			{
				Short_T vt = lp.GetVarType( j );

				if( vt == VT_LO && IsZero( lp.GetL( j ) ) )
					continue;

				const char *lab = MakeLabel( j, scen );

				//
				//	Detect FREE variables.
				//
				if( !( vt & ( VT_LO | VT_UP ) ) )
				{
					if( fprintf( fp, " %2s %-8s  %-8s\n", "FR", "BND", lab )
						== EOF )
						return False;
				}
				//
				//	Detect FIXED variables.
				//
				else if( vt == VT_FIXED )
				{
					if( fprintf( fp, " %2s %-8s  %-8s  %12G\n", "FX", "BND",
						lab, lp.GetL( j ) ) == EOF )
						return False;
				}
				//
				//	Detect non-positive (MI) variables.
				//
				else if( ( vt & VT_UP ) && !( vt & VT_LO ) )
				{
					if( fprintf( fp, " %2s %-8s  %-8s\n", "MI", "BND", lab )
						== EOF )
						return False;
					if( IsNonZero( lp.GetU( j ) ) )
						if( fprintf( fp, " %2s %-8s  %-8s  %12G\n", "UP", "BND",
							lab, lp.GetU( j ) ) == EOF )
							return False;
				}
				//--------------------------------------------------------------
				//	For all other variables print out non-zero lower bounds and
				//	finite upper bounds.
				//
				else
				{
					if( vt & VT_LO && IsNonZero( lp.GetL( j ) ) )
						if( fprintf( fp, " %2s %-8s  %-8s  %12G\n", "LO", "BND",
							lab, lp.GetL( j ) ) == EOF )
							return False;

					if( vt & VT_UP )
						if( fprintf( fp, " %2s %-8s  %-8s  %12G\n", "UP", "BND",
							lab, lp.GetU( j ) ) == EOF )
							return False;
				}
			} // End of loop on columns in the given scenario/period pair.

		} // End of loop on scenarios.
	} // End of loop on periods

	return True;
}


Bool_T ScenarioTree::OutputStochVectors( FILE *fp, Int_T p0, Int_T s0, // )
	Int_T s1 )
const
{
	const Int_T mod = s1 - s0;

	for( Int_T p = p0; p < PeriodNum; p++ )
		for( Int_T i = s0, scen = -1; i < s1; i++ )
		{
			if( Structure[i][0][p].scen == scen ) continue;
			scen = Structure[i][0][p].scen;

			const Real_T prob = TotalProbability( scen, p, p0 );
			ExtendArray<ChangeLP *> *clp = NULL;
			
			clp = StochasticData[ Structure[i][0][p].block ];

			if( clp == NULL ) continue;

			//------------------------------------------------------------------
			//	Now output the random RHS and objective coefficients.
			//
			for( Int_T j = 0, len = clp->Len(); j < len; j++ )
				switch( (*clp)[j]->GetType() )
				{
				case ChangeLP::RHS:
					(*clp)[j]->Write( fp, scen % mod, prob );
					break;

				case ChangeLP::COST:
					(*clp)[j]->Write( fp, scen % mod, prob );
					break;

				case ChangeLP::MATRIX:
					break;

				default:
					abort();
				}
		}

	return True;
}


Bool_T ScenarioTree::OutputStochMatrix( FILE *fp, SolvableLP &lp, // )
	const Periods &periods, Int_T p0, Int_T s0, Int_T s1 )
const
{
	const Int_T mod = s1 - s0;

	for( Int_T p = 1; p < p0; p++ )
		for( Int_T i = s0, scen = -1; i < s1; i++ )
		{
			if( Structure[i][0][p].scen == scen ) continue;
			scen = Structure[i][0][p].scen;
			ApplyScenario( lp, scen, p );

			Int_T ColStart, ColEnd, RowStart;

			periods.ColumnRange( p, ColStart, ColEnd );
			RowStart = periods.FirstRow( p0 );

			for( Int_T j = ColStart; j < ColEnd; j++ )
			{
				Ptr<Real_T> a;
				Ptr<Int_T> row;
				Int_T len;
				char lab[9];
				
				strcpy( lab, MakeLabel( j, scen ) );

				lp.MPS_LP::GetColumn( j, a, row, len );
				for( ; len; --len, ++a, ++row )
					if( *row >= RowStart )
						if( fprintf( fp, "    %-8s  %-8s  %12G\n", lab,
							MakeLabel( *row, scen % mod ), *a ) == EOF )
							return False;
			}
		}

	return True;
}


/*------------------------------------------------------------------------------

	Real_T ScenarioTree::TotalProbability( Int_T scen, Int_T period,
		Int_T CutAt ) const

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

Real_T ScenarioTree::TotalProbability( Int_T scen, Int_T period, Int_T CutAt )
const
{
	assert( scen >= 0 && scen < ScenNum );
	assert( period >= 0 && period < PeriodNum );

	Real_T prob = 1.0;

	for( Int_T p = Int_T( (period < CutAt) ? 1 : CutAt+1); p <= period; p++ )
		prob *= Structure[scen][0][p].prob;

	assert( prob >= 0.0 && prob <= 1.0 );
	return prob;
}


Int_T ScenarioTree::HowManyNodeDuplicates( Int_T period )
const
{
	assert( period >= 0 && period <= PeriodNum );

	if( period == PeriodNum ) return 0;

	Int_T scen = Structure[0][0][period].scen,
		i;

	for( i = 1; i < ScenNum && Structure[i][0][period].scen == scen; i++ )
		;

	return i;
}
