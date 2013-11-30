/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	scenario.cpp
CREATED:			1994.07.27
LAST MODIFIED:		1996.02.22

DEPENDENCIES:		stdtype.h, smartptr.h, scenario.h, error.h, std_math.h,
					lexer.h
					<assert.h>, <stdlib.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __DETERMLP_H__
#	include "determlp.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PARSSTOC_H__
#	include "parsstoc.h"
#endif
#ifndef __LEXER_H__
#	include "lexer.h"
#endif
#ifndef __RAND01_H__
#	include "rand01.h"
#endif


#ifdef COMPILE_DISTANCE_MEASUREMENT
//------------------------------------------------------------------------------
//	Static data.
//
Scenario::DistanceMode Scenario::DistMode	= Scenario::NONE;
Long_T Scenario::DistanceCalculationCnt		= 0;
//
//	End of the static data.
//------------------------------------------------------------------------------
#endif


//@BEGIN--------------------------
extern bool sadsam; 
//@END----------------------------



/*------------------------------------------------------------------------------

	IndepStochVarBlock::~IndepStochVarBlock( void )

PURPOSE:
	Destructor. Frees memory occupied by the stochastic variable block. A
block consists of a number of Delta structures.

PARAMETERS:
	None.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

IndepStochVarBlock::~IndepStochVarBlock( void )
{
	Int_T i;
	for( i = 0; i < len; i++ )
	{
		assert( block[i] != NULL );
		delete block[i];
	}

#ifndef NDEBUG
	for( ; i < maxLen; i++ )
		assert( block[i] == NULL );
#endif
}


/*------------------------------------------------------------------------------

	void IndepStochVarBlock::Append( Delta *d )

PURPOSE:
	Appends a single object change information to the stochastic bata block
structure.

PARAMETERS:
	Delta *d
		An object change information.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void IndepStochVarBlock::Append( Delta *d )
{
	assert( d != NULL );
	assert( len <= maxLen );

	if( len == maxLen )
	{
		Int_T NewMaxLen = Max( Int_T( 3 * maxLen / 2 ), Int_T( maxLen + 10 ) );

		block.Resize( NewMaxLen );
		block.Fill( (Delta *)NULL, NewMaxLen, maxLen );
		maxLen = NewMaxLen;
	}

	block[len++] = d;
}


/*------------------------------------------------------------------------------

	Distribution::~Distribution( void )

PURPOSE:
	Destructor of an object representing a single distribution.

PARAMETERS:
	None.

RETURN VALUE:
	Non applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Distribution::~Distribution( void )
{
	for( Int_T i = 0; i < len; i++ )
	{
		assert( block[i] != NULL );
		delete block[i];
	}
}


/*------------------------------------------------------------------------------

	void Distribution::Append( StochasticDataBlock *bl )

PURPOSE:
	Adds a new block to a stochastic data distribution stucture.

PARAMETERS:
	StochasticDataBlock *bl
		A block to be appended to the distribution.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Distribution::Append( StochasticDataBlock *bl )
{
	assert( bl != NULL );
	assert( len <= maxLen );

	if( len == maxLen )
	{
		Int_T NewMaxLen = Max( Int_T( 3 * maxLen / 2 ), Int_T( maxLen + 10 ) );

		block.Resize( NewMaxLen );
		block.Fill( (StochasticDataBlock *)NULL, NewMaxLen, maxLen );
		maxLen = NewMaxLen;
	}

	block[len++] = bl;
}


/*------------------------------------------------------------------------------

	Scenarios::Scenarios( void )

PURPOSE:
	Scenario repository constructor. Initializes an empty object.

PARAMETERS:
	None.

RETURN VALUE:
	Non applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Scenarios::Scenarios( void )
	: ScenNum( 0 ), ArrayOfScenarios(), status( EMPTY ),

	RowLabels( NULL ),		ColumnLabels( NULL ),
	ProblemName( NULL ),	RHS_Name( NULL ),
	ObjName( NULL ),		BoundsName( NULL ),
	RangesName( NULL ),
	ProblemNameLen( -1 ),	RHS_NameLen( -1 ),
	ObjNameLen( -1 ),		BoundsNameLen( -1 ),
	RangesNameLen( -1 ),
	SemanticError( False ),

	len( 0 ), maxLen( INIT_DIST_NUM ),
	dist( maxLen, (Distribution *)NULL ),
	MaxScen( 0.0 ),

	LastStochBlockType( NONE ), LastIndepDiscrete(),
	LastBlock( NULL )
{
	LastBlockName[0] = '\0';
}


/*------------------------------------------------------------------------------

	void Scenarios::Cleanup( void )
	void Scenarios::CleanupDistributions( void )

PURPOSE:
	The first function removes all data of the scenario repository object. It is
called by the destructor. The other one removes just the distributions. It is
called by the "Cleanup".

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	The scenario object becomes empty after the first function is called.

------------------------------------------------------------------------------*/

void Scenarios::Cleanup( void )
{
	if( status == EMPTY ) return;
	status = EMPTY;

	//--------------------------------------------------------------------------
	//	Clear all data regarding the deterministic linear problem.
	//
	RowLabels		= ColumnLabels	= NULL;
	ProblemName		= RHS_Name		= ObjName		= BoundsName	=
		RangesName		= NULL;
	ProblemNameLen	= RHS_NameLen	= ObjNameLen	= BoundsNameLen	=
		RangesNameLen	= -1;
	SemanticError	= False;

	//--------------------------------------------------------------------------
	//	Free the array of scenarios.
	//
	for( Int_T i = 0; i < ScenNum; i++ )
		if( ArrayOfScenarios[i] != NULL ) delete ArrayOfScenarios[i];
	ArrayOfScenarios.Resize( 0 );

	//--------------------------------------------------------------------------
	//	Re-initialize data describing the dimension of scenario space.
	//
	ScenNum = 0;

	CleanupDistributions();
}


void Scenarios::CleanupDistributions( void )
{
	for( Int_T i = 0; i < len; i++ )
	{
		assert( dist[i] != NULL );
		delete dist[i];
	}

	dist.Resize( 0 );
	len = maxLen = 0;
}


/*------------------------------------------------------------------------------

	void Scenarios::AppendNewDistribution( void )

PURPOSE:
	Appends a new empty distribution to the distribution array "dist" of the
"Scenarios" object. This function is called during reading of the distributions.
After it's call, the distribution is filled with stochastic data read from the
file.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenarios::AppendNewDistribution( void )
{
	assert( len <= maxLen );

	if( len == maxLen )
	{
		Int_T NewMaxLen = Max( Int_T( 3 * maxLen / 2 ), Int_T( maxLen + 10 ) );

		dist.Resize( NewMaxLen );
		dist.Fill( (Distribution *)NULL, NewMaxLen, maxLen );
		maxLen = NewMaxLen;
	}

	Distribution *d = new Distribution;

	if( d == NULL ) FatalError( "Out of memory." );
	dist[len++] = d;
}


/*------------------------------------------------------------------------------

	void Scenarios::AppendLastBlock( void )

PURPOSE:
	This function append the last stochastic data block from the BLOCKS
DISCRETE section to the list of blocks in the current distribution.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenarios::AppendLastBlock( void )
{
	if( LastBlock == NULL ) return;

	assert( dist[len-1] != NULL );

	dist[len-1]->Append( LastBlock );
	LastBlock = NULL;
}


/*------------------------------------------------------------------------------

	Bool_T Scenarios::ReadScenarioFile( const char *StochFileName,
		FILE *StochFilePtr, DeterministicLP &DeterministicMatrix )

PURPOSE:
	This function reads and parses a scenario file. After a discrete
distribution of every random data object is read, it is stored in the
"dist" array.

PARAMETERS:
	const char *StochFileName, FILE *StochFilePtr,
		Scenario file name and file pointer.
	
	DeterministicLP &DeterministicMatrix
		Deterministic linear problem.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	File pointer is moved until a parsing error is encountered or "ENDATA"
statement is found.

------------------------------------------------------------------------------*/

Bool_T Scenarios::ReadScenarioFile( const char *StochFileName, // )
	FILE *StochFilePtr, DeterministicLP &DeterministicMatrix )
{
	status = READING_FILE;

	//--------------------------------------------------------------------------
	//	Initialize object data concerning the deterministic linear problem.
	//
	RowLabels = &( DeterministicMatrix.RevealRowLabels() );
	ColumnLabels = &( DeterministicMatrix.RevealColumnLabels() );

	ProblemName	= DeterministicMatrix.RevealProblemName();
	RHS_Name	= DeterministicMatrix.RevealRHS_Name();
	ObjName		= DeterministicMatrix.RevealObjectiveName();
	BoundsName	= DeterministicMatrix.RevealBoundsName();
	RangesName	= DeterministicMatrix.RevealRangesName();

	assert( ProblemName != NULL );
	assert( RHS_Name != NULL );
	assert( ObjName != NULL );
	assert( BoundsName != NULL );
	assert( RangesName != NULL );

	ProblemNameLen	= strlen( ProblemName );
	RHS_NameLen		= strlen( RHS_Name );
	ObjNameLen		= strlen( ObjName );
	BoundsNameLen	= strlen( BoundsName );
	RangesNameLen	= strlen( RangesName );

	assert( ProblemNameLen >= 0 );
	assert( RHS_NameLen >= 0 );
	assert( BoundsNameLen >= 0 );
	assert( ObjNameLen >= 0 );
	assert( RangesNameLen >= 0 );

	SemanticError = False;

	//--------------------------------------------------------------------------
	//	Read in the stoch file.
	//
	if( !GetStochFile( StochFileName, StochFilePtr, this ) || SemanticError )
	{
#ifndef NDEBUG
		CheckArrayOfDistributions();
#endif
		Cleanup();
		status = EMPTY;
		return False;
	}
	else
	{
		status = DISTRIBUTIONS_READY;
		return True;
	}
}


//==============================================================================
//
//  ##  #### #   #  ##  #   # ##### #  ##      ##   ##  ##### #  ##  #   #  ##
// #    #    ## ## #  # ##  #   #   # #  #    #  # #  #   #   # #  # ##  # # 
//  ##  ###  # # # #  # # # #   #   # #       #  # #      #   # #  # # # #  ##
//    # #    #   # #### # # #   #   # #       #### #      #   # #  # # # #    #
// #  # #    #   # #  # #  ##   #   # #  #    #  # #  #   #   # #  # #  ## #  #
//  ##  #### #   # #  # #   #   #   #  ##     #  #  ##    #   #  ##  #   #  ##
//
//==============================================================================



/*------------------------------------------------------------------------------

	void Scenarios::NewIndepDiscreteEntry( const char *Row, const char *Col,
		const char *PeriodLabel, double Val, double Prob )
	void Scenarios::NewBlocksDiscrete( const char *BlockName,
		const char *PeriodLabel, Real_T Prob )
	void Scenarios::NewBlocksDiscreteEntry( const char *Row, const char *Col,
		Real_T Val )
	void Scenarios::DeltaEntry( const char *Row, const char *Col,
		Delta::DeltaType &type, Int_T &row, Int_T &col )
	void Scenarios::EndData( void )

PURPOSE:
	The first function processes data from a single line of the "INDEP DISCRETE"
section of the stochastic data file at a time. It produces a probability
distribution for one stochastic datum.
	The second one is called when the begining of a new block is encountered in
while reading the BLOCKS DISCRETE section of the stochastic data file. The third
one processes one entry of one block at a time.
	The fourth one is used by two of the wprevoiusly mentioned ones to convert
the row/column labels into indice in the deterministic LP.
	The last function is called when the ENDATA statement is found during file
reading.

PARAMETERS:
	const char *Row, const char *Col
		Row and column labels of the datum. Sometimes either "Row" or "Column"
		may actually be the names of some vectors describing the linear problem.

	const char *PeriodLabel
		Perid label. Ignored (so far at least).

	const char *BlockLabel
		Current block label in the BLOCKS DISCRETE section.

	double Val, double Prob
		Value its probability.

	Delta::DeltaType &type, Int_T &row, Int_T &col
		Reference arguments used to pass back the indice.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenarios::DeltaEntry( const char *Row, const char *Col, // )
	Delta::DeltaType &type, Int_T &row, Int_T &col )
{
	//--------------------------------------------------------------------------
	//	First detect the stochastic entry type.
	//
	row	= col = -1;
	type = Delta::NONE;

	if( strcmp( Col, RHS_Name ) == 0 && (int)strlen( Col ) == RHS_NameLen )
		type = Delta::RHS;
	else if( strcmp( Row, ObjName ) == 0 && (int)strlen( Row ) == ObjNameLen )
		type = Delta::COST;
	else
		type = Delta::MATRIX;

	//--------------------------------------------------------------------------
	//	Then decode the necessary values (row/column indice). Indicate error on
	//	failure.
	//
	if( type == Delta::RHS || type == Delta::MATRIX )
		if( ( row = (Int_T) RowLabels->FindLabel( Row ) ) < 0 )
		{
			Error( "File %s, line %d: Row label %s not found in the 2nd stage "
				"problem.", Lexer::FileName, Lexer::LineNumber, Row );
			SemanticError = True;
		}

	if( type == Delta::COST || type == Delta::MATRIX )
		if( ( col = (Int_T) ColumnLabels->FindLabel( Col ) ) < 0 )
		{
			Error( "File %s, line %d: Column label %s not found in the "
				"2nd stage problem.", Lexer::FileName, Lexer::LineNumber,
				Col );
			SemanticError = True;
		}
}


void Scenarios::NewIndepDiscreteEntry( const char *Row, const char *Col, // )
	const char *, double Value, double Prob )
{
	//--------------------------------------------------------------------------
	//	If this is the beginning of a new section, make sure the last block
	//	is stored on the distribution list.
	//
	if( LastBlock != NULL ) AppendLastBlock();

	Int_T row, col;
	Delta::DeltaType type;
	DeltaEntry( Row, Col, type, row, col );

	//--------------------------------------------------------------------------
	//	See if the probability makes (some) sense.
	//
	if( Prob < 0.0 || Prob > 1.0 )
	{
		Error( "File %s, line %d: Invalid value %f for a probability.",
				Lexer::FileName, Lexer::LineNumber, Prob );
		SemanticError = True;
		Prob = 0.0;
	}

	//--------------------------------------------------------------------------
	//	Now is the time to decide whether the distribution of the variable
	//	that was just presented is complete or not. If it is - we
	//	generate appropriate entries for the scenarios.
	//

#ifndef NDEBUG
	if( SemanticError )
		assert( type != Delta::NONE );
#endif

	//--------------------------------------------------------------------------
	//	If the same distribution continues, add a new entry to the current
	//	distribution. When a new distribution starts create a new distribution
	//	and store the entry there.
	//
	IndepStochVar *isv = new IndepStochVar( Prob, type, row, col, Value );
	if( isv == NULL ) FatalError( "Out of memory." );

	if( !( LastStochBlockType == INDEP_DISCRETE && row == LastIndepDiscrete.row
		&& col == LastIndepDiscrete.col && type == LastIndepDiscrete.type ) )
	{
		LastStochBlockType		= INDEP_DISCRETE;
		LastIndepDiscrete.row	= row;
		LastIndepDiscrete.col	= col;
		LastIndepDiscrete.type	= type;

		AppendNewDistribution();
	}

	dist[len-1]->Append( isv );
}


void Scenarios::NewBlocksDiscrete( const char *BlockName, // )
	const char * /* PeriodLabel */, Real_T Prob )
{
	AppendLastBlock();
	LastStochBlockType = BLOCKS_DISCRETE;

	//--------------------------------------------------------------------------
	//	See if the probability makes (some) sense.
	//
	if( Prob < 0.0 || Prob > 1.0 )
	{
		Error( "File %s, line %d: Invalid value %f for a probability.",
				Lexer::FileName, Lexer::LineNumber, Prob );
		SemanticError = True;
		Prob = 0.0;
	}

	//--------------------------------------------------------------------------
	//	A new name starts a new distribution.
	//
	if( strcmp( LastBlockName, BlockName ) )
	{
		strncpy( LastBlockName, BlockName, 20 );
		LastBlockName[19] = '\0';

		AppendNewDistribution();
	}

	//--------------------------------------------------------------------------
	//	Create a new block in the current distribution (whichever this may be).
	//
	LastBlock = new IndepStochVarBlock( Prob );
}


void Scenarios::NewBlocksDiscreteEntry( const char *Row, const char *Col, // )
	Real_T Val )
{
	Int_T row, col;
	Delta::DeltaType type;
	DeltaEntry( Row, Col, type, row, col );

	Delta *d = new Delta( type, row, col, Val );
	if( d == NULL ) FatalError( "Out of memory." );

	LastBlock->Append( d );
}


void Scenarios::EndData( void )
{
	switch( LastStochBlockType )
	{
	case INDEP_DISCRETE:
		//
		//	There is nothing to do - the INDEP DISCRETE entries are added to
		//	the array of distributions immediately after they were detected.
		//
		LastIndepDiscrete.row	= -1;
		LastIndepDiscrete.col	= -1;
		LastIndepDiscrete.value	= 0.0;
		break;

	case BLOCKS_DISCRETE:
		//
		//	The last block has to be added to the array of distributions as
		//	it was being collected before the ENDATA was found.
		//
		AppendLastBlock();
		LastBlockName[0] = '\0';
		break;

	case NONE:
		break;
	}

	LastStochBlockType  = NONE;
}


//==============================================================================
//
//  M   M  OOO  RRRR  EEEE     FFFF U   U N   N  CCC  TTTTT I  OOO  N   N  SSS
//  MM MM O   O R   R E        F    U   U NN  N C   C   T   I O   O NN  N S
//  M M M O   O RRRR  EEE      FFF  U   U N N N C       T   I O   O N N N  SSS
//  M   M O   O R R   E        F    U   U N  NN C       T   I O   O N  NN     S
//  M   M O   O R  R  E        F    U   U N   N C   C   T   I O   O N   N S   S
//  M   M  OOO  R   R EEEE     F     UUU  N   N  CCC    T   I  OOO  N   N  SSS
//
//==============================================================================


#ifdef COMPILE_DISTANCE_MEASUREMENT

/*------------------------------------------------------------------------------

	void SetDistanceMode( Scenario::DistanceMode dm )

	Real_T GetDistance( const Scenario &s1, const Scenario &s2 )

	Real_T GetStdNormDistance( const Scenario &s1, const Scenario &s2 )
	Real_T GetWeightedNormDistance( const Scenario &s1, const Scenario &s2 )
	Real_T GetDiscreteMeasureDistance( const Scenario &s1, const Scenario &s2 )

PURPOSE:
	Measures a "distance" between two scenarios. The actual type of distance
measurement depends on a static member datum of the "Scenario" class
"Scenario::DistMode". The "distance" is supposed to be a measure of
dissimilarity: the smaller the "distance" the more similar the scenarios.

PARAMETERS:
	Scenario::DistanceMode dm
		Distance measurement method.

	const Scenario &s1, const Scenario &s2
		The scanerios to be compared.

RETURN VALUE:
	A nonnegative number.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenario::SetDistanceMode( Scenario::DistanceMode dm )
{
	assert( dm != NONE );
	DistMode = dm;
}

// @BEGIN Rebecca's code to get samples
Real_T Scenario::GetScenarioComponents(const Scenario &s1)
{
	int value = 0;
	char scenarios[30] = {0};	
	char temp[30];
		
	for( Int_T l = s1.GetLength(), i = 0; i < l; i++)
	{
		const Delta &d1	= s1[i][0];
		value = (int) d1.value;
		sprintf(temp, "%d ", value);
		strcat(scenarios, temp);
	}
	Print(scenarios);
	
	return 0;
}
// @END Rebecca's code

Real_T Scenario::GetDistance( const Scenario &s1, const Scenario &s2 )
{
	DistanceCalculationCnt++;

	switch( DistMode )
	{
	case NONE:
#ifndef NDEBUG
	abort();
#endif
	case STD_NORM:			return GetStdNormDistance( s1, s2 );
	case WEIGHTED_NORM:		return GetWeightedNormDistance( s1, s2 );
	case DISCRETE_MEASURE:	return GetDiscreteMeasureDistance( s1, s2 );
	default:
#ifndef NDEBUG
		abort();
#endif
		return 0.0;
	}
}


Real_T Scenario::GetStdNormDistance( const Scenario &s1, const Scenario &s2 )
{
	assert( DistMode == STD_NORM );
	assert( s1.GetLength() == s2.GetLength() );

	Real_T norm = 0.0;

	for( Int_T l = s1.GetLength(), i = 0; i < l; i++ )
	{
		assert( s1[i].Len() == s2[i].Len() );

		for( Int_T j = 0, bl = s1[i].Len(); j < bl; j++ )
		{
			const Delta &d1	= s1[i][j],
				&d2			= s2[i][j];

			assert( d1.type == d2.type );
			assert( d1.row == d2.row );
			assert( d1.col == d2.col );

			norm += fabs( d1.value - d2.value );
		}
	}

	return norm;
}


Real_T Scenario::GetWeightedNormDistance( const Scenario &s1, // )
	const Scenario &s2 )
{
	assert( DistMode == WEIGHTED_NORM );
	assert( s1.GetLength() == s2.GetLength() );

	Real_T norm = 0.0;
	const Real_T factor = 2.0;

	for( Int_T l = s1.GetLength(), i = 0; i < l; i++ )
	{
		assert( s1[i].Len() == s2[i].Len() );

		for( Int_T j = 0, bl = s1[i].Len(); j < bl; j++ )
		{
			const Delta &d1	= s1[i][j],
				&d2			= s2[i][j];

			assert( d1.type == d2.type );
			assert( d1.row == d2.row );
			assert( d1.col == d2.col );

			norm += ( ( d1.type == Delta::COST ) ? 1.0 : factor ) *
					fabs( d1.value - d2.value );
		}
	}

	return norm;
}


Real_T Scenario::GetDiscreteMeasureDistance( const Scenario &s1, // )
	const Scenario &s2 )
{
	assert( DistMode == DISCRETE_MEASURE );
	assert( s1.GetLength() == s2.GetLength() );

	Real_T norm = 0.0;

	for( Int_T l = s1.GetLength(), i = 0; i < l; i++ )
	{
		assert( s1[i].Len() == s2[i].Len() );

		for( Int_T j = 0, bl = s1[i].Len(); j < bl; j++ )
		{
			const Delta &d1	= s1[i][j],
				&d2			= s2[i][j];

			assert( d1.type == d2.type );
			assert( d1.row == d2.row );
			assert( d1.col == d2.col );

			if( IsEqual( d1.value, d2.value ) ) norm += 1.0;
		}
	}

	return (Real_T)norm;
}


/*------------------------------------------------------------------------------

	Real_T Scenarios::GetDistance( Int_T s1, Int_T s2 )

PURPOSE:
	Duplicates the functionality of the "Scenario::GetDistance" function, but
uses scenario numbers on the scenario list rather than the references to the
particular scenarios.

PARAMETERS:
	Int_T s1, Int_T s2
		Numbers of scenarios to be compared.

RETURN VALUE:
	The distance (see the description of "Scenario::GetDistance()" for details.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T Scenarios::GetDistance( Int_T s1, Int_T s2 )
{
	assert( s1 >= 0 && s1 < ScenNum );
	assert( s2 >= 0 && s2 < ScenNum );

	//--------------------------------------------------------------------------
	//	A distance between a scenario and the same scenario is always zero.
	//	Non-trivial distance computation may be performed.
	//
	return ( s1 == s2 ) ? 0.0 : Scenario::GetDistance( *ArrayOfScenarios[s1],
		*ArrayOfScenarios[s2] );
}

// @BEGIN Rebecca's code to get sample data
Real_T Scenarios::GetScenarioComponents( Int_T s1 )
{
	return Scenario::GetScenarioComponents( *ArrayOfScenarios[s1] );

}
// @END Rebecca's code
#endif


/*------------------------------------------------------------------------------

	void Scenarios::CheckArrayOfDistributions( void )

PURPOSE:
	Checks the existence and some basic properties of the distributions
collected during the scenario file reading. Compiled only when self debugging
is turned on (macro "NDEBUG" not defined).

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

#ifndef NDEBUG

void Scenarios::CheckArrayOfDistributions( void )
{
	assert( status == DISTRIBUTIONS_READY );

	//--------------------------------------------------------------------------
	//	Loop on all distributions.
	//
	for( Int_T d = 0; d < len; d++ )
	{
		//
		//	Loop on blocks of the distribution.
		//
		for( Int_T i = 0, l = dist[d]->Len(); i < l; i++ )
		{
			StochasticDataBlock &sdb = (*dist[d])[i];

			for( Int_T j = 0, bl = sdb.Len(); j < bl; j++ )
			{
				assert( sdb[j].type != Delta::NONE );
				assert( sdb[j].col >= 0 );
				assert( sdb[j].row >= 0 );
			}
		}
	}
}

#endif


/*------------------------------------------------------------------------------

	Bool_T Scenarios::GenerateScenarios( Int_T &num )

PURPOSE:
	This routine uses the distributions to create scenarios. The number of
scenarios is given by the "num" argument. If it is negative, it is assumed that
we ought to generate all possible scenarios. If it is positive, then it
represents exactly the number of requested scenarios.
	The total number of scenarios to be generated may not exceed a preset
maximum (MAX_SCEN_NUM, defined in "scenario.h"). If it exceeds the total number
of possible scenarios, it is adjusted as well.
	Either all possible scenarios are generated, or just a sample of them.

PARAMETERS:
	Int_T &num
		Nober of scenarios (see above).

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Scenarios::GenerateScenarios( Int_T &num )
{
	
	//I am commenting out the first assertion as it is much easier 
	//this way -- GB. 
	//assert( status == DISTRIBUTIONS_READY );



	//--------------------------------------------------------------------------
	//	Calculate the actual number of scenarios to be generated. This number
	//	should not exceed
	//	a)	the number of all possible scenarios or
	//	b)	the scenario number limit MAX_SCEN_NUM.
	//
	//	"num == -1" means that we ought to generate all scenarios. In such case
	//	"num" is set to Min( MaxScen, MAX_SCEN_NUM );
	//
	MaxScen = GetNumberOfPossibleScenarios();
	Print( "\n\tThe total number of scenarios is %G.\n", MaxScen );

	Bool_T GenerateAll = False;

	/*if( num > MaxScen )
	{
		assert( (Int_T) MaxScen > 0 );
		Print( "\tSpecified number of scenarios exceeds the number of all "
			"scenarios.\n\tAdjusting to %d.\n", int( num = (Int_T) MaxScen ) );
		GenerateAll = True;
	}
        */
	if( num > MAX_SCEN_NUM )
	{
		Print( "\tSpecified number of scenarios exceeds the maximum number "
			"allowed.\n\tAdjusting to %d.\n",
			int( num = (Int_T) MAX_SCEN_NUM ) );
		GenerateAll = False;
	}
	else if( num < 0 )
	{
		if( MaxScen <= (double) MAX_SCEN_NUM )
		{
			assert( ( (Int_T) MaxScen ) > 0 );
			num = Int_T( MaxScen );
			GenerateAll = True;
		}
		else
		{
			Print( "\tToo many possible scenarios. Sampling only %d.\n",
				int( MAX_SCEN_NUM ) );
			num = (Int_T) MAX_SCEN_NUM;
			GenerateAll = False;
		}
	}
	assert( num > 0 && num <= MAX_SCEN_NUM );
	ScenNum = num;

	//@BEGIN-------------------------------------------------------------
	//If running sadsam, scenarios need to be generated without
	//cumulative probabilities...

	Int_T AdjustedProbabilities; 

	if( sadsam )
		AdjustedProbabilities = CalculateAccumulatedProbabilities(False); 
	else
	//@END---------------------------------------------------------------
		AdjustedProbabilities = CalculateAccumulatedProbabilities(
		GenerateAll ? False : True );

	if( AdjustedProbabilities )
	{
		Error( "%d probabilities had to be adjusted! Check the "
			"distributions.", (int) AdjustedProbabilities );
		ScenNum = 0;
		return False;
	}

	//--------------------------------------------------------------------------
	//	Call the appropriate scenario generation routine (generating all
	//	scenarios is a combinatorial task, while generating a sample requires
	//	a random choice procedure).
	//
	Print( "\tGenerating %d scenarios now...\n", num );

	ArrayOfScenarios.Resize( num );
	ArrayOfScenarios.Fill( NULL, num );

	if( GenerateAll )
		GenerateAllScenarios( num );
	else
		GenerateScenarioSample( num );

	status = SCENARIOS_READY;
	return True;
}


/*------------------------------------------------------------------------------

	Real_T Scenarios::GetNumberOfPossibleScenarios( void ) const

PURPOSE:
	Calculates the total number of all possible distinct scenarios.

PARAMETERS:
	None.

RETURN VALUE:
	The number (as a floating point number). It may easily reach 1E+100!

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T Scenarios::GetNumberOfPossibleScenarios( void )
	const
{
	Real_T max = 1.0;

	//--------------------------------------------------------------------------
	//	Loop on all distributions.
	//	Compute the total number of all possible scenarios.
	//
	for( Int_T d = 0; d < len; d++ )
		max *= dist[d]->Len();
	return max;
}


/*------------------------------------------------------------------------------

	void Scenarios::GenerateAllScenarios( Int_T num )

PURPOSE:
	In cases when the distributions are not to numerous and not too long, we
may attempt to generate all possible scenarios. This procedure is only called
when it is discovered (by the "GenerateScenarios()" procedure), that it
is actually possible.

PARAMETERS:
	Int_T num
		The (prevoiusly calculated) number of scenarios to generate.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenarios::GenerateAllScenarios( Int_T num )
{
	assert( num > 0 && num <= MAX_SCEN_NUM );

	//--------------------------------------------------------------------------
	//	Loop on scenarios. Create one at a time.
	//
	Int_T i = 0;
	for( Int_T seclen; i < num; i++ )
	{
		seclen = num;		//	Section length (interval [0,...,num-1] is
							//	divided into sections).

		//------------------------------------------------------------------
		//	Allocate a new scenario and put it on the scenario list.
		//
		Scenario *scen = new Scenario( len, 0.0 );
		Real_T prob = 1.0;

		if( !scen ) FatalError( "Out of memory." );
		ArrayOfScenarios[i] = scen;

		for( Int_T j = 0; j < len; j++ )	// "len" denotes scenario legth.
		{
			//------------------------------------------------------------------
			//	Add a new stochastic data block to the scenario.
			//
			Int_T distlen = dist[j]->Len();		//	Distribution length.

			seclen /= distlen;
			assert( ( j < Int_T( len - 1 ) && seclen > 1 ) || seclen == 1 );

			StochasticDataBlock *sdb = &( (*dist[j])[ i/seclen%distlen ] );

			scen->Set( j, sdb );

			//------------------------------------------------------------------
			//	Update the probability of the scenario.
			//
			prob *= sdb->GetProbability();
		}

		scen->SetProbability( prob );
	}

	//--------------------------------------------------------------------------
	//	Now simple check of the scenarios: do the probabilities sum up to one.
	//
	Real_T TotalProbability = 0.0;

	for( i = 0; i < ScenNum; i++ )
	{
		assert( ArrayOfScenarios[i] != NULL );
		TotalProbability += ArrayOfScenarios[i]->GetProbability();
	}

	if( TotalProbability - 1.0 > PROBABILITY_TOL )
		Warning( "Scenario generation problem: total probability = %G.\n",
			TotalProbability );
}


/*------------------------------------------------------------------------------

	void Scenarios::GenerateScenarioSample( Int_T num )

PURPOSE:
	When it is impossible or undesireable to generate all possible scenarios,
then this procedure produces a random sample of "num" scenarios.

PARAMETERS:
	Int_T num
		The size of the sample (number of scenarios to be generated).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Scenarios::GenerateScenarioSample( Int_T num )
{
	assert( num > 0 && num <= MAX_SCEN_NUM );

	//--------------------------------------------------------------------------
	//	Loop on all scenarios. Create one at a time.
	//
	for( Int_T s = 0; s < num; s++ )
	{
		Scenario *scen = new Scenario( len, 1.0/(double)num );

		if( !scen ) FatalError( "Out of memory." );
		ArrayOfScenarios[s] = scen;

		//----------------------------------------------------------------------
		//	Loop on distributions. Choose one block from each distribution and
		//	put it in the scenario.
		//
		for( Int_T i = 0; i < len; i++ )
		{
			Real_T prob = Random01::Next();

			Int_T j, l = dist[i]->Len();
			for( j = 0; j < l; j++ )
				if( (*dist[i])[j].GetProbability() > prob ) break;

			assert( j < l );

			scen->Set( i, &( (*dist[i])[j] ) );
		}
	}

	//--------------------------------------------------------------------------
	//	Check if all scenarios have been created.
	//
#ifndef NDEBUG
	{
		for( Int_T s = 0; s < ScenNum; s++ )
			assert( ArrayOfScenarios[s] != NULL );
	}
#endif
}


/*------------------------------------------------------------------------------

	Int_T Scenarios::CalculateAccumulatedProbabilities( Bool_T Store )

PURPOSE:
	During scenario generation it is always needed that a sum of probabilities 
of all possible stochastic data realizations of each random variable equals
one. Additionally, when only a random sample of scenarios is generated, a so 
called accumulated probability is calculated for each random variable
distribution. It allows one random real number in the range [0,1] inclusive to
choose unambigously a random variable realization.

PARAMETERS:
	Bool_T Store
		Says whether or not to store the accumulated probabilities.

RETURN VALUE:
	A non-negative number of probabilities (not distributions) that had to be
adjusted because they did not sum up to one. So far the adjustment is PERFORMED
(and not only determined to be necessary) only when "Store" is "True".

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T Scenarios::CalculateAccumulatedProbabilities( Bool_T Store )
{
	Int_T AdjustedProbabilities	= 0;
	Bool_T AdjustmentWarning	= False;

	//--------------------------------------------------------------------------
	//	Loop on all distributions. Calculate accumulated probabilities,
	//	see if they add up to 1. Compute the total number of all possible
	//	scenarios.
	//
	for( Int_T d = 0; d < len; d++ )
	{
		//----------------------------------------------------------------------
		//	All stochastic data blocks in a distribution should have a total
		//	probability of 1.
		//
		Real_T TotalProb = 0.0;

		//	Loop on blocks in the distribution.
		//
		for( Int_T i = 0, l = dist[d]->Len(); i < l; i++ )
		{
			StochasticDataBlock &sdb = (*dist[d])[i];

			TotalProb += sdb.GetProbability();

			if( IsPositive( TotalProb - 1.0 ) && !AdjustmentWarning )
			{
				Warning( "Some probabilities in distributions do not sum up to "
					"unity." );
				AdjustmentWarning = True;
			}
			else if( TotalProb - 1.0 > PROBABILITY_TOL )
			{
				AdjustedProbabilities++;
				TotalProb = 1.0;
			}

			if( Store ) sdb.SetProbability( TotalProb );
		}
	}
	//	End of loop on distributions.
	//--------------------------------------------------------------------------

	return AdjustedProbabilities;
}


/*------------------------------------------------------------------------------

	void Scenarios::RenumberIndiceInScenarios( Array<Int_T> &NewRowNumber,
		Int_T rLen, Array<Int_T> &NewColNumber, Int_T cLen,
		Int_T ActualStage2Rows, Int_T Stage2Col )

PURPOSE:
	Scenarios are generated when only a deterministic LP formulation is
available. Therefore the row and column numbers stored in "Delta" structures
correspond to the row and column numbers in the deterministic LP.
	But, as we know, we ought to store the indice relative to the second stage
problem. This means, that new indice have to be calculated for the random data
in the right hand side, the cost vector, the technology matrix, the bound
vectors and so on.
	The appropriate permutation tables are calculated by the "DeterministicLP"
object after it is divided into stages. It then calls this function so that
the scenario indice are adjusted.

PARAMETERS:
	Array<Int_T> &NewRowNumber, Int_T rLen
	Array<Int_T> &NewColNumber, Int_T cLen
		'Permutation' tables and their respective lengths. Note that both per-
		mutation tables actually store two permutations: one for the first and
		one for the second stage.

		Table lengths are used only for self-debugging. If "NDEBUG" macro is
		defined, they are ignored.

	Int_T ActualStage2Rows, Int_T Stage2Col
		More data used for self-debugging only. They hold the total number of
		stage two rows and columns, resp.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

#ifndef NDEBUG
void Scenarios::RenumberIndiceInScenarios( Array<Int_T> &NewRowNumber, // )
	Int_T rLen, Array<Int_T> &NewColNumber, Int_T cLen, Int_T ActualStage2Rows,
	Int_T Stage2Col )
#else
void Scenarios::RenumberIndiceInScenarios( Array<Int_T> &NewRowNumber,
	Int_T, Array<Int_T> &NewColNumber, Int_T, Int_T, Int_T )
#endif
{
	//--------------------------------------------------------------------------
	//	Loop on all distributions.
	//
	for( Int_T d = 0; d < len; d++ )
	{
		//	Loop on blocks in distribution "*dist[d]".
		//
		for( Int_T i = 0, l = dist[d]->Len(); i < l; i++ )
		{
			//	Loop on delta's in block "(*dist[d])[i]".
			//
			for( Int_T j = 0, bl = (*dist[d])[i].Len(); j < bl; j++ )
			{
				Delta &delta = (*dist[d])[i][j];

				Int_T &row				= delta.row,
					&col				= delta.col;
				Delta::DeltaType &type	= delta.type;

				switch( type )
				{
				case Delta::RHS:
					assert( row >= 0 && row < rLen );

					row = NewRowNumber[ row ];

					assert( row >= 0 && row < ActualStage2Rows );
					break;

				case Delta::MATRIX:
					assert( col >= 0 && col < cLen );
					assert( row >= 0 && row < rLen );

					row = NewRowNumber[ row ];
					col = NewColNumber[ col ];

					assert( row >= 0 && row < ActualStage2Rows );
					assert( col >= 0 && col < Stage2Col );
					break;

				case Delta::COST:
					assert( col >= Stage2Col && col < cLen );

					col = NewColNumber[ col ];

					assert( col >= 0 && col < cLen - Stage2Col );
					break;

				default:
#ifndef NDEBUG
					abort();
#endif
					break;
				}
			} // End of loop on deltas.
		} // End of loop on stoch. data blocks.
	} // End of loop on all distributions.
}

void Scenarios::ScaleObjective( Real_T div )
{
	//--------------------------------------------------------------------------
	//	Loop on all distributions.
	//
	for( Int_T d = 0; d < len; d++ )
		//	Loop on blocks in distribution "*dist[d]".
		//
		for( Int_T i = 0, l = dist[d]->Len(); i < l; i++ )
			//	Loop on delta's in block "(*dist[d])[i]".
			//
			for( Int_T j = 0, bl = (*dist[d])[i].Len(); j < bl; j++ )
			{
				Delta &delta = (*dist[d])[i][j];

				if( delta.type == Delta::COST )
					delta.value *= div;
			}
}




//@BEGIN----------------------------------------------------------------------------
//  This section contains a number of scenario related functions 
//  to be used in SADSAM and in SINGLE and TWO-REPLICATION procedures
//  for testing solution quality.
//

/*------------------------------------------------------------------------------
	Real_T Distribution::CalculateExpectedValue( void )

PURPOSE:
	To calculate expected value of a distribution
	(Used in SADSAM implementation)

    Works only for Stochastic Independent Variables right now
	  and only when probabilities are "not" stored as cumulative probabilities

------------------------------------------------------------------------------*/

Real_T Distribution::CalculateExpectedValue( void )
{
	assert( sadsam == true );
	Real_T ExpVal = 0.0; 
	Int_T i; 

	for(i=0; i<len; i++){
		ExpVal += block[i]->GetProbability() * 5; 
										//just for the time being
	}

	return ExpVal; 
}




/*------------------------------------------------------------------------------
	  void Scenarios::AppendScenario( Int_T &num )

PURPOSE:
	To append more scenarios to an already existing scenario repository
	This function adds random scenarios.
USE:
	Used generically and in SADSAM... eventually...

PARAMETERS:
	Int_T &num	is the total number of scenarios there'll be in the scenario 
				repository after it is being appended.

RETURN VALUE:	None
SIDE EFFECTS:	None.
------------------------------------------------------------------------------*/


void Scenarios::AppendScenario( Int_T &num )
{  

	int NumAppend, InitScenNum = ScenNum;
		
	NumAppend = num - InitScenNum; 
	
	assert( NumAppend >= 0 ); 
	
	//if there is at least one scenario to be appended...
	if( NumAppend >0){

		//First, a bunch of assertions and checks 
		
		MaxScen = GetNumberOfPossibleScenarios();
/*		if( num > MaxScen )
		{
			Print( "\nTOTAL NUMBER OF SCENARIOS EXCEEDS THE NUMBER OF ALL "
				"SCENARIOS FOR THIS PROBLEM.\n"); 
			NumAppend = (Int_T) MaxScen - InitScenNum; 
			
			assert( NumAppend >0 ); 
			//I am assuming this check would be done before the code enters AppendScenario...
			
			Print("INSTEAD OF %d, ADDING ONLY %d SCENARIO(S)\n", num - InitScenNum, NumAppend );
			num = (Int_T) MaxScen;
		} 
*/
		if( num > MAX_SCEN_NUM )
		{
			Print( "\tSpecified number of scenarios exceeds the maximum number of scenarios"
				"allowed.\n\tAdjusting to %d.\n",
				int( num = (Int_T) MAX_SCEN_NUM ) );
		}
		
/*		if( MaxScen > (double) MAX_SCEN_NUM )
		{
			Print( "\tToo many possible scenarios. Sampling only %d.\n", int( MAX_SCEN_NUM ) );
			num = (Int_T) MAX_SCEN_NUM;
		}
*/
		assert( num > 0 && num <= MAX_SCEN_NUM );
		
		//after the assertions, now resize the scenario repository
		ScenNum = num;
	
		Print( "\nAPPENDING %d SCENARIO(S)...", NumAppend );
		
		this ->ArrayOfScenarios.Resize( ScenNum );
		
		#ifndef NDEBUG
		{
			for( int s = 0; s < ScenNum; s++ )
				assert( ArrayOfScenarios[s] != NULL );
		}
		#endif


		//Adjust probabilities of the scenarios (e.g. 2 scen->0.5, 4 scen->0.25) 
		int s = 0; 
		for( s = 0; s < InitScenNum; s++ )
			ArrayOfScenarios[s]->SetProbability ( 1.0/(double)ScenNum );


		//Loop on all scenarios to be appended
		for( s = 0; s < NumAppend; s++ )
		{
			//create a new scenario; resize the array and add the new scenario to the array
			Scenario *scen_a = new Scenario( len, 1.0/(double)num );
			if( !scen_a ) FatalError( "Out of memory." );
		
			ArrayOfScenarios[InitScenNum + s] = scen_a;

			// Form scen_a:
			// Loop on distributions. Choose one block from each distribution and put it in scen_a
			for( int i = 0; i < len; i++ )
			{
				Real_T prob_a = Random01::Next();
				int j, l = dist[i]->Len();
			
				for( j = 0; j < l; j++ )
					if( (*dist[i])[j].GetProbability() > prob_a ) break;
			
				scen_a->Set( i, &( (*dist[i])[j] ) );
			}
		}
		Print(" DONE.\n"); 
	}

} // end of AppendScenario




/*------------------------------------------------------------------------------
	  void Scenarios::AppendScenarioRep( const Scenarios *Sc )

PURPOSE:
	This version of AppendScenario appends the scenarios of Sc to an already 
	existing scenario repository.
USE:
	Used in Two-Replication Procedures

PARAMETERS:
	const Scenarios *Sc	:  contains the scenarios to be appended.Passed by
		reference for convenience (large data) but const so that it is not 
		changed.

RETURN VALUE:	None	
SIDE EFFECTS:	
	Changes the scenario probabilities of Sc! Be careful!
------------------------------------------------------------------------------*/


void Scenarios::AppendScenarioRep( const Scenarios *Sc )
{  

	int s; 
	int NumAppend, TotalScen, InitScenNum = ScenNum;
		
	NumAppend = Sc->GetNumberOfScenarios(); 
	
	assert( NumAppend >= 0 ); 
	
	TotalScen = InitScenNum + NumAppend; 


	//if there is at least one scenario to be appended...
	if( NumAppend >0){

		//First, a bunch of assertions and checks 

		MaxScen = GetNumberOfPossibleScenarios();
		if( TotalScen > MaxScen )
		{
			Print( "\nTOTAL NUMBER OF SCENARIOS EXCEEDS THE NUMBER OF ALL "
				"SCENARIOS FOR THIS PROBLEM. ADJUSTING...\n"); 
			NumAppend = (Int_T) MaxScen - InitScenNum; 
			
			assert( NumAppend >0 ); 
			
			Print("INSTEAD OF %d, ADDING ONLY %d SCENARIO(S)\n", TotalScen - InitScenNum, NumAppend );
			TotalScen = (Int_T) MaxScen;
		} 

		if( TotalScen > MAX_SCEN_NUM )
		{
			Print( "\tSpecified number of scenarios exceeds the maximum number of scenarios"
				"allowed.\n\tAdjusting to %d.\n",
				int( TotalScen = (Int_T) MAX_SCEN_NUM ) );
		}
		
		if( MaxScen > (double) MAX_SCEN_NUM )
		{
			Print( "\tToo many possible scenarios. Sampling only %d.\n", int( MAX_SCEN_NUM ) );
			TotalScen = (Int_T) MAX_SCEN_NUM;
		}

		assert( TotalScen > 0 && TotalScen <= MAX_SCEN_NUM );
		
		//After assertions, now resize the scenario repository
		ScenNum = TotalScen;
	
		Print( "\nAPPENDING SCENARIO REPOSITORY... " );
		
		this ->ArrayOfScenarios.Resize( TotalScen );
		
		#ifndef NDEBUG
		{
			for( int s = 0; s < ScenNum; s++ )
				assert( ArrayOfScenarios[s] != NULL );
		}
		#endif


		//Adjust probabilities of the existing scenarios (e.g. 2 scen->0.5, 4 scen->0.25) 
		for( s = 0; s < InitScenNum; s++ )
			ArrayOfScenarios[s]->SetProbability ( 1.0/(double)ScenNum );


		//Loop on all scenarios to be appended and point to Sc's scenarios
		for( s = 0; s < NumAppend; s++ )
		{
			ArrayOfScenarios[InitScenNum + s] = Sc->ArrayOfScenarios[s];			
		}

		//Adjust probabilities of new scenarios (e.g. 2 scen->0.5, 4 scen->0.25) 
		//Note that since it points to Sc's scenarios, it changes probabilities of 
		//Sc's scenarios as well...
		for( s = 0; s < NumAppend; s++ )
			ArrayOfScenarios[InitScenNum + s]->SetProbability ( 1.0/(double)ScenNum );

		Print(" DONE.\n"); 
	}

} // end of AppendScenarioRep




/*------------------------------------------------------------------------------
	  void Scenarios::ReGenerateScenarios( Int_T scennum )

PURPOSE:
	 This function is used to repeat the sampling process a number of times, 
	 when a different sample need to be used each time... 
USE:	
	 SRP... testing solution quality

PARAMETERS:
	Int_T scennum	# of scenarios to generate

RETURN VALUE:	None	
SIDE EFFECTS:	(Hopefully...) None.
------------------------------------------------------------------------------*/

void Scenarios::ReGenerateScenarios( Int_T scennum, int gamma )
{
	Int_T s, i, j, l; 
	Real_T prob; 
        // David Love ---  Counts the number of batches that we use.
        // static int batchNumber = 0;
        // David Love -- hold seeds

	assert( scennum > 0 && scennum <= MAX_SCEN_NUM );

	assert( scennum <= ScenNum ); 

	// David Love -- Sets the same random seed for every set of scenarios
        // Random01::Seed( 7625, 3293, 41);

        static MTRand::uint32 batchState[ MTRand::SAVE ];
        static bool firstTime = true;
        if( firstTime ) 
        {
           Random01::GetSeed( batchState );
           firstTime = false;
        }
	//--------------------------------------------------------------------------
	//	If needed, resize. Also, Re-assign the scenario probabilities.
	//
	if( scennum < ScenNum ){
		ScenNum = scennum; 
		this ->ArrayOfScenarios.Resize( scennum );		
	}

	for(s = 0; s < scennum; s++ )
			ArrayOfScenarios[s]->SetProbability ( 1.0/(double)ScenNum );

	//(Asumes that the probabilities of distributions have been filled by an  
	// earlier call of GenerateScenarios...)

	//--------------------------------------------------------------------------
	//	Loop on all scenarios. Then, loop on distributions and re-form scenarios
	//
 
        // David Love -- Set up the random number seed again
        // Print( "Seed Values: %d %d %d\n", seedx, seedy, seedz );
        Random01::ReSeed( batchState );

        // David Love -- Print out the distribution locations, check that overlapping is working properly
        // Print( "Batch random variables\n" );
	for(s = 0; s < scennum; s++ )
	{
		for(i = 0; i < len; i++ )
		{
			prob = Random01::Next();
                        // David Love -- quick debugging statement for probability
                        // printf( "prob = %lf\n", prob );

			l = dist[i]->Len();
			for(j = 0; j < l; j++ )
                        {
                                // David Love -- more debugging check for probability break for large values
                                // printf( "j = %d, probability = %lf\n", j ,(*dist[i])[j].GetProbability() );
				if( (*dist[i])[j].GetProbability() > prob ) break;
                        }

			assert( j < l );

			ArrayOfScenarios[s]->SetAgain( i, &( (*dist[i])[j] ) );

                        // David Love -- Print out the distribution locations, check that overlapping is working properly
                        //Print( "j[%d][%d] = %3d\t", s, i, j );
		}
                // David Love -- Get the seed values to use for the overlapping
                if( s + 1 == gamma )
                {
                   Random01::GetSeed( batchState );
                }
                // David Love -- Print out the distribution locations, check that overlapping is working properly
                //Print("\n");
	}

        // David Love -- A quick ckech of how many batches we are using.
        // batchNumber++;
        // fprintf( stderr, "Batch # %d, len %d\n", batchNumber, len );


}  //end of ReGenerateScenarios



//@END------------------------------------------------------------------------------------------

