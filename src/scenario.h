/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	scenario.h
CREATED:			1994.07.27
LAST MODIFIED:		1996.02.07

DEPENDENCIES:		stdtype.h, smartptr.h,
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SCENARIO_H__
#define __SCENARIO_H__

#include <stdio.h>
#include <assert.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif
#ifndef __MY_DEFS_H__
#	include "my_defs.h"
#endif

#define COMPILE_DISTANCE_MEASUREMENT
#define MAX_SCEN_NUM			(5000000000)
#define PROBABILITY_TOL			(1.0e-4)

//==============================================================================
//	Type definitions for scenario list elements.
//
struct Delta
{
	enum DeltaType { NONE, RHS, MATRIX, COST /*, LO_BOUND, UP_BOUND */ };

	DeltaType type;
	Int_T row, col;
	Real_T value;

	Delta( void );
	Delta( const Delta &d );
	Delta( DeltaType t, Int_T r, Int_T c, Real_T v );
};

inline
Delta::Delta( void )
	: type( NONE ), row( -1 ), col( -1 ), value( 0.0 )
{}

inline
Delta::Delta( const Delta &d )
	: type( d.type ), row( d.row ), col( d.col ), value( d.value )
{}

inline
Delta::Delta( DeltaType t, Int_T r, Int_T c, Real_T v )
	: type( t ), row( r ), col( c ), value( v )
{
	assert( row >= 0 || type == COST );
	assert( col >= 0 || type == RHS );
	assert( value >= -INFINITY && value <= +INFINITY );
}
//
//	End of type definitions for scenario list elements.
//==============================================================================


//==============================================================================
//	An abstract class representing a block of random data (a part of a
//	scenario).
//
abstract class StochasticDataBlock
{
private:
	Real_T Prob;

public:
	StochasticDataBlock( Real_T prob );
	virtual ~StochasticDataBlock( void );

	Real_T GetProbability( void ) const;
	void SetProbability( Real_T p );

	virtual Int_T Len( void ) const pure;
	virtual const Delta & operator []( Int_T i ) const pure;
	virtual Delta & operator []( Int_T i ) pure;
};


inline
StochasticDataBlock::StochasticDataBlock( Real_T prob )
	: Prob( prob )
{ assert( prob >= 0.0 && prob <= 1.0 + PROBABILITY_TOL ); }


inline
StochasticDataBlock::~StochasticDataBlock( void )
{}


inline
Real_T StochasticDataBlock::GetProbability( void )
const
{ return Prob; }


inline
void StochasticDataBlock::SetProbability( Real_T p )
{
	assert( p >= 0.0 && p <= 1.0 + PROBABILITY_TOL );
	Prob = p;
}
//
//	The end of the abstract class representing a block of random data.
//==============================================================================


//==============================================================================
//	Two classes for handling two actual types of data blocks:
//	-	a class for handing just one independent random variable called
//		"IndepStochVar" and
//	-	a class for handling a true stochastic data block called
//		"IndepStochVarBlock".
//
class IndepStochVar : public StochasticDataBlock, public Delta
{
public:
	IndepStochVar( Real_T prob, Delta::DeltaType t, Int_T r, Int_T c,
		Real_T v );
	virtual ~IndepStochVar( void );

	virtual Int_T Len( void ) const;
	virtual const Delta & operator []( Int_T i ) const;
	virtual Delta & operator []( Int_T i );
};


inline
IndepStochVar::IndepStochVar( Real_T prob, Delta::DeltaType t, Int_T r, // )
	Int_T c, Real_T v )
	: StochasticDataBlock( prob ), Delta( t, r, c, v )
{}


inline
IndepStochVar::~IndepStochVar( void )
{}


inline
Int_T IndepStochVar::Len( void )
const
{ return 1; }


inline
const Delta & IndepStochVar::operator[]( Int_T // )
#ifndef NDEBUG
	i
#endif
	)
	const
{
	assert( i == 0 );
	return (const Delta &)*this;
}


inline
Delta & IndepStochVar::operator[]( Int_T // )
#ifndef NDEBUG
	i
#endif
	)
{
	assert( i == 0 );
	return (Delta &)*this;
}


//
//	And now the second class.
//


class IndepStochVarBlock : public StochasticDataBlock
{
private:
	Int_T len, maxLen;
	Array<Delta *> block;
	enum { INIT_DATA_BLOCK_SIZE = 10 };

public:
	IndepStochVarBlock( Real_T prob, Int_T InitLen = INIT_DATA_BLOCK_SIZE );
	virtual ~IndepStochVarBlock( void );

	virtual Int_T Len( void ) const;
	virtual const Delta & operator []( Int_T i ) const;
	virtual Delta & operator []( Int_T i );

	void Append( Delta *d );
};


inline
IndepStochVarBlock::IndepStochVarBlock( Real_T prob, Int_T InitLen )
	: StochasticDataBlock( prob ), len( 0 ), maxLen( InitLen ),
	block( InitLen, (Delta *)NULL )
{ assert( InitLen >= 2 ); }


inline
Int_T IndepStochVarBlock::Len( void )
const
{ return len; }


inline
const Delta & IndepStochVarBlock::operator[]( Int_T i )
	const
{
	assert( i >= 0 && i < len );
	return *(block[i]);
}


inline
Delta & IndepStochVarBlock::operator[]( Int_T i )
{
	assert( i >= 0 && i < len );
	return *(block[i]);
}
//
//	The end of the two classes handling the data blocks.
//==============================================================================


//==============================================================================
//	A class representing a single scenario.
//
class Scenario
{
private:
	Int_T ScenLen;
	Array<StochasticDataBlock *> Scen;
	Real_T Probability;

#ifdef COMPILE_DISTANCE_MEASUREMENT

public:
	enum DistanceMode { NONE, STD_NORM, WEIGHTED_NORM, DISCRETE_MEASURE };

private:
	static DistanceMode DistMode;
	static Long_T DistanceCalculationCnt;
#endif

public:
	Scenario( Int_T Len, Real_T Prob = 0.0 );
	~Scenario( void );

	StochasticDataBlock & operator[] ( int i );
	const StochasticDataBlock & operator[] ( int i ) const;

	Int_T GetLength( void ) const;

	void Set( Int_T pos, StochasticDataBlock *s );

	//@BEGIN----------------------------------------------
	void SetAgain( Int_T pos, StochasticDataBlock *d );
	//@END------------------------------------------------

	void SetProbability( Real_T Prob );
	Real_T GetProbability( void ) const;

#ifdef COMPILE_DISTANCE_MEASUREMENT
	//--------------------------------------------------------------------------
	//	Data and functions used for computing distances between the scenarios.
	//
public:
	static void SetDistanceMode( DistanceMode dm );

        //Rebecca's code to get scenario information
	static Real_T GetScenarioComponents( const Scenario &s1 );
	static Real_T GetDistance( const Scenario &s1, const Scenario &s2 );

	static void ResetDistanceCalculationCnt( void );
	static Long_T GetDistanceCalculationCnt( void );

private:
	static Real_T GetStdNormDistance( const Scenario &s1, const Scenario &s2 );
	static Real_T GetWeightedNormDistance( const Scenario &s1,
		const Scenario &s2 );
	static Real_T GetDiscreteMeasureDistance( const Scenario &s1,
		const Scenario &s2 );
#endif
};
//
//  End of a single scenario class.
//==============================================================================

//==============================================================================
//	A class reprsenting a single (discrete) distribution of a random variable.
//
class Distribution
{
private:
	enum { INIT_DIST_LEN = 10 };

	Int_T len, maxLen;
	Array<StochasticDataBlock *> block;

public:
	Distribution( void );
	~Distribution( void );

	Int_T Len( void );
	const StochasticDataBlock &operator []( Int_T i ) const;
	StochasticDataBlock &operator []( Int_T i );

	void Append( StochasticDataBlock *bl );

	//@BEGIN-------------------------------------------------
	Real_T CalculateExpectedValue ( void ); 
	//@END---------------------------------------------------

};
//
//	End of the class representing a single random distribution.
//==============================================================================

class DeterministicLP;

//==============================================================================
//	A class representing a scenario repository.
//
class Scenarios
{
protected:
	//--------------------------------------------------------------------------
	//	Scenarios as an array of arrays of simple LP changes (Deltas).
	//
	Int_T ScenNum;
	Array<Scenario *> ArrayOfScenarios;

	enum Status { EMPTY, READING_FILE, DISTRIBUTIONS_READY, SCENARIOS_READY };
	Status status;

	//--------------------------------------------------------------------------
	//	Data used for scenario file parsing (and only during parsing!).
	//
private:
	SortedArrayOfLabels *RowLabels, *ColumnLabels;
	const char *ProblemName, *RHS_Name, *ObjName, *BoundsName, *RangesName;
	int ProblemNameLen, RHS_NameLen, ObjNameLen, BoundsNameLen, RangesNameLen;

	Bool_T SemanticError;

	//--------------------------------------------------------------------------
	//	More data used by semantic actions.
	//
	//	The data collected during the scenario file reading consists of an array
	//	of distributions of all random stochastic data blocks (incl. one-datum
	//	blocks). An array of distributions is kept here.
	//
	enum { INIT_DIST_NUM = 10 };
	Int_T len, maxLen;
	Array<Distribution *> dist;
	Real_T MaxScen;				// The total number of possible scenarios.

	//
	//	What was the last item that was processed?
	//
	enum StochBlockType { NONE, INDEP_DISCRETE, BLOCKS_DISCRETE };

	StochBlockType LastStochBlockType;
	Delta LastIndepDiscrete;
	char LastBlockName[20];
	IndepStochVarBlock *LastBlock;
	//
	//  End of data used by the semantic actions.
	//--------------------------------------------------------------------------

public:
	//--------------------------------------------------------------------------
	//  Constructor, destructor.
	//
	Scenarios( void );
	virtual ~Scenarios( void );

private:
	void Cleanup( void );
	void CleanupDistributions( void );

	void AppendNewDistribution( void );
	void AppendLastBlock( void );

	//--------------------------------------------------------------------------
	//	Procedure for reading a scenario file and generating all scenarios.
	//
public:
	virtual Bool_T ReadScenarioFile( const char *StochFileName,
		FILE *StochFilePtr, DeterministicLP &DeterministicMatrix );
	virtual Bool_T GenerateScenarios( Int_T &num );
	void RenumberIndiceInScenarios( Array<Int_T> &NewRowNumber, Int_T rLen,
		Array<Int_T> &NewColNumber, Int_T cLen, Int_T ActualStage2Rows,
		Int_T Stage2Col );
	void ScaleObjective( Real_T div );

	//@BEGIN-----------------------------------
	
	void AppendScenario( Int_T &num );
	void AppendScenarioRep( const Scenarios *Sc );
	void ReGenerateScenarios( Int_T scennum, int gamma ); 

	//@END-------------------------------------

private:
	void GenerateAllScenarios( Int_T num );
	void GenerateScenarioSample( Int_T num );
	Real_T GetNumberOfPossibleScenarios( void ) const;
	Int_T CalculateAccumulatedProbabilities( Bool_T Store = True );

#ifndef NDEBUG
private:
	void CheckArrayOfDistributions( void );
#endif

	//--------------------------------------------------------------------------
	//	Data access functions: allow reading the scenarios once they are
	//	generated.
	//
public:
	Int_T GetNumberOfScenarios( void ) const;
	virtual Scenario & operator[]( int ScenarioNumber );
	virtual const Scenario & operator[]( int ScenarioNumber ) const;

	virtual Int_T PreviousBlockNumber( Int_T CurrentBlock ) const;

public:
	//--------------------------------------------------------------------------
	//	Auxiliary functions for scenario file reading and scenario generation
	//	(semantic actions, called from outside the class).
	//
	void DeltaEntry( const char *Row, const char *Col,
		Delta::DeltaType &type, Int_T &row, Int_T &col );
	void NewIndepDiscreteEntry( const char *Col, const char *Row,
		const char *PeriodLab, double Val, double Prob );
	void NewBlocksDiscrete( const char *BlockName, const char *PeriodLabel,
		Real_T Prob );
	void NewBlocksDiscreteEntry( const char *Row, const char *Col, Real_T Val );
	void EndData( void );

#ifdef COMPILE_DISTANCE_MEASUREMENT
	//--------------------------------------------------------------------------
	//	Functions used for computing distances between the scenarios.
	//
public:
	Real_T GetDistance( Int_T s1, Int_T s2 );
	Real_T GetScenarioComponents( Int_T s1 );
#endif
};
//
//	End of the scenario repository class declaration.
//
//==============================================================================


//==============================================================================
//
//	Some inline function's implementations.
//
//==============================================================================

//------------------------------------------------------------------------------
//	Class Scenario
//
inline
Scenario::Scenario( Int_T Len, Real_T Prob )
	: ScenLen( Len ), Scen( ScenLen, (StochasticDataBlock *)NULL ),
	Probability( Prob )
{ assert( Len > 0 ); }


//
//	Note: a scenario does not hold it's own data, just the pointers into the
//	data stored in the distribution array.
//
inline
Scenario::~Scenario( void )
{}


//
//	Access to scenario data.
//
inline
StochasticDataBlock & Scenario::operator[]( int i )
{
	assert( i >= 0 && i < ScenLen );
	assert( Scen[i] != NULL );

	return *(Scen[i]);
}


inline
const StochasticDataBlock & Scenario::operator[]( int i )
const
{
	assert( i >= 0 && i < ScenLen );
	assert( Scen[i] != NULL );

	return *(Scen[i]);
}


inline
Int_T Scenario::GetLength( void )
const
{ return ScenLen; }


//
//	Scenario creation.
//
inline
void Scenario::Set( Int_T pos, StochasticDataBlock *d )
{
	assert( pos >= 0 && pos < ScenLen );
	assert( Scen[pos] == NULL );
	assert( d != NULL );

	Scen[pos] = d;
}

//@BEGIN------------------------------------------------------
// The difference between the above is that Scen is not null
// it exists, but we are changing its values... 

inline 
void Scenario::SetAgain(Int_T pos, StochasticDataBlock *d )
{
	assert( pos >= 0 && pos < ScenLen );
	assert( d != NULL );

	Scen[pos] = d;
}
//@END--------------------------------------------------------


inline
void Scenario::SetProbability( Real_T prob )
{
	assert( prob >= 0.0 && prob <= 1.0 + PROBABILITY_TOL );
	Probability = prob;
}


//
//	More access functions.
//
inline
Real_T Scenario::GetProbability( void )
const
{ return Probability; }


//
//	Handling the counter of scenario distance measurements.
//
#ifdef COMPILE_DISTANCE_MEASUREMENT
inline
void Scenario::ResetDistanceCalculationCnt( void )
{ DistanceCalculationCnt = 0; }


inline
Long_T Scenario::GetDistanceCalculationCnt( void )
{ return DistanceCalculationCnt; }
#endif


//------------------------------------------------------------------------------
//	Class Scenarios
//
inline
Scenarios::~Scenarios( void )
{ Cleanup(); }


inline
Int_T Scenarios::GetNumberOfScenarios( void )
const
{
	assert( status == SCENARIOS_READY );
	return ScenNum;
}


inline
Scenario & Scenarios::operator[]( int ScenarioNumber )
{
	assert( status == SCENARIOS_READY );
	assert( ScenNum > 0 );
	assert( ScenarioNumber >= 0 && ScenarioNumber < ScenNum );

	return *ArrayOfScenarios[ ScenarioNumber ];
}


inline
const Scenario & Scenarios::operator[]( int ScenarioNumber )
const
{
	assert( status == SCENARIOS_READY );
	assert( ScenNum > 0 );
	assert( ScenarioNumber >= 0 && ScenarioNumber < ScenNum );

	return *ArrayOfScenarios[ ScenarioNumber ];
}


inline
Int_T Scenarios::PreviousBlockNumber( Int_T CurrentBlock )
const
{
	assert( ScenNum > 0 );
	assert( CurrentBlock >= 0 && CurrentBlock < ScenNum );

	return Int_T( (CurrentBlock == Int_T( ScenNum - 1 )) ? 0 : --CurrentBlock );
}



//------------------------------------------------------------------------------
//	Class Distribution
//
inline
Distribution::Distribution( void )
	: len( 0 ), maxLen( INIT_DIST_LEN ),
	block( maxLen, (StochasticDataBlock *)NULL )
{}


inline
Int_T Distribution::Len( void )
{ return len; }


inline
const StochasticDataBlock & Distribution::operator[]( Int_T i )
	const
{
	assert( i >= 0 && i < len );
	return *block[i];
}


inline
StochasticDataBlock & Distribution::operator[]( Int_T i )
{
	assert( i >= 0 && i < len );
	return *block[i];
}


#endif
