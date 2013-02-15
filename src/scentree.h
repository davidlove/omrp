/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	scentree.h
CREATED:			1994.08.07
LAST MODIFIED:		1995.09.28

DEPENDENCIES:		stdtype.h, smartptr.h, scenario.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SCENTREE_H__
#define __SCENTREE_H__

#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif

#define PRIM_DIJKSTRA_ALGORITHM

//==============================================================================
//
//	A class representing a scenario repository with a tree-like structure.
//	Scenarios are initially treated as a complete graph. A spanning tree is then
//	found. It is formed so as to ensure that the distance between a parent node
//	and its children be as small as possible.
//
class TreeOfScenarios : public Scenarios
{
private:
	Array<Int_T> Order;			// Scenario order generated via a spanning tree
								// construction.
	Array<Int_T> Predecessor;	// Array of predecessor numbers: for each entry
								// of the "Order" array we store the index of
								// it's predecessor in graph relative to the
								// "Order" array.

public:
	//--------------------------------------------------------------------------
	//  Constructor, destructor.
	//
	TreeOfScenarios( void );
	virtual ~TreeOfScenarios( void );

	//--------------------------------------------------------------------------
	//  Procedure for reading a scenario file and generating all scenarios.
	//
	virtual Bool_T GenerateScenarios( Int_T &num );

	//--------------------------------------------------------------------------
	//  Data access functions: allow reading the scenarios once they are
	//  generated.
	//
	virtual Scenario & operator[]( int ScenarioNumber );
	virtual const Scenario & operator[]( int ScenarioNumber ) const;

	virtual Int_T PreviousBlockNumber( Int_T CurrentBlock ) const;

	//--------------------------------------------------------------------------
	//	Auxiliary functions needed for scenario ordering.
	//
private:

#ifdef PRIM_DIJKSTRA_ALGORITHM
	void MakeTreeByPrimDijkstraAlgorithm( Array<Int_T> &Graph );
#else
	void MakeInitialGraph( Array<Int_T> &Graph );
	void CutCyclesInGraph( Array<Int_T> &Graph, Array<Int_T> &ChainArray,
		Int_T &Chain );
	void MarkSubTrees( Array<Int_T> &Graph, Array<Int_T> &ChainArray,
		Int_T &Chain );
	Int_T CountSubTrees( Array<Int_T> &ChainArray, Int_T &Chain );
	void ConnectTrees( Array<Int_T> &Graph, Array<Int_T> &ChainArray,
		Int_T &Chain );
#endif

	void DetermineNodeOrder( Array<Int_T> &Graph );

	Real_T TotalStraightCost( void );
	Real_T TotalTreeCost( Array<Int_T> &Graph );

#ifndef NDEBUG
	Bool_T HasCycles( Array<Int_T> &Graph );
#endif
};
//
//	The end of the class representing the scenario tree.
//
//==============================================================================

//==============================================================================
//
//  Some inline function's implementations.
//
//==============================================================================

inline
TreeOfScenarios::TreeOfScenarios( void )
	: Scenarios(), Order(), Predecessor()
{}


inline
TreeOfScenarios::~TreeOfScenarios( void )
{}


inline
Scenario & TreeOfScenarios::operator[]( int ScenarioNumber )
{
	assert( ScenNum > 0 );
	assert( Order[ScenarioNumber] >= 0 && Order[ScenarioNumber] < ScenNum );

	return *ArrayOfScenarios[ Order[ScenarioNumber] ];
}


inline
const Scenario & TreeOfScenarios::operator[]( int ScenarioNumber )
const
{
	assert( ScenNum > 0 );
	assert( Order[ScenarioNumber] >= 0 && Order[ScenarioNumber] < ScenNum );

	return *ArrayOfScenarios[ ScenarioNumber ];
}


#endif
