/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	scentree.cpp
CREATED:			1994.08.07
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		scentree.h, scenario.h, error.h, print.h
		            <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __SCENTREE_H__
#   include "scentree.h"
#endif
#ifndef __ERROR_H__
#   include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T TreeOfScenarios::GenerateScenarios( Int_T &num );

PURPOSE:
	Generates a set of scenarios using the underlying class "Scenarios" and its
"GenerateScenarios()" function. Then proceeds to create a tree of scenarios. It
is supposed to be the "cheapest" spanning tree with respect to the sum of
distances between scenarios connected by arcs.
	For this purpose a smart (well, at least I believe it's smart) heuristic
strategy is developed.

PARAMETERS:
	Int_T &num
		Reference parameter passed down to the "Scenarios::GenerateScenarios()"
		routine, which uses it to determine the actual number of scenarios to be		generated. Then it puts the result in the same argument.

RETURN VALUE:
	Boolean success indicator.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T TreeOfScenarios::GenerateScenarios( Int_T &num )
{
	//--------------------------------------------------------------------------
	//	First call the actual scenario generation procedure.
	//
	if( !Scenarios::GenerateScenarios( num ) )
		return False;

	assert( ScenNum > 0 );

	Order.Resize( ScenNum );		Order.Fill( -1, ScenNum );
	Predecessor.Resize( ScenNum );	Predecessor.Fill( -1, ScenNum );

	//--------------------------------------------------------------------------
	//	Calculate the cost of traversing the scenarios in the order in which
	//	they are generated.
	//
	Scenario::SetDistanceMode( Scenario::DISCRETE_MEASURE );
	Scenario::ResetDistanceCalculationCnt();

	Print( "\n\tSCENARIO ORDERING ROUTINE: SPANNING TREE GENERATION.\n"
		"\t%-40s%12.2E.\n",
		"The cost in scenario generation order:", TotalStraightCost() );

#ifdef PRIM_DIJKSTRA_ALGORITHM
	//==========================================================================
	//
	//	Cheapest spanning tree generation by Prim-Dijkstra algorithm.
	//
	//==========================================================================

	//
	//	Some data needed for generation of the "Order" array.
	//
	Array<Int_T> Graph( ScenNum );		//	Graph of the nodes.

	MakeTreeByPrimDijkstraAlgorithm( Graph );
	Print( "\t%-40s%12ld\n", "The number of scenario comparisons:", 
		(long) Scenario::GetDistanceCalculationCnt() );

#else
	//==========================================================================
	//
	//	Cheapest spanning tree generation by my own heuristic.
	//
	//==========================================================================

	//
	//	Some data needed for generation of the "Order" array.
	//
	Array<Int_T> Graph( ScenNum ),		//	Graph of the nodes.
		ChainArray( ScenNum );			//	A work array needed for spanning
										//	tree generation.
	Int_T Chain = -1;					//	The number of chains.

	//--------------------------------------------------------------------------
	//	Make graph in which the arcs will join the closest nodes (corresponding
	//	to scenarios).
	//	
	MakeInitialGraph( Graph );
	Print( "\t%-40s%12ld\n", "The number of scenario comparisons:", 
		(long) Scenario::GetDistanceCalculationCnt() );

	//--------------------------------------------------------------------------
	//	Cut cycles and thus create a graph that would consist of a number of
	//	trees.
	//
	CutCyclesInGraph( Graph, ChainArray, Chain );
	Print( "\t%-40s%12d\n", "Chains found:", (int)Chain );

	//--------------------------------------------------------------------------
	//	Connect the trees into one spanning tree of the graph.
	//
	MarkSubTrees( Graph, ChainArray, Chain );
	ConnectTrees( Graph, ChainArray, Chain );
#endif

	Print( "\t%-40s%12.2E.\n", "The cost in tree order:",
		TotalTreeCost( Graph ) );

	//--------------------------------------------------------------------------
	//	Fill the "Order" array.
	//
	DetermineNodeOrder( Graph );

	//--------------------------------------------------------------------------
	Print( "\tSCENARIO ORDERING ROUTINE FINISHED SUCCESSFULLY.\n" );

	return True;
}


#ifdef PRIM_DIJKSTRA_ALGORITHM

/*------------------------------------------------------------------------------

	void TreeOfScenarios::MakeTreeByPrimDijkstraAlgorithm( Array<Int_T> &Graph )

PURPOSE:
	Creates a cheapest spanning tree in the complete graph of "ScenNum" nodes.
Uses algorithm of R.C. Prim [1] and E. W. Dijkstra [2].

[1]	Prim R. C. "Shortest connection networks and some generalizations" Bell
	Syst. Tech J. 1957 Vol. 36
[2]	Dijkstra E. W. "A note on two problems in connection with graphs" Numerische
Mathematik 1959 Vol. 1

PARAMETERS:
	Array<Int_T> &Graph
		The array that will store the graph of the spanning tree.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::MakeTreeByPrimDijkstraAlgorithm( Array<Int_T> &Graph )
{
	const Int_T n = ScenNum;

	Array<Bool_T> mark( n, False );		// Marks nodes belonging to the tree.

	//	For nodes not in the tree: the shotrest distance to the tree current
	//	and the node in the tree to which i-th node would have to be attached
	//	to be just that shortest distance away.
	//
	Array<Real_T> bestCost( n, +INFINITY );
	Array<Int_T> bestInd( n, 0 );

	Graph.Fill( -1, ScenNum );			// Empty tree (no arcs).

	//--------------------------------------------------------------------------
	//	Creae the initial tree consisting of node 0 only. Initialize the
	//	"bestCost" array ("bestInd" is already initialize to point to node 0).
	//
	Int_T i;

	mark[0] = True;
	for( i = 1; i < n; i++ )
		bestCost[i] = GetDistance( i, 0 );

	//--------------------------------------------------------------------------
	//	Loop to find the rest of the connections.
	//
	for( Int_T iter = 1; iter < n; iter++ )
	{
		//----------------------------------------------------------------------
		//	What shall we add to the tree?
		//
		Int_T ind	= -1,
			i;
		Real_T cost	= +INFINITY;

		for( i = 0; i < n; i++ )
			if( !mark[i] && bestCost[i] < cost )
			{
				cost	= bestCost[i];
				ind		= i;
			}

		//----------------------------------------------------------------------
		//	Now append the new node to the tree.
		//
		assert( ind >= 0 && ind < n );
		assert( cost >= 0.0 && cost < +INFINITY );
		assert( mark[ bestInd[ ind ] ] );

		mark[ind]	= True;
		Graph[ind]	= bestInd[ind];

		//----------------------------------------------------------------------
		//	Update "bestCost" and "bestInd" arrays.
		//
		for( i = 0; i < n; i++ )
			if( !mark[i] )
			{
				cost = GetDistance( i, ind );
				if( cost < bestCost[i] )
				{
					bestCost[i]	= cost;
					bestInd[i]	= ind;
				}
			}
	}

#ifndef NDEBUG
	//--------------------------------------------------------------------------
	//	We won't believe it till we see it.
	//
	for( i = 0; i < n; i++ )
	{
		assert( mark[i] );
		assert( i == 0 || Graph[i] >= 0 );
	}
#endif
}


#else

/*------------------------------------------------------------------------------

	void TreeOfScenarios::MakeInitialGraph( void )

PURPOSE:
	Creates a graph in which scenarios (in the scenario repository) are the
nodes. An arc connects each scenario with another scenario, which is the most
similar to the latter. The graph is guaranteed to have cycles (unfortunately).

%	An option for future:
%		To reduce the number of scenario-to-scenario distance computations we
%		may only scan some candidates.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::MakeInitialGraph( Array<Int_T> &Graph )
{
	Graph.Fill( -1, ScenNum );

	Array<Real_T> MinDist( ScenNum, +INFINITY );

	//--------------------------------------------------------------------------
	//	For each scenario "i" (which is not connected to the graph yet) ...
	//
	for( Int_T i = 0; i < ScenNum; i++ )
	{
		assert( Graph[i] < 0 );

		//
		//	Look for the most similar scenario "j" (exclude cycles of length 1).
		//
		for( Int_T j = i; j < ScenNum; j++ )
			if( j != i /* && Graph[j] != i */ )
			{
				Real_T d = GetDistance( i, j );

				if( d < MinDist[i] )
				{
					Graph[i]	= j;
					MinDist[i]	= d;
				}
				if( d < MinDist[j] )
				{
					Graph[j]	= i;
					MinDist[j]	= d;
				}
			}
	}

	//--------------------------------------------------------------------------
	//	Now let's check if all nodes are connected to something.
	//
#ifndef NDEBUG
	{
		for( Int_T i = 0; i < ScenNum; i++ )
			assert( Graph[i] >= 0 && Graph[i] < ScenNum );
	}
#endif
}


/*------------------------------------------------------------------------------

	void TreeOfScenarios::CutCyclesInGraph( void )

PURPOSE:
	The cycles are detected and cut here. Additionally, "ChainConversio" array
is filled with data needed in order to mark sub-trees. Tis data will be later
processed in procedure "TreeOfScenarios::ConnectTrees()".

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::CutCyclesInGraph( Array<Int_T> &Graph, // )
	Array<Int_T> &ChainArray, Int_T &Chain )
{
	ChainArray.Fill( -1, ScenNum );

	//--------------------------------------------------------------------------
	//	Scan all nodes of the graph.
	//
	Chain = 0;
	for( Int_T i = 0; i < ScenNum; i++ )
		//
		//	If a node does not belong to any chain we found so far ...
		//
		if( ChainArray[i] < 0 )
		{
			//
			//	... scan the chain it belongs to ...
			//
			for( Int_T j = i, next; ; j = next )
			{
				//
				//	... and mark all the nodes you find.
				//
				ChainArray[j] = Chain;

				Int_T nextChain = ChainArray[ next = Graph[j] ];

				//
				//	If you find a cycle, cut it. This ends the scan of the
				//	current chain.
				//
				if( nextChain == Chain )
				{
					Graph[j] = -1;
					break;
				}
				//
				//	If you run into another chain, we know it was already
				//	checked, so it does not contin any cycles.
				//
				else if( nextChain >= 0 )
					break;
			}

			//
			//	Increment the chain number.
			//
			Chain++;
		}

	assert( ! HasCycles( Graph ) );
}


/*------------------------------------------------------------------------------

	void TreeOfScenarios::MarkSubTrees( Array<Int_T> &Graph,
		Array<Int_T> &ChainArray, Int_T &Chain )

PURPOSE:
	x

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::MarkSubTrees( Array<Int_T> &Graph, // )
	Array<Int_T> &ChainArray,
#ifndef NDEBUG
	Int_T &Chain )
#else
	Int_T & )
#endif
{
	assert( Chain > 0 );

	//--------------------------------------------------------------------------
	//	Loop on all entries of the chain number array.
	//
	Int_T i;
	for( i = 0; i < ScenNum; i++ )
	{
		//
		//	Skip roots and sub-trees that were scanned previously (they are
		//	marked by negative entries in the "ChainArray" - see below).
		//
		if( ChainArray[i] < 0 ) continue;

		//
		//	Follow the conversion chain (till the end).
		//
		Int_T ChangeTo = -1, j;
		for( j = i; j >= 0 ; j = Graph[j] )
		{
			ChangeTo = ChainArray[j];

			//
			//	If we hit a sub-tree branch that was scanned before, we end
			//	the search.
			//
			if( ChangeTo <= -2 )
			{
				ChangeTo = Int_T( -2 - ChangeTo );
				break;
			}
		}

		assert( ChangeTo >= 0 && ChangeTo < Chain );

		//
		//	If a meaningful conversion was found, perform it on the chain that
		//	was just searched. Store a negative number so as to be able to skip
		//	some searches (see above).
		//
		for( j = i; j >= 0; j = Graph[j] )
			ChainArray[j] = Int_T( -2 - ChangeTo );
	}

	//--------------------------------------------------------------------------
	//	Reverse all "ChainArray" entries when required.
	//
	for( i = 0; i < ScenNum; i++ )
	{
		Int_T &c = ChainArray[i];

		if( c <= -2 )
			c = Int_T( -2 - c );
	}

	//--------------------------------------------------------------------------
	//	Some simple checking: check the range of the "ChainArray" entries.
	//
#ifndef NDEBUG
	for( i = 0; i < ScenNum; i++ )
		assert( ChainArray[i] >= -1 && ChainArray[i] < Chain );
#endif
}


/*------------------------------------------------------------------------------

	void TreeOfScenarios::ConnectTrees( void )

PURPOSE:
	This procedure does a few things:
1)	first it makes sure each tree has exactly one number in the "ChainArray",
2)	then for each sub-tree root it searches for the closest node between the
	nodes of the other sub-trees.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::ConnectTrees( Array<Int_T> &Graph, // )
	Array<Int_T> &ChainArray, Int_T &Chain )
{
	assert( Chain > 0 );

	//--------------------------------------------------------------------------
	//	Count the sub-trees. Exit if there is only one sub-tree.
	//
	Int_T SubTreeCnt = CountSubTrees( ChainArray, Chain );

	//
	//	See if the number calculated above matches the number of sub-tree roots.
	//
#ifndef NDEBUG
	{
		Int_T cnt = 0;
		for( Int_T i = 0; i < ScenNum; i++ )
			if( Graph[i] < 0 )
				cnt++;
		assert( cnt == SubTreeCnt );
	}

#endif
/**/
	Print( "\t%-40s%12d\n", "Sub-trees found:", (int) SubTreeCnt );
/**/

	if( SubTreeCnt == 1 ) return;

	//--------------------------------------------------------------------------
	//	Find which sub-tree root should be connected where.
	//
	for( Int_T i = 0; i < ScenNum; i++ )
	{
		//----------------------------------------------------------------------
		//	Look only for sub-tree roots.
		//
		if( Graph[i] >= 0 ) continue;

		Int_T MinInd	= -1;
		Real_T MinDist	= +INFINITY;

		//----------------------------------------------------------------------
		//	Scan all other sub-trees for the cheapest connection point.
		//
		for( Int_T j = 0; j < ScenNum; j++ )
			if( ChainArray[j] != ChainArray[i] )
			{
				assert( i != j );

				Real_T d = GetDistance( i, j );

				if( d < MinDist )
				{
					MinInd	= j;
					MinDist	= d;
				}
			}

		//----------------------------------------------------------------------
		//	If a connection point was found ...
		//
		if( MinInd >= 0 )
		{
			//
			//	... connect ...
			//
			Graph[i] = MinInd;

			//
			//	... and mark as belonging to the same sub-tree.
			//
			Int_T ParentTree	= ChainArray[ MinInd ],
				CurrentTree		= ChainArray[i];

			assert( ParentTree != CurrentTree );

			for( Int_T j = 0; j < ScenNum; j++ )
				if( ChainArray[j] == CurrentTree )
					ChainArray[j] = ParentTree;
		}
		else
		//----------------------------------------------------------------------
		//	If a connection point was not found - it means we have only one
		//	sub-tree by now.
		//
		{
#ifndef NDEBUG
			Int_T SubTrees = 0;
			for( Int_T j = 0; j < ScenNum; j++ )
				if( Graph[j] < 0 )
					SubTrees++;
			assert( SubTrees == 1 );
#endif
			break;
		}
	}

	//--------------------------------------------------------------------------
	//	The tree should not have any cycles.
	//
	assert( ! HasCycles( Graph ) );
}


/*------------------------------------------------------------------------------

	void TreeOfScenarios::CountSubTrees( void )

PURPOSE:
	Of no real importance. Counts sub-trees.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T TreeOfScenarios::CountSubTrees( Array<Int_T> &ChainArray, // )
	Int_T &Chain )
{
	//--------------------------------------------------------------------------
	//	We may calculate how many sub-trees there are.
	//
	Array<Bool_T> b( Chain );
	b.Fill( False, Chain );

	Int_T i;
	for( i = 0; i < ScenNum; i++ )
	{
		Int_T ch = ChainArray[i];

		assert( ch < Chain );

		if( ch >= 0 ) b[ch] = True;
	}

	Int_T SubTreeCnt;
	for( SubTreeCnt = i = 0; i < Chain; i++ )
		if( b[i] )
			SubTreeCnt++;

	return SubTreeCnt;
}

#endif // ifdef PRIM_DIJKSTRA_ALGORITHM


/*------------------------------------------------------------------------------

	double TreeOfScenarios::TotalStraightCost( void )
	double TreeOfScenarios::TotalTreeCost( void )

PURPOSE:
	Those two functions calculate the cost of traversing the set scenarios
a)	in the order in which they were generated,
b)	using the tree.

PARAMETERS:
	None.

RETURN VALUE:
	Non-negative (and practically always positive) number.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Real_T TreeOfScenarios::TotalStraightCost( void )
{
	Real_T cost = 0.0;

	for( Int_T i = 1; i < ScenNum; i++ )
		cost += GetDistance( Int_T( i - 1 ), i );

	return cost;
}


Real_T TreeOfScenarios::TotalTreeCost( Array<Int_T> &Graph )
{
	Real_T cost = 0.0;

	for( Int_T i = 0; i < ScenNum; i++ )
		if( Graph[i] >= 0 )
			cost += GetDistance( i, Graph[i] );

	return cost;
}


/*------------------------------------------------------------------------------

	void TreeOfScenarios::DetermineNodeOrder( void )

PURPOSE:
	x

PARAMETERS:
	Int_T Len
		x

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void TreeOfScenarios::DetermineNodeOrder( Array<Int_T> &Graph )
{
	Order.Fill( -1, ScenNum );

	//--------------------------------------------------------------------------
	//	Loop on all nodes to construct the "Order" array.
	//
	Int_T i = 0;
	for( Int_T OrdStart = 0; i < ScenNum; i++ )
	{
		//----------------------------------------------------------------------
		//	Skip nodes that were already placed in the "Order" array.
		//
		if( Graph[i] < -1 ) continue;

		//----------------------------------------------------------------------
		//	Compute the length of the chain from the current node to the
		//	last node that was not already placed in the "Order" array.
		//
		//	Note: the first search will find the tree root.
		//
		Int_T node, len;
		for( len = 1, node = i; Graph[node] >= 0; node = Graph[node] )
		{
			len++;
			assert( len <= ScenNum );
		}

		if( Graph[node] < -1 )	// Correction if we did not finish at the root
			len--;				// or if the root was already visited.

		//----------------------------------------------------------------------
		//	Scan the chain once again. This time store the entries in the
		//	"Order" array and mark the entries as done with by assigning them
		//	negative numbers:
		//		(-2)	for the root,
		//		(-3-k)	for any other node which is a successor of node "k".
		//
		Int_T OrdEnd = Int_T( OrdStart + len ),
			next;

		for( node = i; --len >= 0; node = next )
		{
			assert( node >= 0 );
			next = Graph[node];
			assert( next >= -1 );

			assert( ( OrdStart + len == ScenNum ) ||
				Order[ OrdStart + len ] == -1 );
			Order[ OrdStart + len ] = node;
			Graph[node] =
				Int_T( ( next == -1 || next == -2 ) ? -2 : -3 - next );
		}

		assert( Order[OrdStart] >= 0 );
		assert( Order[OrdEnd-1] >= 0 );
		OrdStart = OrdEnd;
	}

	//--------------------------------------------------------------------------
	//	See if the "Order" array has been filled correctly.
	//
#ifndef NDEBUG
	//	The root should not have predecessors. "Order" array should begin with
	//	the root node. The root node is (just now) marked by (-2).
	//
	assert( Graph[ Order[0] ] == -2 );

	//	All entries should be non-negative and each graph node should occur
	//	exactly once.
	//
	Array<Bool_T> Check( ScenNum, False );

	for( i = 0; i < ScenNum; i++ )
	{
		assert( Order[i] >= 0 );

		assert( ! Check[Order[i]] );
		Check[Order[i]] = True;
	}

	for( i = 0; i < ScenNum; i++ )
		assert( Check[Order[i]] );
#endif

	//--------------------------------------------------------------------------
	//	Restore the contents of the "Graph" array.
	//
	for( i = 0; i < ScenNum; i++ )
		if( Graph[i] == -2 )
			Graph[i] = -1;
		else if( Graph[i] < -2 )
			Graph[i] = Int_T( -3 - Graph[i] );
#ifndef NDEBUG
		else
			abort();
#endif

	//--------------------------------------------------------------------------
	//	Now is the time to construct the "Predecessor" array.
	//
	Predecessor.Fill( -1, ScenNum );

	//
	//	Construct a permutation reverse to "Order" array. It will be needed for
	//	the "Predecessor" array construction.
	//
	Array<Int_T> BackOrd( ScenNum, -1 );
	for( i = 0; i < ScenNum; i++ )
	{
		assert( BackOrd[Order[i]] == -1 );
		BackOrd[Order[i]] = i;
	}

	//
	//	"Graph" describes the tree by pointing from each node to its parent
	//	node. Thus it is sufficient to permute the "Graph" by "Order" to obtain
	//	what we need.
	//
	for( i = 0; i < ScenNum; i++ )
		Predecessor[i] = Int_T( ( Graph[ Order[i] ] == -1 ) ? -1 :
			BackOrd[ Graph[ Order[i] ] ] );

	assert( Predecessor[0] == -1 );
}


#ifndef NDEBUG

/*------------------------------------------------------------------------------

	Bool_T TreeOfScenarios::HasCycles( Array<Int_T> &Graph )

PURPOSE:
	x

PARAMETERS:
	Int_T Len
		x

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T TreeOfScenarios::HasCycles( Array<Int_T> &Graph )
{
	//--------------------------------------------------------------------------
	//	Loop on all nodes. For each node scan the entire path to the root and
	//	see if the path does not go back to the node you started from.
	//
	for( Int_T i = 0; i < ScenNum; i++ )
		for( Int_T j = Graph[i], len = 0; j >= 0; j = Graph[j], len++ )
			if( j == i || len > ScenNum )
				return True;

	return False;
}

#endif


/*------------------------------------------------------------------------------

	Int_T TreeOfScenarios::PreviousBlockNumber( Int_T CurrentBlock ) const

PURPOSE:
	When solving a set of scenarios (nodes) in the tree order, a solution
procedure has to know from which prevoiusly solved subproblem the solution of
the current one has to be restarted. But the solution procedure obly knows the
tree order of the scenarios. Thus it is needed to give the scenario number
back-permuted by the "Order" array. This information is stored in the
"Predecessor" array.

PARAMETERS:
	Int_T CurrentBlock
		Current block number.

RETURN VALUE:
	The index into the "order" array which points at the current node's
predecessor in the graph.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Int_T TreeOfScenarios::PreviousBlockNumber( Int_T CurrentBlock )
const
{
	assert( CurrentBlock >= 0 && CurrentBlock < ScenNum );

	Int_T Pred = Predecessor[CurrentBlock];

	assert( Pred >= 0 && Pred < ScenNum );

	return Pred;
}
