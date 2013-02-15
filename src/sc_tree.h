/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	sc_tree.h
CREATED:			1996.01.20
LAST MODIFIED:		1996.04.06

DEPENDENCIES:		extendar.h, error.h, smartptr.h, smartdcl.h, sptr_ndb.h,
					myalloc.h, std_tmpl.h, sort_lab.h, stdtype.h, simplex.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SC_TREE_H__
#define __SC_TREE_H__

#include <stdio.h>
#include <assert.h>

#ifndef __EXTENDAR_H__
#	include "extendar.h"
#endif
#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif


class Periods;
class ChangeLP;
class SolvableLP;


class ScenarioTree
{
private:
	struct ScStruct
	{
		Int_T scen, block;
		Real_T prob;
	};

	Int_T ScenNum, BlockNum, PeriodNum;
	SortedArrayOfLabels ScenLab;
	ExtendArray<ExtendArray<ScStruct> *> Structure;

	ExtendArray<ExtendArray<ChangeLP *> *> StochasticData;

private:
	void Clean( void );
	void PrintInfo( void );

	Bool_T CheckTree( void );

public:
	ScenarioTree( void );
	~ScenarioTree( void );

	Bool_T ReadScenarioFile( FILE *fp, const char *fname,
		const SolvableLP &lp, const Periods &per, VerbLevel Verbosity );

	Bool_T ConvertToTwoStage( SolvableLP &lp, const Periods &periods,
		Int_T CutAt, FILE *fpCOR, FILE *fpSTO, FILE *fpTIM,
		VerbLevel Verbosity ) const;
	void ApplyScenario( SolvableLP &lp, Int_T scen, Int_T period ) const;

private:
	//--------------------------------------------------------------------------
	//	Functions for problem output.
	//
	Bool_T OutputCORE( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T CutAt, VerbLevel Verbosity ) const;
	Bool_T OutputTIME( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T CutAt, VerbLevel Verbosity ) const;
	Bool_T OutputSTOCH( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T CutAt, VerbLevel Verbosity ) const;

	Bool_T OutputROWS( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T ScenMax ) const;
	Bool_T OutputCOLUMNS( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T ScenMax ) const;
	Bool_T OutputCOLUMNSFromSubTree( FILE *fp, SolvableLP &lp,
		const Periods &periods, Int_T p0, Int_T p1, Int_T s0, Int_T s1,
		Int_T CutAt ) const;
	Bool_T OutputRHS( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T ScenMax ) const;
	Bool_T OutputBOUNDS( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T ScenMax ) const;

	Bool_T OutputStochVectors( FILE *fp, Int_T p0, Int_T s0, Int_T s1 ) const;
	Bool_T OutputStochMatrix( FILE *fp, SolvableLP &lp, const Periods &periods,
		Int_T p0, Int_T s0, Int_T s1 ) const;

	Real_T TotalProbability( Int_T scen, Int_T period, Int_T CutAt ) const;
	Int_T HowManyNodeDuplicates( Int_T period ) const;

private:
	//--------------------------------------------------------------------------
	//	Static data and member functions used for scenario file reading.
	//
	static Bool_T ParseError;
	static const SolvableLP *LP;
	static const Periods *PER;

	//
	//	Scenario file parsing.
	//
	Bool_T ReadFile( const SolvableLP &lp, const Periods &per );
	Bool_T GetSections( void );
	Bool_T GetScenarioSection( void );
	Bool_T GetScenarioIndicator( Int_T &AttachAtPeriod );
	Bool_T GetScenarioEntry( Int_T AttachAtPeriod );

	void IgnoreTillEndOfSection( void );

	//
	//	Semantic actions.
	//
	void AddNewScenario( const char *ScLab, const char *RootLab,
		const char *PeriodLab, Real_T prob, Int_T &AttachAtPeriod );
	void AddEntry( const char *Column, const char *Row, Real_T Val,
		Int_T AttachAtPeriod );

	//
	//	Computing conditional probabilities in the tree.
	//
	void CalculateProbabilities( void );
};


/*------------------------------------------------------------------------------

	ScenarioTree::~ScenarioTree( void )

PURPOSE:
	Class "ScenarioTree" destructor. Calls the "Clean" function, which
deallocates all memory and resets the object to its initial (empty) state.

PARAMETERS:
	None.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ScenarioTree::~ScenarioTree( void )
{ Clean(); }


#endif
