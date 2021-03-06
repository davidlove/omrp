/*------------------------------------------------------------------------------
MODULE TYPE:		Main file header
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	main.h
CREATED:			1994.08.17
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		stdtype.h, simplex.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __MAIN_H__
#define __MAIN_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif


//==============================================================================
//
//	Function prototypes.
//

class Solver;
class TimeInfo;
class StochSolution;

void run( int argc, char *argv[] );

void PrintHelpScreen( const char *ProgName );
void PrintCopyright( void );
void PrintSubproblemStatistics( const Solver &s );
void PrintSolution( const char *file, const StochSolution *sol );
void PrintTimings( TimeInfo &TI, bool first );

struct DecompOptions;
Bool_T ParseArguments( int argc, char *argv[], DecompOptions &DecOpt );

//
//	End of function prototypes.
//
//==============================================================================


struct DecompOptions
{
public:
	enum { FILE_NAME_LEN = 512 };
	typedef char FileName[FILE_NAME_LEN+1];

	Int_T ScenNum;
	Bool_T AllScen;

	FileName CoreFile, TimeFile, StochFile, SolutionFile;

	RD_SubproblemManager::RestartMode Restart;

	VerbLevel Verbosity;

	PricingScheme Pricing;

	Bool_T DoCrash;
	Real_T InitPen;

        // David Love -- Added the parameter gamma
        Int_T NonOverlap;
        // David Love -- Added to read replication from command line
        Int_T Replications;
        // David Love -- Added to read number nonoverlapping batches from command line
        Int_T NumNonOverBatches;

public:
	DecompOptions( void );	// Sets default values for all options.
};


inline
DecompOptions::DecompOptions( void )
	: ScenNum( 0 ), AllScen( False ), Restart( RD_SubproblemManager::SELF ),
	Verbosity( V_LOW ), Pricing( PRS_ASE ), DoCrash( True ), InitPen( 1 ),
        NonOverlap( 0 )
{ *CoreFile = *StochFile = *TimeFile = *SolutionFile = '\0'; }

//For L-shaped, set InitPen (1e-6)

#endif
