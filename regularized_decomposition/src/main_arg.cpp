/*------------------------------------------------------------------------------
MODULE TYPE:		Argument parsing for the main module.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	main_arg.cpp
CREATED:			1995.08.25
LAST MODIFIED:		1995.11.15

DEPENDENCIES:		error.h, print.h, stdtype.h, simplex.h, std_tmpl.h, 
					<stdio.h>, <stdlib.h>, <ctype.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/


#include <assert.h>
#include <ctype.h>


#ifndef __MAIN_H__
#	include "main.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __MAST_SOL_H__
#	include "mast_sol.h"
#endif
#ifndef __CONFIG_H__
#	include "config.h"
#endif


//==============================================================================
//	Static function prototypes.
//
static void SetScenarioNumber( const char *argument );
static void SetCoreFile( const char *argument );
static void SetTimeFile( const char *argument );
static void SetStochFile( const char *argument );
static void SetSolutionFile( const char *argument );
static void SetRestartMode( const char *argument );
static void SetProblem( const char *argument );
static void SetVerbosity( const char *argument );
static void SetPricingMode( const char *argument );
static void SetCrashType( const char *argument );
static void SetInitialPenalty( const char *argument );
// David Love -- get gamma from the input
static void SetNonOverlap( const char *argument );
// David Love -- get number of replications from input
static void SetReplications( const char *argument );
// David Love -- Get number of nonoverlapping batches
static void SetNonOverBatches( const char *argument );

//
//	Static data.
//
static Bool_T valid = True;
static DecompOptions *DecOpt = NULL;
//
//	End of the static function prototypes.
//==============================================================================


Bool_T ParseArguments( int argc, char *argv[], DecompOptions &opt )
{
	DecOpt = &opt;

	Config Cfg;

	valid = True;

	Cfg.AddOption( new OptionWithArgument(	"cor",		SetCoreFile ) );
	Cfg.AddOption( new OptionWithArgument(	"tim",		SetTimeFile ) );
	Cfg.AddOption( new OptionWithArgument(	"sto",		SetStochFile ) );
	Cfg.AddOption( new OptionWithArgument(	"txt_sol",	SetSolutionFile ) );
	Cfg.AddOption( new OptionWithArgument(	"problem",	SetProblem ) );
	Cfg.AddOption( new Argument(						SetProblem ) );
	Cfg.AddOption( new OptionWithArgument(	"s",		SetScenarioNumber ) );
	Cfg.AddOption( new OptionWithArgument(	"scen",		SetScenarioNumber ) );
	Cfg.AddOption( new OptionWithArgument(	"restart",	SetRestartMode ) );
	Cfg.AddOption( new OptionWithArgument(	"v",		SetVerbosity ) );
	Cfg.AddOption( new OptionWithArgument(	"pric",		SetPricingMode ) );
	Cfg.AddOption( new OptionWithArgument(	"crash",	SetCrashType ) );
	Cfg.AddOption( new OptionWithArgument(	"penalty",	SetInitialPenalty ) );
        // David Love -- Get the value of gamma NonOverlap
	Cfg.AddOption( new OptionWithArgument(	"g",	        SetNonOverlap ) );
        // David Love -- Get the number of replications
	Cfg.AddOption( new OptionWithArgument(	"r",	        SetReplications ) );
        // David Love -- Get the number of non-overlapping batches
	Cfg.AddOption( new OptionWithArgument(	"k",	        SetNonOverBatches ) );
        

	if( !Cfg.ReadArguments( argc, argv ) )
		valid = False;

	if( !*(DecOpt->CoreFile) || !*(DecOpt->TimeFile) || !*(DecOpt->StochFile) )
	{
		valid = False;
		Error( "Three input files required." );
	}

	if( DecOpt->ScenNum == 0 )
	{
		valid = False;
		Error( "Number of scenarios to generate not specified." );
	}

        // David Love -- If no gamma is provided, set gamma = m
        if( DecOpt->NonOverlap == 0 )
        {
           DecOpt->NonOverlap = DecOpt->ScenNum;
        }

        // David Love -- Check if value of gamma makes sense
        if( DecOpt->NonOverlap < 0 || DecOpt->NonOverlap > DecOpt->ScenNum )
        {
           valid = False;
           Error( "NonOverlap, gamma, has incorrect value." );
        }

        // David Love -- If no replication count is provided, rep = 1000
        if( DecOpt->Replications == 0 )
        {
           DecOpt->Replications = 1000;
        }
        
        // David Love -- If number nonoverlapping batches not specified, set to 30
        if( DecOpt->NumNonOverBatches == 0 )
        {
           DecOpt->NumNonOverBatches = 30;
        }

        // David Love -- Check if value of replication makes sense
        if( DecOpt->Replications < 0 )
        {
           valid = False;
           Error( "Number of repitions must be positive" );
        }

	assert( !valid || ( DecOpt->ScenNum == -1 && DecOpt->AllScen ) ||
		( DecOpt->ScenNum > 0 && !DecOpt->AllScen ) );

	return valid;
}


static void SetScenarioNumber( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( strcmp( argument, "all" ) == 0 )
	{
		DecOpt->ScenNum = -1;
		DecOpt->AllScen = True;
	}
	else
	{
		DecOpt->ScenNum = (Int_T) atoi( argument );
		DecOpt->AllScen = False;

		if( DecOpt->ScenNum <= 0 )
		{
			Error( "Invalid number of scenarios (non-positive or too large)." );
			valid = False;
		}
	}
}


static void SetCoreFile( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( DecOpt->CoreFile[0] == '\0' )
	{
		strncpy( DecOpt->CoreFile, argument, DecompOptions::FILE_NAME_LEN );
		DecOpt->CoreFile[DecompOptions::FILE_NAME_LEN] = '\0';
	}
	else
	{
		Error( "More than one core file specified." );
		valid = False;
	}
}


static void SetTimeFile( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( DecOpt->TimeFile[0] == '\0' )
	{
		strncpy( DecOpt->TimeFile, argument, DecompOptions::FILE_NAME_LEN );
		DecOpt->TimeFile[DecompOptions::FILE_NAME_LEN] = '\0';
	}
	else
	{
		Error( "More than one time file specified." );
		valid = False;
	}
}


static void SetStochFile( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	strncpy( DecOpt->StochFile, argument, DecompOptions::FILE_NAME_LEN );
	DecOpt->StochFile[DecompOptions::FILE_NAME_LEN] = '\0';
}


static void SetSolutionFile( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( DecOpt->SolutionFile[0] == '\0' )
	{
		strncpy( DecOpt->SolutionFile, argument, DecompOptions::FILE_NAME_LEN );
		DecOpt->SolutionFile[DecompOptions::FILE_NAME_LEN] = '\0';
	}
	else
	{
		Error( "More than one solution file specified." );
		valid = False;
	}
}


static void SetRestartMode( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( strcmp( argument, "tree" ) == 0 )
		DecOpt->Restart = RD_SubproblemManager::TREE;
	else if( strcmp( argument, "random" ) == 0 )
		DecOpt->Restart = RD_SubproblemManager::RANDOM;
	else if( strcmp( argument, "self" ) == 0 )
		DecOpt->Restart = RD_SubproblemManager::SELF;
	else
		Warning( "Unrecognized subproblem restart mode. Ignoring." );
}


static void SetProblem( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( DecOpt->StochFile[0] )
	{
		Error( "More than one stochastic data file specified." );
		valid = False;
	}
	if( DecOpt->TimeFile[0] )
	{
		Error( "More than one time file specified." );
		valid = False;
	}
	if( DecOpt->CoreFile[0] )
	{
		Error( "More than one core file specified." );
		valid = False;
	}

	//--------------------------------------------------------------------------
	//	Strip possible .TIM, .STO or .COR extension.
	//
	Int_T len = strlen( argument );
	if( len > 4 )
	{
		char ext[5];

		strcpy( ext, argument + len - 4 );
		for( int i=0; i < 4; i++ )
			if( islower( ext[i] ) ) ext[i] = toupper( ext[i] );

		if( strcmp( ext, ".TIM" ) == 0 ||
			strcmp( ext, ".STO" ) == 0 ||
			strcmp( ext, ".COR" ) == 0 )
			((char *)argument)[len-4] = '\0';
	}

	if( valid )
	{
		strncpy( DecOpt->CoreFile, argument, DecompOptions::FILE_NAME_LEN-4 );
		strncpy( DecOpt->TimeFile, argument, DecompOptions::FILE_NAME_LEN-4 );
		strncpy( DecOpt->StochFile, argument, DecompOptions::FILE_NAME_LEN-4 );
		DecOpt->CoreFile[DecompOptions::FILE_NAME_LEN-4] = '\0';
		strcat( DecOpt->CoreFile,	".cor" );
		strcat( DecOpt->TimeFile,	".tim" );
		strcat( DecOpt->StochFile, 	".sto" );
		DecOpt->CoreFile[DecompOptions::FILE_NAME_LEN] = '\0';
	}
}


static void SetVerbosity( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( strcmp( argument, "none" ) == 0 )
		DecOpt->Verbosity = V_NONE;
	else if( strcmp( argument, "low" ) == 0 )
		DecOpt->Verbosity = V_LOW;
	else if( strcmp( argument, "high" ) == 0 )
		DecOpt->Verbosity = V_HIGH;
	else
	{
		Error( "Unrecognized verbosity level: %s.", argument );
		valid = False;
	}
}


static void SetPricingMode( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( strcmp( argument, "rc" ) == 0 )
		DecOpt->Pricing = PRS_RC;
	else if( strcmp( argument, "se" ) == 0 )
		DecOpt->Pricing = PRS_SE;
	else if( strcmp( argument, "ase" ) == 0 )
		DecOpt->Pricing = PRS_ASE;
	else
	{
		Error( "Unrecognized pricing mode: %s.", argument );
		valid = False;
	}
}


static void SetCrashType( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	if( strcmp( argument, "on" ) == 0 )
	{
		DecOpt->DoCrash = True;
		DecOpt->InitPen = 1.0;
	}
	else if( strcmp( argument, "off" ) == 0 )
	{
		DecOpt->DoCrash = False;
		DecOpt->InitPen = 0.01;
	}
	else
	{
		Error( "Unrecognized crash type: %s.", argument );
		valid = False;
	}
}


static void SetInitialPenalty( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

	Real_T pen = atof( argument );

	if( pen <= 0.0 )
		Warning( "Non-positive penalty ignored." );
	else if( pen <= PENALTY_LO )
		Warning( "Penalty value to small; adjusted to %g.",
			(double)PENALTY_LO );
	else if( pen >= PENALTY_HI )
		Warning( "Penalty value to large; adjusted to %g.",
			(double)PENALTY_HI );
	else
		DecOpt->InitPen = pen;
}

// David Love -- Set the value of gamma
static void SetNonOverlap( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

        DecOpt->NonOverlap = (Int_T) atoi( argument );

        if( DecOpt->NonOverlap < 0 || DecOpt->NonOverlap > DecOpt->ScenNum )
        {
                Error( "Invalid number of nonoverlap." );
                valid = False;
        }
}

// David Love -- Set the number of replications
static void SetReplications( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

        DecOpt->Replications = (Int_T) atoi( argument );

        if( DecOpt->Replications < 0 )
        {
                Error( "Invalid number of replications." );
                valid = False;
        }
}

// David Love -- Set number of nonoverlapping batches
static void SetNonOverBatches( const char *argument )
{
	assert( DecOpt != NULL );
	assert( argument != NULL );

        DecOpt->NumNonOverBatches = (Int_T) atoi( argument );

        if( DecOpt->NumNonOverBatches < 0 )
        {
                Error( "Invalid number of nonoverlapping batches." );
                valid = False;
        }
}
