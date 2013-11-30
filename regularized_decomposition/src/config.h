/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	config.h
CREATED:			1994.12.28
LAST MODIFIED:		1995.01.07

DEPENDENCIES:		option.h, smartptr.h, stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/


#ifndef __CONFIG_H__
#define __CONFIG_H__

#ifndef __OPTION_H__
#	include "option.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

//==============================================================================
//
//	Class "Config" declaration.
//

class Config
{
private:
	Bool_T valid;
	Int_T len, maxlen;
	Array<OptionBase *> options;

public:
	Config( void );
	~Config( void );

	Bool_T ReadArguments( int argc, char *argv[] );

	void AddOption( OptionBase *o );
};

//------------------------------------------------------------------------------
//
//	Inline functions' definitions.
//

inline
Config::Config( void )
	: valid( True ), len( 0 ), maxlen( 10 ), options( maxlen )
{
	options.Fill( NULL, maxlen );
}


inline
Config::~Config( void )
{}


//
//	End of class "Config" declaration.
//
//==============================================================================

#endif
