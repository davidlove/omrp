/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	config.cpp
CREATED:			1994.12.28
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		config.h, option.h, smartptr.h, stdtype.h, error.h
					<assert.h>, <stdlib.h>, <string.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

------------------------------------------------------------------------------*/


#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifndef __CONFIG_H__
#	include "config.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif

#if defined( explicit_templates )
	template class SmartPointerBase<OptionBase *>;
	template class Array<OptionBase *>;
	template class Ptr<OptionBase *>;

#ifdef NDEBUG
	template OptionBase **MALLOC( OptionBase **& Table, size_t len );
	template OptionBase **REALLOC( OptionBase **& Table, size_t len );
	template void FREE( OptionBase **& Table );
#endif
#endif


/*------------------------------------------------------------------------------

	Bool_T Config::ReadArguments( int argc, char *argv[] )

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

Bool_T Config::ReadArguments( int argc, char *argv[] )
{
	assert( argv != NULL );
	assert( argc >= 0 );
	assert( len > 0 );

	valid = True;

	while( argc > 0 )
	{
		Bool_T Recognized = False;

		for( Int_T i = 0; i < len; i++ )
			if( argc <= 0 )
				break;
			else if( options[i]->Try( argc, argv ) )
			{
				Recognized = True;
				break;
			}

		if( !Recognized )
		{
			valid = False;

			if( argc > 0 )
			{
				Error( "Unrecognized argument or option '%s'.", *argv );

				argc--;
				argv++;
			}
		}
	}

	return valid;
}


/*------------------------------------------------------------------------------

	void Config::AddOption( OptionBase *o )

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

void Config::AddOption( OptionBase *o )
{
	assert( o != NULL );

	if( len >= maxlen )
	{
		Int_T newmaxlen = Int_T( maxlen + 10 );

		options.Resize( newmaxlen );
		options.Fill( NULL, newmaxlen, maxlen );
		maxlen = newmaxlen;
	}

	options[ len++ ] = o;
}
