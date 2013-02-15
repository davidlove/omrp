/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	option.cpp
CREATED:			1994.12.28
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		option.h, stdtype.h, error.h, 
					<stdlib.h>, <assert.h>, <string.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

------------------------------------------------------------------------------*/

#include <string.h>
#include <assert.h>

#ifndef __OPTION_H__
#	include "option.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __STRDUP_H__
#	include "strdup.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif


#if defined( explicit_templates )
	template class SmartPointerBase<char *>;
	template class Array<char *>;
	template class Ptr<char *>;

	template char **MALLOC( char **& Table, size_t len );
	template char **REALLOC( char **& Table, size_t len );
	template void FREE( char **& Table );
#endif

//------------------------------------------------------------------------------
//
//	Class "SimpleOption" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

/*------------------------------------------------------------------------------

	SimpleOption::SimpleOption( const char *label, void (*action)( void ),
		int cnt )

PURPOSE:
	Class "SimpleOption" constructor.

PARAMETERS:
	const char *label
		Option label; e.g. for option "-fexternal-templates" label should be
		"fexternal-templates".

	void (*action)( void )
		Pointer of a function to be called upon recognition of this option.
		If NULL, the functio will not be called.

	int cnt
		How many times the option may be repeated; 0 if there is no limit.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SimpleOption::SimpleOption( const char *label, void (*action)( void ), // )
	int cnt )
	: RepeatCount( 0 ), MaxRepeatCount( cnt ), Label( NULL ), Action( action )
{
	assert( cnt >= 0 );
	assert( label != NULL );

	Label = DuplicateString( label );

	if( !Label )
		FatalError( "Out of memory." );
}


/*------------------------------------------------------------------------------

	virtual
	Bool_T SimpleOption::Try( int &argc, char **&argv )

PURPOSE:
	If the current option passed via "argc"/"argv" pair of arguments is
compatiblie with option "*this" return "True", otherwise return "False".

	If option is recognized, its corresponding "action" function is called.

PARAMETERS:
	int &argc
	char **&argv
		This argument pair is interpreted as if it was passed to the "main"
		function. It's passed by reference so that modifications are seen after
		the function returns. "argc" is the numer of arguments left to be
		processed; "argv" is an array of strings holding the actual program
		parameters.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	If option was succesfully recognized, its match count is incremented and
"argc"/"argv" are advanced.

------------------------------------------------------------------------------*/

Bool_T SimpleOption::Try( int &argc, char **&argv )
{
	assert( argc > 0 );
	assert( argv != NULL );
	assert( *argv != NULL );

	if( **argv == '-' && strcmp( (*argv)+1, Label ) == 0 )
	{
		RepeatCount++;
		if( MaxRepeatCount > 0 && RepeatCount > MaxRepeatCount )
		{
			Error( "Option '%s' repeated too many times.", *argv );

			argc--;
			argv++;

			return False;
		}
		else
		{
			if( Action ) Action();

			argc--;
			argv++;

			return True;
		}
	}
	else
		return False;
}


//------------------------------------------------------------------------------
//
//	Class "Argument" non-inline function implementation.
//
//------------------------------------------------------------------------------


/*------------------------------------------------------------------------------

	Bool_T Argument::Try( int &argc, char **&argv )

PURPOSE:
	Accepts a non-emty argument and passes it to the ArgumentAction()
function.

PARAMETERS:
	int &argc
	char **&argv
		This argument pair is interpreted as if it was passed to the "main"
		function. It's passed by reference so that modifications are seen after
		the function returns. "argc" is the numer of arguments left to be
		processed; "argv" is an array of strings holding the actual program
		parameters.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	If option was succesfully recognized, its match count is incremented and
	"argc"/"argv" are advanced.

------------------------------------------------------------------------------*/

Bool_T Argument::Try( int &argc, char **&argv )
{
	assert( argc > 0 );
	assert( *argv != NULL );

	if( **argv == '-' ) return False;

	if( MaxRepeatCount > 0 && ++RepeatCount > MaxRepeatCount )
	{
		Error( "Too many options without arguments: stopped at %s.", *argv );

		argc--;
		argv++;

		return False;
	}

	ArgumentAction( *argv );
	argc--;
	argv++;

	return True;
}


//------------------------------------------------------------------------------
//
//	Class "OptionWithArgument" non-inline functions' implementations.
//
//------------------------------------------------------------------------------


/*------------------------------------------------------------------------------

	OptionWithArgument::OptionWithArgument( const char *label,
		void (*arg_action)( const char *arg ), int cnt )

PURPOSE:
	Class "OptionWithArgument" constructor.

PARAMETERS:
	const char *label
		Option label; e.g. for option "-mps_in" label should be "mps_in".

	void (*arg_action)( const char *arg )
		Pointer of a function to be called upon recognition of this option.
		If NULL, the functio will not be called. Function will be called with
		the option argument.

	int cnt
		How many times the option may be repeated; 0 if there is no limit.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

OptionWithArgument::OptionWithArgument( const char *label, // )
	void (*arg_action)( const char *arg ), int cnt )
	: SimpleOption( label, NULL, cnt ), ArgumentAction( arg_action )
{}


/*------------------------------------------------------------------------------

	Bool_T OptionWithArgument::Try( int &argc, char **&argv )

PURPOSE:
	If the current option passed via "argc"/"argv" pair of arguments is
	compatiblie with option "*this" return "True", otherwise return "False".

	If option is recognized, its corresponding "arg_action" function is called.

PARAMETERS:
	int &argc
	char **&argv
		This argument pair is interpreted as if it was passed to the "main"
		function. It's passed by reference so that modifications are seen after
		the function returns. "argc" is the numer of arguments left to be
		processed; "argv" is an array of strings holding the actual program
		parameters.

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	If option was succesfully recognized, its match count is incremented and
	"argc"/"argv" are advanced.

------------------------------------------------------------------------------*/

Bool_T OptionWithArgument::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;
	
	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	if( ArgumentAction ) ArgumentAction( *argv );

	argc--;
	argv++;

	return True;
}


//------------------------------------------------------------------------------
//
//	Class "StoredArgument" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

void StoredArgument::SetDefault( const char *c )
{
	if( arg ) delete (char *)arg;

	if( c != NULL )
	{
		char *tmp = new char[ strlen( c ) + 1 ];
		strcpy( tmp, c );
		arg = c;
	}
}


Bool_T StoredArgument::Try( int &argc, char **&argv )
{
	assert( argc > 0 );
	assert( *argv != NULL );

	if( **argv == '-' ) return False;

	if( defined )
	{
		Error( "Too many arguments: stopped at %s.", *argv );
		return False;
	}

	char *tmp = new char[ strlen( *argv ) + 1 ];

	if( tmp == NULL ) FatalError( "Out of memory." );

	strcpy( tmp, *argv );
	arg = tmp;
	defined = True;

	argc--;
	argv++;

	return True;
}


//------------------------------------------------------------------------------
//
//	Class "OptionWithStoredArgument" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

Bool_T OptionWithStoredArgument::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;

	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	return StoredArgument::Try( argc, argv );
}


//------------------------------------------------------------------------------
//
//	Class "MultiStateOption" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

MultiStateOption::MultiStateOption( const char *label, int NumberOfStates, // )
	const char **labels, const int *states, Bool_T repeat )
	: SimpleOption( label, (voidF_void)NULL, repeat ? 0 : 1 ),
	n( NumberOfStates ),
	value( NumberOfStates, 0 ), str( NumberOfStates, (char *)NULL )
{
	assert( n > 0 );
	assert( labels != NULL );
	assert( states != NULL );

	for( int i = 0; i < n; i++ )
	{
		assert( *labels[i] != '\0' );

		value[i]		= states[i];
		str[i]			= new char[ strlen( labels[i] ) + 1 ];

		if( str[i] == NULL ) FatalError( "Out of memory." );
		strcpy( str[i], labels[i] );
	}
}


MultiStateOption::~MultiStateOption( void )
{
	for( int i = 0; i < n; i++ )
		if( str[i] ) delete str[i];
}


Bool_T MultiStateOption::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;
	
	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	Bool_T found = False;
	for( int i = 0; i < n; i++ )
		if( strcmp( *argv, str[i] ) == 0 )
		{
			state = value[i];
			found = True;
			break;
		}

	if( found )
	{
		argc--;
		argv++;

		return True;
	}
	else
	{
		Error( "Unrecognized argument of option '-%s': '%s'.", Label, *argv );
		return False;
	}
}


//------------------------------------------------------------------------------
//
//	Class "FloatValueOption" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

Bool_T FloatValueOption::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;
	
	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	value = atof( *argv );

	if( Limits && ( value < l || value > u ) )
	{
		Error( "Invalid numeric value '%g' (should be in range <%g, %g>).",
			value, l, u );
		return False;
	}

	argv++;
	argc--;

	return True;
}


//------------------------------------------------------------------------------
//
//	Class "BitmapOption" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

Bool_T BitmapOption::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;
	
	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	char *arg = *argv;

	argc--;
	argv++;

	//--------------------------------------------------------------------------
	//	Loop on the composite option argument. The syntax is:
	//	<str_i1> [ + <str_i2> [ + .... ] ]
	//	(without any seperating spaces). <str_i1>, <str_i2> are supposed to be
	//	the strings defined when the object is constructed.
	//
	while( *arg )
	{
		Bool_T found = False;
		for( int i = 0; i < n; i++ )
			if( strncmp( arg, str[i], strlen( str[i] ) ) == 0 )
			{
				state |= value[i];
				arg += strlen( str[i] );
				found = True;
				break;
			}

		if( !found )
		{
			char *end = NULL;
			for( end = arg; *end && *end != '+'; end++ )
				;

			Error( "Unrecognized flag '%.*s' in option '%s'.", end-arg, arg,
				Label );
			return False;
		}

		if( *arg )
		{
			if( *arg != '+' )
			{
				Error( "Invalid character in composite option; '+' expected." );
				return False;
			}

			arg++;
			if( *arg == '\0' )
				Warning( "Possible error in composite option: "
					"a space after '+'." );
		}
	}

	return True;
}


//------------------------------------------------------------------------------
//
//	Class "ON_OFF_Option" non-inline functions' implementations.
//
//------------------------------------------------------------------------------

Bool_T ON_OFF_Option::Try( int &argc, char **&argv )
{
	if( !SimpleOption::Try( argc, argv ) )
		return False;
	
	if( argc <= 0 )
	{
		Error( "Unexpected end of argument list after option '%s'.", argv[-1] );
		return False;
	}

	assert( *argv != NULL );

	if( strcmp( *argv, "on" ) == 0 )
		on = True;
	else if( strcmp( *argv, "off" ) == 0 )
		on = False;
	else
	{
		Error( "Option '%s' must be followed by 'on' or 'off'.", argv[-1] );
		return False;
	}

	argc--;
	argv++;

	return True;
}
