/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	option.h
CREATED:			1994.12.28
LAST MODIFIED:		1995.10.28

DEPENDENCIES:		stdtype.h, smartptr.h, my_defs.h
					<stdlib.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A number of classes representing different types of command line options.

------------------------------------------------------------------------------*/

#ifndef __OPTION_H__
#define __OPTION_H__

#include <stdlib.h>
#include <assert.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __MY_DEFS_H__
#	include "my_defs.h"
#endif


typedef void (*voidF_char_ptr)( const char *arg );
typedef void (*voidF_void)( void );

//------------------------------------------------------------------------------
//	Abstract class "OptionBase" declaration and inline functions'
//	implementations.
//
abstract class OptionBase
{
public:
	OptionBase( void );
	virtual ~OptionBase( void );

	virtual Bool_T Try( int &argc, char **&argv ) pure;
	virtual Bool_T Defined( void ) pure;
};


inline
OptionBase::OptionBase( void )
{}


inline
OptionBase::~OptionBase( void )
{}


//------------------------------------------------------------------------------
//	Class "SimpleOption" declaration and inline functions' implementations.
//
class SimpleOption : public virtual OptionBase
{
private:
	int RepeatCount, MaxRepeatCount;

protected:
	char *Label;

private:
	voidF_void Action;

public:
	SimpleOption( const char *label, voidF_void action = NULL, int cnt = 1 );
	virtual ~SimpleOption( void );

	virtual Bool_T Try( int &argc, char **&argv );
	virtual Bool_T Defined( void );
};


inline
SimpleOption::~SimpleOption( void )
	{ if( Label ) free( Label ); }


inline
Bool_T SimpleOption::Defined( void )
	{ return RepeatCount > 0 ? True : False; }


//------------------------------------------------------------------------------
//	Class "Argument" declaration and inline functions' implementations.
//
class Argument : public virtual OptionBase
{
private:
	int RepeatCount, MaxRepeatCount;
	voidF_char_ptr ArgumentAction;

public:
	Argument( voidF_char_ptr arg_action, int cnt = 1 );
	virtual ~Argument( void );

	virtual Bool_T Try( int &argc, char **&argv );
	virtual Bool_T Defined( void );
};


inline
Argument::Argument( voidF_char_ptr arg_action, int cnt )
	: RepeatCount( 0 ), MaxRepeatCount( cnt ),
	ArgumentAction( arg_action )
{
	assert( cnt >= 0 );
	assert( arg_action != NULL );
}


inline
Argument::~Argument( void )
{}


inline
Bool_T Argument::Defined( void )
	{ return RepeatCount > 0 ? True : False; }

//------------------------------------------------------------------------------
//	Class "StoredArgument" declaration and inline functions' implementations.
//
class StoredArgument : public virtual OptionBase
{
private:
	Bool_T defined;

public:
	const char *arg;

public:
	StoredArgument( void );
	virtual ~StoredArgument( void );

	virtual void SetDefault( const char *c );
	virtual Bool_T Try( int &argc, char **&argv );
	virtual Bool_T Defined( void );
};


inline
StoredArgument::StoredArgument( void )
	: defined( False ), arg( NULL )
	{}


inline
StoredArgument::~StoredArgument( void )
	{ if( arg ) delete (char *)arg; }


inline
Bool_T StoredArgument::Defined( void )
	{ return defined; }


//------------------------------------------------------------------------------
//	Class "OptionWithArgument" declaration and inline functions'
//	implementations.
//
class OptionWithArgument : public SimpleOption
{
private:
	voidF_char_ptr ArgumentAction;

public:
	OptionWithArgument( const char *label, voidF_char_ptr arg_action,
		int cnt = 1 );
	virtual ~OptionWithArgument( void );

	virtual Bool_T Try( int &argc, char **&argv );
};


inline
OptionWithArgument::~OptionWithArgument( void )
{}


//------------------------------------------------------------------------------
//	Class "OptionWithStoredArgument" declaration and inline functions'
//	implementations.
//
class OptionWithStoredArgument : public SimpleOption, public StoredArgument
{
public:
	OptionWithStoredArgument( const char *label ); 
	virtual ~OptionWithStoredArgument( void );

	virtual Bool_T Try( int &argc, char **&argv );
	virtual Bool_T Defined( void );
};


inline
OptionWithStoredArgument::OptionWithStoredArgument( const char *label )
	: SimpleOption( label ), StoredArgument()
	{}


inline
OptionWithStoredArgument::~OptionWithStoredArgument( void )
	{}


inline
Bool_T OptionWithStoredArgument::Defined( void )
	{ return SimpleOption::Defined(); }


//------------------------------------------------------------------------------
//	Class 'MultiStateOption' declaration and inline function's definitions.
//
class MultiStateOption : public SimpleOption
{
protected:
	int n;
	Array<int> value;
	Array<char *> str;

public:
	int state;

public:
	MultiStateOption( const char *label, int NumberOfStates,
		const char **labels, const int *states, Bool_T repeat = False );
	virtual ~MultiStateOption( void );

	virtual void SetDefault( int v );
	virtual Bool_T Try( int &argc, char **&argv );
};


inline
void MultiStateOption::SetDefault( int v )
	{ state = v; }


//------------------------------------------------------------------------------
//	Class 'FloatValueOption' declaration and inline function's definitions.
//
class FloatValueOption : public SimpleOption
{
public:
	double value;

private:
	Bool_T Limits;
	double l, u;

public:
	FloatValueOption( const char *label, double def = 0.0 );
	void SetLimits( double lo, double up );
	virtual ~FloatValueOption( void );

	virtual Bool_T Try( int &argc, char **&argv );
};


inline
FloatValueOption::FloatValueOption( const char *label, double def )
	: SimpleOption( label ), value( def ), Limits( False ), l( 0.0 ), u( 0.0 )
	{}


inline
FloatValueOption::~FloatValueOption( void )
	{}


inline
void FloatValueOption::SetLimits( double lo, double up )
{
	assert( lo <= up );
	l = lo; u = up;
	Limits = True;
}

//------------------------------------------------------------------------------
//	Class 'BitmapOption' declaration and inline function's definitions.
//
class BitmapOption : public MultiStateOption
{
public:
	BitmapOption( const char *label, int NumberOfStates,
		const char **labels, const int *states );
	virtual ~BitmapOption( void );

	virtual Bool_T Try( int &argc, char **&argv );
};


inline
BitmapOption::BitmapOption( const char *label, int NumberOfStates,  // )
	const char **labels, const int *states )
	: MultiStateOption( label, NumberOfStates, labels, states, True )
{}


inline
BitmapOption::~BitmapOption( void )
	{}


//------------------------------------------------------------------------------
//	Class 'ON_OFF_Option' declaration and inline function's definitions.
//
class ON_OFF_Option : public SimpleOption
{
public:
	Bool_T on;

public:
	ON_OFF_Option( const char *label, Bool_T def = False );
	virtual ~ON_OFF_Option( void );

	virtual Bool_T Try( int &argc, char **&argv );
};


inline
ON_OFF_Option::ON_OFF_Option( const char *label, Bool_T def )
	: SimpleOption( label ), on( def )
{}


inline
ON_OFF_Option::~ON_OFF_Option( void )
{}


#endif
