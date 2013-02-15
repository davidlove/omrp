/*------------------------------------------------------------------------------
MODULE TYPE:		Set of functions.
PROJECT CODE:		General purpose module.
PROJECT FULL NAME:	-----------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	error.cpp
CREATED:			1993.09.11
LAST MODIFIED:		1995.10.29

DEPENDENCIES:		error.h,
					<stdarg.h>, <stdio.h>, <stdlib.h>, <string.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	The module provides some simple error message output procedures. There are 3
distinct categories of errors:
*	warnings,
*	ordinary errors and
*	fatal errors.
	Warnings and errors result in displaying appropriate error messages and
updating their respective counters. There are only two differences between
warnings and errors:
*	message format (see below) and
*	possession of separate counters.
The counters may be read and reset; they are automatically incremented
when appropriate error message output functions are called.
	Fatal error function call displays an error message and immediately
terminates the program with a call to ANSI C function "abort" (semantics of
"abort" may be found in your local ANSI C reference manual). Its purpose is to
inform the user of some severe error that makes further execution of the
program impossible (or pointless).
	The error messages have the following format:

	error_type [(error_code)]: custom_message

Error_type may be any of: WARNING, ERROR, FATAL ERROR. Optional error_code is
a non-negative integer value. Custom_message is defined by the user.
	There are two ways of displaying error messages:
*	by a variable length argument list with syntax identical to that of "printf"
	ANSI C function,
*	by a non-negative integer number interpreted as an index into previously
	declared array of message format strings followed by a variable number of
	parameters.
	When the latter method is used, the argument passed as an index may be
checked (if macro NDEBUG is not defined). The messages by default are output to
standard error output stream (stderr). They may be redirected (at any time
during program execution) to some other stream, or suppressed.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void DeclareErrorMessages( int num, const char ** messages )
	void SetErrorOutputStream( FILE *fp = stderr )
	int Warning( int error_code, ... )
	int Warning( const char *format, ...)
	int Error( int error_code, ... )
	int Error( const char *format, ...)
	void FatalError( int error_code, ... )
	void FatalError( const char *format, ...)
	int WarningCount( void )
	int ErrorCount( void )
	void ResetWarningCount( void )
	void ResetErrorCount( void )
	void MaxWarnCount( int num = 0)
	void MaxErrCount( int num = 0 )

--------------------------------------------------------------------------------

STATIC TYPE DEFINITIONS:
	enum ErrorCode { EC_WARNING, EC_ERROR, EC_FATAL_ERROR }

STATIC DATA:
	static const char **MessageTable
		A table of error message format strings declared by the program.

	static int MessageTableLen
		Number of strings in the above table.
	
	static FILE *ErrorOutputStream
		Stream descriptor for message output, or NULL when message
		output is to be supressed.

	static int _WarningCount = 0
	static int _ErrorCount = 0
		Error and warning counters.

	static int _MaxWarningCount = 0
	static int _MaxErrorCount = 0
		Maximum allowed numbers of warnings and errors respectively.

STATIC FUNCTIONS:
	void ErrorMessageOutput( ErrorCode ec, const char *format, va_list params,
		int error_code = -1 )
	const char *GetFormatString( int msg_num )
	
------------------------------------------------------------------------------*/

//#define ERROR_ABORT

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>

#include "error.h"

//==============================================================================
//
//	Static type definitions.
//
enum ErrorCode { EC_WARNING, EC_ERROR, EC_FATAL_ERROR };
//
//	End of static type definitions.
//
//==============================================================================

//==============================================================================
//
//	Static data.
//

//------------------------------------------------------------------------------
//	Stream descriptor for error message output (may be NULL when
//	output is to supressed).
//
static FILE *ErrorOutputStream = stderr;

//------------------------------------------------------------------------------
//	Error message format string table (a pointer and table length). Each message
//	format string is assumed to be a zero terminated character string. Table is
//	initially empty.
//
static const char **MessageTable = NULL;
static int MessageTableLen = 0;

//------------------------------------------------------------------------------
//	Warning and error counters are held here. They are incremented by
//	respectively "Warning" and "Error" calls. They may be reset by
//	"ResetWarningCount" and "ResetErrorCount" functions. Their values are
//	returned by "Warning", "Error", "_WarningCount" and "_ErrorCount" calls.
//	Also "MaxWarnCount" / "MaxErrCount" calls reset their values.
//
static int _WarningCount = 0,
	_ErrorCount = 0;

//------------------------------------------------------------------------------
//	Maximum allowed numbers of warning / error messages are held here. If set to
//	zero - checking for maximum is disabled.
//
static int _MaxWarningCount = 0,
	_MaxErrorCount = 0;
//
//	End of static data.
//
//==============================================================================

//==============================================================================
//
//	Static functions' prototypes.
//
void ErrorMessageOutput( ErrorCode ec, const char *format, va_list params,
	int error_code = -1 );
const char *GetFormatString( int msg_num );
//
//	End of static functions' prototypes.
//
//==============================================================================

//==============================================================================
//
//	Public functions' definitions.
//

/*------------------------------------------------------------------------------

	void DeclareErrorMessages( int num, const char ** messages )

PURPOSE:
	This function may be called at any time during the program execution. It
declares (or redeclares) a table of string literals representing error /
warning messages. From the moment of this call and until the next call any
error message output function with a parameter specifying the message format
string number will result in displaying a message from this table (or a fatal
runtime error if no message with corresponding number exists).
	It is recommended that this function is called only once during the program
execution. Frequent use of this function may cause confusion and thus become a
source of programming errors.
	Before the first call with non-zero values of both num and messages the
table is empty and functions with single integral argument may not be used.

PARAMETERS:

	int num
		Specifies the number of messages (length of messages table). It should
		be greater than or equal to zero. Sub-zero value will generate a fatal
		runtime error; num == 0 will be interpreted as a request to forget the
		previously stored table and substitute it with an empty one (which will
		effectively disable all corresponding error message output functions
		using an integral argument). If num == 0, the other argument also has to
		be set to zero (as it is a pointer NULL symbolic value may be more
		appropriate). If one of these parameters is equal to zero and the other
		is not, a fatal runtime error is generated and the program aborts.

	const char ** messages
		Passes the pointer to a (preferably static and constant) table of
		strings representing error message format strings, which will later be
		interpreted as "printf"-style format strings. The table is assumed to
		contain at least num valid string entries. The table is not copied!  Thus
		any change done to the strings between declaration and error message
		output calls will be reflected in the message output. Taking advantage of
		this fact is strongly discouraged.

		All message format string pointers are required to be non-NULL.  If
		NDEBUG macro is not declared, this property is ensured during table
		declaration.

		If there is a need for more flexibility than this function together with
		appropriate error message output functions provide, the alternative error
		output syntax should be used (see below: message output functions with
		explicitly given format string).

		NULL pointer may only be passed together with num == 0. For the meaning
		of such combination of argument values - see description of the num
		parameter above.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Stores the "message" pointer together with the number of messages "num", or
declares the message format string list empty.

------------------------------------------------------------------------------*/

void DeclareErrorMessages( int num, const char ** messages )
{
	assert( num > 0 );
	assert( messages != NULL );

	if( !num && !messages )
	{
		MessageTable	= NULL;
		MessageTableLen	= 0;
	}
	else
	{
		//----------------------------------------------------------------------
		//	Check (optionally) if message pointers are non-NULL.
		//
#ifndef NDEBUG
		for( int i = 0; i < num; i++ )
			assert( messages[ num ] != NULL );
#endif
		//----------------------------------------------------------------------
		//	Store the table address and its length.
		//
		MessageTable	= messages;
		MessageTableLen	= num;
	}
}


/*------------------------------------------------------------------------------

	void SetErrorOutputStream( FILE *fp )

PURPOSE:
	This function redirects error message output to specified output stream fp.
The default stream corresponds to ANSI C's standard error (stderr). If fp ==
NULL, message output will be suppressed (but all the side effects of error
message output functions as well as parameter checking inside these functions
will remain unaffected).  "FatalError" functions are not affected by this
function call with fp == NULL. Consequently any erroneous calls to functions
from the module that would normally cause the program to be aborted with some
fatal error message, continue to stop the program with the same error message
output to the last defined stream or stderr if all other error / warning
message output was suppressed.
	The stream descriptor pointer is only stored, and not duplicated.
Therefore it is possible to pass a pointer to some stream and later close it.
If this should be detected by the message output functions, a fatal error is
generated. The responsibility to keep the stream open remains with the user.

PARAMETERS:
	FILE *fp = stderr
		Passes a pointer to C-type output stream fp. If omitted, stderr is
		assumed.  Any valid stream may be used here. NULL suppresses message
		output until a later call with another value of the file descriptor
		pointer.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Stores the stream descriptor pointer.

------------------------------------------------------------------------------*/

void SetErrorOutputStream( FILE *fp )
{
	ErrorOutputStream = fp;
}


/*------------------------------------------------------------------------------

	int Warning( int error_code, ... )

PURPOSE:
	Outputs a predefined warning message to a stream. The message format is:

	WARNING (error_code): custom_message

The message custom_message is constructed by applying "printf"-style formatting
to the format string (read from the previously declared message table) and the
variable length argument list (see: function "DeclareErrorMessages"). If the
error_code number is below zero or greater than the predeclared number of
messages in the table, a fatal error is generated . Due to use of ellipsis,
only limited argument checking is available.
	The message output may be suppressed by a call to "SetErrorOutputStream"
with NULL stream descriptor pointer. The side effects and error checking inside
this function remain unaffected.

PARAMETERS:
	error_code
	Passes the message format string number. It is interpreted as an index into
messages array - see function "DeclareErrorMessages" description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	Returns current warning counter value (updated by this call).

SIDE EFFECTS:
	Increments warning counter.

------------------------------------------------------------------------------*/

int Warning( int error_code, ... )
{
	va_list params;

	va_start( params, error_code );
	ErrorMessageOutput( EC_WARNING, GetFormatString( error_code ), params,
		error_code );
	va_end( params );

	_WarningCount++;
	if( _MaxWarningCount && _WarningCount >= _MaxWarningCount )
		FatalError( "Too many warning messages." );
	return _WarningCount;
}


/*------------------------------------------------------------------------------

	int Warning( const char *format, ...)

PURPOSE:
	Outputs a custom warning message specified by a "printf"-style format
string format and a variable length argument list referred to with ellipsis.
The message format is:

	WARNING: custom_message

	Message output may be suppressed by a call to "SetErrorOutputStream" with
NULL stream descriptor pointer. Due to use of ellipsis, only limited argument
checking is available.

PARAMETERS:
	const char *format
	This is a format string identical to that of "printf" functions' family. In
fact one of these functions is used to actually output the message to the
output stream (predefined - stderr - or specified with a call of
"SetErrorOutputStream" function). For details of interpretation of this string
- see your C or C++ reference manual, "printf" function description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	Returns current warning counter value (updated by this call).

SIDE EFFECTS:
	Increments warning counter.

------------------------------------------------------------------------------*/

int Warning( const char *format, ...)
{
	va_list params;

	va_start( params, format );
	ErrorMessageOutput( EC_WARNING, format, params );
	va_end( params );

	_WarningCount++;
	if( _MaxWarningCount && _WarningCount >= _MaxWarningCount )
		FatalError( "Too many warning messages." );
	return _WarningCount;
}


/*------------------------------------------------------------------------------

	int Error( int error_code, ... )

PURPOSE:
	Outputs a predefined error message to a stream. The message format is:

	ERROR (error_code): custom_message

The message custom_message is constructed by applying "printf"-style formatting
to the format string read from the previously declared message table and the
variable length argument list (see:  "DeclareErrorMessages"). If the error_code
number is below zero or greater than the predeclared number of messages in the
table, a fatal error is generated. Due to use of ellipsis, only limited
argument checking is available.
	The message output may be suppressed by a call to "SetErrorOutputStream"
with NULL stream descriptor pointer. The side effects and error checking inside
this function remain unaffected.

PARAMETERS:
	error_code
	Passes the message format string number. It is interpreted as an index into
messages array - see function "DeclareErrorMessages" description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	Returns current error counter value (updated by this call).

SIDE EFFECTS:
	Increments error counter.

------------------------------------------------------------------------------*/

int Error( int error_code, ... )
{
	va_list params;

	va_start( params, error_code );
	ErrorMessageOutput( EC_ERROR, GetFormatString( error_code ), params,
		error_code );
	va_end( params );

	_ErrorCount++;
	if( _MaxErrorCount && _ErrorCount >= _MaxErrorCount )
		FatalError( "Too many error messages." );
	return _ErrorCount;
}


/*------------------------------------------------------------------------------

	int Error( const char *format, ...)

PURPOSE:
	Outputs a custom error message specified by a "printf"-style format string
format and a variable length argument list referred to with ellipsis. The
message format is:

	ERROR: custom_message

	Message output may be suppressed by a call to "SetErrorOutputStream" with
NULL stream descriptor pointer. Due to use of ellipsis, only limited argument
checking is available.

PARAMETERS:
	const char *format
	This is a format string identical to that of "printf" functions' family. In
fact one of these functions is used to actually output the message to the
output stream (predefined - stderr - or specified with a call of
"SetErrorOutputStream" function). For details of interpretation of this string
- see your C or C++ reference manual, "printf" function description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	Returns current error counter value (updated by this call).

SIDE EFFECTS:
	Increments error counter.

------------------------------------------------------------------------------*/

int Error( const char *format, ...)
{
	va_list params;

	va_start( params, format );
	ErrorMessageOutput( EC_ERROR, format, params );
	va_end( params );

	_ErrorCount++;
	if( _MaxErrorCount && _ErrorCount >= _MaxErrorCount )
		FatalError( "Too many error messages." );
	return _ErrorCount;
}


/*------------------------------------------------------------------------------

	void FatalError( int error_code, ... )

PURPOSE:
	Outputs a predefined fatal error message to previously declared stream (or
if error output was suppressed - to an ANSI C's standard error stream stderr),
then terminates the program with a call to ANSI C function "abort". On some
systems this will make post-mortem debugging possible. The message format is:

	FATAL ERROR (error_code): custom_message

The message custom_message is constructed by applying "printf"-style formatting
to the format string read from the previously declared message table and the
variable length argument list (see:  "DeclareErrorMessages"). On some systems
output produced by "abort" function may follow the custom message (see
description of "abort" function in your reference manual). If the error_code
number is below zero or greater than the predeclared number of messages in the
table, a fatal error message is substituted by "undefined internal error".  Due
to use of ellipsis, only limited argument checking is available.
	The message output is never suppressed by a call to "SetErrorOutputStream"
with NULL stream descriptor pointer. If the stream pointer is NULL, the message
is output to the standard error stream stderr.

PARAMETERS:
	error_code
	Passes the message format string number. It is interpreted as an index into
messages array - see function "DeclareErrorMessages" description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	None (this function never returns).

SIDE EFFECTS:
	None (one would not be able to observe them anyway).

------------------------------------------------------------------------------*/

void FatalError( int error_code, ... )
{
	va_list params;

	va_start( params, error_code );
	ErrorMessageOutput( EC_FATAL_ERROR, GetFormatString( error_code ), params,
		error_code );
	va_end( params );
}


/*------------------------------------------------------------------------------

	void FatalError( const char *format, ...)

PURPOSE:
	Outputs a custom fatal error message specified by a "printf"-style format
string format and a variable length argument list referred to with ellipsis to
previously declared stream (or if error output was suppressed - to an ANSI C's
standard error stream stderr), then terminates the program with a call to ANSI
C function "abort". On some systems this will make post-mortem debugging
possible. The message format is:

	FATAL ERROR: custom_message

	On some systems output produced by "abort" function may follow the custom
message (see description of "abort" function in your reference manual). The
message output is never suppressed by a call to "SetErrorOutputStream" with
NULL stream descriptor pointer. If the stream pointer is NULL, the message is
output to the standard error output stream stderr. Due to use of ellipsis, only
limited argument checking is available.

PARAMETERS:
	const char *format
	This is a format string identical to that of "printf" functions' family. In
fact one of these functions is used to actually output the message to the
output stream (predefined - stderr - or specified with a call of
"SetErrorOutputStream" function). For details of interpretation of this string
- see your C or C++ reference manual, "printf" function description.

	... (ellipsis)
	The meaning of the remaining actual parameters depends on the contents of
the format string (see above). For details - see your C or C++ reference
manual, "printf" function description.

RETURN VALUE:
	None (this function never returns).

SIDE EFFECTS:
	None (one would not be able to observe them anyway).

------------------------------------------------------------------------------*/

void FatalError( const char *format, ...)
{
	va_list params;

	va_start( params, format );
	ErrorMessageOutput( EC_FATAL_ERROR, format, params );
	va_end( params );
}


/*------------------------------------------------------------------------------

	int WarningCount( void )

PURPOSE:
	This function reveals current value of the warning counter. The value
returned is identical to that returned by the last call of one of overloaded
"Warning" functions, or is zero if:
*	none of the "Warning" functions was ever called or
*	the warning counter was reset by "ResetWarningCount" call and was not
	incremented by subsequent "Warning" calls.

PARAMETERS:
	None.

RETURN VALUE:
	Current warning counter value.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int WarningCount( void )
{
	return _WarningCount;
}


/*------------------------------------------------------------------------------

	int ErrorCount( void )

PURPOSE:
	This function reveals current value of the error counter. The value
returned is identical to that returned by the last call of one of overloaded
"Error" functions, or is zero if:
*	none of the "Error" functions was ever called or
*	the error counter was reset by "ResetErrorCount" call and was not
	incremented by subsequent "Error" calls.

PARAMETERS:
	None.

RETURN VALUE:
	Current error counter value.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int ErrorCount( void )
{
	return _ErrorCount;
}


/*------------------------------------------------------------------------------

	void ResetWarningCount( void )

PURPOSE:
	This function sets the warning counter to zero (regardless of its value).

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Sets warning counter to zero.

------------------------------------------------------------------------------*/

void ResetWarningCount( void )
{
	_WarningCount = 0;
}


/*------------------------------------------------------------------------------

	void ResetErrorCount( void )

PURPOSE:
	This function sets the error counter to zero (regardless of its value).

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Sets error counter to zero.

------------------------------------------------------------------------------*/

void ResetErrorCount( void )
{
	_ErrorCount = 0;
}

/*------------------------------------------------------------------------------

	void MaxWarnCount( int num )

PURPOSE:
	This function declares the maximum number of warning messages that may be
issued. It automatically resets the warning counter to zero. When the limit is
reached, an appropriate fatal error message is forced by subsequent "Warning"
call.
	The number may not be redeclared, that is you have to reset the maximum
counter before you can set it again.

PARAMETERS:
	int num
	When greater than zero, indicates the number of warning messages allowed.
When zero, clears the maximum warning number and resets the warning counter.
When less than zero, generates a runtime fatal error.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Sets maximum warning count to num and sets warning counter to zero, or (when
num is equal to zero) resets both counters.

------------------------------------------------------------------------------*/

void MaxWarnCount( int num )
{
	assert( num >= 0 );

	if( num > 0 && _MaxWarningCount == 0 )
	{
		_MaxWarningCount	= num;
		_WarningCount		= 0;
	}
	else if( num == 0 )
		_MaxWarningCount = _WarningCount = 0;
}

/*------------------------------------------------------------------------------

	void MaxErrCount( int num )

PURPOSE:
	This function declares the maximum number of error messages that may be
issued. It automatically resets error counter to zero. When the limit is
reached, an appropriate fatal error message is forced by subsequent "Warning"
call.
	The number may not be redeclared, that is you have to reset the maximum
counter before you can set it again.

PARAMETERS:
	int num
	When greater than zero, indicates the number of error messages allowed.
When zero, clears the maximum error number and resets the error counter.  When
less than zero, generates a runtime fatal error.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Sets maximum warning count to num and sets error counter to zero, or (when
num is equal to zero) resets both counters.

------------------------------------------------------------------------------*/

void MaxErrCount( int num )
{
	assert( num >= 0 );

	if( num > 0 && _MaxErrorCount == 0 )
	{
		_MaxErrorCount	= num;
		_ErrorCount		= 0;
	}
	else if( num == 0 )
		_MaxErrorCount = _ErrorCount = 0;
}

//
//	End of public functions' definitions.
//
//==============================================================================

//==============================================================================
//
//	Static functions' definitions.
//

/*------------------------------------------------------------------------------

	void ErrorMessageOutput( ErrorCode ec, const char *format, va_list params,
		int error_code )

PURPOSE:
	This function actually diplays the error message with appropriate format.
It is called by all "Warning", "Error" and "FatalError" functions. If
NDEBUG is not defined it checks if the file output is successful. It also
ensures that FatalError messages are displayed even when output stream pointer
"ErrorOutputStream" is NULL.
	The error messages have the following format:

	error_type [(error_code)]: custom_message

"Error_type" may be any of: WARNING, ERROR, FATAL ERROR. Optional "error_code"
is a non-negative integer value passed as an optional parameter.
"Custom_message" is made up of the format string and the variable length
parameter list - see parameters' descriptions below. Newline is automatically
appended to the message.

PARAMETERS:
	ErrorCode ec
		Error code (one of EC_WARNING, EC_ERROR, EC_FATAL_ERROR) using
		enumerated type. Compiler ensures that this parameter always takes one
		of these values. "ErrorCode" is a locally defined (file scope) type.

	const char *format
		Error message format string (later passed to the "vfprintf" function).
		Its' syntax is identical tho that of "printf". 

	va_list params
		Optional extra paramaters. Their number and types depend on the contents
		of the fomat string. They are passed to "vfprintf" function without any
		checking.

	int error_code
		It is an optional parameter. If it is greater than or equal to zero,
		then it is considered to be an error message code number and is
		displayed (changing the message format). If it is equal to -1, then we
		assume there is no error code to display. If it is less than -1, we
		issue a message "Theoretically impossible error" with all properties of
		a fatal error type message.

RETURN VALUE:
	None.

SIDE EFFECTS:
	May terminate the program (if parameter "ec" is equal to EC_FATAL_ERROR or
if an error occurs when the function is executed).

------------------------------------------------------------------------------*/

void ErrorMessageOutput( ErrorCode ec, const char *format, va_list params, // )
	int error_code )
{
	FILE *fp;
	int IO_status = 0;

Beginning:

	fp = ErrorOutputStream;

	//--------------------------------------------------------------------------
	//	See if error_code has correct value. Generate an error if it hasn't.
	//
	if( error_code < -1 )
	{
		ec			= EC_FATAL_ERROR;
		format		= "Theoretically impossible error";
		error_code	= -1;
	}

	//--------------------------------------------------------------------------
	//	Fatal error messages are always displayed.
	//
	if( ec == EC_FATAL_ERROR || fp == NULL )
		fp = stderr;

	//--------------------------------------------------------------------------
	//	First display the error type if fp != NULL.
	//
	switch( ec )
	{
	case EC_WARNING:
		if( error_code == -1 && fp )
			IO_status = fprintf( fp, "WARNING: " );
		else
			IO_status = fprintf( fp, "WARNING (%d): ", error_code );
		break;

	case EC_ERROR:
		if( error_code == -1 && fp )
			IO_status = fprintf( fp, "ERROR: " );
		else
			IO_status = fprintf( fp, "ERROR (%d): ", error_code );
		break;

	case EC_FATAL_ERROR:
		if( error_code == -1 )
			IO_status = fprintf( fp, "\nFATAL ERROR: " );
		else
			IO_status = fprintf( fp, "\nFATAL ERROR (%d): ", error_code );
		break;
	}

	//--------------------------------------------------------------------------
	//	See if output succeeded.
	//
	if( IO_status == EOF ) 
		if( ec != EC_FATAL_ERROR )
		{
			ec			= EC_FATAL_ERROR;
			format		= "Error message output to defined stream failed";
			error_code	= -1;
			goto Beginning;
		}
		else
			goto Abort;

	//--------------------------------------------------------------------------
	//	Then display the custom message (if fp != NULL).
	//
	if( fp ) IO_status = vfprintf( fp, format, params );

	//--------------------------------------------------------------------------
	//	Again see if output succeeded.
	//
	if( IO_status == EOF ) 
		if( ec != EC_FATAL_ERROR )
		{
			ec			= EC_FATAL_ERROR;
			format		= "Error message output to defined stream failed";
			error_code	= -1;
			goto Beginning;
		}
		else
			goto Abort;

	//--------------------------------------------------------------------------
	//	Finally display a newline (if fp != NULL).
	//
	if( fp ) IO_status = vfprintf( fp, "\n", NULL );

	//--------------------------------------------------------------------------
	//	See if output succeeded.
	//
	if( IO_status == EOF )
		if( ec != EC_FATAL_ERROR )
		{
			ec			= EC_FATAL_ERROR;
			format		= "Error message output to defined stream failed";
			error_code	= -1;
			goto Beginning;
		}
		else
			goto Abort;

Abort:
	//--------------------------------------------------------------------------
	//	If it is a fatal error - terminate the program.
	//
	if( ec == EC_FATAL_ERROR ) exit( -1 );
}

/*------------------------------------------------------------------------------

	const char *GetFormatString( int MsgNum )

PURPOSE:
	This function attempts to find in static table of error message format
strings "MessageTable" string number "MsgNum". "MsgNum" is required to denote a
valid "MessageTable" table index. If macro NDEBUG is not defined, "MsgNum"
is checked.

PARAMETERS:
	int MsgNum
		Error message format string number - index into "MessageTable" table.
		Has to be non-negative and less than "MessageTableLen". If
		"MessageTableLen" is zero and "MessageTable" is NULL, any call to this
		function will generate a runtime fatal error.

RETURN VALUE:
	Pointer to format string corresponding to parameter "MsgNum".

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

const char *GetFormatString( int MsgNum )
{
	assert( MessageTableLen != 0 && MessageTable != NULL &&
		MsgNum < MessageTableLen && MsgNum >= 0 );

	return MessageTable[ MsgNum ];
}

//
//	End of static functions' definitions.
//
//==============================================================================
