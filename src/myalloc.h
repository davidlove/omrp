/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose memory allocation routines.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	myalloc.h
CREATED:			1993.09.03
LAST MODIFIED:		1994.12.16

DEPENDENCIES:		error.h,
					<stdlib.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains template memory allocation / reallocation /
release routines, which do some simple checking and possibly report errors (via
"FatalError" global function).
	This file makes use of templates to force the compiler to generate
appropriate functions for all possible types of allocated tables.

------------------------------------------------------------------------------*/

#ifndef __MYALLOC_H__
#define __MYALLOC_H__


#include <stdlib.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif

/*------------------------------------------------------------------------------

	template <class T, class L> inline T *MALLOC( T *&table, L len )
	template <class T, class L> inline T *CALLOC( T *&table, L len )
	template <class T, class L> inline T *REALLOC( T *&table, L len )

PURPOSE:
	These three functions provide a convenient interface for three basic memory
allocation functions defined in ANSI C standard. Instead of lengthy statements
with plenty of type casting and size computing, we may now issue a brief call
of the form

	xALLOC( table, length )

where "table" is a pointer to the type of table's elements and "length" may be
converted to "size_t" integral type, and expect the compiler to find out what
type of data we want to store and how much space we will need. The type of the
pointer is used for computing the amount of memory needed for a single object
of type "T".
	Each of the above functions will call its ANSI C lowercase counterpart in
order to actually allocate memory.
	In case memory cannot be allocated, "FatalError" message is displayed and
program terminates. "REALLOC" does not require "table" pointer to be non-NULL.
If "table" is NULL, memory is simply allocated. Otherwise standard reallocation
proceeds.
	It is recommended that these function be used only for memory allocation
for simple data types, because constructors / destructors will not be called
for class objects allocated that way. If you want a variable length table of
class objects, you may allocate (and then reallocate) a table of pointers to
this class objects and then individually allocate / delete the objects using
"operator new". In that case you are responsible for keeping score of the
number of objects you need to create or destroy.

PARAMETERS:
	T *&table
		Reference to pointer in which the address of allocated block of memory
		will be stored.

	L len
		Number of "T" type objects we want to store at location "table"

RETURN VALUE:
	Address of allocated block is not only stored in "table" pointer, but is
also returned by the functions.

SIDE EFFECTS:
	If case of allocation failure, each of the xALLOC functions will terminate
the program with a call to "FatalError".

------------------------------------------------------------------------------*/

template <class T, class L>
inline
T *MALLOC( T *&table, L len )
{
	table = (T *)malloc( (size_t)len * sizeof( T ) );
	if( !table )
		FatalError( "Cannot allocate space" );
	return table;
}

template <class T, class L>
inline
T *CALLOC( T *&table, L len )
{
	table = (T *)calloc( (size_t)len, sizeof( T ) );
	if( !table )
		FatalError( "Cannot allocate space" );
	return table;
}

template <class T, class L>
inline
T *REALLOC( T *&table, L len )
{
	if( !table )
		table = (T *)malloc( (size_t)len * sizeof( T ) );
	else
		table = (T *)realloc( (void *)table, (size_t)len * sizeof( T ) );
	if( !table )
		FatalError( "Cannot allocate space" );
	return table;
}

/*------------------------------------------------------------------------------

	template <class T> inline void FREE( T *&table )

PURPOSE:
	This function deallocates memory previously allocated with one of the above
mentioned xALLOC functions (of their ANSI C counterparts). After deallocation,
"table" pointer is set to NULL to prevent future use.

PARAMETERS:
	T *&table
		A pointer to a block of memory to be released. If on entry "table" ==
NULL, then no action is performed. Otherwise memory id freed and "table" is set
to NULL.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

template <class T>
inline
void FREE( T *&table )
{
	if( !table )
		return;
	free( table );
	table = NULL;
}

#endif
