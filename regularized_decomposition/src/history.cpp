/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	history.cpp
CREATED:			1993.10.24
LAST MODIFIED:		1994.08.16

DEPENDENCIES:		history.cpp, stdtype.h, history.h, solvcode.h, error.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __HISTORY_H__
#	include "history.h"
#endif
#ifndef __SOLVCODE_H__
#	include "solvcode.h"
#endif


/*------------------------------------------------------------------------------

	Auxiliary class: Class HistoryEvent
		A single object of this class stores information sufficient to go back
		one step in history.

------------------------------------------------------------------------------*/

class HistoryEvent
{
	Int_T in, out;

	HistoryEvent( void ) : in( SLV_HIST_LIST_EMPTY ),
		out( SLV_HIST_LIST_EMPTY ) {};
	friend class History;
};


/*------------------------------------------------------------------------------

	History::~History( void )

PURPOSE:
	Class 'History' destructor. Deallocates storage.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

History::~History( void )
	{ if( h ) delete [] h; }


/*------------------------------------------------------------------------------

	void History::ResetHistory( void )

PURPOSE:
	Clear the history list.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void History::ResetHistory( void )
	{ start = end = 0; empty = True; }


/*------------------------------------------------------------------------------

	History::History( Int_T hLen )

PURPOSE:
	Class 'History' constructor. Allocates storage necessary during operation.
Initializes the class to represent an empty list.

PARAMETERS:
	Int_T hLen
		Maximum number of history list entries to be stored. If the limit is
		exceeded, the oldest entries are removed and replaced with new ones.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

History::History( Int_T hLen )
	: start( 0 ), end( 0 ), empty( True )
{
	assert( hLen > 0 );

	h = new HistoryEvent[ len = hLen ];

	if( !h ) FatalError( "history.cpp: History: Out of memory" );
}


/*------------------------------------------------------------------------------

	void History::PushStep( Int_T in, Int_T out )

PURPOSE:
	Stores pair (in, out) on the list. If necessary removes the oldest entry
from the list. The list is accessed as a stack would be (FIFO type access), but
is actually a fixed maximum length cyclical list. Therefore when an attempt is
made to store more then 'len' entries on the list, an oldest entry is replaced
with the new one and the 'start' and 'end' markers are shifted one position.

PARAMETERS:
	Int_T in, Int_T out
		Data to be stored.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void History::PushStep( Int_T in, Int_T out )
{
	if( start == end && !empty )
		( start < len - 1 ) ? start++ : (start = 0);

	empty			= False;
	h[ end ].in		= in;
	h[ end ].out	= out;

	( end < len - 1 ) ? end++ : (end = 0);
}


/*------------------------------------------------------------------------------

	Bool_T History::PopStep( Int_T &in, Int_T &out )

PURPOSE:
	This function retrieves the latest information stored on the history list
and removes it from the list.

PARAMETERS:
	Int_T &in, Int_T &out
		These two parameters are used only to return the information to the
		caller. Their values on entry is not important. On failure they are
		both assigned -1.

RETURN VALUE:
	Boolean succes status (True - if an entry was successfully retrieved, False
if the list is empty).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T History::PopStep( Int_T &in, Int_T &out )
{
	if( empty )
	{
		in = out = SLV_HIST_LIST_EMPTY;
		return False;
	}

	if( start == end )
		empty = True;
	else 
		( end > 0 ) ? end-- : (end = Int_T( len - 1 ) );

	in	= h[ end ].in;
	out	= h[ end ].out;

	return True;
}
