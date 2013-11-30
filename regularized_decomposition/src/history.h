/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

HEADER FILE NAME:	history.h
CREATED:			1993.10.30
LAST MODIFIED:		1994.08.16

DEPENDENCIES:		stdtype.h

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __HISTORY_H__
#define __HISTORY_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

//==============================================================================
//
//	Class "History" declaration.
//
//==============================================================================

class HistoryEvent;
class History
{
	HistoryEvent *h;
	Int_T len, start, end;
	Bool_T empty;

public:
	History( Int_T HistoryLength );
	~History( void );

	void PushStep( Int_T in, Int_T out );
	Bool_T PopStep( Int_T &in, Int_T &out );
	void ResetHistory( void );
};

//==============================================================================
//
//	End of class "History" declaration.
//
//==============================================================================

#endif
