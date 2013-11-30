/*------------------------------------------------------------------------------
MODULE TYPE:		Specialized list.
PROJECT CODE:		-----------------
PROJECT FULL NAME:	-----------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	cl_list.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1995.10.05

DEPENDENCIES:		cl_list.h, stdtype.h, smartptr.h, error.h, memblock.h,
					smartdcl.h, sptr_deb.h, sptr_ndb.h
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

------------------------------------------------------------------------------*/


#include <assert.h>

#ifndef __CL_LIST_H__
#	include "cl_list.h"
#endif


/*------------------------------------------------------------------------------

	ListByLength::ListByLength( Int_T VecNum, Int_T LenRng, Ptr<Int_T> LenArr )

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

ListByLength::ListByLength( Int_T VecNum, Int_T LenRng, Ptr<Int_T> LenArr )
	: VecNumber( VecNum ), LengthRange( LenRng ),
	Start( LenRng + 1, -1 ), LinkFor( VecNum, -1 ),
	LinkBack( VecNum, -1 ), End( LenRng + 1, -1 ),
	Next( -1 ), CurLen( -1 ), LenArray( LenArr )
{
	for( Int_T i = 0; i < VecNumber; i++ )
		AddToList( i );
}


/*------------------------------------------------------------------------------

	void ListByLength::Initialize( Int_T VecNum, Int_T LenRng,
		Ptr<Int_T> LenArr )

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

void ListByLength::Initialize( Int_T VecNum, Int_T LenRng, Ptr<Int_T> LenArr )
{
	assert( VecNumber == -1 );
	assert( LengthRange == -1 );

	VecNumber	= VecNum;
	LengthRange	= LenRng;

	Start.Resize( LenRng + 1 );			Start.Fill(		-1,	LenRng + 1 );
	LinkFor.Resize( VecNum );			LinkFor.Fill(	-1,	VecNum );
	LinkBack.Resize( VecNum );			LinkBack.Fill(	-1,	VecNum );
	End.Resize( LenRng + 1 );			End.Fill(		-1,	LenRng + 1 );

	Next		= -1;
	CurLen		= -1;
	LenArray	= LenArr;

	for( Int_T j = 0; j < VecNum; j++ )
		AddToList( j );
}


/*------------------------------------------------------------------------------

	Int_T ListByLength::GetFirst( Int_T j )

PURPOSE:
	Returns the index of the first column of length "j".

PARAMETERS:
	Int_T j
		Column length (in range <0,m>).

RETURN VALUE:
	Column number (or -1 if none was found).

SIDE EFFECTS:
	The requested column length is stored as "CurLen". Also the next column
of the same length is remebered as "Next".

------------------------------------------------------------------------------*/

Int_T ListByLength::GetFirst( Int_T j )
{
	assert( j >= 0 && j <= LengthRange );

	j = Start[ CurLen = j ];
	Next = Int_T( ( j >= 0 ) ? LinkFor[j] : -1 );
	return j;
}


/*------------------------------------------------------------------------------

	Int_T ListByLength::GetNext( void )

PURPOSE:
	Retrieves the number of the next column of the same length, as was requested
before. Reads the number from the "Next" pointer.

PARAMETERS:
	None.

RETURN VALUE:
	Column number.

SIDE EFFECTS:
	Moves the "Next" pointer to the next column of the same length.

------------------------------------------------------------------------------*/

Int_T ListByLength::GetNext( void )
{
	Int_T j = Next;
	Next = Int_T( ( Next >= 0 ) ? LinkFor[ Next ] : -1 );
	return j;
}


/*------------------------------------------------------------------------------

	void ListByLength::AddToList( Int_T j )

PURPOSE:
	Adds a column number "j" to a list corresponding to the column's current
length. If "Next" pointer for columns of that length (which is read from the
"len" table) indicates, that there are no more such columns, then "Next"
is updated to point at the column added.

PARAMETERS:
	Int_T j
		Number of column to be added to some list.

RETURN VALUE:
	None.

SIDE EFFECTS:
	If necessary, updates "Next" (see explanation above).

------------------------------------------------------------------------------*/

void ListByLength::AddToList( Int_T j )
{
	assert( j >= 0 && j < VecNumber );

	Int_T len	= LenArray[j],
		&head	= Start[len],
		&end	= End[len];

	if( head < 0 && end < 0 )
	{
		LinkBack[j]		= -1;
		head				= j;
	}
	else
	{
		LinkFor[ end ]	= j;
		LinkBack[j]		= end;
	}

	LinkFor[j]			= -1;
	end					= j;

	if( CurLen == len && Next <= 0 ) Next = j;
}


/*------------------------------------------------------------------------------

	void ListByLength::RmFromList( Int_T j )

PURPOSE:
	Removes column "j" from the list of columns. If the column being removed is
the one stored in "Next" then "Next" is updated - set to point to the
next entry on the appropriate list.

PARAMETERS:
	Int_T j
		Number of column to be removed from some list.

RETURN VALUE:
	None.

SIDE EFFECTS:
	If necessary, "Next" is updated (see explanation above).

------------------------------------------------------------------------------*/

void ListByLength::RmFromList( Int_T j )
{
	assert( j >= 0 && j < VecNumber );

	Int_T next	= LinkFor[j],
		prev	= LinkBack[j],
		len		= LenArray[j];

	if( next >= 0 )
		LinkBack[ next ]	= prev;
	else
		End[ len ] 		= prev;

	if( prev >= 0 )
		LinkFor[ prev ]		= next;
	else
		Start[ len ]		= next;

	if( CurLen == len && Next == j ) Next = LinkFor[j];
}


/*------------------------------------------------------------------------------

	void ListByLength::TruncateLength( Int_T j, Int_T l )

PURPOSE:
	This function is an interface to "RmFromList" and "AddToList" pair of
functions. It performs action corresponding to shortening column "j" by "l"
entries.

PARAMETERS:
	Int_T j
		Number of column, which length is to be updated.

	Int_T l
		Number of non-zero entries to be substracted form the column's current
		length.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None (except the side effects of the called functions, namely possible
update to "Next").

------------------------------------------------------------------------------*/

void ListByLength::TruncateLength( Int_T j, Int_T l )
{
	assert( j >= 0 && j < VecNumber );
	assert( LenArray[j] >= l );

	RmFromList( j );
	LenArray[j] -= l;
	AddToList( j );
}
