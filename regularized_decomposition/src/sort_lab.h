/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose implementation of a sorted array of labels.
PROJECT CODE:		-------------------
PROJECT FULL NAME:	-------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	sort_lab.h
CREATED:			1994.04.28
LAST MODIFIED:		1996.09.20

DEPENDENCIES:		myalloc.h, std_tmpl.h, smartptr.h, stdtype.h,
					<stdlib.h>, <string.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SORT_LAB_H__
#define __SORT_LAB_H__

#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif


//==============================================================================
//
//	Declaration of "SortedArrayOfLabels" class.
//
//==============================================================================

enum LabelOrder { LO_SORTED, LO_UNSORTED };
typedef char charLEN[LAB_LEN+1];
struct HTable
{
	char *LPtr;
	int LNr;
};


class SortedArrayOfLabels
{
private:
	enum LabelStatus { EMPTY, FILLED, SORTED };

	static int CmpLabels( const void *lab1, const void *lab2 );

private:
	charLEN *Label;
	HTable *h;
	size_t len, maxlen;
	LabelStatus Status;

#ifndef NDEBUG
	static HTable *hstart, *hend, *hentry;
#endif

public:
	SortedArrayOfLabels( size_t InitLen );
	~SortedArrayOfLabels( void );

	void FreeMemory( void );

public:
	int AddLabel( const char *lab );

	const char *FindLabel( size_t i, LabelOrder lo = LO_UNSORTED ) const;
	int FindLabel( const char *lab, int last_n = -1 ) const;

	int Duplicates( void ) const;

	int SortLabels( void );

	void DiscardLabels( void );

	int NumberOfLabels( void ) const;

	int RemoveLabels( const Array<Bool_T> &Exclude );

	int IsSorted( void ) const;
	int IsFilled( void ) const;
};


//==============================================================================
//
//	End of declaration of "SortedArrayOfLabels" class.
//
//==============================================================================


//==============================================================================
//
//	Definitions of "SortedArrayOfLabels" class member functions.
//
//==============================================================================


//------------------------------------------------------------------------------
//
//	Object creation / destruction and data removal.
//

inline
SortedArrayOfLabels::~SortedArrayOfLabels( void )
{
	FREE( Label );
	FREE( h );
}


inline
void SortedArrayOfLabels::FreeMemory( void )
{
	FREE( Label );
	FREE( h );
	
	maxlen = len = 0;
	Status = EMPTY;
}


//------------------------------------------------------------------------------
//
//	Data manipulation - storage, access, sorting.
//


inline
int SortedArrayOfLabels::NumberOfLabels( void )
const
{ return len; }


inline
int SortedArrayOfLabels::IsSorted( void )
const
{ return Status == SORTED; }


inline
int SortedArrayOfLabels::IsFilled( void )
const
{ return Status == FILLED; }


#endif
