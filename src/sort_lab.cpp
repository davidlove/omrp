/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose implementation of a sorted array of labels.
PROJECT CODE:		-------------------
PROJECT FULL NAME:	-------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	sort_lab.cpp
CREATED:			1994.04.28
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		myalloc.h, std_tmpl.h, smartptr.h, stdtype.h, sort_lab.h
					<stdlib.h>, <string.h>, <assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif

#include <string.h>
#include <assert.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif


//------------------------------------------------------------------------------
//	Implicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )
	template charLEN * MALLOC( charLEN *& Table, size_t Len );
	template charLEN * REALLOC( charLEN *& Table, size_t Len );
	template void FREE( charLEN *& Table );

	template HTable * REALLOC( HTable *& Table, size_t Len );
	template void FREE( HTable *& Table );
#endif
//
//------------------------------------------------------------------------------


#ifndef NDEBUG

HTable *SortedArrayOfLabels::hstart		= NULL,
	*SortedArrayOfLabels::hend			= NULL,
	*SortedArrayOfLabels::hentry		= NULL;

#endif


//------------------------------------------------------------------------------
//
//	Object creation / destruction and data removal.
//

SortedArrayOfLabels::SortedArrayOfLabels( size_t InitLen )
	: Label( NULL ), h( NULL ), len( 0 ), maxlen( InitLen ), Status( EMPTY )
{
	MALLOC( Label, maxlen );

	for( int i = 0; i < (int)maxlen; i++ )
		Label[i][0] = '\0';
}


void SortedArrayOfLabels::DiscardLabels( void )
{
	if( Status == EMPTY ) return;

	FREE( h );

	for( int i = 0; i < (int)len; i++ )
		Label[i][0] = '\0';

	len = 0;
	Status = EMPTY;
}


//------------------------------------------------------------------------------
//
//	Static private function for comparing two labels (called by 'qsort' and
//	'bsearch' standard C library functions).
//

int SortedArrayOfLabels::CmpLabels( const void *lab1, const void *lab2 )
{
	assert( ( hentry && (HTable *)lab1 == hentry ) ||
		( (HTable *)lab1 >= hstart && (HTable *)lab1 < hend ) );
	assert( ( hentry && (HTable *)lab1 == hentry ) ||
		( (HTable *)lab2 >= hstart && (HTable *)lab2 < hend ) );

	return strncmp( ((HTable *)lab1)->LPtr, ((HTable *)lab2)->LPtr, LAB_LEN );
}


//------------------------------------------------------------------------------
//
//	Data manipulation - storage, access, sorting.
//

int SortedArrayOfLabels::AddLabel( const char *lab )
{
	if( Status == SORTED )
		FREE( h );
	Status = FILLED;

	if( len >= maxlen )
	{
		size_t maxlenNew = maxlen + Max( unsigned(maxlen/2), 50U );
		REALLOC( Label,	maxlenNew );

		while( maxlen < maxlenNew )
			Label[maxlen++][0] = '\0';
	}

	strncpy( Label[len], lab, LAB_LEN );
	Label[len][LAB_LEN]	= '\0';
	len++;

	return 0;
}


int SortedArrayOfLabels::SortLabels( void )
{
	if( Status == FILLED )
	{
//		maxlen = len;

//		REALLOC( Label,	maxlen );
		REALLOC( h,		len );

		for( int i = 0; i < (int)len; i++ )
		{
			h[i].LPtr	= Label[i];
			h[i].LNr	= i;
		}

#ifndef NDEBUG
		hstart	= h;
		hend	= h+len;
		hentry	= NULL;
#endif
		qsort( h, len, sizeof( HTable ), SortedArrayOfLabels::CmpLabels );
		Status = SORTED;

#ifndef NDEBUG
		hstart	= NULL;
		hend	= NULL;
#endif
		return 0;
	}
	else if( Status == SORTED )
		return 0;

	return -1;
}


const char *SortedArrayOfLabels::FindLabel( size_t i, LabelOrder lo )
	const
{
	assert( i < len );

	switch( lo )
	{
	case LO_SORTED:
		assert( h != NULL && Status == SORTED );
		return h[i].LPtr;

	case LO_UNSORTED:
		return Label[i];

#ifndef NDEBUG
	default:
		abort();
#endif
	}
	return NULL;
}


int SortedArrayOfLabels::FindLabel( const char *lab, int last_n )
	const
{
	switch( Status )
	{
	case EMPTY:
		break;

	case FILLED:
		{
			for( int i = (last_n > 0) ? len - last_n : 0; i < (int)len; i++ )
				if( strncmp( Label[i], lab, LAB_LEN ) == 0 )
					return i;
		}
		break;

	case SORTED:
		{
			assert( h != NULL );

			HTable hEntry = { (char *) lab, 0 };

#ifndef NDEBUG
			hstart	= h;
			hend	= h+len;
			hentry	= &hEntry;
#endif
			HTable *item = (HTable *) bsearch( &hEntry, h, len, sizeof( HTable ),
				SortedArrayOfLabels::CmpLabels );

#ifndef NDEBUG
			hstart	= NULL;
			hend	= NULL;
			hentry	= NULL;
#endif
			if( item ) return item->LNr;
		}
	}

	return -1;
}



int SortedArrayOfLabels::Duplicates( void )
	const
{
	assert( Status == SORTED );

	int count = 0;

	for( int i = 1; i < (int)len; i++ )
		if( strncmp( h[i-1].LPtr, h[i].LPtr, LAB_LEN ) == 0 )
			count++;

	return count;
}


int SortedArrayOfLabels::RemoveLabels( const Array<Bool_T> &Exclude )
{
	if( Status == EMPTY )
		return -1;

	FREE( h );

	int i, j;
	for( j = i = 0; i < (int)len && !Exclude[i]; i++, j++ )
		;

	for( ; i < (int)len; i++ )
		if( !Exclude[i] )
			strncpy( Label[j++], Label[i], LAB_LEN );

//	maxlen =
	len = j;

//	REALLOC( Label,	maxlen );
	Status = FILLED;

	return 0;
}
