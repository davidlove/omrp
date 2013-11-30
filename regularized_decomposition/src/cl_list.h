/*------------------------------------------------------------------------------
MODULE TYPE:		Specialized list.
PROJECT CODE:		-----------------
PROJECT FULL NAME:	-----------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	cl_list.h
CREATED:			1995.03.09
LAST MODIFIED:		1995.10.05

DEPENDENCIES:		stdtype.h, smartptr.h, error.h, memblock.h, smartdcl.h,
					sptr_deb.h, sptr_ndb.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/


#ifndef __CL_LIST_H__
#define __CL_LIST_H__

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif


class ListByLength
{
private:
	Int_T VecNumber, LengthRange;
	Array<Int_T> Start, LinkFor, LinkBack, End;
	Int_T Next, CurLen;

	Ptr<Int_T> LenArray;

public:
	ListByLength( void );
	ListByLength( Int_T VecNum, Int_T LenRng, Ptr<Int_T> LenArr );
	~ListByLength( void );

	void Initialize( Int_T VecNum, Int_T LenRng, Ptr<Int_T> LenArr );

	Int_T GetNext( void );
	Int_T GetFirst( Int_T len );
	void AddToList( Int_T j );
	void RmFromList( Int_T j );
	void TruncateLength( Int_T i, Int_T l );
};


inline
ListByLength::ListByLength( void )
	: VecNumber( -1 ), LengthRange( -1 ), Next( -1 ), CurLen( -1 )
{}

inline
ListByLength::~ListByLength( void )
{}

#endif
