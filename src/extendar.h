/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose header file with templates.
PROJECT CODE:		SMART PTR
PROJECT FULL NAME:	Implementation of smart pointers that would guarantee full
					control of memory array and pointer operations.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	extendar.h
CREATED:			1995.12.30
LAST MODIFIED:		1996.02.15

DEPENDENCIES:		error.h, memblock.h, smartptr.h, smartdcl.h, sptr_deb.h,
					sptr_ndb.h
					<alloc.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __EXTENDAR_H__
#define __EXTENDAR_H__


#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif



template <class Obj>
class ExtendArray : public Array<Obj>
{
protected:
	size_t len, maxlen;

public:
	//--------------------------------------------------------------------------
	//	Constructors and destructors.
	//
	ExtendArray( void );
	ExtendArray( size_t l, size_t maxl = 0 );
	ExtendArray( size_t l, const Obj &o, size_t maxl = 0 );
	ExtendArray( const SmartPointerBase<Obj> &spb, size_t len );
	ExtendArray( const ExtendArray<Obj>& ear );

	virtual ~ExtendArray( void );

	//--------------------------------------------------------------------------
	//	Some operations that were not possible without explicit knowledge of
	//	array length.
	//
	void FillFragment( const Obj &o, size_t start, size_t end );

	//--------------------------------------------------------------------------
	//	Resize operator.
	//
	void Resize( size_t l, double grow = 1.5 );

	size_t Len( void ) const;
};


template < class Obj >
inline
ExtendArray<Obj>::ExtendArray( void )
	: Array<Obj>(), len( 0 ), maxlen( 0 )
{}


template < class Obj >
inline
ExtendArray<Obj>::ExtendArray( size_t l, size_t maxl )
	: Array<Obj>( maxl > l ? maxl : l ),
	len( l ), maxlen( maxl > l ? maxl : l )
{}


template < class Obj >
inline
ExtendArray<Obj>::ExtendArray( size_t l, const Obj &o, size_t maxl )
	: Array<Obj>( maxl > l ? maxl : l, o ),
	len( l ), maxlen( maxl > l ? maxl : l )
{}


template < class Obj >
inline
ExtendArray<Obj>::ExtendArray( const SmartPointerBase<Obj> &spb, size_t len )
	: Array<Obj>( spb, len ),
	len( len ), maxlen( len )
{}


template < class Obj >
inline
ExtendArray<Obj>::ExtendArray( const ExtendArray<Obj> &ear )
	: Array<Obj>( ear, ear.maxlen ),
	len( ear.len ), maxlen( ear.maxlen )
{}


template < class Obj >
inline
ExtendArray<Obj>::~ExtendArray( void )
{}


template < class Obj >
void ExtendArray<Obj>::FillFragment( const Obj &o, size_t start, size_t end )
{
	assert( end <= len );

	for( size_t i = start; i < end; i++ )
		(*this)[i] = o;
}


//------------------------------------------------------------------------------
//	Object size change procedure.
//

template < class Obj >
void ExtendArray<Obj>::Resize( size_t l, double grow )
{
	assert( grow > 1.0 && grow < 100 );

	if( l <= maxlen )
		len = l;
	else
	{
		maxlen	= Max( size_t( maxlen * grow ), l );
		len		= l;

		Array<Obj>::Resize( maxlen );
	}
}


template < class Obj >
inline
size_t ExtendArray<Obj>::Len( void )
const
{ return len; }


#endif
