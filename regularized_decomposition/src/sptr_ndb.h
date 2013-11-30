/*------------------------------------------------------------------------------
MODULE TYPE:		Internal header file with templates.
PROJECT CODE:		SMART PTR
PROJECT FULL NAME:	Implementation of smart pointers that would guarantee full
					control of memory array and pointer operations.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	sptr_ndb.h
CREATED:			1994.08.20
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		error.h, memblock.h, smartdcl.h, smartptr.h
					<alloc.h>, <assert.h>, <string.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A set of class templates implementing smart pointers and arrays is proposed.
Both the smart arrays and the smart pointers have the ability to check array
index before even attempting actual memory access. This control is enabled when
macro "NDEBUG" is not defined.

------------------------------------------------------------------------------*/

#ifndef __SPTR_NDB_H__
#define __SPTR_NDB_H__

#if defined( xlc )
#	include <memory.h>
#endif

#ifdef watcom
#	include <string.h>
#endif

#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif

//==============================================================================
//
//	Class template "SmartPointerBase<class Obj>" functions.
//

template < class Obj >
inline
SmartPointerBase<Obj>::SmartPointerBase( void )
	: ptr( NULL )
{}


template < class Obj >
inline
SmartPointerBase<Obj>::SmartPointerBase( const SmartPointerBase<Obj> &spb )
	: ptr( spb.ptr )
{}


template < class Obj >
inline
SmartPointerBase<Obj>::~SmartPointerBase( void )
{
	ptr = NULL;
}


template < class Obj >
inline
const SmartPointerBase<Obj> & SmartPointerBase<Obj>::operator=( // )
	const SmartPointerBase<Obj> &spb )
{
	ptr = spb.ptr;

	return *this;
}


//==============================================================================
//
//	INLINE AND NON-INLINE IMPLEMENTATIONS OF THE METHODS OF THE CLASS
//	TEMPLATE "Array".
//
//==============================================================================


//------------------------------------------------------------------------------
//	Object creation and destruction methods.
//

template < class Obj >
inline
Array<Obj>::Array( void )
	: SmartPointerBase<Obj>()
{}


template < class Obj >
Array<Obj>::Array( size_t l )
	: SmartPointerBase<Obj>()
{
	if( l != 0 )
		MALLOC( ptr, l );
}


template < class Obj >
Array<Obj>::Array( const SmartPointerBase<Obj> &spb, size_t len )
	: SmartPointerBase<Obj>()
{
	if( len > 0 )
	{
		MALLOC( ptr, len );
		for( Obj *src = spb.ptr, *dst = ptr; len; len--, src++, dst++ )
			*dst = *src;
	}
}


template < class Obj >
inline
Array<Obj>::Array( size_t l, const Obj &o )
	: SmartPointerBase<Obj>()
{
	if( l > 0 )
	{
		MALLOC( ptr, l );
		for( Obj *p = ptr; l; l--, p++ ) *p = o;
	}
}


template < class Obj >
inline
Array<Obj>::~Array( void )
	{ if( ptr ) FREE( ptr ); }

//------------------------------------------------------------------------------
//	Access methods.
//

template < class Obj >
inline
Obj & Array<Obj>::operator[]( int i )
	{ return ptr[i]; }


template < class Obj >
inline
const Obj & Array<Obj>::operator[]( int i )
	const
	{ return ptr[i]; }


template < class Obj >
#ifdef gnucc
void Array<Obj>::Fill( const Obj o, size_t len, size_t Start )
#else
void Array<Obj>::Fill( const Obj &o, size_t len, size_t Start )
#endif
{
	for( Obj *p = ptr + Start, *end = ptr + len; p < end; p++ )
		*p = o;
}


//------------------------------------------------------------------------------
//	Pointer arithmetics - operators "+" and "-".
//

template < class Obj >
inline
Ptr<Obj> Array<Obj>::operator+( int i )
{
	ptr += i;
	Ptr<Obj> p( *this );
	ptr -= i;
	return p;
}


template < class Obj >
inline
Ptr<Obj> Array<Obj>::operator-( int i )
{
	ptr -= i;
	Ptr<Obj> p( *this );
	ptr += i;
	return p;
}


//------------------------------------------------------------------------------
//	Object size change procedure.
//

template < class Obj >
inline
void Array<Obj>::Resize( size_t l )
{
	if( l == 0 )
		FREE( ptr );
	else
		REALLOC( ptr, l );
}


//------------------------------------------------------------------------------
//	Two versions (member and friend) of array copying.
//

template < class Obj >
void Copy( Array<Obj> &Dst, const Array<Obj> &Src, size_t SrcLen, size_t DstLen,
	int nelem, int first )
{
	if( nelem < 0 )
		nelem = SrcLen - first;

	if( nelem > DstLen ) Dst.Resize( nelem );

	memcpy( Dst.ptr, Src.ptr + first, nelem * sizeof( Obj ) );
}

template < class Obj >
void Array<Obj>::Copy( const Array<Obj> &Src, size_t SrcLen, size_t DstLen,
	int nelem, int first )
{
	if( nelem < 0 )
		nelem = SrcLen - first;

	if( (int)DstLen < nelem ) Resize( nelem );

	memcpy( ptr, Src.ptr + first, nelem * sizeof( Obj ) );
}


//==============================================================================
//
//	INLINE AND NON-INLINE IMPLEMENTATIONS OF THE METHODS OF THE CLASS
//	TEMPLATE "Ptr".
//
//==============================================================================

//------------------------------------------------------------------------------
//	Object creation methods.
//

template < class Obj >
inline
Ptr<Obj>::Ptr( void )
	: SmartPointerBase<Obj>()
{}


template < class Obj >
inline
Ptr<Obj>::Ptr( const Ptr<Obj> &p )
	: SmartPointerBase<Obj>( p )
{}


template < class Obj >
inline
Ptr<Obj>::Ptr( const Array<Obj> &a )
	: SmartPointerBase<Obj>( a )
{}


template < class Obj >
inline
Ptr<Obj>::Ptr( const Array<Obj> &a, const size_t Start, const size_t )
	: SmartPointerBase<Obj>( a )
	{ ptr += Start; }


template < class Obj >
inline
void Ptr<Obj>::ExtractFragment( const SmartPointerBase<Obj> &spb, // )
	const size_t Start, const size_t )
	{ ptr = spb.ptr + Start; }


template < class Obj >
inline
Ptr<Obj>::~Ptr( void )
{}


//------------------------------------------------------------------------------
//	Assignment operators.
//

template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator=( const Array<Obj> &a )
{
	ptr = a.ptr;
	return *this;
}


template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator=( const Ptr<Obj> &p )
{
	ptr = p.ptr;
	return *this;
}


//------------------------------------------------------------------------------
//	Access methods.
//

template < class Obj >
inline
Obj& Ptr<Obj>::operator[]( int i )
	{ return ptr[i]; }


template < class Obj >
inline
const Obj& Ptr<Obj>::operator[]( int i )
	const
	{ return ptr[i]; }


template < class Obj >
inline
Obj& Ptr<Obj>::operator*( void )
	{ return *ptr; }


template < class Obj >
inline
const Obj& Ptr<Obj>::operator*( void )
	const
	{ return *ptr; }


//------------------------------------------------------------------------------
//	Increment and decrement operators (prefix).
//

template < class Obj >
inline
Ptr<Obj>& Ptr<Obj>::operator++( void )
{
	ptr++;
	return *this;
}


template < class Obj >
inline
Ptr<Obj>& Ptr<Obj>::operator--( void )
{
	ptr--;
	return *this;
}


//------------------------------------------------------------------------------
//	Pointer arithmetics.
//

template < class Obj >
inline
const Ptr<Obj>& Ptr<Obj>::operator+=( int i )
{
	ptr += i;
	return *this;
}


template < class Obj >
inline
const Ptr<Obj>& Ptr<Obj>::operator-=( int i )
{
	ptr -= i;
	return *this;
}


template < class Obj >
Ptr<Obj> Ptr<Obj>::operator+( int i )
	const
{
	Ptr<Obj> p( *this );
	p.ptr += i;

	return p;
}


template < class Obj >
Ptr<Obj> Ptr<Obj>::operator-( int i )
	const
{
	Ptr<Obj> p( *this );
	p.ptr -= i;

	return p;
}


#endif
