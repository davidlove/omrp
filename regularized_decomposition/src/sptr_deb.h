/*------------------------------------------------------------------------------
MODULE TYPE:		Internal header file with templates.
PROJECT CODE:		SMART PTR
PROJECT FULL NAME:	Implementation of smart pointers that would guarantee full
					control of memory array and pointer operations.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	smartptr.h
CREATED:			1994.08.20
LAST MODIFIED:		1995.09.12

DEPENDENCIES:		error.h, memblock.h, smartdcl.h, smartptr.h
					<alloc.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A set of class templates implementing smart pointers and arrays is proposed.
Both the smart arrays and the smart pointers have the ability to check array
index before even attempting actual memory access. This control is enabled when
macro "NDEBUG" is not defined.

WARNING: Never include this file directly. Include "smartptr.h" only!

------------------------------------------------------------------------------*/

#ifndef __SPTR_DEB_H__
#define __SPTR_DEB_H__

//==============================================================================
//
//	Class template "SmartPointerBase<class Obj>" functions.
//

template < class Obj >
SmartPointerBase<Obj>::SmartPointerBase( void )
	: start( NULL ), end( NULL ), mem( NULL )
{}


template < class Obj >
SmartPointerBase<Obj>::SmartPointerBase( const SmartPointerBase<Obj> &spb )
	: start( spb.start ), end( spb.end ), mem( spb.mem )
{
	Link();
}


template < class Obj >
SmartPointerBase<Obj>::~SmartPointerBase( void )
{
	UnLink();
	start = end = NULL;
}


template < class Obj >
void SmartPointerBase<Obj>::UnLink( void )
{
	if( mem && mem->UnLink() == 0 ) delete mem;
	mem = NULL;
}


template < class Obj >
void SmartPointerBase<Obj>::Link( void )
{
	if( mem ) mem->Link();
}


template < class Obj >
const SmartPointerBase<Obj> & SmartPointerBase<Obj>::operator=( // )
	const SmartPointerBase<Obj> &spb )
{
	UnLink();
	start = spb.start;
	end = spb.end;
	mem = spb.mem;
	Link();

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
Array<Obj>::Array( void )
	: SmartPointerBase<Obj>()
{}


template < class Obj >
Array<Obj>::Array( size_t l )
	: SmartPointerBase<Obj>()
{
	if( l != 0 )
	{
		SmartPointerBase<Obj>::mem = new MemoryBlock( sizeof( Obj ), l );
		if( !SmartPointerBase<Obj>::mem ) FatalError( "Not enough memory." );
		SmartPointerBase<Obj>::start = (Obj *) (void *) *SmartPointerBase<Obj>::mem;
		SmartPointerBase<Obj>::end = SmartPointerBase<Obj>::start + l;
	}
}


template < class Obj >
Array<Obj>::Array( const SmartPointerBase<Obj> &spb, size_t len )
	: SmartPointerBase<Obj>( spb )
{
	assert( (int)len == spb.SmartPointerBase<Obj>::end - spb.SmartPointerBase<Obj>::start );

	if( SmartPointerBase<Obj>::start )
	{
		SmartPointerBase<Obj>::mem = new MemoryBlock( *SmartPointerBase<Obj>::mem );
		if( !SmartPointerBase<Obj>::mem ) FatalError( "Not enough memory." );

		Obj *NewStart = (Obj *) (void *) *SmartPointerBase<Obj>::mem;

		SmartPointerBase<Obj>::end = NewStart + size_t( SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );
		SmartPointerBase<Obj>::start = NewStart;
	}
}


template < class Obj >
Array<Obj>::Array( size_t l, const Obj &o )
	: SmartPointerBase<Obj>()
{
	if( l != 0 )
	{
		SmartPointerBase<Obj>::mem = new MemoryBlock( sizeof( Obj ), l );
		if( !SmartPointerBase<Obj>::mem ) FatalError( "Not enough memory." );

		SmartPointerBase<Obj>::start = (Obj *) (void *) *SmartPointerBase<Obj>::mem;
		SmartPointerBase<Obj>::end = SmartPointerBase<Obj>::start + l;
		for( Obj *ptr = SmartPointerBase<Obj>::start; l; l--, ptr++ ) *ptr = o;
	}
}


template < class Obj >
Array<Obj>::~Array( void )
{}

//------------------------------------------------------------------------------
//	Access methods.
//

template < class Obj >
Obj & Array<Obj>::operator[]( int i )
{
	assert( SmartPointerBase<Obj>::start != NULL );
	assert( i >= 0 && i < SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	return SmartPointerBase<Obj>::start[i];
}


template < class Obj >
const Obj & Array<Obj>::operator[]( int i )
	const
{
	assert( SmartPointerBase<Obj>::start != NULL );
	assert( i >= 0 && i < SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	return SmartPointerBase<Obj>::start[i];
}


template < class Obj >
#ifdef gnucc
void Array<Obj>::Fill( const Obj o, size_t len, size_t Start )
#else
void Array<Obj>::Fill( const Obj &o, size_t len, size_t Start )
#endif
{
	assert( (int)len == SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );
	assert( SmartPointerBase<Obj>::start != NULL );
	assert( Start == 0 || Start < len );

	for( Obj *ptr = SmartPointerBase<Obj>::start + Start; ptr < SmartPointerBase<Obj>::end; ptr++ )
		*ptr = o;
}


//@BEGIN------------------------------------------------------------------------
// This template object is written to fill only one part of the Array
// (more specifically, this is used in resizing the SolverState Array when a 
// new scenario is added. 

template < class Obj >
void Array<Obj>::FillOne( const Obj &o, size_t theOne )
{
	assert( (int)theOne <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );
	
	Obj *ptr = SmartPointerBase<Obj>::start + theOne -1; 
	*ptr = o;
}

//@END--------------------------------------------------------------------------



//------------------------------------------------------------------------------
//	Pointer arithmetics - operators "+" and "-".
//

template < class Obj >
Ptr<Obj> Array<Obj>::operator+( int i )
{
	assert( SmartPointerBase<Obj>::start != NULL );
	assert( i >= 0 && i < SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	SmartPointerBase<Obj>::start += i;
	Ptr<Obj> p( *this );
	SmartPointerBase<Obj>::start -= i;
	p.ToStart();
	return p;
}


template < class Obj >
Ptr<Obj> Array<Obj>::operator-( int i )
{
	assert( SmartPointerBase<Obj>::start != NULL );
	assert( -i >= 0 && -i < SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	SmartPointerBase<Obj>::start -= i;
	Ptr<Obj> p( *this );
	p.ToStart();
	SmartPointerBase<Obj>::start += i;
	return p;
}


//------------------------------------------------------------------------------
//	Object size change procedure.
//

template < class Obj >
void Array<Obj>::Resize( size_t l )
{
	if( l == 0 )
	{
		if( SmartPointerBase<Obj>::start == NULL )
			return;
		else
		{
			SmartPointerBase<Obj>::UnLink();
			SmartPointerBase<Obj>::start = SmartPointerBase<Obj>::end = NULL;
		}
	}
	else if( !SmartPointerBase<Obj>::start || ( SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start != (int)l ) )
	{
		if( SmartPointerBase<Obj>::start == NULL )
		{
			SmartPointerBase<Obj>::mem = new MemoryBlock( sizeof( Obj ), l );
			if( !SmartPointerBase<Obj>::mem ) FatalError( "Not enough memory." );
		}
		else
			SmartPointerBase<Obj>::mem->Resize( sizeof( Obj ), l );
		SmartPointerBase<Obj>::start = (Obj *) (void *) *SmartPointerBase<Obj>::mem;
		SmartPointerBase<Obj>::end = SmartPointerBase<Obj>::start + l;
	}
}


//------------------------------------------------------------------------------
//	Two versions (member and friend) of array copying.
//

template < class Obj >
void Copy( Array<Obj> &Dst, const Array<Obj> &Src, size_t SrcLen, // )
	size_t DstLen, int nelem, int first )
{
	assert( SrcLen == Src.SmartPointerBase<Obj>::end - Src.SmartPointerBase<Obj>::start );
	assert( DstLen <= Dst.SmartPointerBase<Obj>::end - Dst.SmartPointerBase<Obj>::start );
	assert( first >= 0 && first < SrcLen );

	if( nelem < 0 )
		nelem = SrcLen - first;

	assert( nelem >= 0 && nelem <= SrcLen - first );

	if( DstLen < nelem ) Dst.Resize( nelem );

	memcpy( Dst.SmartPointerBase<Obj>::start, Src.SmartPointerBase<Obj>::start + first, nelem * sizeof( Obj ) );
}

template < class Obj >
void Array<Obj>::Copy( const Array<Obj> &Src, size_t SrcLen, size_t DstLen, // )
	int nelem, int first )
{
	assert( (int)SrcLen == Src.SmartPointerBase<Obj>::end - Src.SmartPointerBase<Obj>::start );
	assert( (int)DstLen <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );
	assert( first >= 0 && first < (int)SrcLen );

	if( nelem < 0 )
		nelem = SrcLen - first;

	assert( nelem >= 0 && nelem <= (int)SrcLen - first );

	if( SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start < nelem ) Resize( nelem );

	memcpy( SmartPointerBase<Obj>::start, Src.SmartPointerBase<Obj>::start + first, nelem * sizeof( Obj ) );
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
Ptr<Obj>::Ptr( void )
	: SmartPointerBase<Obj>(), ptr( NULL )
{}


template < class Obj >
Ptr<Obj>::Ptr( const Ptr<Obj> &p )
	: SmartPointerBase<Obj>( p ), ptr( p.ptr )
{}


template < class Obj >
Ptr<Obj>::Ptr( const Array<Obj> &a )
	: SmartPointerBase<Obj>( a )
{
	ToStart();
}


template < class Obj >
Ptr<Obj>::Ptr( const Array<Obj> &a, const size_t Start, const size_t End )
	: SmartPointerBase<Obj>( a )
{
	assert( Start <= End );
	assert( (int)Start <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start && (int)End <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	SmartPointerBase<Obj>::end = SmartPointerBase<Obj>::start + End;
	ptr = SmartPointerBase<Obj>::start += Start;
}


template < class Obj >
void Ptr<Obj>::ExtractFragment( const SmartPointerBase<Obj> &spb, // )
	const size_t Start, const size_t End )
{
	//--------------------------------------------------------------------------
	//	If "spb" is the same object, to which "this" points, then we may skip
	//	the assignment of the base object (and "Link"/"UnLink" call pair).
	//
	if( (SmartPointerBase<Obj> *)this != &spb )
	{
		SmartPointerBase<Obj>::UnLink();
		SmartPointerBase<Obj>::start	= spb.SmartPointerBase<Obj>::start;
		SmartPointerBase<Obj>::end	= spb.SmartPointerBase<Obj>::end;
		SmartPointerBase<Obj>::mem	= spb.SmartPointerBase<Obj>::mem;
		SmartPointerBase<Obj>::Link();
	}

	assert( Start <= End );
	assert( (int)Start <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start && (int)End <= SmartPointerBase<Obj>::end - SmartPointerBase<Obj>::start );

	SmartPointerBase<Obj>::end = SmartPointerBase<Obj>::start + End;
	ptr = SmartPointerBase<Obj>::start += Start;
}


template < class Obj >
Ptr<Obj>::~Ptr( void )
{
	ptr = NULL;
}


//------------------------------------------------------------------------------
//	Assignment operators.
//

template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator=( const Array<Obj> &a )
{
	SmartPointerBase<Obj>::operator =( a );
    ToStart();

	return *this;
}


template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator=( const Ptr<Obj> &p )
{
	if( this == &p ) return *this;

	SmartPointerBase<Obj>::operator =( p );
	ptr = p.ptr;

	return *this;
}


//------------------------------------------------------------------------------
//	Access methods.
//

template < class Obj >
Obj& Ptr<Obj>::operator[]( int i )
{
	assert( ptr != NULL );
	assert( i >= 0 && i < SmartPointerBase<Obj>::end - ptr );

	return ptr[i];
}


template < class Obj >
const Obj& Ptr<Obj>::operator[]( int i )
	const
{
	assert( ptr != NULL );
	assert( i >= 0 && i < SmartPointerBase<Obj>::end - ptr );

	return ptr[i];
}


template < class Obj >
Obj& Ptr<Obj>::operator*( void )
{
	assert( ptr != NULL );

	return *ptr;
}


template < class Obj >
const Obj& Ptr<Obj>::operator*( void )
	const
{
	assert( ptr != NULL );

	return *ptr;
}


//------------------------------------------------------------------------------
//	Increment and decrement operators (prefix).
//

template < class Obj >
Ptr<Obj>& Ptr<Obj>::operator++( void )
{
	assert( ptr != NULL );
	assert( ptr < SmartPointerBase<Obj>::end );

	ptr++;
	return *this;
}


template < class Obj >
Ptr<Obj>& Ptr<Obj>::operator--( void )
{
	assert( ptr != NULL );
	assert( ptr > SmartPointerBase<Obj>::start );

	ptr--;
	return *this;
}


//------------------------------------------------------------------------------
//	Pointer arithmetics.
//

template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator+=( int i )
{
	assert( ptr != NULL );
	assert( ptr + i >= SmartPointerBase<Obj>::start && ptr + i <= SmartPointerBase<Obj>::end );

	ptr += i;
	return *this;
}


template < class Obj >
const Ptr<Obj>& Ptr<Obj>::operator-=( int i )
{
	assert( ptr != NULL );
	assert( ptr - i >= SmartPointerBase<Obj>::start && ptr - i <= SmartPointerBase<Obj>::end );

	ptr -= i;
	return *this;
}


template < class Obj >
Ptr<Obj> Ptr<Obj>::operator+( int i )
	const
{
	assert( ptr != NULL );
	assert( ptr + i >= SmartPointerBase<Obj>::start && ptr + i <= SmartPointerBase<Obj>::end );

	Ptr<Obj> p( *this );
	p.ptr += i;

	return p;
}


template < class Obj >
Ptr<Obj> Ptr<Obj>::operator-( int i )
	const
{
	assert( ptr != NULL );
	assert( ptr - i >= SmartPointerBase<Obj>::start && ptr - i <= SmartPointerBase<Obj>::end );

	Ptr<Obj> p( *this );
	p.ptr -= i;

	return p;
}


//------------------------------------------------------------------------------
//	Other pointer manipulation methods.
//
template < class Obj >
void Ptr<Obj>::ToStart( void )
{
	ptr = SmartPointerBase<Obj>::start;
}


template < class Obj >
void Ptr<Obj>::ToEnd( void )
{
	ptr = SmartPointerBase<Obj>::end;
}

#endif
