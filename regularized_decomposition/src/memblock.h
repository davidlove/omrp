/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose block defining a memory management class.
PROJECT CODE:		MEMORY BLOCKS
PROJECT FULL NAME:	Implementation of memory block manager.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	memblock.h
CREATED:			1994.04.11
LAST MODIFIED:		1996.09.20

DEPENDENCIES:		error.h,
					<string.h>, <stdlib.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains a declaration and some definitions of class
"MemoryBlock". An object of this class holds information about one block of
allocated memory. The object knows the address of the allocated memory (member
"void *mem"), its size (member "size_t len") and the number of existing
(registered) references to the block (member "long link").
	Additionally when self-debugging is turned on (macro "NDEBUG" is not
defined) static class member "long count" stores the total number of allocated
memory blocks.

	The objects of the class should be used as follows:
1.	When a block of memory of a specific size is needed, a new object of class
	"MemoryBlock" should be created:

	MemoryBlock *mb = new MemoryBlock( sizeof( T ), N );

	where "T" denotes the name of the SIMPLE type of data to be stored in the
	memory block and "N" specifies the maximum number of "T" objects to be
	stored in the block.

2.	When we want to access the allocated memory, we call member function
	"operator void *":

	MemoryBlock *mb = new MemoryBlock( sizeof( int ), 100 );
	int *i = (int *)((void *)*mb);

3.	When we want to create an additional reference to the same memory block, we
	call member function "Link()":

	int *j = (int *)((void *)*mb);
	mb->Link();

4.	When one of the referenes is not needed any more, you should decrement the
	link count. Only when the link count reaches zero should the "MemoryBlock"
	be destroyed.

	i = NULL;
	if( mb->Unlink() == 0 ) { delete mb; mb = NULL; }

	Class "MemoryBlock" is a relatively low level tool. It was created as a
building block of smart pointer class templates, in which all those usage
details were hidden from the "end user".

WARNING: Class "MemoryBlock" may only hold objects not requiring neither
constructors nor destructors. E.g. classes which have virtual functions or
virtual base classes (and thus use virtual tables) may not be managed by a
"MemoryBlock" object. Class "MemoryBlock" was meant to be used as a manager of
simple types' arrays and pointers.

------------------------------------------------------------------------------*/

#ifndef __MEMBLOCK_H__
#define __MEMBLOCK_H__

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif


//==============================================================================
//
//	Class "MemoryBlock" declaration.
//

class MemoryBlock
{
private:
	void *mem;
	size_t len;
	long link;

#ifndef NDEBUG
	static long count;
#endif

public:
	MemoryBlock( const MemoryBlock &mb );
	MemoryBlock( size_t ObjLen, size_t ObjNum = 1 );

	~MemoryBlock( void );

	//--------------------------------------------------------------------------
	//	This function is declared but is never defined. No reasonable semantics
	//	for "MemoryBlock::operator =" seems to exist, so this declaration is
	//	here only in order to prevent the compiler from generating the default 
	//	"operator =" for this class.
	//
private:
	const MemoryBlock & operator =( const MemoryBlock &mb );

public:
	operator void *( void );

	void Resize( size_t ObjNum, size_t ObjSize = 1 );

	long Link( void );
	long UnLink( void );

	size_t Len( void );

#ifndef NDEBUG
	static long Count( void );
#endif
};

//
//	End of class "MemoryBlock" declaration.
//
//==============================================================================


//==============================================================================
//
//	Definitions of class "MemoryBlock" inline member functions.
//

//------------------------------------------------------------------------------
//
//	Object creation methods.
//

inline
MemoryBlock::MemoryBlock( size_t ObjLen, size_t ObjNum )
	: len( ObjLen * ObjNum ), link( 1L )
{
	assert( len > 0 );

	mem = malloc( len );
	if( ! mem )
		FatalError( "Not enough memory." );

#ifndef NDEBUG
	count++;
#endif
}


inline
MemoryBlock::MemoryBlock( const MemoryBlock &mb )
	: len( mb.len ), link( 1L )
{
	mem = malloc( len );
	if( ! mem )
		FatalError( "Not enough memory." );
	memcpy( mem, mb.mem, len );

#ifndef NDEBUG
	count++;
#endif
}


//------------------------------------------------------------------------------
//
//	Object destruction methods.
//

inline
MemoryBlock::~MemoryBlock( void )
{
	assert( --count >= 0 );
	assert( link == 0 );

	free( mem );
}


inline
long MemoryBlock::UnLink( void )
{
	assert( --link >= 0 );

	return link;
}


//------------------------------------------------------------------------------
//
//	Access method(s).
//

inline
MemoryBlock::operator void *( void )
{
	return mem;
}


//------------------------------------------------------------------------------
//
//	Manipulation method(s).
//

inline
long MemoryBlock::Link( void )
{
	return ++link;
}


inline
size_t MemoryBlock::Len( void )
{
	return len;
}


#ifndef NDEBUG
inline
long MemoryBlock::Count( void )
{
	return count;
}
#endif


//
//	End of definitions of class "MemoryBlock" inline member functions.
//
//==============================================================================

#endif
