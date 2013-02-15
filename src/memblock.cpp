/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose block defining a memory management class.
PROJECT CODE:		MEMORY BLOCKS
PROJECT FULL NAME:	Implementation of memory block manager.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	memblock.cpp
CREATED:			1994.04.11
LAST MODIFIED:		1996.09.20

DEPENDENCIES:		memblock.h, error.h,
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This source file contains non-inline definitions of some of the
"MemoryBlock" member functions. For explanation of class "MemoryBlock" purpose
and description of its usage - see header file "memblock.h".

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void MemoryBlock::Resize( size_t ObjLen, size_t ObjNum );

STATIC FUNCTIONS:
	None.

STATIC DATA (INCL. STATIC CLASS MEMBERS):
	long int MemoryBlock::count
		A count of allocated memory blocks. It is initialised at zero. Should
		never be decremented below zero. Only when NDEBUG macro is not defined
		this static data member of class "MemoryBlock" is actually created.

------------------------------------------------------------------------------*/


#include <assert.h>

#include "memblock.h"


#ifndef NDEBUG
	long int MemoryBlock::count = 0L;
#endif


/*------------------------------------------------------------------------------

	void MemoryBlock::Resize( size_t ObjLen, size_t ObjNum )

PURPOSE:
	This function reallocates an existing block of memory. The new size of block
will be "ObjLen * ObjNum".

PARAMETERS:
	size_t ObjLen
		Size of a single object to be held in the table (obtained e.g. by
		"sizeof()" operator.
	
	size_t ObjNum
		Number of objects to be held in the table.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void MemoryBlock::Resize( size_t ObjNum, size_t ObjSize )
{
	len = ObjSize * ObjNum;

	assert( link == 1 && len > 0 );

	mem = realloc( mem, len );
	if( ! mem )
		FatalError( "Not enough memory." );
}
