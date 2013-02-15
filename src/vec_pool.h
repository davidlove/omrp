/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose - work vector manager.
PROJECT CODE:		--------------------------
PROJECT FULL NAME:	--------------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	vec_pool.h
CREATED:			1994.04.20
LAST MODIFIED:		1996.09.20

DEPENDENCIES:		smartptr.h
					<stdlib.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __VEC_POOL_H__
#define __VEC_POOL_H__


#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __MEMBLOCK_H__
#	include "memblock.h"
#endif


//==============================================================================
//
//	Definition of class "WorkVectorPool".
//

class WorkVectorPool
{
	//--------------------------------------------------------------------------
	//	Publicly available enumerations.
	//
public:
	enum
	{
		NO_VECTOR = -1,
		LABEL_LEN = 30
	};

	enum VectorStatus { NOT_ALLOC, USED, UNUSED, LOCKED };
	enum { NO_TYPE = 0 };

	//--------------------------------------------------------------------------
	//	Private data and type definitions: 
	//	-	a work vector handle data structure,
	//	-	array of poiters to vector handles and
	//	-	it's maximum (allocated) length and current (filled) length.
	//
private:
	struct VectorHandle
	{
		VectorStatus status;
		MemoryBlock *mb;
		char label[ LABEL_LEN + 1 ];
		int type;

		VectorHandle( void );
		~VectorHandle( void );
	};

	static Array<VectorHandle *> VectorHandles;
	static size_t len, maxlen;


	//--------------------------------------------------------------------------
	//	Functions used by the higher level template-based interface to work
	//	vectors' pool.
	//
protected:
	static int Create( size_t len, int type, MemoryBlock *&mb,
		const char *label = NULL );
	static void Detach( int id );
	static void Destroy( int id );


	//--------------------------------------------------------------------------
	//	Public available functions.
	//	
public:
	static int Attach( size_t len, int type, MemoryBlock *&mb,
		const char *label );
	static int Find( const char *label, int type = NO_TYPE );
	static void Detach( const char *label );
	static void Destroy( const char *label );

	static void Initialize( size_t NumberOfVectors = 10 );
	static void CleanUp( void );

	//--------------------------------------------------------------------------
	//	Each vector is recognized by it's label and type. The types of work
	//	vectors are registered as they (the vectors) are instantiated. A static
	//	variable holds the number of the last registered type. The associated
	//	function generates a new type id (a number).
	//
protected:
	static int GetNewTypeID( void );
};

//
//	End of definition of class "WorkVectorPool".
//
//==============================================================================


//==============================================================================
//
//	Inline definition of class "WorkVectorPool::VectorHandle" constructor.
//

inline
WorkVectorPool::VectorHandle::VectorHandle( void )
	: status( NOT_ALLOC ), mb( NULL ), type( NO_TYPE )
{
	*label = '\0';
}

#endif
