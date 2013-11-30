/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose - work vector manager
PROJECT CODE:		--------------------------
PROJECT FULL NAME:	--------------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	vec_pool.cpp
CREATED:			1994.04.20
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		smartptr.h, sptr_ndb.h, sptr_deb.h, error.h, vec_pool.h,
					std_tmpl.h
					<stdlib.h>, <string.h>, <assert.h>

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

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/

#include <string.h>
#include <assert.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __VEC_POOL_H__
#	include "vec_pool.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif


//------------------------------------------------------------------------------
//	Implicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )
	template class SmartPointerBase<WorkVectorPool::VectorHandle *>;
	template class Array<WorkVectorPool::VectorHandle *>;
	template class Ptr<WorkVectorPool::VectorHandle *>;

	template WorkVectorPool::VectorHandle **MALLOC(
		WorkVectorPool::VectorHandle **& Table, size_t len );
	template WorkVectorPool::VectorHandle **REALLOC(
		WorkVectorPool::VectorHandle **& Table, size_t len );
	template void FREE( WorkVectorPool::VectorHandle **& Table );
#endif
//
//------------------------------------------------------------------------------

Array<WorkVectorPool::VectorHandle *> WorkVectorPool::VectorHandles;
size_t WorkVectorPool::len	= 0,
	WorkVectorPool::maxlen	= 0;

/*------------------------------------------------------------------------------

	WorkVectorPool::VectorHandle::~VectorHandle( void )

PURPOSE:
	Destructor of class "WorkVectorPool::VectorHandle" object. Disconnects the
vector handle from the associated memory block if the vector is marked as
unused. If it was not allocated --- does nothing. All other situations must be a
result of an error.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

WorkVectorPool::VectorHandle::~VectorHandle( void )
{
	switch( status )
	{
	case NOT_ALLOC:
		break;

	case USED:
	case LOCKED:
		assert( ! "Attempt to destruct a work vector that is being used" );

	case UNUSED:
		if( mb->UnLink() == 0 ) delete mb;
		break;
	}
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::Initialize( size_t NumberOfVectors )

PURPOSE:
	This function initializes the vector pool manager. It allocates an array of
vector handle poiters. Vectors are not allocated.

PARAMETERS:
	size_t NumberOfVectors
		The initial number of work vector poiters to store.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::Initialize( size_t NumberOfVectors )
{
	assert( NumberOfVectors > 0 && len == 0 && maxlen == 0 );

	len = 0;
	maxlen = Max( unsigned(NumberOfVectors), 2U );

	VectorHandles.Resize( maxlen );
	VectorHandles.Fill( (VectorHandle *)NULL, maxlen );
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::CleanUp( void )

PURPOSE:
	Deallocates all work vectors and zeros the lengths and pointer to array of
vector handles. It essentially works as a destructor would, but of course class
"WorkVectorPool" has only static member data. It is meant to be called at the
end of a program that uses the vector pool manager.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::CleanUp( void )
{
	if( len == 0 && maxlen == 0 ) return;

	for( int i = 0; i < (int)len; i++ )
		delete VectorHandles[i];

	VectorHandles.Resize( 0 );

	len = maxlen = 0;
}


/*------------------------------------------------------------------------------

	int WorkVectorPool::Create( size_t l, int type, MemoryBlock *&mb,
		const char *label )

PURPOSE:
	This function "creates" a work vector of length "l" bytes, of type "type"
and with a label "label". If a vector with the same combination of a nonempty
"label" and "type" already exists, the function exits reporting failure (setting
"mb" to NULL and returning "NO_VECTOR" instead of a vector ID). If no label is
given the vector is assumed to be a temporary object to which it will not be
possible to refer later.
	If possible memory belonging to another vector that was created previously
but is currently not used any more is used. Otherwise a new block of memory is
allocated. In general: we try to keep as few vectors as possible and waste as
little memory as possible.

PARAMETERS:
	size_t l
		Number of bytes required.

	const char *label
	int type
		Label (optional) and type (mandatory, although

	MemoryBlock *&mb

RETURN VALUE:
	Vector ID (non-negative number) on success, "NO_VECTOR" on failure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int WorkVectorPool::Create( size_t l, int type, MemoryBlock *&mb, // )
	const char *label )
{
	int i;

	//--------------------------------------------------------------------------
	//	Phase one: Check if no vector with identical label exists. NULL pointer
	//	in 'label' means no label.
	//
	if( label && *label && ( Find( label, type ) != NO_VECTOR ) )
	{
		mb = NULL;
		return NO_VECTOR;
	}

	//--------------------------------------------------------------------------
	//	Phase two: Try to find an unused but allocated vector of matching size.
	//
	//	"Best match" criterion is employed if possible (to avoid allocations and
	//	avoid wasting memory resources at the same time).
	//

	//
	//	Pass 1:	Look for the smallest memory block that can be used "as is"
	//			(i.e. without reallocation).
	//
	int BestMatch = -1;
	size_t BestSize = 0;
	for( i = 0; ( i < (int)len ) && ( BestSize != l ) ; i++ )
	{
		VectorHandle &vh = *VectorHandles[i];
		size_t size = vh.mb->Len();

		if( vh.status == UNUSED && size >= l )
			if( !BestSize || size < BestSize )
			{
				BestSize = size;
				BestMatch = i;
			}
	}
			
	//
	//	Pass 2:	If pass 1 failed, look for an allocated, but unused block of any
	//			size. Expand the memory block to match the needed size.
	//
	if( BestMatch == -1 )
		for( i = 0; i < (int)len; i++ )
		{
			VectorHandle &vh = *VectorHandles[i];

			if( vh.status == UNUSED )
			{
				BestMatch = i;
				vh.mb->Resize( l );
				break;
			}
		}

	//
	//	Finally: If a block was found in one of the passes, return it.
	//
	if( BestMatch != -1 )
	{
		VectorHandle &vh = *VectorHandles[ BestMatch ];

		if( label )
		{
			strncpy( vh.label, label, LABEL_LEN );
			vh.label[ LABEL_LEN ] = '\0';
		}
		else
			vh.label[0] = '\0';

		vh.status = USED;
		mb = vh.mb;
		vh.type = type;
		return BestMatch;
	}

	//--------------------------------------------------------------------------
	//	Phase three (if phase two failed): create a new work vector and fill
	//	it's handle's data fields.
	//
	if( len >= maxlen )
	{
		size_t newlen = maxlen + Max( int( maxlen/2 ), 10 );

		VectorHandles.Resize( newlen );
		VectorHandles.Fill( (VectorHandle *)NULL, newlen, maxlen );

		maxlen = newlen;
	}

	if( ( VectorHandles[len] = new VectorHandle ) == NULL )
		FatalError( "Not enough memory." );

	VectorHandle &vh = *VectorHandles[len];

	if( label )
	{
		strncpy( vh.label, label, LABEL_LEN );
		vh.label[ LABEL_LEN ] = '\0';
	}
	else
		vh.label[0] = '\0';

	vh.mb = new MemoryBlock( l );

	if( !vh.mb )
		FatalError( "Not enough memory." );

	vh.type = type;
	vh.status = USED;
	mb = vh.mb;
	return len++;
}


/*------------------------------------------------------------------------------

	int WorkVectorPool::Attach( size_t l, int type, MemoryBlock *&mb,
		const char *label )

PURPOSE:
	Attaches to a vecor identified by "type" and "label" (if such already
exists). 

PARAMETERS:
	int type
	const char *label
		Type and label of the vector. Label has to be non-empty.

	size_t l
		Minimum actual length of the vector. Checked only in debugging mode
		(i.e. when macro "NDEBUG" is not defined).

	MemoryBlock *&mb
		Actually one of the return values and not an argument. On return stores
		the pointer to a memory block (if a vector was found) or "NULL".

RETURN VALUE:
	Vector ID on success, otherwise "NO_VECTOR".

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int WorkVectorPool::Attach( // )
#ifdef NDEBUG
	size_t,
#else
	size_t l,
#endif
	int type, MemoryBlock *&mb, const char *label )
{
	int i = ( label && *label ) ? Find( label, type ) : (int)NO_VECTOR;

	mb = NULL;
	if( i == NO_VECTOR ) return NO_VECTOR;

	VectorHandle &vh = *VectorHandles[i];

	if( vh.status != LOCKED && vh.status != USED ) return NO_VECTOR;

	assert( vh.mb->Len() >= l );

	vh.status = USED;
	vh.type = type;
	mb = vh.mb;
	return i;
}


/*------------------------------------------------------------------------------

	int WorkVectorPool::Find( const char *label, int type )

PURPOSE:
	Looks for a work vector with a given combination of type and label. If type
is not given, it returns the first vector with the given label.

PARAMETERS:
	const char *label
	int type
		Type and label of a vector which is to be found in the pool.
		Label has to be non-empty.

RETURN VALUE:
	Vector ID on success, "NO_VECTOR" on failure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int WorkVectorPool::Find( const char *label, int type )
{
	assert( label != NULL && *label != '\0' );

	if( type == NO_TYPE )
		for( int i = 0; i < (int)len; i++ )
		{
			VectorHandle &vh = *VectorHandles[i];

			if( strncmp( vh.label, label, LABEL_LEN ) == 0 )
				return i;
		}
	else
		for( int i = 0; i < (int)len; i++ )
		{
			VectorHandle &vh = *VectorHandles[i];

			if( vh.type == type && strncmp( vh.label, label, LABEL_LEN ) == 0 )
				return i;
		}

	return NO_VECTOR;
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::Detach( const char *label )

PURPOSE:
	Detaching a vector does not destroy it or its data: the vector is preserved
until it is either explicitly destroyed or attached to.
	Public interface to vector detachment: uses the label (which is publicly
available) instead of vector ID. Detaches all vectors with a given label using
"WorkVectorPool::Detach( int id )" as lower level interface.

PARAMETERS:
	const char *label
		Vector label.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::Detach( const char *label )
{
	assert( label != NULL && *label != '\0' );

	int id;

	while( ( id = Find( label ) ) != NO_VECTOR )
		Detach( id );
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::Detach( int id )

PURPOSE:
	Internal vector detachment routine.

PARAMETERS:
	int id
		Vector ID of the vector to be detached.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::Detach( int id )
{
	assert( id >= 0 && id < (int)len &&
		VectorHandles[id]->status != NOT_ALLOC );

	VectorHandles[id]->status = LOCKED;
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::Destroy( const char *label )

PURPOSE:
	Public interface to explicit vector destruction. When a vector is destroyed,
it's data may be overwritten or the memory may be re- or deallocated. All
vectors with a given label are destroyed.

PARAMETERS:
	const char *label
		Vector label.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::Destroy( const char *label )
{
	assert( label != NULL && *label != '\0' );

	int id;

	while( ( id = Find( label ) ) != NO_VECTOR )
		Destroy( id );
}


/*------------------------------------------------------------------------------

	void WorkVectorPool::Destroy( int id )

PURPOSE:
	Internal vector destruction routine. Marks the given vector as unused (and
thus subject to use by other vector, reallocation or deallocation).

PARAMETERS:
	int id
		Vector ID.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void WorkVectorPool::Destroy( int id )
{
	assert( id >= 0 && id < (int)len &&
		VectorHandles[id]->status != NOT_ALLOC );

	VectorHandles[id]->status = UNUSED;
	VectorHandles[id]->type = NO_TYPE;
	VectorHandles[id]->label[0] = '\0';
}


/*------------------------------------------------------------------------------

	int WorkVectorPool::GetNewTypeID( void )

PURPOSE:
	Returns a new (unique) non-negative type identificator.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

int WorkVectorPool::GetNewTypeID( void )
{
	static int LastRegType = NO_TYPE;

	assert( LastRegType >= 0 );
	assert( LastRegType != NO_TYPE || len == 0 );

	return ++LastRegType;
}
