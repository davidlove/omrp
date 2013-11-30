/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose header file with templates.
PROJECT CODE:		SMART PTR
PROJECT FULL NAME:	Implementation of smart pointers that would guarantee full
					control of memory array and pointer operations.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	smartptr.h
CREATED:			1994.04.11
LAST MODIFIED:		1996.09.21

DEPENDENCIES:		error.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h
					<alloc.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A set of class templates implementing smart pointers and arrays is proposed.
Both the smart arrays and the smart pointers have the ability to check array
index before even attempting actual memory access. This control is enabled when
macro "NDEBUG" is not defined.

------------------------------------------------------------------------------*/

#ifndef __SMARTPTR_H__
#define __SMARTPTR_H__


#include <assert.h>

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef NDEBUG
#	ifndef	__MEMBLOCK_H__
#		include "memblock.h"
#	endif
#endif
#ifndef __SMARTDCL_H__
#	include "smartdcl.h"
#endif


#ifdef NDEBUG
#	define SMART_PTR_INLINE inline
#else
#	define SMART_PTR_INLINE
#endif


//
//  DDDD  EEEEE  CCC  L      AAA  RRRR   AAA  TTTTT III  OOO  NN   N  SSS
//  D   D E     C   C L     A   A R   R A   A   T    I  O   O N N  N S   S
//  D   D E     C     L     A   A R   R A   A   T    I  O   O N N  N S
//  D   D EEE   C     L     AAAAA RRRR  AAAAA   T    I  O   O N  N N  SSS
//  D   D E     C     L     A   A R R   A   A   T    I  O   O N  N N     S
//  D   D E     C   C L     A   A R  R  A   A   T    I  O   O N   NN S   S
//  DDDD  EEEEE  CCC  LLLLL A   A R   R A   A   T   III  OOO  N    N  SSS
//


//==============================================================================
//
//	Class template "SmartPointerBase<class Obj>" declaration.
//

template < class Obj >
class SmartPointerBase
{
//------------------------------------------------------------------------------
//	Data members.
//
	public:

#ifndef NDEBUG
	Obj *start, *end;		// Pointer range ( start <= ptr <= end )
	MemoryBlock *mem;		// Associated memory block.
#else
	Obj *ptr;				// Raw pointer.
#endif

//------------------------------------------------------------------------------
//	Methods.
//

//--------------------------------------
//	Creation and destruction methods.
//
protected:
	SmartPointerBase( void );
	SmartPointerBase( const SmartPointerBase<Obj> &spb );

public:
	virtual ~SmartPointerBase( void );

//--------------------------------------
//	Create / delete a link to memory block "*mem".
//
#ifndef NDEBUG
protected:
	void UnLink( void );
	void Link( void );
#endif

//--------------------------------------
//	Assignment operator.
//
public:
	const SmartPointerBase<Obj>& operator =( const SmartPointerBase<Obj> &spb );
};

//
//	End of class template "SmartPointerBase" declaration.
//
//==============================================================================


//==============================================================================
//
//	Class template "Array<Obj>" declaration.
//

template < class Obj >
class Array : public SmartPointerBase<Obj>
{
public:
	Array( void );
	Array( size_t l );
	Array( size_t l, const Obj &o );
	Array( const SmartPointerBase<Obj> &spb, size_t len );

	//
	//	Copy construction is not allowed! Object length is (in general) not
	//	known when copying, thus copy construction may not take place.
	//
	Array( const SmartPointerBase<Obj> & /* spb */ )
		{ assert( ! "Array copy constructor called" ); }

	virtual ~Array( void );

	Obj & operator []( int i );
	const Obj & operator []( int i ) const;

	Ptr<Obj> operator +( int i );
	Ptr<Obj> operator -( int i );

#ifdef gnucc
	void Fill( const Obj o, size_t len, size_t start = 0 );
#else
	void Fill( const Obj &o, size_t len, size_t start = 0 );
#endif

	//@BEGIN-------------------------------------------------------

	void FillOne (const Obj &o, size_t theOne );

	//END----------------------------------------------------------


	void Resize( size_t l );

	friend void Copy( Array<Obj> &Dst, const Array<Obj> &Src,
		size_t SrcLen, size_t DstLen, int nelem = -1, int start = 0 );
	void Copy( const Array<Obj> &Src, size_t SrcLen, size_t dst_len,
		int nelem = -1, int start = 0 );
};

//
//	End of class template "Array<Obj>" declaration.
//
//==============================================================================


//==============================================================================
//
//	Class template "Ptr<Obj>" declaration.
//
//==============================================================================


template < class Obj >
class Ptr : public SmartPointerBase<Obj>
{
#ifndef NDEBUG
private:
	Obj *ptr;
#endif

public:
	Ptr( void );
	Ptr( const Ptr<Obj> &p );
	Ptr( const Array<Obj> &a );
	Ptr( const Array<Obj> &a, const size_t Start, const size_t End );

	virtual ~Ptr( void );

	void ExtractFragment( const SmartPointerBase<Obj> &spb, const size_t Start,
		const size_t End );

	const Ptr<Obj> & operator =( const Array<Obj> &a );
	const Ptr<Obj> & operator =( const Ptr<Obj> &p );

	Obj & operator []( int i );
	const Obj & operator []( int i ) const;

	Obj & operator *( void );
	const Obj & operator *( void ) const;

	Ptr<Obj> & operator ++( void );		//	Only prefix version defined!
	Ptr<Obj> & operator --( void );		//

	const Ptr<Obj> & operator +=( int i );
	const Ptr<Obj> & operator -=( int i );

	Ptr<Obj> operator +( int i ) const;
	Ptr<Obj> operator -( int i ) const;

#ifndef NDEBUG
	void ToStart( void );
	void ToEnd( void );
#endif
};

//
//	End of class template "Ptr<Obj>" declaration.
//
//==============================================================================



//                                                                            //
//       DDDD  EEEEE FFFFF III NN   N III TTTTT III  OOO  NN   N  SSS         //
//       D   D E     F      I  N N  N  I    T    I  O   O N N  N S   S        //
//       D   D E     F      I  N N  N  I    T    I  O   O N N  N S            //
//       D   D EEE   FFF    I  N  N N  I    T    I  O   O N  N N  SSS         //
//       D   D E     F      I  N  N N  I    T    I  O   O N  N N     S        //
//       D   D E     F      I  N   NN  I    T    I  O   O N   NN S   S        //
//       DDDD  EEEEE F     III N    N III   T   III  OOO  N    N  SSS         //
//                                                                            //


#	ifndef NDEBUG
#		include "sptr_deb.h"
#	else
#		include "sptr_ndb.h"
#	endif


#endif
