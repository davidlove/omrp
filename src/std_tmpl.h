/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose template definitions.
PROJECT CODE:		-------------------
PROJECT FULL NAME:	-------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	std_tmpl.h
CREATED:			1993.09.26
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		<string.h>, <assert.h>, <math.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains some simple general purpose templates, namely
"Min", "Max" and "Swap". All are implemented inline.
	"Min" and "Max" are designed to work well on simple types (they take "T"
not "&T" as arguments).  They rely on existence of "operator <" and "operator
>" respectively.
	"Swap" exchanges two objects to which references are passed. It relies on
existence of a copy constructor and "operator =" for the class and on their
"compatibility" (i.e. that both initialization of one class object "t1" with
another object "t2" and assignment "t1 = t2" would yield the same value of
"t1").

------------------------------------------------------------------------------*/

#ifndef __STD_TMPL_H__
#define __STD_TMPL_H__

#include <string.h>
#include <assert.h>
#include <math.h>

//==============================================================================
//
//	Function templates.
//
//==============================================================================

//------------------------------------------------------------------------------
//	'Min' / 'Max' templates and 'Swap' operation.
//	WARNING: 'Min' and 'Max' (as they are provided now) can be effectively used
//	only for simple data types. Also 'Swap' will be inefficient for large
//	structures - for such structures a type-specific 'Swap' should be added.
//
template <class T>
inline
T Min( T t1, T t2 )
	{ return ( t1 < t2 ) ? t1 : t2; }


template <class T>
inline
T Max( T t1, T t2 )
	{ return ( t1 > t2 ) ? t1 : t2; }


template <class T>
inline
void Swap( T &t1, T &t2 )
	{ T tmp = t1; t1 = t2; t2 = tmp; }


template <class T>
inline
T Abs( T t )
	{ return ( t >= 0 ) ? t : -t; }


//------------------------------------------------------------------------------
//	Simple vector operations.
//
template <class T, class L>
inline
void VecCopy( T *dst, const T *src, L len )
{
	assert( dst != NULL || src != NULL );

	memcpy( (void *)dst, (void *)src, len * sizeof( T ) );
}


//------------------------------------------------------------------------------
//	Arithmetic functions.
//
template <class T>
inline
T Sqr( const T x )
	{ return x*x; }

//------------------------------------------------------------------------------
//	Simple array templates.
//
template <class T, class L>
#if defined( watcom )
void Fill( T *array, L len, const T val )
#else
void Fill( T *array, L len, const T &val )
#endif
{
	for( size_t i = 0, l = (size_t) len; i < l; i++, array++ )
		*array = val;
}


//==============================================================================
//
//	End of function templates.
//
//==============================================================================

#endif
