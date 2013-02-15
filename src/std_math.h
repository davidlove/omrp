/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose maths functions.
PROJECT CODE:		------------
PROJECT FULL NAME:	------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	std_math.h
CREATED:			1994.01.12
LAST MODIFIED:		1995.01.22

DEPENDENCIES:		std_tmpl.h,
					<math.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	The file contains definitions of some simple inline functions designed to
provide functionality of arithmetic comparisons with a zero tolerance. The zero
tolerance 'tol' is such a positive, small number, that if

	-tol <= x <= tol

then we consider 'x' to be equal to zero. From this basic definition other
definitions may be derived:
*	definition of a positive / negative / non-positive / non-negative number,
*	definition of equal / unequal numbers.

	The functions in this file provide all such comparison predicates.
Additionally, a (logically) global variable which stores the current value of
the zero tolerance is added. Naturally, methods for setting and inspecting this
tolerance are also provided.

	A class is created in order to make it possible to hide the "SMALL_ELEM"
zero tolerance and thus protect it from possible abuse.

------------------------------------------------------------------------------*/

#ifndef __STD_MATH_H__
#define __STD_MATH_H__

#include <math.h>

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif


//==============================================================================
//
//	Class "StdMath" declaration
//
//==============================================================================

class StdMath
{
private:
	static double SMALL_ELEM;

public:
	static void SetZeroTolerance( double tol );
	static double GetZeroTolerance( void );

	static int IsZero( double x );
	static int IsNonZero( double x );
	static int IsPositive( double x );
	static int IsNegative( double x );
	static int IsEqual( double x, double y );
	static int IsNotEqual( double x, double y );
};

//==============================================================================
//
//	INLINE IMPLEMENTATIONS OF ALL OF THE FUNCTIONS OF THE CLASS
//
//==============================================================================

inline
double StdMath::GetZeroTolerance( void )
	{ return SMALL_ELEM; }

inline
void StdMath::SetZeroTolerance( double tol )
{
	tol = fabs( tol );
	if( tol > 0.0 && tol < 1.0 ) SMALL_ELEM = tol;
}

inline
int StdMath::IsZero( double x )
	{ return fabs( x ) <= SMALL_ELEM; }

inline
int StdMath::IsNonZero( double x )
	{ return fabs( x ) > SMALL_ELEM; }

inline
int StdMath::IsPositive( double x )
	{ return x > SMALL_ELEM; }

inline
int StdMath::IsNegative( double x )
	{ return x < - SMALL_ELEM; }

inline
int StdMath::IsEqual( double x, double y )
{
	double m = Max( fabs( x ), fabs( y ) );

	return IsZero( m ) || IsZero( ( x - y ) / m );
}

inline
int StdMath::IsNotEqual( double x, double y )
{
	double m = Max( fabs( x ), fabs( y ) );

	return IsNonZero( m ) && IsNonZero( fabs( x - y ) / m );
}

//==============================================================================
//
//	INLINE IMPLEMENTATIONS OF THE CLASS METHODS' INTERFACE.
//
//==============================================================================

inline
double GetZeroTolerance( void )
	{ return StdMath::GetZeroTolerance(); }

inline
void SetZeroTolerance( double tol )
	{ StdMath::SetZeroTolerance( tol ); }

inline
int IsZero( double x )
	{ return StdMath::IsZero( x ); }

inline
int IsNonZero( double x )
	{ return StdMath::IsNonZero( x ); }

inline
int IsPositive( double x )
	{ return StdMath::IsPositive( x ); }

inline
int IsNegative( double x )
	{ return StdMath::IsNegative( x ); }

inline
int IsEqual( double x, double y )
	{ return StdMath::IsEqual( x, y ); }

inline
int IsNotEqual( double x, double y )
	{ return StdMath::IsNotEqual( x, y ); }


//==============================================================================
//
//	END OF THE INLINE IMPLEMENTATIONS
//
//==============================================================================

#endif
