/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose mathematical module.
PROJECT CODE:		REGULARIZED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Amdrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	rand01.cpp
CREATED:			1994.08.13
LAST MODIFIED:		1996.09.14

DEPENDENCIES:		rand01.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This source file contains function definitions and static data of class
"Random01". Class "Random01" generates pseudo random numbers rectangularly
distributed between 0 and 1. Integer arithmetic up to 30323 is required.
	Algorithm from 183 Applied Statistics (1982) vol. 31, No. 2

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void Random01::Seed( int a, int b, int c )
	double Random01::Next( void )

STATIC FUNCTIONS:
	None.

STATIC DATA:
	int Random01::ix,		// The seed of the random number generator.
		Random01::iy,		//
		Random01::iz		//

------------------------------------------------------------------------------*/

#include <assert.h>

/* RS addendum */
#include "mersenne_twister.h"
/* end */

#ifndef __RAND01_H__
#	include "rand01_twister.h"
#endif

//------------------------------------------------------------------------------
//	Static data for class 'Random01'.
//

//	The seed for the random number generator.
//	You can use arbitrary numbers between 1 and 29999.
//	Four alternative seeds are given (this was just needed for numerical tests
//	using different samples. Any set is good.

int Random01::ix	= 2907,
	Random01::iy	= 1951,
	Random01::iz	= 706;

/* Rs addendum */
MTRand Random01::mtrand1;
/* end */

//int Random01::ix	= 12,
//	Random01::iy	= 26333,
//	Random01::iz	= 18976;

//int Random01::ix	= 8794,
//	Random01::iy	= 22101,
//	Random01::iz	= 22300;

//int Random01::ix	= 28765,
//	Random01::iy	= 12345,
//	Random01::iz	= 3765;

//
//	End of static data for class 'Random01'.
//------------------------------------------------------------------------------


/*------------------------------------------------------------------------------

	void Random01::Seed( int a, int b, int c )

PURPOSE:
	Sets the seed for the random number generator.

PARAMETERS:
	int a, int b, int c
		Integer numbers in range 1 - 29999. They will be used as random
		generator's seed.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Random01::Seed( int a, int b, int c )
{
	assert( a >= 1 && a <= 29999 );		ix = a;
	assert( b >= 1 && b <= 29999 );		iy = b;
	assert( c >= 1 && c <= 29999 );		iz = c;
	
	
	 //MTRand::uint32 seed[ MTRand::N ];
	 //for( int n = 0; n < MTRand::N; ++n )
            //seed[n] = 23 * n;  // fill with anything
	 //MTRand mtrand1( seed );
	
}


/*------------------------------------------------------------------------------

	double Random01::Next( void )

PURPOSE:
	Generate next pseudo-random number in the range 0-1.

PARAMETERS:
	None.

RETURN VALUE:
	A random nuber in the range [0;1] inclusive.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

double Random01::Next( void )
{
	double r;

	/*ix = 171 * ( ix % 177 ) -  2 * ( ix / 177 );
	iy = 172 * ( ix % 176 ) - 35 * ( iy / 176 );
	iz = 170 * ( ix % 178 ) - 63 * ( iz / 178 );

	if( ix < 0 ) ix += 30269;
	if( iy < 0 ) iy += 30307;
	if( iz < 0 ) iz += 30323;

	r = ix / 30269.0 + iy / 30307.0 + iz / 30323.0;
	r -= (int)r;

	assert( r >= 0.0 && r <= 1.0 );

 	return r; */
	
	r = mtrand1.randDblExc();
	
	return r; 
	
}

// David Love -- Functions below used to implement overlapping batches.

void Random01::GetSeed( MTRand::uint32* seedArray )
{
   mtrand1.save( seedArray );
   return;
}

int Random01::GetSeedx( void )
{
   return ix;
}

int Random01::GetSeedy( void )
{
   return iy;
}

int Random01::GetSeedz( void )
{
   return iz;
}

void Random01::ReSeed( MTRand::uint32* seedArray )
{
   mtrand1.load( seedArray );
   return;
}

void Random01::ReSeed( int a, int b, int c )
{
   ix = a;
   iy = b;
   iz = c;
}
