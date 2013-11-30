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
#	include "rand01.h"
#endif

//------------------------------------------------------------------------------
//	Static data for class 'Random01'.
//

//	The seed for the random number generator.
//	You can use arbitrary numbers between 1 and 29999.
//	Four alternative seeds are given (this was just needed for numerical tests
//	using different samples. Any set is good.

/* Rs addendum */
MTRand Random01::mtrand1;
/* end */

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
	
	r = mtrand1.randDblExc();
	
	return r; 
	
}

// David Love -- Functions below used to implement overlapping batches.

void Random01::GetSeed( MTRand::uint32* seedArray )
{
   mtrand1.save( seedArray );
   return;
}

void Random01::ReSeed( MTRand::uint32* seedArray )
{
   mtrand1.load( seedArray );
   return;
}
