/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose mathematical module.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	rand01.h
CREATED:			1994.08.13
LAST MODIFIED:		1995.08.23

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	Class "Random01" declaration. The class is a random number generator, which
generates numbers in range 0 to 1. It has only static members. Thus the only
reason why it is a class is data protection.

------------------------------------------------------------------------------*/

#ifndef __RAND01_H__
#define __RAND01_H__

#include "mersenne_twister.h"
//==============================================================================
//
//	Class "Random01" declaration.
//
//==============================================================================

class Random01
{
private:
        static MTRand mtrand1;

	Random01( void );	// In order to make it impossible to create objects
						// of this type.

public:
	static void Seed( int s1, int s2, int s3 );
	static double Next( void );
        // David Love -- Functions below used to implement overlapping batches
        static void GetSeed( MTRand::uint32* seedArray );
        static void ReSeed( MTRand::uint32* seedArray );
};

//==============================================================================
//
//	End of class "Random01" declaration.
//
//==============================================================================

#endif
