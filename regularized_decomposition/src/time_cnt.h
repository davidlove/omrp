/*------------------------------------------------------------------------------
	PROJECT CODE:		Simplex
	PROJECT FULL NAME:	Advanced implementation of revised simplex method
						for large scale linear problems.

	AUTHOR:				Artur Swietanowski (solver driver),
						prof. Andrzej Ruszczynski (factorization routines).

	PROJECT SUPERVISOR:	dr Jacek Gondzio.

--------------------------------------------------------------------------------

	HEADER FILE NAME:	time_cnt.h
	CREATED:			1993.06.13
	LAST MODIFIED:		1996.10.07

	DEPENDENCIES:		smartptr.h
						<sys/times.h> (in UNIX System V or later)
						<time.h> (in DOS v. 3.30 or later)
						<stdlib.h>

--------------------------------------------------------------------------------

	HEADER CONTENTS:
		This header contains declarations and definitions needed when UNIX
		System V time measuring functions (namely "times()") are to be used.
		The "time_cnt.cc" module defines functions needed for measuring 
		performance of simplex algorithm.

------------------------------------------------------------------------------*/
#ifndef __TIME_CNT_H__
#define __TIME_CNT_H__


#if defined( solaris ) && defined( gnucc )
	extern "C" {
#	include <sys/times.h>
	}
#else
#	include <time.h>
#endif

#include <stdlib.h>

#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif

//------------------------------------------------------------------------------
//
//	Class TimeInfo
//		Stores information on execution times of consequent phases of program's
//		work.
//
class TimeInfo
{
#if defined( solaris ) && defined( gnucc )
	Array<tms>
#else
//	Array<time_t>
//@BEGIN---------

	Array<clock_t>

//@END-----------

#endif
	Mark;

	size_t Events;

public:
	TimeInfo( size_t NumberOfEvents );

	void MarkTime( int WhatTime );
	double TimeDifference( int From, int To );
};

#endif
