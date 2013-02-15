/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose - time counter.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

AUTHOR:				Artur Swietanowski (solver driver),
					prof. Andrzej Ruszczynski (factorization routines).

PROJECT SUPERVISOR:	dr Jacek Gondzio.

--------------------------------------------------------------------------------

MODULE NAME:		time_cnt.cpp
CREATED:			1993.06.02
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		time_cnt.h, smartptr.h

--------------------------------------------------------------------------------

MODULE FUNCTIONS:
	PUBLIC:
		TimeInfo::TimeInfo()
		TimeInfo::MarkTime()
		TimeInfo::TimeDifference()

--------------------------------------------------------------------------------

MODULE PURPOSE:
	This file contains to simple time measuring functions (based on UNIX
	System V 'times' function and DOS 'time' function). 'MarkTime' stores
	current time (in clock ticks) in one of the fields of the TimeInfo class
	object. 'TimeDifference' computes the time between two events (stored
	in fields of class object) and returns result in seconds.

	The module is compiled only when TIMER_ON macro is declared.

------------------------------------------------------------------------------*/

#ifndef __TIME_CNT_H__
#	include "time_cnt.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif


#if defined( solaris ) && defined( gnucc )
	extern "C"
	{
#	include <unistd.h>
	}

	extern "C"
	{
		clock_t times( tms *tmsPtr );
		long sysconf( int name );
	}
#endif

//------------------------------------------------------------------------------
//	Implicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )
#	if defined( solaris )
		template class SmartPointerBase<tms>;
		template class Array<tms>;
		template class Ptr<tms>;

		template tms *MALLOC( tms *& Table, size_t len );
		template tms *REALLOC( tms *& Table, size_t len );
		template void FREE( tms *& Table );
#	elif !defined( ibmpc386 )
		template class SmartPointerBase<time_t>;
		template class Array<time_t>;
		template class Ptr<time_t>;

		template time_t *MALLOC( time_t *& Table, size_t len );
		template time_t *REALLOC( time_t *& Table, size_t len );
		template void FREE( time_t *& Table );
#	endif
#endif
//
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//
//	TimeInfo::TimeInfo( size_t NumberOfEvents )
//
//		Class "TimeInfo" construtor. Takes the maximum number of remembered
//		events as an argument. The number of events has to be positive.
//
TimeInfo::TimeInfo( size_t NumberOfEvents )
	: Mark( NumberOfEvents > 0 ? NumberOfEvents : 1 ), Events( NumberOfEvents )
{
	assert( Events > 1 );
}


//------------------------------------------------------------------------------
//
//	TimeInfo::MarkTime()
//		Stores the current time in one of the fields of the TimeInfo object.
//		TimeInfo object can store times of six predefined events. The event
//		number is passed in 'WhatTime' parameter. See "time_cnt.h" for
//		definitions of 'TimeInfo' class and meaning and codes of events.
//
void TimeInfo::MarkTime( int WhatTime )
{
	assert( WhatTime >= 0 && WhatTime < (int)Events );

#if defined( solaris ) && defined( gnucc )
	times( &(Mark[WhatTime]) );
#else
	//@BEGIN----------
	//	Mark[WhatTime] = time( NULL );  -->This was original

	Mark[WhatTime] = clock(); 

	//@END-------------
#endif
}


//------------------------------------------------------------------------------
//
//	TimeInfo::TimeDifference
//		Computes time in seconds that elapsed between to specific events. Times
//		are assumed to have been stored with 'MarkTime' calls. No checking
//		of times are performed. They may either be uninitialized, or erroneous.
//		In such case the behaviour of this function may be undetermined.
//
double TimeInfo::TimeDifference( int From, int To )
{
	assert( From >= 0 && From < (int)Events && To > From && To < (int)Events );

#if defined( solaris ) && defined( gnucc )
	clock_t tStart	= Mark[From].tms_utime,
		tEnd		= Mark[To].tms_utime;
	long clockTicks	= sysconf( _SC_CLK_TCK );

	assert( clockTicks != -1 );

	return ( (double)tEnd - (double)tStart ) / (double)(clockTicks);
#else
	//@BEGIN----------
	
	//return difftime( Mark[To], Mark[From] );  -->This was original

	return ( (double)Mark[To] - (double)Mark[From] ) / (double)(CLOCKS_PER_SEC);

	//@END------------
#endif
}
