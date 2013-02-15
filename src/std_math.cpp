/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose maths functions.
PROJECT CODE:		------------
PROJECT FULL NAME:	------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	std_math.cpp
CREATED:			1994.03.30
LAST MODIFIED:		1994.03.30

DEPENDENCIES:		std_tmpl.h, std_math.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	The file contains only a definition of a static member variable of the
"StdMath" class. See "std_math.h" file for more information on the class as well
as the static variable declared here.
	

------------------------------------------------------------------------------*/

#include "std_math.h"

double StdMath::SMALL_ELEM = 1.0e-10;
