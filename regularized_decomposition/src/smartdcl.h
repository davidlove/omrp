/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose header file with templates.
PROJECT CODE:		SMART PTR
PROJECT FULL NAME:	Implementation of smart pointers that would guarantee full
					control of memory array and pointer operations.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	-----------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	smartdcl.h
CREATED:			1994.04.11
LAST MODIFIED:		1995.05.14

DEPENDENCIES:		none

--------------------------------------------------------------------------------

HEADER CONTENTS:
	Names of class templates implementing smart pointers and arrays are
declared. Declarations of classes themselves are given in "smartptr.h" header
file.

------------------------------------------------------------------------------*/

#ifndef __SMARTDCL_H__
#define __SMARTDCL_H__

template < class Obj >
class Array;

template < class Obj >
class Ptr;


#endif
