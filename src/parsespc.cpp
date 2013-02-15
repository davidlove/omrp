/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - program configuration file reading.
PROJECT CODE:		Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	parsespc.cpp
CREATED:			1992.09.17
LAST MODIFIED:		1995.10.27

DEPENDENCIES:		compile.h, stdtype.h, std_tmpl.h, error.h, parsespc.h,
					lexer.h, simplex.h, solvtol.h
					<stdio.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This module implements reading a specification file that configures a
linear programming solver. 

	This is an example of what the specification file looks like:

= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
$ Recognized verbosity levels: NONE, LINE, LOW, HIGH
$
VERBOSITY LEVEL				HIGH

SMALL ELEMENTS				1.0E-10
LU PIVOT TOLERANCE			1.0E-1
FEASIBILITY TOLERANCE		1.0E-8
OPTIMALITY TOLERANCE		1.0E-8
PIVOT TOLERANCE				1.0E-6

GOOD RESIDUALS				1.0E-10
SATISFACTORY RESIDUALS		1.0E-10
POOR RESIDUALS				1.0E-6
ALARMING RESIDUALS			1.0E-2

GROWTH FACTOR				1.0E4
LENGTH FACTOR				3

$	Integer values expected here (may be expressed in floating point format
$	though). Non-integral value will cause a warning to be issued.
$
REFACTORIZATION FREQENCY	50
RESIDUALS CHECK FREQUENCY	15
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	Some comments on ".spc" file syntax:
*	empty lines and lines beginning with a dollar character are ignored,
*	line ordering is irrelevant,
*	spaces are irrelevant (except for separating spaces),
*	exactly one statement should be placed in one line,
*	you may choose a shorthand form for some or all of the options (see below -
	a list of abbreviations (note that the file with full parameter names is
	much easier to read and understand),
*	if "COMP_LEXER_UPPERCASE" is defined (in "compile.h") all keywords have to
	be uppercase.

	The parameters that are not set in the specification file assume default
values. You may include some, all or none of the above mentioned options in
your specification file.

Table: Abbreviations for parameter names
		+-------------------------------+-----------------------+
		|	Full name					|	Abbreviation		|
		+-------------------------------+-----------------------+
		|	VERBOSITY LEVEL				|	VERB				|
		|								|						|
		|	SMALL ELEMENTS				|	SMALL				|
		|	LU PIVOT TOLERANCE			|	LUPIV				|
		|	FEASIBILITY TOLERANCE		|	FTOL				|
		|	OPTIMALITY TOLERANCE		|	OPTOL				|
		|	PIVOT TOLERANCE				|	PVTOL				|
		|								|						|
		|	GOOD RESIDUALS				|	GRES				|
		|	SATISFACTORY RESIDUALS		|	SRES				|
		|	POOR RESIDUALS				|	PRES				|
		|	ALARMING RESIDUALS			|	ARES				|
		|								|						|
		|	GROWTH FACTOR				|	GRWF				|
		|	LENGTH FACTOR				|	LNGF				|
		|								|						|
		|	REFACTORIZATION FREQENCY	|	RFCT				|
		|	RESIDUALS CHECK FREQUENCY	|	RESCH				|
		+-------------------------------+-----------------------+

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	Spc::Spc()
	Spc::Read()

STATIC FUNCTIONS:
	x

STATIC DATA:
	None.

--------------------------------------------------------------------------------

USED MACROS FROM COMPILE.H AND THEIR MEANING:
	COMP_XXXXX			- x

------------------------------------------------------------------------------*/


#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif

#ifndef __COMPILE_H__
#	include "compile.h"
#endif
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
#ifndef __SOLVTOL_H__
#	include "solvtol.h"
#endif




/*------------------------------------------------------------------------------

	Spc::Spc( void )

PURPOSE:
	Constructor initializes all numerical tolerances to default (sane) values.
Strings (file names) are initialized to NULL pointers.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Spc::Spc( void )
	: MPS_File			( NULL ),
	OutputFile			( NULL ),
	SolutionFile		( NULL ),

	Verbosity			( V_LOW ),

	LU_PIVOT_TOL		( LU_PIVOT_TOL_DEF ),
	FEASIBILITY_TOL		( FEASIBILITY_TOL_DEF ),
	OPTIMALITY_TOL		( OPTIMALITY_TOL_DEF ),
	MIN_STEP_LENGTH		( MIN_STEP_LENGTH_DEF ),
	PIVOT_TOL			( PIVOT_TOL_DEF ),

	GROWTH_FACTOR		( GROWTH_FACTOR_DEF ),
	LENGTH_FACTOR		( LENGTH_FACTOR_DEF ),

	GOOD_RESID			( GOOD_RESID_DEF ),
	SATISF_RESID		( SATISF_RESID_DEF ),
	POOR_RESID			( POOR_RESID_DEF ),
	ALARM_RESID			( ALARM_RESID_DEF ),

	REFACT_FREQ			( REFACT_FREQ_DEF ),
	RESID_CHECK_FREQ	( RESID_CHECK_FREQ_DEF ),

	Pricing 			( PRS_ASE )
{}


/*------------------------------------------------------------------------------

	Spc::~Spc()

PURPOSE:
	Destructor deallocates the only memory dynamically allocated for the class
object: the copies of three file names.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Spc::~Spc()
{
	FREE( MPS_File );
	FREE( OutputFile );
	FREE( SolutionFile );
}
