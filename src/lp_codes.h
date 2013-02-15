/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	lp_codes.h
CREATED:			1993.09.16
LAST MODIFIED:		1995.05.08

DEPENDENCIES:		none,
					<none>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	This header file contains symbolic dafinitions (enumeratios) of variable
and constraint types in a general LP problem (read from an MPS file and then
processed into a solvable form).

------------------------------------------------------------------------------*/

#ifndef __LP_CODES_H__
#define __LP_CODES_H__

//==============================================================================
//
//	Variable and constraint type symbolic constants' definitions.
//
//==============================================================================

//------------------------------------------------------------------------------
//	Row types (as read from MPS file).
//
enum
{
	RT_LE			= 0x01,	// Less-than non-equality row.
	RT_GE			= 0x02,	// Greater-than non-equality row.
	RT_EQ			= 0x04,	// Equality row.
	RT_FR			= 0x08,	// Free (non-binding) row.
	RT_RNG			= 0x10,	// Range on the row (Bit mask to be OR'ed with
							// ordinary row types - LE, GE, EQ, FR).
	RT_UNDEFINED	= 0x00,	// Used to mark unknown row types (during reading
							// the MPS file.
	RT_TYPE			= RT_LE | RT_GE | RT_EQ | RT_FR
};

//------------------------------------------------------------------------------
//	Variable types (as read from MPS file). Stored as a bit mask.
//
enum
{
	VTM_FR			= 0x0001,	// Free variable.
	VTM_MI			= 0x0002,	// Non-positivity declared.
	VTM_PL			= 0x0004,	// Non-negatitivity declared.

	VTM_LO			= 0x0008,	// Finite lower bound declared.
	VTM_UP			= 0x0010,	// Finite upper bound declared.
	VTM_FX			= 0x0020,	// Variable declared fixed.
	VTM_BV			= 0x0040,	// Binary variable (used in mixed integer
								// problems). Lower bound of 0 and upper of 1
								// assumed.

	VTM_UNDEFINED	= 0x0080	// Variable whose name was not mentioned in
								// BOUNDS section is initially given this code.
};

//------------------------------------------------------------------------------
//	Variable types (after preprocessing into standard form). Stored as a bit
//	mask.
//
enum
{
	//--------------------------------------------------------------------------
	//	Bit primitives. The more elaborate variable descriptions enumerated
	//	below (in the next section) are constructed of those primitives.
	//
	VT_LO			= 0x01,	// Variable posessing finite lower bound.
	VT_UP			= 0x02,	// Variable posessing finite upper bound.
	VT_FX			= 0x04,	// Fixed variable.

	VT_EXTRA_1		= 0x10,	// These two bits may be used for other
	VT_EXTRA_2		= 0x20,	// purposes.

	VT_ARTIF		= 0x40,	// Artificial variable.
	VT_SLACK		= 0x80,	// Slack variable (also artificial slack added in
							// the crash procedure).

	//--------------------------------------------------------------------------
	//	These are variable types usable for e.g. simplex or interior point
	//	solver. They are compared to using bitwise logical AND.
	//
	VT_FREE			= 0x00,
	VT_NORM			= VT_LO,
	VT_MI			= VT_UP,
	VT_FIXED		= VT_FX | VT_LO | VT_UP,
	VT_BOUNDED		= VT_LO | VT_UP,

	VT_HAS_LO_BND	= VT_LO,
	VT_HAS_UP_BND	= VT_UP,

	//--------------------------------------------------------------------------
	// A pair of mark flag and unmark pattern used by simplex solver.
	//
	VT_MARKED		= VT_EXTRA_1,
	VT_UNMARK		= ~VT_MARKED,

	//--------------------------------------------------------------------------
	//	Bit mask that - after bitwise logical multiplication - leaves only code
	//	meaningful for an LP solver.
	//
	VT_TYPE			= VT_LO | VT_UP | VT_FX
};

//==============================================================================
//
//	End of variable and constraint type symbolic constants definitions
//
//==============================================================================

#endif
