/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose code
PROJECT CODE:		--------------------
PROJECT FULL NAME:	--------------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	--------------------

--------------------------------------------------------------------------------

HEADER FILE NAME:	cust_opt.h
CREATED:			1995.10.27
LAST MODIFIED:		1995.11.01

DEPENDENCIES:		stdtype.h, smartptr.h, my_defs.h, option.h
					<stdlib.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	A number of classes representing different types of command line options
used in the context of the linear (incl. stochastic) programming. Options for
a simplex optimizer, a standalone presolver/postsolver package, a regularized
decomposition package are all placed here.

------------------------------------------------------------------------------*/

#ifndef __CUST_OPT_H__
#define __CUST_OPT_H__

#ifndef __OPTION_H__
#	include "option.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __PRE_CODE_H__
#	include "pre_code.h"
#endif


//------------------------------------------------------------------------------
//	File input/output options. Types of files are:
//	-	MPS (I/O) for a simplex optimizer or a modular presolver,
//	-	DIT linear problem or LP solution (I/O) for a simplex optimizer or a
//		modular presolver and postsolver,
//	-	text file containing a solution,
//	-	text file storing all presolve actions (for communication between a
//		presolver and a postsolver),
//	-	text file storing the optimal basis of a linear program (used by the
//		simplex optimizer),
//
struct MPS_IN_Option : public OptionWithStoredArgument
	{ MPS_IN_Option( void ) : OptionWithStoredArgument( "mps_in" ) {} };

struct MPS_LP_Option : public OptionWithStoredArgument
	{ MPS_LP_Option( void ) : OptionWithStoredArgument( "mps_lp" ) {} };

struct MPS_OUT_Option : public OptionWithStoredArgument
	{ MPS_OUT_Option( void ) : OptionWithStoredArgument( "mps_out" ) {} };

struct ACTION_Option : public OptionWithStoredArgument
	{ ACTION_Option( void ) : OptionWithStoredArgument( "action" ) {} };

struct DIT_IN_Option : public OptionWithStoredArgument
	{ DIT_IN_Option( void ) : OptionWithStoredArgument( "dit_in" ) {} };

struct DIT_LP_Option : public OptionWithStoredArgument
	{ DIT_LP_Option( void ) : OptionWithStoredArgument( "dit_lp" ) {} };

struct DIT_OUT_Option : public OptionWithStoredArgument
	{ DIT_OUT_Option( void ) : OptionWithStoredArgument( "dit_out" ) {} };

struct DIT_SOL_Option : public OptionWithStoredArgument
	{ DIT_SOL_Option( void ) : OptionWithStoredArgument( "dit_sol" ) {} };

struct TXT_OUT_Option : public OptionWithStoredArgument
	{ TXT_OUT_Option( void ) : OptionWithStoredArgument( "txt_out" ) {} };

struct TXT_SOL_Option : public OptionWithStoredArgument
	{ TXT_SOL_Option( void ) : OptionWithStoredArgument( "txt_sol" ) {} };

struct ERR_Option : public OptionWithStoredArgument
	{ ERR_Option( void ) : OptionWithStoredArgument( "err" ) {} };

struct LOG_Option : public OptionWithStoredArgument
	{ LOG_Option( void ) : OptionWithStoredArgument( "log" ) {} };

struct BASIS_Option : public OptionWithStoredArgument
	{ BASIS_Option( void ) : OptionWithStoredArgument( "basis" ) {} };


//------------------------------------------------------------------------------
//	Verbosity level option used by all programs of the LP package. There are
//	four verbosity level options: V_NONE, V_LINE, V_LOW, V_HIGH. The V_LINE
//	verbosity level option is not always available.
//
struct VERBOSITY_Option : public MultiStateOption
{
private:
	static const char *labels[];
	static const int states[];

public:
	VERBOSITY_Option( Bool_T line = False )
		: MultiStateOption( "v", line ? 5 : 3, labels, states )
		{ SetDefault( V_LOW ); }
};


//------------------------------------------------------------------------------
//	Presolve mode option is used by the simplex optimizer as well as by the
//	standalone presolver. When called from the presolver, the option is
//	preceeded by '-mode', while the simplex optimizer uses '-presolve' keyword.
//
struct PRESOLVE_MODE_Options : public BitmapOption
{
private:
	static const char *labels[];
	static const int states[];

public:
	PRESOLVE_MODE_Options( const char *keyword  )
		: BitmapOption( keyword, 14, labels, states )
		{ SetDefault( LPR_NONE ); }
};


//------------------------------------------------------------------------------
//	Pricoing mode option is used by the simpelx optimizer and possibly by
//	other programs that use the simplex algorightm. It allows the user to
//	specify the preferred pricing method.
//
struct PRICING_MODE_Option : public MultiStateOption
{
private:
	static const char *labels[];
	static const int states[];

public:
	PRICING_MODE_Option( int def = PRS_ASE )
		: MultiStateOption( "pricing", 3, labels, states )
		{ SetDefault( def ); }
};

//------------------------------------------------------------------------------
//	The following two options are used for specyfying numeric parameters for
//	the simplex optimizer as well as possibly for other optimization algorithms.
//
struct FTOL_Option : public FloatValueOption
{
	FTOL_Option( void ) : FloatValueOption( "ftol", 1e-8 )
		{ SetLimits( 1e-16, 1e-2 ); }
};

struct OPTOL_Option : public FloatValueOption
{
	OPTOL_Option( void ) : FloatValueOption( "optol", 1e-8 )
		{ SetLimits( 1e-16, 1e-2 ); }
};

struct PENALTY_Option : public FloatValueOption
{
	PENALTY_Option( void ) : FloatValueOption( "penalty", 1e-2 )
		{ SetLimits( 0.0, 1e+10 ); }
};


//------------------------------------------------------------------------------
//
//

//------------------------------------------------------------------------------
//
//

//------------------------------------------------------------------------------
//
//

//------------------------------------------------------------------------------
//
//

#endif
