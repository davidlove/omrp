/*------------------------------------------------------------------------------
MODULE TYPE:		Linear optimization supporting header
PROJECT CODE:		Simplex.
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	lp_sol.cc
CREATED:			1995.03.03
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		lp_sol.h, stdtype.h, smartptr.h, msp_lp.h, sort_lab.h,
					lp_codes.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	XXXXX			- x

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __LP_SOL_H__
#	include "lp_sol.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif


/*------------------------------------------------------------------------------

	LP_Solution::LP_Solution( MPS_LP *lp, Int_T mm, Int_T nn, int cntnts )

PURPOSE:
	Constructor. Stores the data passed by arguments.

PARAMETERS:
	MPS_LP *lp
		Associated linear problem.
	
	int cntnts
		Contents of the solution vector.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

LP_Solution::LP_Solution( MPS_LP *lp, int cntnts )
	:
	SolutionWithLabels(
		( lp != NULL ) ? &lp->RevealColumnLabels(): (SortedArrayOfLabels *)NULL,
		( lp != NULL ) ? &lp->RevealRowLabels(): (SortedArrayOfLabels *)NULL,
		Int_T( ( lp != NULL ) ? lp->GetM() : 1 ),
		Int_T( ( lp != NULL ) ? lp->GetN() : 1 ),
		cntnts
	), LP( lp )
{
	assert( LP != NULL );
}


/*------------------------------------------------------------------------------

	void LP_Solution::SetLP( MPS_LP &lp )

PURPOSE:
	Passes a new LP to the object.

PARAMETERS:
	MPS_LP &lp
		The linear problem.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void LP_Solution::SetLP( MPS_LP &lp )
{
	LP = &lp;

	SolutionWithLabels::AddLabels( &lp.RevealColumnLabels(),
		&lp.RevealRowLabels() );

	SetContents( lp.GetM(), lp.GetN(), contents );
}


/*------------------------------------------------------------------------------

	void LP_Solution::WriteText( FILE *fp, int cntnts )
		const

PURPOSE:
	Write selected data from the solution vector to a text file.

PARAMETERS:
	FILE *fp
		Output file.

	int cntnts
		Data selection.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void LP_Solution::WriteText( FILE *fp, int cntnts )
	const
{
	assert( fp != NULL );
	assert( LP->GetM() == m );
	assert( LP->GetN() <= n );	// Slacks are not included in an MPS_LP object.

	if( contents == Empty ) return;
	if( cntnts == Empty ) cntnts = contents;

	fprintf( fp, "LP PROBLEM NAME:  %-8s\n\n", LP->RevealProblemName() );

	if( contents & ( RowAct | Slack | Dual ) )
	{
		assert( m > 0 );

		fprintf( fp,
			"EQUATIONS\n"
			"Row No  Label    Type   Lower limit    Upper limit  %s%s%s\n"
			"------  -------- ---- -------------  -------------  %s%s%s\n",

			( contents & RowAct )	? "Row activity" : "",
			( contents & Slack )	? "  Slack" : "",
			( contents & Dual )		? "  Dual var." : "",

			( contents & RowAct )	? "-------------" : "",
			( contents & Slack )	? "  -------------" : "",
			( contents & Dual )		? "  -------------" : ""
		);

		SortedArrayOfLabels &lab = LP->RevealRowLabels();

		for( Int_T i = 0; i < m; i++ )
		{
			Short_T	rt = Short_T( LP->GetRowType( i ) & RT_TYPE );

			fprintf( fp, "%6d  %-8s  %-2s  ",
				(int) i,
				lab.FindLabel( i ),
				( rt == RT_LE ) ? "LE" : ( rt == RT_GE ) ? "GE" :
					( rt == RT_EQ ) ? "EQ" : ( rt == RT_FR ) ? "FR" : "??"
				);

			Real_T	rl, ru;

			LP->GetConstraintRange( i, rl, ru );
			if( rl > -INFINITY )	fprintf( fp, "%13.6E  ", rl );
			else					fprintf( fp, "%13s  ", "-inf" );
			if( ru < +INFINITY )	fprintf( fp, "%13.6E  ", ru );
			else					fprintf( fp, "%13s  ", "+inf" );

			if( contents & RowAct )	fprintf( fp, "%13.6E", ra[i] );
			if( contents & Slack )	fprintf( fp, "  %13.6E", s[i] );
			if( contents & Dual )	fprintf( fp, "  %13.6E", y[i] );

			fprintf( fp, "\n" );
		}

		fprintf(  fp, "\n" );
	}

	if( contents & ( Primal | RC ) )
	{
		assert( n > 0 );

		fprintf( fp,
			"COLUMNS SECTION\n"
			"Var No  Label    Type   Lower bound    Upper bound  %s%s\n"
			"------  -------- ---- -------------  -------------  %s%s\n",

			( contents & Primal )	? " Primal value  " : "",
			( contents & RC )		? " Reduced cost" : "",

			( contents & Primal )	? "-------------  " : "",
			( contents & RC )		? "-------------" : ""
		);

		SortedArrayOfLabels &lab = LP->RevealColumnLabels();

		for( Int_T i = 0; i < n; i++ )
		{
			Short_T vt = Short_T( LP->GetVarType( i ) & VT_TYPE );

			fprintf( fp,
				"%6d  %-8s  %2s  ",

				(int) i,
				lab.FindLabel( i ),
				( vt == VT_FREE ) ? "FR" :
					( vt == VT_NORM ) ? "PL" :
					( vt == VT_MI ) ? "MI" :
					( vt == VT_FIXED ) ? "FX" :
					( vt == VT_BOUNDED ) ? "UP" : "??"
				);

			if( vt & VT_LO )		fprintf( fp, "%13.6E  ", LP->GetL( i ) );
			else					fprintf( fp, "%13s  ", "-inf" );
			if( vt & VT_UP )		fprintf( fp, "%13.6E  ", LP->GetU( i ) );
			else					fprintf( fp, "%13s  ", "+inf" );

			if( contents & Primal )	fprintf( fp, "%13.6E  ", x[i] );
			if( contents & RC )		fprintf( fp, "%13.6E", z[i] );

			fprintf( fp, "\n" );
		}

		fprintf(  fp, "\n" );
	}

	if( contents & Res )
		fprintf( fp, "%-25s%20.12E\n", "OBJECTIVE:", result );
}
