/*------------------------------------------------------------------------------
MODULE TYPE:		General purpose templates
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	stochsol.cpp
CREATED:			1996.06.21
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		solution.h, stochsol.h
					<stdio.h>

------------------------------------------------------------------------------*/

#ifndef __STOCHSOL_H__
#	include "stochsol.h"
#endif


void StochSolution::WriteText( FILE *fp, int cntnts )
	const
{
	assert( fp != NULL );

	if( cntnts & Res )
		fprintf( fp, "OPTIMAL VALUE: %20.10E\n", result );

	Int_T i;

	if( cntnts & Primal )
	{
		fprintf( fp, "\nSOLUTION PRINTOUT:\n" );

		for( i = 0; i < GetN(); i++ )
			fprintf( fp, "\tx[%6d] = %20.10E\n", (int) i, (double) x[i] );
	}

	if( cntnts & Weights )
	{
		fprintf( fp, "\nWEIGHTS AND RECOURSE COSTS FOR REALIZATIONS:\n" );
	
		fprintf( fp, "\t %10s    %10s    %10s\n",
			"NUMBER", "PROBABIL.", "COST" );

		for( i = 0; i < Scenarios; i++ )
			fprintf( fp, "\t#%10d    %10.4E    %10.4E\n", (int) i,
				weights[i], f[i] );
	}
}
