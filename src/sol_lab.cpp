/*------------------------------------------------------------------------------
MODULE TYPE:		Linear optimization supporting header
PROJECT CODE:		Simplex.
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	sol_lab.cpp
CREATED:			1995.03.14
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		sol_lab.h, stdtype.h, smartptr.h, sort_lab.h
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

#ifndef __SOL_LAB_H__
#	include "sol_lab.h"
#endif


/*------------------------------------------------------------------------------

	SolutionWithLabels::SolutionWithLabels( SortedArrayOfLabels *cols,
		SortedArrayOfLabels *rows, Int_T mm, Int_T nn, int cntnts )

PURPOSE:
	Constructor. Stores the data passed with arguments.

PARAMETERS:
	SortedArrayOfLabels *cols, SortedArrayOfLabels *rows
		Column and row labels.

	Int_T mm, Int_T nn
		Solution vectors' dimensions.
	
	int cntnts
		Solution contents (bit mask --- see "solution.h" for possible values).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

SolutionWithLabels::SolutionWithLabels( SortedArrayOfLabels *cols, // )
	SortedArrayOfLabels *rows, Int_T mm, Int_T nn, int cntnts )
	: ColumnLabels( cols ), RowLabels( rows ), OwnLabels( False ),
	Solution( mm, nn, cntnts )
{}


/*------------------------------------------------------------------------------

	(virtual)
	void SolutionWithLabels::SetN( Int_T nn )
	void SolutionWithLabels::SetM( Int_T mm )

PURPOSE:
	Change the dimension of groups of solution vectors. Calls the same function
of the underlying class "Solution", then performs an appropriate adjustment of
the row/column label data.

PARAMETERS:
	Int_T nn, Int_T mm
		The new dimensions >= 0.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SolutionWithLabels::SetN( Int_T nn )
{
	Solution::SetN( nn );
	if( nn == 0 && !OwnLabels ) ColumnLabels = NULL;
}


void SolutionWithLabels::SetM( Int_T mm )
{
	Solution::SetM( mm );
	if( mm == 0 && !OwnLabels ) RowLabels = NULL;
}


/*------------------------------------------------------------------------------

	void SolutionWithLabels::FreeSpace( void )

PURPOSE:
	Free the memory allocated by the object.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	Object becomes unusable (until new data is stored in it).

------------------------------------------------------------------------------*/

void SolutionWithLabels::FreeSpace( void )
{
	if( OwnLabels )
	{	
		if( ColumnLabels )	delete ColumnLabels;
		if( RowLabels )		delete RowLabels;
	}

	ColumnLabels = RowLabels = NULL;
	Solution::FreeSpace();
}


/*------------------------------------------------------------------------------

	void SolutionWithLabels::AddLabels( SortedArrayOfLabels *cols,
		SortedArrayOfLabels *rows )

PURPOSE:
	Adds labels for the primal and dual variables.

PARAMETERS:
	SortedArrayOfLabels *cols, SortedArrayOfLabels *rows
		The labels.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SolutionWithLabels::AddLabels( SortedArrayOfLabels *cols, // )
	SortedArrayOfLabels *rows )
{
	if( OwnLabels )
	{	
		if( ColumnLabels )	delete ColumnLabels;
		if( RowLabels )		delete RowLabels;
	}

	ColumnLabels = cols;
	RowLabels = rows;

	OwnLabels = False;
}


/*------------------------------------------------------------------------------

	void SolutionWithLabels::WriteText( FILE *fp, int cntnts )
		const

PURPOSE:
	Writes the solution to a text file.

PARAMETERS:
	FILE *fp
		Output stream.
	
	int cntnts
		What parts of the current solution data are to be written.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void SolutionWithLabels::WriteText( FILE *fp, int cntnts )
	const
{
	assert( fp != NULL );
	assert( ColumnLabels == NULL || ColumnLabels->NumberOfLabels() >= m );
	assert( RowLabels == NULL || RowLabels->NumberOfLabels() >= m );

	if( ColumnLabels == NULL && RowLabels == NULL )
	{
		Solution::WriteText( fp, cntnts );
		return;
	}

	if( contents == Empty ) return;
	if( cntnts == Empty ) cntnts = contents;

	fprintf( fp, "LP PROBLEM SOLUTION.\n\n" );

	if( contents & ( RowAct | Slack | Dual ) )
	{
		assert( m > 0 );

		fprintf( fp,
			"EQUATIONS\n"
			"Row No  Label    %s%s%s\n"
			"------  -------- %s%s%s\n",

			( contents & RowAct )	? "  Row activity" : "",
			( contents & Slack )	? "  Slack" : "",
			( contents & Dual )		? "  Dual var." : "",

			( contents & RowAct )	? "  -------------" : "",
			( contents & Slack )	? "  -------------" : "",
			( contents & Dual )		? "  -------------" : ""
		);

		SortedArrayOfLabels *lab = RowLabels;

		for( Int_T i = 0; i < m; i++ )
		{
			fprintf( fp, "%6d  %-8s  ", (int) i,
				( lab ) ? lab->FindLabel( i ) : "   ??   " );

			if( contents & RowAct )	fprintf( fp, "  %13.6E", ra[i] );
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
			"Var No  Label    %s%s\n"
			"------  -------- %s%s\n",

			( contents & Primal )	? " Primal value  " : "",
			( contents & RC )		? " Reduced cost" : "",

			( contents & Primal )	? "-------------  " : "",
			( contents & RC )		? "-------------" : ""
		);

		SortedArrayOfLabels *lab = ColumnLabels;

		for( Int_T i = 0; i < n; i++ )
		{
			fprintf( fp,
				"%6d  %-8s  ",
				(int) i, ( lab ) ? lab->FindLabel( i ) : "   ??   " );

			if( contents & Primal )	fprintf( fp, "%13.6E  ", x[i] );
			if( contents & RC )		fprintf( fp, "%13.6E", z[i] );

			fprintf( fp, "\n" );
		}

		fprintf(  fp, "\n" );
	}

	if( contents & Res )
		fprintf( fp, "%-25s%20.12E", "OBJECTIVE:", result );

	fprintf(  fp, "\n" );
}


/*------------------------------------------------------------------------------

	f()

PURPOSE:
	x

PARAMETERS:
	int x
	xxx

RETURN VALUE:
	x

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/
