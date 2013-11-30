/*------------------------------------------------------------------------------
MODULE TYPE:		Linear optimization supporting header
PROJECT CODE:		Simplex.
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

HEADER FILE NAME:	sol_lab.h
CREATED:			1995.03.14
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		solution.h, stdtype.h, smartptr.h, sort_lab.h
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	The header contains a declaration of a structure used for holding a linear
program solution with variable and constraint labels.

------------------------------------------------------------------------------*/

#ifndef __SOL_LAB_H__
#define __SOL_LAB_H__

#ifndef __SOLUTION_H__
#	include "solution.h"
#endif
#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif
#ifdef SUPPORT_LP_DIT
#	include "lp_dit.h"
#endif


//==============================================================================
//
//	Structure "SolutionWithLabels" declaration.
//

class SolutionWithLabels : public Solution
{
public:
	SortedArrayOfLabels *ColumnLabels,
		*RowLabels;

protected:
	Bool_T OwnLabels;

public:
	SolutionWithLabels( SortedArrayOfLabels *cols = NULL,
		SortedArrayOfLabels *rows = NULL,
		Int_T mm = 0, Int_T nn = 0, int cntnts = Empty );
	virtual ~SolutionWithLabels( void );

	virtual void FreeSpace( void );

	virtual void SetN( Int_T nn );
	virtual void SetM( Int_T mm );

	void AddLabels( SortedArrayOfLabels *cols = NULL,
		SortedArrayOfLabels *rows = NULL );

	int AddRowLabel( const char *label );
	int AddColumnLabel( const char *label );

	virtual void WriteText( FILE *fp, int cntnts ) const;
};

//
//	End of structure "Solution" declaration.
//
//==============================================================================


inline
SolutionWithLabels::~SolutionWithLabels( void )
{
	FreeSpace();
}


inline
int SolutionWithLabels::AddRowLabel( const char *label )
{
	return ( RowLabels ) ? RowLabels->AddLabel( label ) : -1;
}


inline
int SolutionWithLabels::AddColumnLabel( const char *label )
{
	return ( ColumnLabels ) ? ColumnLabels->AddLabel( label ) : -1;
}


#endif
