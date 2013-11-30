/*------------------------------------------------------------------------------
MODULE TYPE:		File input routine.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	periods.h
CREATED:			1995.12.27
LAST MODIFIED:		1996.02.08

DEPENDENCIES:		mps_lp.h, stdtype.h,
					<stdio.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	Public interface to time file reader, which does both parsing and performs
the semantic actions (reads the numbers rows and columns beloging to the
consecutive stages).

------------------------------------------------------------------------------*/

#ifndef __PERIODS_H__
#define __PERIODS_H__

#include <stdio.h>

#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif


//==============================================================================
//
//	Class "Periods" declaration.
//

class Periods
{
private:
	Int_T maxLen, len;
	SortedArrayOfLabels PeriodLab;
	Array<Int_T> Row,
		Col;
	char ProblemName[LAB_LEN+1];
	Int_T TotalRows, TotalCols;

public:
	Periods( void );
	~Periods( void );

	Bool_T ReadPeriodFile( FILE *fp, MPS_LP &LP, VerbLevel Verbosity );

	Int_T NumberOfPeriods( void ) const;
	Int_T PeriodNumber( const char *lab ) const;
	const char *PeriodLabel( Int_T i ) const;
	Int_T FirstRow( Int_T i ) const;
	Int_T FirstColumn( Int_T i ) const;

	Int_T ColumnInPeriod( Int_T col ) const;
	Int_T RowInPeriod( Int_T row ) const;

	void RowRange( Int_T period, Int_T &start, Int_T &end ) const;
	void ColumnRange( Int_T period, Int_T &start, Int_T &end ) const;

private:
	Bool_T AddPeriod( MPS_LP &LP, const char *label, const char *row,
		const char *col );
	Bool_T Check( MPS_LP &LP );
	void PrintPeriods( MPS_LP &LP );
};

//
//	End of class "Periods" declaration.
//
//==============================================================================


//------------------------------------------------------------------------------
//
//	Inline functions' definitions.
//

/*------------------------------------------------------------------------------

	Periods::~Periods( void )

PURPOSE:
	Object destructor. Empty.

PARAMETERS:
	None.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Periods::~Periods( void )
	{}


/*------------------------------------------------------------------------------

	Int_T Periods::NumberOfPeriods( void ) const

PURPOSE:
	Data access routine.

PARAMETERS:
	None.

RETURN VALUE:
	Returns the number of defined time periods.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Int_T Periods::NumberOfPeriods( void )
const
	{ return len; }


/*------------------------------------------------------------------------------

	const char *Periods::PeriodLabel( Int_T i ) const

PURPOSE:
	Data access routine. Returns period label of a given period.

PARAMETERS:
	Int_T i
		Period number.

RETURN VALUE:
	Period label or "NULL" pointer if the object is empty.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
const char *Periods::PeriodLabel( Int_T i )
const
	{ return (len <= 0) ? (const char *)NULL : PeriodLab.FindLabel( i ); }


/*------------------------------------------------------------------------------

	Int_T Periods::FirstRow( Int_T i ) const
	Int_T Periods::FirstColumn( Int_T i ) const

PURPOSE:
	Calculate the number of the first row/column belonging to the given period.

PARAMETERS:
	Int_T i
		Period number.

RETURN VALUE:
	First row/column number or (-1), if the object is empty.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Int_T Periods::FirstRow( Int_T i )
const
	{ return Int_T( (len <= 0) ? -1 : Row[i] ); }


inline
Int_T Periods::FirstColumn( Int_T i )
const
	{ return Int_T( (len <= 0) ? -1 : Col[i] ); }


/*------------------------------------------------------------------------------

	Int_T Periods::PeriodNumber( const char *lab ) const

PURPOSE:
	Calculates the period number corresponding to a given label.

PARAMETERS:
	const char *lab
		Period label.

RETURN VALUE:
	Period number.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Int_T Periods::PeriodNumber( const char *lab )
const
{
	assert( lab != NULL );
	return (Int_T) PeriodLab.FindLabel( lab );
}


#endif
