/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		----------------
PROJECT FULL NAME:	----------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	----------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	writ_mps.cpp
CREATED:			1993.09.27
LAST MODIFIED:		1996.04.16

DEPENDENCIES:		compile.h, mps_lp.h, lp_codes.h, print.h,
					simplex.h, solv_lp.h, std_math.h, sort_lab.h
					<stdio.h>, <math.h>, <assert.h>

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

------------------------------------------------------------------------------*/


#include <math.h>
#include <stdio.h>
#include <assert.h>

#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __COMPILE_H__
#	include "compile.h"
#endif
#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif


/*------------------------------------------------------------------------------

	Bool_T MPS_LP::WriteMPS( FILE *MPS_File, VerbLevel Verbosity ) const

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

Bool_T MPS_LP::WriteMPS( const char *MPS_File, VerbLevel Verbosity )
	const
{
	FILE *fp = fopen( MPS_File, "wt" );
	if( fp == NULL )
		return False;

	//
	//	NOTE: this is a virtual call - it may call
	//	SolvableLP::WriteMPS( FILE *, VerbLevel )
	//
	Bool_T result = WriteMPS( fp, Verbosity );

	fclose( fp );
	return result;
}


Bool_T MPS_LP::WriteMPS( FILE *MPS_File, VerbLevel Verbosity )
	const
{
	Int_T i;
	Int_T cs, ce;
	char c = '\0';
	Bool_T RHS_Present, RangesPresent, BoundsPresent;

	assert( MPS_File != NULL && OS == OS_FULL );

	//--------------------------------------------------------------------------
	//	NAME line.
	//
	if( Verbosity >= V_HIGH )
		Print( "\nWriting output file.\n" );
	if( fprintf( MPS_File, "%-14s%-8s\n", "NAME", Name ) == EOF )
		goto error;

	//--------------------------------------------------------------------------
	//	ROWS section.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ROWS\n" );
	if( fprintf( MPS_File, "%s\n", "ROWS" ) == EOF )
		goto error;
	for( i = 0; i < m; i++ )
	{
		switch( RowType[i] & RT_TYPE )
		{
		case RT_FR:	c = 'N'; break;
		case RT_EQ:	c = 'E'; break;
		case RT_GE:	c = 'G'; break;
		case RT_LE:	c = 'L'; break;
		default:	abort();
		}
		if( fprintf( MPS_File, " %c  %-8s\n", c,
			RowLabels.FindLabel( i ) ) == EOF )
			goto error;
	}

	//--------------------------------------------------------------------------
	//	COLUMNS section.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting COLUMNS\n" );
	if( fprintf( MPS_File, "%s\n", "COLUMNS" ) == EOF )
		goto error;
	for( i = 0; i < n; i++ )
	{
		cs = ColStart[i];
		ce = ColStart[i+1];

		for( register Int_T j = cs; j < ce; j++ )
		{
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G",
				ColLabels.FindLabel( i ),
				RowLabels.FindLabel( Row[ (size_t) j ] ),
				ac[ (size_t) j ] ) == EOF )
				goto error;
			if( j != ce - 1 )
			{
				j++;
				if( fprintf( MPS_File, "   %-8s  %12G",
					RowLabels.FindLabel( Row[ (size_t) j ] ),
					ac[ (size_t) j ] )
					== EOF )
					goto error;
			}
			if( fprintf( MPS_File, "\n" ) == EOF )
				goto error;
		}
	}

	//--------------------------------------------------------------------------
	//	RHS section.
	//
	RHS_Present = False;
	if( Verbosity >= V_HIGH )
		Print( "\tWriting RHS\n" );
	if( fprintf( MPS_File, "%s\n", "RHS" ) == EOF )
		goto error;
	for( i = 0; i < m; i++ )
		if( IsNonZero( b[i] ) && !( RowType[i] & RT_FR ) )
		{
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
				RHS_Name, RowLabels.FindLabel( i ), b[i] ) == EOF )
				goto error;
			RHS_Present = True;
		}
	if( !RHS_Present )
		if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
			RHS_Name, RowLabels.FindLabel( 0, LO_UNSORTED ), 0.0 ) == EOF )
			goto error;

	//--------------------------------------------------------------------------
	//	RANGES section.
	//
	RangesPresent = False;
	for( i = 0; i < m; i++ )
		if( RowType[ i ] & RT_RNG )
		{
			if( !RangesPresent )
			{
				RangesPresent = True;
				if( Verbosity >= V_HIGH )
					Print( "\tWriting RANGES\n" );
				if( fprintf( MPS_File, "%s\n", "RANGES" ) == EOF )
					goto error;
			}
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
				RangesName, RowLabels.FindLabel( i ), r[i] ) == EOF )
				goto error;
		}

	//--------------------------------------------------------------------------
	//	BOUNDS section.
	//	Entries created only when necessary, i.e. for fixed, bounded and free
	//	variables or for variables with non-zero lower bound.
	//
	BoundsPresent = False;
	for( i = 0; i < n; i++ )
		if( VarType[i] & ( VT_FX | VT_UP ) ||
			( VarType[i] & VT_LO && IsNonZero( l[i] ) ) ||
			!( VarType[i] & ( VT_LO | VT_UP ) ) )
		{
			//------------------------------------------------------------------
			//	If it hasn't been done before - write BONDS section header.
			//
			if( !BoundsPresent )
			{
				BoundsPresent = True;
				if( Verbosity >= V_HIGH )
					Print( "\tWriting BOUNDS\n" );
				if( fprintf( MPS_File, "%s\n", "BOUNDS" ) == EOF )
					goto error;
			}

			//------------------------------------------------------------------
			//	Detect FREE variables.
			//
			if( !( VarType[i] & ( VT_LO | VT_UP ) ) )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s\n",
					"FR", BoundsName, ColLabels.FindLabel( i ) ) == EOF )
					goto error;
			}
			//------------------------------------------------------------------
			//	Detect FIXED variables.
			//
			else if( VarType[i] == VT_FIXED )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
					"FX", BoundsName, ColLabels.FindLabel( i ), l[i] ) == EOF )
					goto error;
			}
			//------------------------------------------------------------------
			//	Detect non-positive (MI) variables.
			//
			else if( ( VarType[i] & VT_UP ) && !( VarType[i] & VT_LO ) )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s\n",
					"MI", BoundsName, ColLabels.FindLabel( i ) ) == EOF )
					goto error;
				if( IsNonZero( u[i] ) )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"UP", BoundsName, ColLabels.FindLabel( i ), u[i] )
						== EOF )
						goto error;
			}
			//------------------------------------------------------------------
			//	For all other variables print out non-zero lower bounds and
			//	finite upper bounds.
			//
			else
			{
				if( VarType[i] & VT_LO && IsNonZero( l[i] )  )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"LO", BoundsName, ColLabels.FindLabel( i ), l[i] )
						== EOF )
						goto error;

				if( VarType[i] & VT_UP )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"UP", BoundsName, ColLabels.FindLabel( i ), u[i] )
						== EOF )
						goto error;
			}
		}

	//--------------------------------------------------------------------------
	//	ENDATA.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ENDATA\n" );
	if( fprintf( MPS_File, "%s\n", "ENDATA" ) == EOF )
		goto error;

	if( Verbosity >= V_HIGH )
		Print( "Finished writing.\n\n" );
	return True;

error:
	return False;
}


/*------------------------------------------------------------------------------

	Bool_T SolvableLP::WriteMPS( FILE *MPS_File, VerbLevel Verbosity ) const

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

Bool_T SolvableLP::WriteMPS( FILE *MPS_File, VerbLevel Verbosity )
	const
{
	Int_T i, j;
	Bool_T RHS_Present, RangesPresent, BoundsPresent;

	assert( MPS_File != NULL && OS == OS_FULL );

	//--------------------------------------------------------------------------
	//	NAME line.
	//
	if( Verbosity >= V_HIGH )
		Print( "\nWriting output file.\n" );
	if( fprintf( MPS_File, "%-14s%-8s\n", "NAME", Name ) == EOF )
		goto error;

	//--------------------------------------------------------------------------
	//	ROWS section (finished by the objective vector name).
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ROWS\n" );

	if( fprintf( MPS_File, "%s\n", "ROWS" ) == EOF )
		goto error;

	for( i = 0; i < m; i++ )
	{
		char C = '\0';

		switch( RowType[i] & RT_TYPE )
		{
		case RT_EQ:	C = 'E'; break;
		case RT_GE:	C = 'G'; break;
		case RT_LE:	C = 'L'; break;
		default:	abort();
		}
		if( fprintf( MPS_File, " %c  %-8s\n", C, RowLabels.FindLabel( i ) )
			== EOF )
			goto error;
	}

	if( fprintf( MPS_File, " N  %-8s\n", Obj ) == EOF )
		goto error;

	//--------------------------------------------------------------------------
	//	COLUMNS section (each column followed by its cost coefficient).
	//	Assumption: no empty columns.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting COLUMNS\n" );

	if( fprintf( MPS_File, "%s\n", "COLUMNS" ) == EOF )
		goto error;

	for( j = 0; j < n; j++ )
	{
		Ptr<Real_T> a;
		Ptr<Int_T> row;
		Int_T len;
		const char *clab = ColLabels.FindLabel( j );

		for( MPS_LP::GetColumn( j, a, row, len ); len; --len, ++a, ++row )
		{
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G",
				clab, RowLabels.FindLabel( *row ), *a ) == EOF )
				goto error;

			if( len > 1 )
			{
				--len; ++a; ++row;
				if( fprintf( MPS_File, "   %-8s  %12G",
					RowLabels.FindLabel( *row ), *a ) == EOF )
					goto error;
			}

			if( fprintf( MPS_File, "\n" ) == EOF )
				goto error;
		}

		if( IsNonZero( c[j] ) )
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
				clab, Obj, c[j] ) == EOF )
				goto error;
	}

	//--------------------------------------------------------------------------
	//	RHS section.
	//
	RHS_Present = False;

	if( Verbosity >= V_HIGH )
		Print( "\tWriting RHS\n" );

	if( fprintf( MPS_File, "%s\n", "RHS" ) == EOF )
		goto error;

	for( i = 0; i < m; i++ )
		if( IsNonZero( b[i] ) )
		{
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
				RHS_Name, RowLabels.FindLabel( i ), b[i] ) == EOF )
				goto error;
			RHS_Present = True;
		}
	if( !RHS_Present )
		if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
			RHS_Name, RowLabels.FindLabel( 0, LO_UNSORTED ), 0.0 ) == EOF )
			goto error;

	//--------------------------------------------------------------------------
	//	RANGES section.
	//
	RangesPresent = False;

	for( i = 0; i < m; i++ )
		if( RowType[i] & RT_RNG )
		{
			if( !RangesPresent )
			{
				RangesPresent = True;
				if( Verbosity >= V_HIGH )
					Print( "\tWriting RANGES\n" );
				if( fprintf( MPS_File, "%s\n", "RANGES" ) == EOF )
					goto error;
			}
			if( fprintf( MPS_File, "    %-8s  %-8s  %12G\n",
				RangesName, RowLabels.FindLabel( i ), r[i] ) == EOF )
				goto error;
		}

	//--------------------------------------------------------------------------
	//	BOUNDS section.
	//	Entries created only when necessary, i.e. for fixed, bounded and free
	//	variables or for variables with non-zero lower bound.
	//
	BoundsPresent = False;
	for( j = 0; j < n; j++ )
	{
		Short_T vt = VarType[j];
		const char *lab = ColLabels.FindLabel( j );

		if( vt & ( VT_FX | VT_UP ) ||
			( vt & VT_LO && IsNonZero( l[j] ) ) ||
			!( vt & ( VT_LO | VT_UP ) ) )
		{
			//------------------------------------------------------------------
			//	If it hasn't been done before - write BOUNDS section header.
			//
			if( !BoundsPresent )
			{
				BoundsPresent = True;

				if( Verbosity >= V_HIGH )
					Print( "\tWriting BOUNDS\n" );

				if( fprintf( MPS_File, "%s\n", "BOUNDS" ) == EOF )
					goto error;
			}

			//------------------------------------------------------------------
			//	Detect FREE variables.
			//
			if( !( vt & ( VT_LO | VT_UP ) ) )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s\n",
					"FR", BoundsName, lab ) == EOF )
					goto error;
			}
			//------------------------------------------------------------------
			//	Detect FIXED variables.
			//
			else if( vt == VT_FIXED )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
					"FX", BoundsName, lab, l[j] ) == EOF )
					goto error;
			}
			//------------------------------------------------------------------
			//	Detect non-positive (MI) variables.
			//
			else if( ( vt & VT_UP ) && !( vt & VT_LO ) )
			{
				if( fprintf( MPS_File, " %2s %-8s  %-8s\n",
					"MI", BoundsName, lab ) == EOF )
					goto error;
				if( IsNonZero( u[j] ) )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"UP", BoundsName, lab, u[j] ) == EOF )
						goto error;
			}
			//------------------------------------------------------------------
			//	For all other variables print out non-zero lower bounds and
			//	finite upper bounds.
			//
			else
			{
				if( vt & VT_LO && IsNonZero( l[j] ) )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"LO", BoundsName, lab, l[j] ) == EOF )
						goto error;

				if( vt & VT_UP )
					if( fprintf( MPS_File, " %2s %-8s  %-8s  %12G\n",
						"UP", BoundsName, lab, u[j] ) == EOF )
						goto error;
			}
		}
	}

	//--------------------------------------------------------------------------
	//	ENDATA.
	//
	if( Verbosity >= V_HIGH )
		Print( "\tWriting ENDATA\n" );

	if( fprintf( MPS_File, "%s\n", "ENDATA" ) == EOF )
		goto error;

	if( Verbosity >= V_HIGH )
		Print( "Finished writing.\n\n" );
	return True;

error:
	return False;
}
