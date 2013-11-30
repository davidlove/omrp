/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	presolve.cpp
CREATED:			1993.10.07
LAST MODIFIED:		1996.02.15

DEPENDENCIES:		myalloc.h, error.h, std_tmpl.h, std_math.h, print.h,
					presolve.h, simplex.h, smartptr.h, memblock.h, smartdcl.h,
					sptr_deb.h, sptr_ndb.h, stdtype.h, lp_codes.h, solv_lp.h,
					compile.h, mps_lp.h, sort_lab.h, cl_list.h, postsolv.h,
					my_defs.h
					<stdio.h>, <assert.h>

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
#include <assert.h>

#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif
#ifndef __POSTSOLV_H__
#	include "postsolv.h"
#endif


/*------------------------------------------------------------------------------

	Presolver::Presolver( SolvableLP &lp, Postsolver *postsolve )

PURPOSE:
	Presolver object constructor. Initializes the object so that:
a)	reference to the linear problem and problem dimensions are stored,
b)	column and row list are made empty,
c)	statistic gathering counters are all set to zero,
d)	a postsolver object reference is stored.

PARAMETERS:
	SolvableLP &lp, Postsolver *postsolve
		Reference of the LP, pointer to postsolver.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Presolver::Presolver( SolvableLP &lp, Postsolver *postsolve )
	:
	//--------------------------------------------------------------------------
	//	Initialize linear problem data.
	//
	LP( lp ), m( LP.GetM() ), n( LP.GetN() ),
	Status( LPS_UNKNOWN ), PrimalInf( 0.0 ), DualInf( 0.0 ),
	ftol( DEFAULT_FEASIBILITY_TOL ),

	//--------------------------------------------------------------------------
	//	Initialize presolver data.
	//
	ExcludeRows( 0 ), ExcludeCols( 0 ), ExplSlackRemoved( 0 ),
	bl( 0 ), bu( 0 ), rt( 0 ), f( 0.0 ), ColList(), RowList(),

	//--------------------------------------------------------------------------
	//	Initialize the statistic counters.
	//
	EliminatedNonZeros( 0 ), EliminatedRows( 0 ), EliminatedCols( 0 ),
	OrigFixedVars( 0 ), OrigFreeSingletonCols( 0 ),
	ImpliedFreeSingletonCols( 0 ),
	RelaxedConstraints( 0 ),
	ForcingRows( 0 ), DominatedRows( 0 ), VariableBoundsTightened( 0 ),

	//--------------------------------------------------------------------------
	//	Link with a post-solver object.
	//
	PostSolve( postsolve )
{
	if( postsolve && lp.GetN() > 0 )
		postsolve->SetSize( lp.GetN() );
}


/*------------------------------------------------------------------------------

	void Presolver::ReleaseWorkMemory( void )

PURPOSE:
	Releases all work memory allocated by the presolver.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::ReleaseWorkMemory( void )
{
	RowLen.Resize( 0 );					ColLen.Resize( 0 );

	ExcludeRows.Resize( 0 );			ExcludeCols.Resize( 0 );
	ExplSlackRemoved.Resize( 0 );

	bl.Resize( 0 );						bu.Resize( 0 );
	rt.Resize( 0 );
}


/*------------------------------------------------------------------------------

	void Presolver::InitializePresolverData( void )

PURPOSE:
	Meant to be invoked before the actual presolve procedure is run. It prepares
all the data structures for presolving the LP:
a)	reads the LP dimensions,
b)	allocates arrays of row and column lenghts,
c)	creates lists of rows and columns of equal lengths,
d)	creates and fills arrays of eliminated rows and columns,
e)	creates and fills with zeros arrays of row activity limits.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::InitializePresolverData( void )
{
	f = 0.0;
	m = LP.GetM();
	n = LP.GetN();

	//--------------------------------------------------------------------------
	//	Allocate arrays of row and column lenghts.
	//
	RowLen.Resize( m );
	ColLen.Resize( n );

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T j, len;

	RowLen.Fill( 0, m );

	for( j = 0; j < n; j++ )
	{
		LP.GetColumn( j, a, row, len );
		ColLen[j] = len;

		for( ; len; --len, ++row )
			RowLen[ *row ]++;
	}

#ifndef NDEBUG
	{
		Int_T i, nzc = 0, nzr = 0;

		for( i = 0; i < m; i++ ) nzr += RowLen[i];
		for( i = 0; i < n; i++ ) nzc += ColLen[i];

		assert( nzc == nzr );
	}
#endif

	//--------------------------------------------------------------------------
	//	Allocate and initialize exclusion masks.
	//
	ExcludeRows.Resize( m + 1 );
	ExcludeCols.Resize( n + 1 );
	ExplSlackRemoved.Resize( m );
	ExcludeRows.Fill( False, m + 1 );
	ExcludeCols.Fill( False, n + 1 );
	ExplSlackRemoved.Fill( False, m );

	ExcludeCols[n] = ExcludeRows[m] = True;

	//--------------------------------------------------------------------------
	//	Initialize and fill lists of rows and columns of the same lengths.
	//
	RowList.Initialize( m, n, RowLen );
	ColList.Initialize( n, m, ColLen );

	//--------------------------------------------------------------------------
	//	Allocate and fill tables of row activity limits.
	//
	bl.Resize( m );
	bu.Resize( m );
	rt.Resize( m );

	RHS2BLU();
}


/*------------------------------------------------------------------------------

	void Presolver::Presolve( const int Mode, VerbLevel Verbosity )

PURPOSE:
	Presolves the linear problem. Calls in a loop theose presolve methods, that
were indicated in the "Mode" argument. Starts with eliminating all original
fixed variables, finishes by removing all empty rows and columns.

PARAMETERS:
	const int Mode
		Bit mask marking which presolve metgods are to be used. For the meanings
		of those bits see "presolve.h", enumeration named "LP_REDUCTIONS".
	
	VerbLevel Verbosity
		Verbosity level (used when generating a report).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Presolver::Presolve( const int Mode, VerbLevel Verbosity )
{
	if( !Mode ) return;

	if( Verbosity >= V_HIGH )
	{
		Print( "\nPRESOLVER INVOKED:\n" );
		if( Mode & LPR_ALL )
		{
			Print( "\tThe following presolve techniques activated:\n" );
			if( Mode & LPR_EMPTY_ROWS )
				Print( "\t\tempty column removal\n" );
			if( Mode & LPR_EMPTY_COLS )
				Print( "\t\tempty row removal\n" );
			if( Mode & LPR_ORIG_FIXED )
				Print( "\t\tfixed variable removal\n" );
			if( Mode & LPR_NUM_ELIM )
				Print( "\t\tnumerical eliminations\n" );
			if( Mode & LPR_SINGL_ROWS )
				Print( "\t\tsingleton row reduction\n" );
			if( Mode & LPR_SINGL_COLS )
				Print( "\t\tfree singleton column reduction\n" );
			if( Mode & LPR_FORC_DOM_CONSTR )
				Print( "\t\tforcing and dominated row detection\n" );
			if( Mode & LPR_DOM_COLS )
				Print( "\t\tdominated and weakly dominated column removal\n" );
			if( Mode & LPR_EXPLICIT_SLACKS )
				Print( "\t\texplicit slack removal\n" );
		}
		else
		{
			Print( "\tNo presolve techniques activated. Exiting." );
			return;
		}
	}

	InitializePresolverData();

	Bool_T Elim;
	Int_T ElNum = 0;

	if( Mode & LPR_ORIG_FIXED )
		if( ( OrigFixedVars = EliminateFixedVariables() ) != 0 )
			if( Verbosity >= V_HIGH )
				Print( "\tFound %d fixed variables in original problem.\n",
					(int) OrigFixedVars );

	do
	{
		Elim = False;

		if( Mode & LPR_SINGL_ROWS )
		{
			if( ( ElNum = EliminateSingletonRows() ) != 0 )
			{
				if( Verbosity >= V_HIGH )
					Print( "\tFound %d singleton rows.\n", (int) ElNum );
				Elim = True;
			}
			if( Status != LPS_UNKNOWN ) break;
		}

		if( Mode & ( LPR_SINGL_COLS | LPR_EXPLICIT_SLACKS ) )
		{
			if( ( ElNum = DealWithSigletonColumns( Mode ) ) != 0 )
			{
				if( Verbosity >= V_HIGH )
					Print( "\tConverted %d singleton columns.\n",
						(int) ElNum );
				Elim = True;
			}
			if( Status != LPS_UNKNOWN ) break;
		}

		if( Mode & LPR_FORC_DOM_CONSTR )
		{
			if( ( ElNum = ForcingAndDominatedRows() ) != 0 )
			{
				if( Verbosity >= V_HIGH )
					Print( "\tFound %d rows forcing or dominated.\n",
						(int) ElNum );
				Elim = True;
			}
			if( Status != LPS_UNKNOWN ) break;
		}

		if( !Elim && ( Mode & LPR_NUM_ELIM ) )
		{
			if( ( ElNum = NumericalEliminations( 100.0 ) ) != 0 )
			{
				if( Verbosity >= V_HIGH )
					Print( "\tEliminated %d non-zeros.\n", (int) ElNum );
				Elim = True;
			}
			if( Status != LPS_UNKNOWN ) break;
		}

	} while( Elim );

	if( Mode & LPR_EMPTY_ROWS && ( ElNum = EliminateEmptyRows() ) != 0 &&
		Verbosity >= V_HIGH )
		Print( "\tFound %d empty rows.\n", (int) ElNum );

	if( Mode & LPR_EMPTY_COLS && ( ElNum = EliminateEmptyColumns() ) != 0 &&
		Verbosity >= V_HIGH )
		Print( "\tFound %d empty columns.\n", (int) ElNum );

	if( PostSolve && IsNonZero( f ) )
		PostSolve->FixedAdjustment( f );

	//--------------------------------------------------------------------------
	{
		Int_T nzc = 0;
		for( Int_T j = 0; j < n; j++ )
			if( !ExcludeCols[j] )
				nzc += ColLen[j];

#ifndef NDEBUG
		Int_T nzr = 0;
		for( Int_T i = 0; i < m; i++ )
			if( !ExcludeRows[i] )
				nzr += RowLen[i];
		assert( nzc == nzr );
#endif

		EliminatedNonZeros = Int_T( LP.GetNZ() - nzc );
	}

	//--------------------------------------------------------------------------
	//	Restore the right-hand-side and range vector from "bl" and "bu" vector
	//	pair.
	//
	BLU2RHS();
	UpdateLP_AfterReductions();

	//--------------------------------------------------------------------------
	//	Set problem status to solved if the dimensions have been reduced to
	//	zero.
	//
	EliminatedRows += ForcingRows;
	if( EliminatedRows == m && EliminatedCols == n &&
		EliminatedNonZeros == LP.GetNZ() )
		Status = LPS_SOLVED;

	//--------------------------------------------------------------------------
	//	Report on the presolver activity. Special (short) note if no reductions
	//	were obtained.
	//
	if( Verbosity >= V_LOW )
	{
		if( EliminatedRows || EliminatedCols || EliminatedNonZeros )
			Print(
				//----------------------------------------------------------
				//	Format string.
				//
				"\nElimination statistics:\n"
				"\tTotals:\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%20.12E\n"

				"\tSingleton column analysis:\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"

				"\tRow analysis:\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"
				"\t\t%-25s%10d\n"

				"\tInfeasibility\n"
				"\t\t%-25s%10G\n"
				"\t\t%-25s%10G\n",

				//----------------------------------------------------------
				//	Data.
				//
				"Original fixed variables",	(int) OrigFixedVars,
				"Eliminated rows:",			(int) EliminatedRows,
				"Eliminated columns:",		(int) EliminatedCols,
				"Eliminated non-zeros:",	(int) EliminatedNonZeros,
				"Fixed adjustment:",		(double) f,

				"Original free singletons:",(int) OrigFreeSingletonCols,
				"Implied free variables:",	(int) ImpliedFreeSingletonCols,
				"Relaxed constraints:",		(int) RelaxedConstraints,

				"Forcing rows:",			(int) ForcingRows,
				"Dominated rows:",			(int) DominatedRows,
				"Tightened bounds:",		(int) VariableBoundsTightened,
				"Primal:",					(double) PrimalInf,
				"Dual:",					(double) DualInf
			);
		else
			Print( "The problem was not reduced.\n" );
	}
	else if( Verbosity >= V_LINE )
	{
		Print(
			//----------------------------------------------------------
			//	Format string.
			//
			" %10d | %10d | %10d | %10d | %10d | %10d | %10d | %10d |",

			//----------------------------------------------------------
			//	Data.
			//
			(int) EliminatedRows,
			(int) EliminatedCols,
			(int) EliminatedNonZeros,

			(int) OrigFreeSingletonCols,
			(int) ImpliedFreeSingletonCols,
			(int) RelaxedConstraints,

			(int) ForcingRows,
			(int) DominatedRows
		);
	}

	//--------------------------------------------------------------------------
	//	If the problem's status changed, display information.
	//
	switch( Status )
	{
	case LPS_UNKNOWN:
		break;

	case LPS_INFEASIBLE:
		if( Verbosity >= V_LOW )
			Print( "Problem is infeasible.\n" );
		break;

	case LPS_UNBOUNDED:
		if( Verbosity >= V_LOW )
			Print( "Problem is unbounded.\n" );
		break;

	case LPS_SOLVED:
		if( Verbosity >= V_LOW )
			Print( "Problem is solved.\n" );
		break;
	}
}


/*------------------------------------------------------------------------------

	void Presolver::UpdateAfterReductions( void )

PURPOSE:
	Convert the matrix - remove the redundant rows and columns (so far they have
only been marked as redundant). Use member data 'ExcludeRows' and 'ExcludeCols'.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	All the non-zeros corresponding to the removed rows, columns and numerically
eliminated nonn-zeros are now actually removed. Row structure is created (or
updated) from scratch. Row and column labels are also updated.

------------------------------------------------------------------------------*/

void Presolver::UpdateLP_AfterReductions( void )
{
	LP.UpdateAfterReduction( &ExcludeRows, &ExcludeCols );
	LP.UpdateRowStructure();
}
