/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem postsolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	postsolv.cpp
CREATED:			1994.05.08
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		std_tmpl.h, solution.h, compile.h, smartptr.h,
					error.h, memblock.h, smartdcl.h, sptr_deb.h, sptr_ndb.h,
					stdtype.h, myalloc.h, simplex.h, postsolv.h, my_defs.h
					<assert.h>, <stdio.h>

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

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SOLUTION_H__
#	include "solution.h"
#endif

#ifndef __POSTSOLV_H__
#	include "postsolv.h"
#endif


#if defined( explicit_templates )
	template class SmartPointerBase<PresolverAction *>;
	template class Array<PresolverAction *>;
	template class Ptr<PresolverAction *>;

#	ifdef NDEBUG
		template PresolverAction **MALLOC( PresolverAction **& Table,
			size_t len );
		template PresolverAction **REALLOC( PresolverAction **& Table,
			size_t len );
		template void FREE( PresolverAction **& Table );
#	endif
#endif


/*------------------------------------------------------------------------------

	int Postsolver::FreeMemory( void )

PURPOSE:
	Frees all memory allocated in the Postsolver object.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Postsolver::FreeMemory( void )
{
	for( Int_T i = 0; i < Len; i++ )
	{
		assert( Actions[i] != NULL );
		delete Actions[i];
	}

	Actions.Resize( 0 );
	Len = MaxLen = 0;
}


void Postsolver::SetSize( Int_T nn )
{
	assert( nn >= 0 );

	FreeMemory();
	n = nn;
	ExcludeCols.Resize( n+1 );
	ExcludeCols.Fill( False, n+1 );
	ExcludeCols[n] = True;
	Len = 0;
	MaxLen = 0;
}

/*------------------------------------------------------------------------------

	Bool_T Postsolver::Undo( Solution &solution )

PURPOSE:
	Restores the original values of all primal variables. The values of all
primal variables that were not removed by presolving are assumed to be stored in
the "Solution" structure. Their labels are passed as the second argument.

PARAMETERS:
	Solution &solution
		Solution structure. Contains a vector of primal variables corresponding
		to the lebels passed as the second argument.

RETURN VALUE:
	Success status. "False" may be returned on error only.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::Undo( Solution &solution )
{
	assert( n >= 0 );

	if( n < solution.GetN() )
		FatalError( "Postsolve information and primal solution mismatch" );

	//--------------------------------------------------------------------------
	//	Expand the solution vector to match the size of the variable array.
	//	Fill the empty places with zeros.
	//
	if( n > 0 )
	{
		Int_T pos = Int_T( solution.GetN() - 1 );

		solution.SetN( n );
		for( Int_T j = Int_T( n - 1 ); j >= 0; j-- )
			if( ExcludeCols[j] )
				solution.x[j] = 0.0;
			else
			{
				assert( pos >= 0 );
				solution.x[j] = solution.x[pos--];
			}

		if( pos != -1 )
			FatalError( "Postsolve information and primal solution mismatch" );
	}
	
	//--------------------------------------------------------------------------
	//	Now undo the presolve actions in a backward scan of the presolve action
	//	list. That should compute the values of all variables that were
	//	eliminated during preprocessing.
	//
	for( Int_T j = Int_T( Len - 1 ); j >= 0; j-- )
		if( Actions[j]->Undo( solution, FTOL, ExcludeCols ) == False )
			return False;

	return True;
}


/*------------------------------------------------------------------------------

	void Postsolver::ExtendActionArray( void )

PURPOSE:
	Internal function needed to lenghten the array of presolver actions. Uses
member data to determine the current length and the desired one. Extends the
array by adding at least 50% of the previous size.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Postsolver::ExtendActionArray( void )
{
	if( Len >= MaxLen )
	{
		Int_T NewMaxLen = Int_T( MaxLen + Max( Int_T( MaxLen/2 ),
			Int_T( POST_SOLV_MIN_LEN ) ) );

		Actions.Resize( NewMaxLen );
		Actions.Fill( (PresolverAction *)NULL, NewMaxLen, MaxLen );
		MaxLen = NewMaxLen;
	}
}


/*------------------------------------------------------------------------------

	void Postsolver::FixedAdjustment( Real_T val )
	void Postsolver::VariableFixing( Int_T ind, Real_T val )
	void Postsolver::ExplicitSlackRemoval( Int_T ind,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
		Short_T vartype, Real_T ls, Real_T us,
		Short_T rowtype, Real_T bl, Real_T bu )
	void Postsolver::FreeSingletonColumnRemoval( Int_T ind,
		const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
		Real_T af, Real_T b, Real_T slOpt )

PURPOSE:
	Postsolver action. Adds a new action object to the action array.

PARAMETERS:
	Passed to the action object constructor.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Postsolver::FixedAdjustment( Real_T val )
{
	assert( n > 0 );
	ExtendActionArray();
	Actions[Len] = new (class FixedAdjustment)( val );
	if( !Actions[Len] ) FatalError( "Not enough memory." );
	Len++;
}


void Postsolver::VariableFixing( Int_T ind, Real_T val )
{
	assert( n > 0 );
	assert( ind >= 0 && ind < n );
	ExcludeCols[ind] = True;
	ExtendActionArray();
	Actions[Len] = new (class VariableFixing)( ind, val );
	if( !Actions[Len] ) FatalError( "Not enough memory." );
	Len++;
}


void Postsolver::ExplicitSlackRemoval( Int_T ind, // )
	const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len, Real_T as,
	Short_T vartype, Real_T ls, Real_T us,
	Short_T rowtype, Real_T bl, Real_T bu )
{
	assert( n > 0 );
	assert( ind >= 0 && ind < n );
	ExcludeCols[ind] = True;
	ExtendActionArray();
	Actions[Len] = new (class ExplicitSlackRemoval)( ind, &ExcludeCols,
		a, col, len, as, vartype, ls, us, rowtype, bl, bu );
	if( !Actions[Len] ) FatalError( "Not enough memory." );
	Len++;
}


void Postsolver::FreeSingletonColumnRemoval( Int_T ind, // )
	const Ptr<Real_T> &a, const Ptr<Int_T> &col, Int_T len,
	Real_T af, Real_T b, Real_T slOpt )
{
	assert( n > 0 );
	assert( ind >= 0 && ind < n );
	ExcludeCols[ind] = True;
	ExtendActionArray();
	Actions[Len] = new (class FreeSingletonColumnRemoval)( ind, &ExcludeCols,
		a, col, len, af, b, slOpt );
	if( !Actions[Len] ) FatalError( "Not enough memory." );
	Len++;
}


/*------------------------------------------------------------------------------

	Bool_T Postsolver::WriteActions( const char *fname )
	Bool_T Postsolver::WriteActions( FILE *fp )

PURPOSE:
	Writes the presolver actions to the stream passed by the only argument.

PARAMETERS:
	FILE *fp
		A stream pointer (ssumed to be opened for writing).

RETURN VALUE:
	Success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T Postsolver::WriteActions( const char *fname )
{
	assert( fname != NULL && *fname != '\0' );

	FILE *fp = fopen( fname, "wt" );
	if( fp == NULL ) return False;

	Bool_T result = WriteActions( fp );
	fclose( fp );
	return result;
}


Bool_T Postsolver::WriteActions( FILE *fp )
{
	assert( fp != NULL );
	assert( n > 0 );

	//--------------------------------------------------------------------------
	//	Write the file header line.
	//
	if( fprintf( fp, "PRESOLVER REPORT\n" ) == EOF ) return False;

	Int_T i;

	//--------------------------------------------------------------------------
	//	Write the indice of the removed rows.
	//
//	if( fprintf( fp, "ROWS    %d\n", m ) == EOF ) return False;
//	for( i = 0; i < m; i++ )
//		if( ExcludeRow[i] && fprintf( fp, "    %-8d\n", (int)i ) == EOF )
//			return False;

	//--------------------------------------------------------------------------
	//	Write the labels of the removed columns.
	//
	if( fprintf( fp, "COLUMNS %d\n", n ) == EOF ) return False;
	for( i = 0; i < n; i++ )
		if( ExcludeCols[i] && fprintf( fp, "    %-d\n", (int)i ) == EOF )
			return False;

	//--------------------------------------------------------------------------
	//	Write the actions.
	//
	for( i = 0; i < Len; i++ )
	{
		assert( Actions[i] != NULL );
		if( (Actions[i])->Write( fp ) == False )
			return False;
	}

	//--------------------------------------------------------------------------
	//	Write the end label.
	//
	if( fprintf( fp, "ENDATA\n" ) == EOF )
		return False;

	return True;
}
