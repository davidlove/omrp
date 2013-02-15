/*------------------------------------------------------------------------------
MODULE TYPE:		Linear optimization supporting header
PROJECT CODE:		Simplex.
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solution.cpp
CREATED:			1995.03.03
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		solution.h, stdtype.h, smartptr.h
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

------------------------------------------------------------------------------*/

#include <assert.h>

#ifndef __SOLUTION_H__
#	include "solution.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif


/*------------------------------------------------------------------------------

	void Solution::SetContents( Int_T mm, Int_T nn, int cntnts )

PURPOSE:
	Resets the solution object by changing the dimensions and contents of the
structure as well as resetting all ADDED numerical data to zeros.

PARAMETERS:
	Int_T mm, Int_T nn
		New dimensions.
	
	int cntnts
		New contents bit mask.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solution::SetContents( Int_T mm, Int_T nn, int cntnts )
{
	assert( mm > 0 || ( mm == 0 && !( cntnts & (Dual|RowAct|Slack) ) ) );
	assert( nn > 0 || ( nn == 0 && !( cntnts & (Primal|RC) ) ) );

	SetM( mm );
	SetN( nn );

	contents = cntnts;
}


/*------------------------------------------------------------------------------

	void Solution::FreeSpace( void )

PURPOSE:
	Resets the solution's dimensions and contents flags. Frees the memory.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solution::FreeSpace( void )
{
	maxM = maxN = n = m = 0;
	contents = Empty;

	x.Resize( 0 );
	y.Resize( 0 );
	z.Resize( 0 );
	ra.Resize( 0 );
	s.Resize( 0 );

	result = 0.0;
}


/*------------------------------------------------------------------------------

	(virtual)
	void Solution::SetN( Int_T nn )
	void Solution::SetM( Int_T mm )

PURPOSE:
	Change the dimension of groups of solution vectors.

PARAMETERS:
	Int_T nn, Int_T mm
		The new dimensions >= 0.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solution::SetN( Int_T nn )
{
	assert( nn >= 0 );

	if( nn == n )
		return;
	
	if( nn == 0 )
	{
		contents &= ~( Primal | RC );
		x.Resize( 0 );
		z.Resize( 0 );

		n = maxN = 0;
	}
	else
	{
		if( nn < n )
			n = nn;
		else
		{
			if( nn > maxN )
			{
				maxN = Max( Int_T( 1.5 * maxN ), nn );
				x.Resize( maxN );
				z.Resize( maxN );
			}

			for( Int_T j = n; j < maxN; j++ )
				x[j] = z[j] = 0.0;

			n = nn;
		}
	}
}


void Solution::SetM( Int_T mm )
{
	assert( mm >= 0 );

	if( mm == m )
		return;
	
	if( mm == 0 )
	{
		contents &= ~( Dual | RowAct | Slack );
		y.Resize( 0 );
		ra.Resize( 0 );
		s.Resize( 0 );

		m = maxM = 0;
	}
	else
	{
		if( mm < m )
			m = mm;
		else
		{
			if( mm > maxM )
			{
				maxM = Max( Int_T( 1.5 * maxM ), mm );
				y.Resize( maxM );
				ra.Resize( maxM );
				s.Resize( maxM );
			}

			for( Int_T j = m; j < maxM; j++ )
				y[j] = ra[j] = s[j] = 0.0;

			m = mm;
		}
	}
}


/*------------------------------------------------------------------------------

	void Solution::WriteText( const char *fname, int cntnts ) const
	virtual void Solution::WriteText( FILE *fp, int cntnts ) const

	void Solution::ReadText( const char *fname, int cntnts )
	virtual void Solution::ReadText( FILE *fp, int cntnts )

PURPOSE:
	Write/read a solution to/from a text file.

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

void Solution::ReadText( const char *fname, int cntnts )
{
	assert( fname != NULL );
	FILE *fp = fopen( fname, "rt" );
	ReadText( fp, cntnts );
	fclose( fp );
}


void Solution::ReadText( FILE * /* fp */, int /* cntnts */ )
{
	FatalError( "Text solution file reading not implemented yet." );
}


void Solution::WriteText( const char *fname, int cntnts )
	const
{
	assert( fname != NULL );
	FILE *fp = fopen( fname, "wt" );
	WriteText( fp, cntnts );
	fclose( fp );
}


void Solution::WriteText( FILE *fp, int cntnts )
	const
{
	assert( fp != NULL );

	if( contents == Empty ) return;
	if( cntnts == Empty ) cntnts = (SolutionMask) contents;

	if( contents & ( RowAct | Slack | Dual ) )
	{
		assert( m > 0 );

		fprintf( fp,
			"EQUATIONS\n"
			"Row No  %s%s%s\n"
			"------  %s%s%s\n",

			( contents & RowAct )	? "  Row activity" : "",
			( contents & Slack )	? "  Slack" : "",
			( contents & Dual )		? "  Dual var." : "",

			( contents & RowAct )	? "  -------------" : "",
			( contents & Slack )	? "  -------------" : "",
			( contents & Dual )		? "  -------------" : ""
		);

		for( Int_T i = 0; i < m; i++ )
		{
			fprintf( fp, "%6d  ", (int) i );

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
			"Var No  %s%s\n"
			"------  %s%s\n",

			( contents & Primal )	? "   Primal value" : "",
			( contents & RC )		? "   Reduced cost" : "",

			( contents & Primal )	? "  -------------" : "",
			( contents & RC )		? "  -------------" : ""
		);

		for( Int_T i = 0; i < n; i++ )
		{
			fprintf( fp, "%6d  ", (int) i );

			if( contents & Primal )	fprintf( fp, "  %13.6E", x[i] );
			if( contents & RC )		fprintf( fp, "  %13.6E", z[i] );

			fprintf( fp, "\n" );
		}

		fprintf(  fp, "\n" );
	}

	if( contents & Res )
		fprintf( fp, "%-25s%20.12E", "OBJECTIVE:", result );

	fprintf( fp, "\n" );
}


#ifdef SUPPORT_LP_DIT

/*------------------------------------------------------------------------------

	void Solution::ReadDIT( const char *fname )
	void Solution::WriteDIT( const char *fname )

PURPOSE:
	Reads/writes a solution in LP-DIT format.

PARAMETERS:
	const char *fname
		Input/output file name.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solution::ReadDIT( const char *fname )
	{ lp_res( 0, 0, (char *)fname, (void *) this ); }


void Solution::WriteDIT( const char *fname )
	{ lp_res( 1, 0, (char *)fname, (void *) this ); }


//
//	These two functions are used to make communication with C library LP-DIT
//	possible. They relay the C calls to C++ code. Incidentally, this makes
//	it possible for us to use e.g. virtual functions. The only problem is that
//	we have to rely on (otherwise dubious) type casts.
//
extern "C"
{
	void lpi_res( LP_SOLUTION *s, void *user )
		{ ((Solution *)user)->lpi_res( s ); }

	void lpo_res( LP_SOLUTION *s, void *user )
		{ ((Solution *)user)->lpo_res( s ); }
}


#if defined( explicit_templates )
	template LP_VECT *MALLOC( LP_VECT *& Table, size_t len );
#endif


/*------------------------------------------------------------------------------

	void Solution::lpi_res( LP_SOLUTION *s )
	void Solution::lpo_res( LP_SOLUTION *s )

PURPOSE:
	Communicates with LP-DIT library (via a C-style wrap-around functions given
	above).

PARAMETERS:
	LP_SOLUTION *s
		Solution structure as seen by LP-DIT.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void Solution::lpi_res( LP_SOLUTION *lps )
{
	if( lps->cols.el == NULL && lps->cols_d.el == NULL &&
		lps->rows.el == NULL && lps->rows_d.el == NULL )
	//
	//	The first call: we are only supposed to store the dimensions of the
	//	problem.
	//
	{
		lps->status = S_UNF;
		lps->comment[0] = '\0';
		lps->time1 = lps->time2 = 0.0;
		lps->cols.elems = lps->cols_d.elems = lps->rows.elems =
		lps->rows_d.elems = 0;
		lps->objv = 0.0;

		if( contents == Empty ) return;

		lps->n = (LP_IND)n;
		lps->m = (LP_IND)m;

		//----------------------------------------------------------------------

		lps->cols.elems		= (LP_IND)( ( contents & Primal )	? n : 0 );
		lps->cols_d.elems	= (LP_IND)( ( contents & RC )		? n : 0 );
		lps->rows.elems		= (LP_IND)( ( contents & RowAct )	? m : 0 );
		lps->rows_d.elems	= (LP_IND)( ( contents & Dual )		? m : 0 );
	}
	else if( ( !( contents & Primal )	|| lps->cols.el != NULL ) &&
		( !( contents & RC )			|| lps->cols_d.el != NULL ) &&
		( !( contents & RowAct )		|| lps->rows.el != NULL ) &&
		( !( contents & Dual )			|| lps->rows_d.el != NULL ) )
	//
	//	The second call: now we store the data.
	//
	{
		if( contents == Empty ) return;

		switch( status )
		{
			case Unknown:		lps->status = S_UNF;		break;
			case Optimal:		lps->status = S_OPT;		break;
			case Unbounded:		lps->status = S_UNB;		break;
			case Infeasible:	lps->status = S_LP_INF;		break;
			case Unsolved:		lps->status = S_UNF;		break;
		}

		lps->objv = ( contents & Res ) ? result : 0.0;

		//----------------------------------------------------------------------
		if( contents & Primal )
		{
			assert( lps->cols.elems == n );
			assert( lps->cols.el != NULL );

			for( Int_T j = 0; j < n; j++ )
			{
				LP_VECT &v = lps->cols.el[j];

				v.index = (LP_IND)j;
				v.value = x[j];
			}
		}
		else
			assert( lps->cols.elems == 0 );

		//----------------------------------------------------------------------
		if( contents & RC )
		{
			assert( lps->cols_d.elems == n );
			assert( lps->cols_d.el != NULL );

			for( Int_T j = 0; j < n; j++ )
			{
				LP_VECT &v = lps->cols_d.el[j];

				v.index = (LP_IND)j;
				v.value = z[j];
			}
		}
		else
			assert( lps->cols_d.elems == 0 );

		//----------------------------------------------------------------------
		if( contents & RowAct )
		{
			assert( lps->rows.elems == m );
			assert( lps->rows.el != NULL );

			for( Int_T i = 0; i < m; i++ )
			{
				LP_VECT &v = lps->rows.el[i];

				v.index = (LP_IND)i;
				v.value = ra[i];
			}
		}
		else
			assert( lps->rows.elems == 0 );

		//----------------------------------------------------------------------
		if( contents & Dual )
		{
			assert( lps->rows_d.elems == m );
			assert( lps->rows_d.el != NULL );

			for( Int_T i = 0; i < m; i++ )
			{
				LP_VECT &v = lps->rows_d.el[i];

				v.index = (LP_IND)i;
				v.value = y[i];
			}
		}
		else
			assert( lps->rows_d.elems == 0 );
	}
	else
		FatalError( "LP-DIT data header invalid." );
}


void Solution::lpo_res( LP_SOLUTION *lps )
{
	maxM = m = lps->m;
	maxN = n = lps->n;

	contents = Empty;

	y.Resize( m );			if( m ) y.Fill( 0.0, m );
	z.Resize( n );			if( n ) z.Fill( 0.0, n );
	ra.Resize( m );			if( m ) ra.Fill( 0.0, m );
	s.Resize( 0 );

	result = 0.0;

	switch( lps->status )
	{
		case S_UNF:		status = Unknown;		break;
		case S_OPT:		status = Optimal;		break;
		case S_UNB:		status = Unbounded;		break;
		case S_LP_INF:	status = Infeasible;	break;
	}

	contents |= Res;
	result = lps->objv;

	//--------------------------------------------------------------------------
	if( lps->cols.elems > 0 )
	{
		contents |= Primal;

		x.Resize( n );
		x.Fill( 0.0, n );

		assert( lps->cols.elems <= n );

		for( Int_T j = 0, nn = lps->cols.elems; j < nn; j++ )
		{
			LP_VECT &v = lps->cols.el[j];

			x[v.index] = v.value;
		}
	}
	else
		x.Resize( 0 );

	//--------------------------------------------------------------------------
	if( lps->cols_d.elems > 0 )
	{
		contents |= RC;

		z.Resize( n );
		z.Fill( 0.0, n );

		assert( lps->cols_d.elems <= n );

		for( Int_T j = 0, nn = lps->cols_d.elems; j < nn; j++ )
		{
			LP_VECT &v = lps->cols_d.el[j];

			z[v.index] = v.value;
		}
	}
	else
		z.Resize( 0 );

	//--------------------------------------------------------------------------
	if( lps->rows.elems > 0 )
	{
		contents |= RowAct;

		ra.Resize( m );
		ra.Fill( 0.0, m );

		assert( lps->rows.elems <= m );

		for( Int_T i = 0, mm = lps->rows.elems; i < mm; i++ )
		{
			LP_VECT &v = lps->rows.el[i];

			ra[v.index] = v.value;
		}
	}
	else
		ra.Resize( 0 );

	//--------------------------------------------------------------------------
	if( lps->rows_d.elems > 0 )
	{
		contents |= Dual;

		y.Resize( m );
		y.Fill( 0.0, m );

		assert( lps->rows_d.elems <= m );

		for( Int_T i = 0, mm = lps->rows_d.elems; i < mm; i++ )
		{
			LP_VECT &v = lps->rows_d.el[i];

			y[v.index] = v.value;
		}
	}
	else
		y.Resize( 0 );
}

#endif
