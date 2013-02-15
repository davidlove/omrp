/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:       Simplex
PROJECT FULL NAME:	Advanced implementation of revised simplex method
					for large scale linear problems.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej P. Wierzbicki, dr Jacek Gondzio.

--------------------------------------------------------------------------------

SOURCE FILE NAME:	solver4.cpp
CREATED:			1993.11.01
LAST MODIFIED:		1996.02.14

DEPENDENCIES:		error.h, stdtype.h, solver.h, solv_lp.h, lp_codes.h,
					mps_lp.h, inverse.h, solvcode.h, simplex.h,

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

#ifndef __ERROR_H__
#	include "error.h"
#endif

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SOLVER_H__
#	include "solver.h"
#endif
#ifndef __SMPLX_LP_H__
#	include "smplx_lp.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __LP_SOL_H__
#	include "lp_sol.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __SOLVTOL_H__
#	include "solvtol.h"
#endif


/*------------------------------------------------------------------------------

	Int_T Solver::BackTrack( Int_T Steps )

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

Int_T Solver::BackTrack( Int_T )
{
	return 0;
}


/*------------------------------------------------------------------------------

	Solution *Solver::GetSolution( int mode, Bool_T LP_Valid )

PURPOSE:
	Creates and fills a Solution structure with the current solution of the
linear problem.

PARAMETERS:
	int mode
		Specifies the information to be stored in the solution structure (see
		header file "solution.h" for enumerations and their meaning).

	Bool_T LP_Valid
		If "True", the linear problem is assumed to be still valid, otherwise
		it is assumed to be out of date. This flag is to be used after
		presolving was applied. Thus a solution that does not refer to the
		LP object has to be created.

RETURN VALUE:
	A pointer to a solution structure created by operator "new". The structure
has to be dealocated by the caling procedure.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Solution *Solver::GetSolution( int mode, Bool_T LP_Valid )
{
	assert( Initialized );

	if( !LP_Valid ) mode &= Solution::Primal | Solution::Res;

	const Int_T n	= LP.GetStructN(),
		m			= Int_T( LP_Valid ? LP.GetM() : 0 );

	//--------------------------------------------------------------------------
	//	Create a solution object.
	//
	Solution *sol = ( LP_Valid ) ?
		new LP_Solution( &LP, mode ) :
		new Solution( m, n, mode );
	
	if( !sol )
		FatalError( "Not enough memory to create solution structure." );

	sol->SetContents( m, n, mode );

	//--------------------------------------------------------------------------
	//	If required, store primal variables.
	//
	if( mode & Solution::Primal )
	{
		for( Int_T i = 0; i < N; i++ )
			if( IsZero( x[i] ) ) x[i] = 0.0;
		sol->x.Copy( x, N, n, n );
	}

	//--------------------------------------------------------------------------
	//	If required, store reduced costs.
	//
	if( mode & Solution::RC )
	{
		for( Int_T i = 0; i < N; i++ )
			if( IsZero( z[i] ) ) z[i] = 0.0;
		sol->z.Copy( z, N, n, n );
	}

	//--------------------------------------------------------------------------
	//	If required, store dual variables.
	//
	if( mode & Solution::Dual )
	{
		for( Int_T i = 0; i < M; i++ )
			if( IsZero( y[i] ) ) y[i] = 0.0;
		sol->y.Copy( y, M, m, m );
	}

	//--------------------------------------------------------------------------
	//	If required, store result. Store solution status.
	//
	if( mode & Solution::Res )
		sol->result = Result;

	sol->SetStatus( Solution::Optimal );

	//--------------------------------------------------------------------------
	//	If required, compute and store row activity.
	//
	if( mode & Solution::RowAct )
	{
		sol->ra.Resize( m );

		for( Int_T i = 0; i < m; i++ )
		{
			Ptr<Real_T> a;
			Ptr<Int_T> col;
			Int_T len;
			Real_T &ra = sol->ra[i];

			ra = 0.0;
			for( LP.SolvableLP::GetRow( i, a, col, len ); len; --len, ++a, ++col )
				ra += x[*col] * *a;
		}
	}

	return sol;
}


/*------------------------------------------------------------------------------

	void Solver::OutputBasisToFile( FILE *fp )

PURPOSE:
	Generates a standard simplex basis (expressed in terms of basic rows and
columns).

PARAMETERS:
	FILE *fp
		A stream to utput the basis file to.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/
static void OutputSection( FILE *fp, Array<Int_T> &Pos, Int_T n,
	const SortedArrayOfLabels &Labels );

void Solver::OutputBasisToFile( FILE *fp )
{
	assert( fp != NULL );

	//--------------------------------------------------------------------------
	//	Output the basis file header if an optimal solution was found.
	//
	assert( Status >= SR_OPTIMUM );

	if( Status != SR_OPTIMUM )
	{
		fprintf( fp, "BASIS         %-8s  ", LP.GetName() );

		if( Status == SR_INFEASIBLE )
			fprintf( fp, "INFEASIBLE\n" );
		else if( Status == SR_UNBOUNDED )
			fprintf( fp, "UNBOUNDED\n" );
		else
			fprintf( fp, "UNSOLVED\n" );

		return;
	}

	//--------------------------------------------------------------------------
	//	If necessary (i.e. if a solution was found) construct the basis.
	//
	Int_T StructN		= LP.GetStructN(),
		i;
	Array<Int_T> &Rows	= w1Ind;
	WorkVector<Int_T> Cols( StructN );

	//--------------------------------------------------------------------------
	//	Find basic columns (a subset of the structural columns).
	//
	Cols.Fill( A2B_UNDEF, StructN );
	Int_T j;
	for( j = 0; j < StructN; j++ )
		if( A2B[j] >= 0 )
			Cols[j] = A2B_BASIC;
		else
			Cols[j] = A2B[j];

	//--------------------------------------------------------------------------
	//	Find basic rows (by searching slack and artificial columns).
	//
	Rows.Fill( A2B_UNDEF, M );
	for( ; j < N; j++ )
	{
		Ptr<Real_T> a;
		Ptr<Int_T> row;
		Int_T len;

		LP.GetColumn( j, a, row, len );
		assert( len == 1 );
		i = *row;

		if( !( VarType[j] & VT_ARTIF ) || ( Rows[i] == A2B_UNDEF ) )
		{
			if( A2B[j] >= 0 )
				Rows[i] = A2B_BASIC;
			else
			{
				switch( A2B[j] )
				{
				case A2B_LO:
					Rows[i] = ( LP.GetRowType(i) & (RT_EQ|RT_GE) )?
						A2B_LO : A2B_UP;
					break;

				case A2B_UP: 
					Rows[i] = ( LP.GetRowType(i) & (RT_EQ|RT_LE) )?
						A2B_LO : A2B_UP;
					break;

				case A2B_IN:
					Rows[i] = A2B[i];
					break;

#ifndef NDEBUG
				default:
					abort();
#endif
				}
			}
		}
	}

#ifndef NDEBUG
	//--------------------------------------------------------------------------
	//	Count the number of the basic variables.
	//
	Int_T InBasis	= 0;

	for( j = 0; j < StructN; j++ )
		if( Cols[j] == A2B_BASIC )
			InBasis++;
	for( i = 0; i < M; i++ )
	{
		assert( Rows[i] != A2B_UNDEF );
		if( Rows[i] == A2B_BASIC )
			InBasis++;
	}
	assert( InBasis == M );
#endif

	//--------------------------------------------------------------------------
	//	Output the optimal basis.
	//
	assert( Status >= SR_OPTIMUM );
	fprintf( fp, "BASIS         %-8s  OPTIMAL\n", LP.GetName() );

	fprintf( fp, "COLUMNS                 %d\n", StructN );
	OutputSection( fp, Cols, StructN, LP.RevealColumnLabels() );

	fprintf( fp, "ROWS                    %d\n", M );
	OutputSection( fp, Rows, M, LP.RevealRowLabels() );

	fprintf( fp, "ENDATA\n" );
}


static void OutputSection( FILE *fp, Array<Int_T> &Pos, Int_T n, // )
	const SortedArrayOfLabels &Labels )
{
	const char *PosCode[] = { "LO", "UP", "IN", "BA", "??" };

	for( Int_T i = 0; i < n; i++ )
	{
		Int_T c = 3;

		switch( Pos[i] )
		{
		case A2B_LO:	c = 0; break;
		case A2B_UP:	c = 1; break;
		case A2B_IN:	c = 2; break;
		case A2B_BASIC:	c = 3; break;
		default:		c = 4; break;
		}

		fprintf( fp, " %2s %-8s\n", PosCode[c], Labels.FindLabel( i ) );
	}
}


/*------------------------------------------------------------------------------

	SolverStateDump::SolverStateDump( Int_T n, const Array<Int_T> &a2b,
		const Array<Real_T> &x, Int_T AllocN )

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

SolverStateDump::SolverStateDump( Int_T n, const Array<Int_T> &a2b, // )
	const Array<Real_T> &x, Int_T AllocN )
	: N( n ), A2B( AllocN, -1 ), X( AllocN, 0.0 )
{
	assert( n > 0 );

	A2B.Copy( a2b, n, n );
	X.Copy( x, n, n );
}


/*------------------------------------------------------------------------------

	SolverStateDump *Solver::GetSolverStateDump( void )

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

SolverStateDump *Solver::GetSolverStateDump( void )
{
	assert( Initialized );

	SolverStateDump *dump = new SolverStateDump( N, A2B, x, AllocN );

	if( dump == NULL ) FatalError( "Out of memory." );
	return dump;
}


/*------------------------------------------------------------------------------

	void Solver::CheckA2B_Consistency( void )

PURPOSE:
	Checks if the "A2B" array is consistent with the values of primal
variables "x" and upper bounds "u".

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

#ifndef NDEBUG
void Solver::CheckA2B_Consistency( void )
{
	assert( Initialized );
	assert( N > 0 );

	Real_T ftol = 10 * FEASIBILITY_TOL_DEF;

	for( Int_T i = 0; i < N; i++ )
		if( A2B[i] >= 0 )
			assert( x[i] >= -ftol && x[i] <= u[i] +ftol );
		else if( A2B[i] == A2B_LO )
			assert( x[i] >= -ftol && x[i] <= ftol );
		else if( A2B[i] == A2B_UP )
			assert( x[i] >= u[i] - ftol && x[i] <= u[i] + ftol );
		else if( A2B[i] == A2B_IN )
			assert( x[i] >= ftol && x[i] <= u[i] -ftol );
#ifndef NDEBUG
		else
			abort();
#endif
}
#endif


Real_T Solver::DualityGap( void )
const
{
	Real_T gap = 0.0;

	for( Int_T j = 0; j < N; j++ )
	{
		switch( A2B[j] )
		{
		case A2B_LO:	gap += fabs( x[j] * z[j] );				break;
		case A2B_UP:	gap += fabs( ( u[j] - x[j] ) * z[j] );	break;
		case A2B_IN:	gap += fabs( x[j] * z[j]);				break;
		default:		assert( A2B[j] >= 0 ); break;
		}
	}

	return fabs( gap );
}
