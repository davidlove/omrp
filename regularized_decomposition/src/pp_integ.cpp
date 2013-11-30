/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	pp_integ.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		presolve.h, simplex.h, error.h, smartptr.h, memblock.h,
					smartdcl.h, sptr_deb.h, sptr_ndb.h, stdtype.h, lp_codes.h,
					solv_lp.h, std_tmpl.h, compile.h, mps_lp.h, sort_lab.h,
					myalloc.h, std_math.h, cl_list.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	One function used (in DEBUG mode only) for checking the integrity of the
linear problem constraint matrix. The matrix is stored both by rows and by
columns. Both representations are compared and any mismatches cause the program
to crash.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	void Presolver::CheckMatrixIntegrity( void )

------------------------------------------------------------------------------*/

#ifndef NDEBUG

#ifndef __PRESOLVE_H__
#	include "presolve.h"
#endif
#ifndef __STD_MATH_H__
#	include "std_math.h"
#endif

void Presolver::CheckMatrixIntegrity( void )
{
	//--------------------------------------------------------------------------
	//	Search on matrix rows. All non-zeros found are then looked up in the
	//	column file. If they're not found, or found different, an alarm is
	//	triggered.
	//
	{
		for( Int_T i = 0; i < m; i++ )
			if( !ExcludeRows[i] )
			{
				Ptr<Real_T> ar;
				Ptr<Int_T> col;
				Int_T rLen;

				for( LP.GetRow( i, ar, col, rLen ); rLen; --rLen, ++ar, ++col )
					if( !ExcludeCols[ *col ] )
					{
						Ptr<Real_T> ac;
						Ptr<Int_T> row;
						Int_T cLen;
						Bool_T Found = False;

						for( LP.GetColumn( *col, ac, row, cLen ); cLen;
							--cLen, ++ac, ++row )
							if( *row == i )
								if( IsEqual( *ar, *ac ) )
								{
									Found = True;
									break;
								}
#ifndef NDEBUG
								else
									abort();
#endif

						assert( Found );
					}
			}
	}

	//--------------------------------------------------------------------------
	//	Search on matrix columns. All non-zeros found are then looked up in the
	//	row file. If they're not found, or found different, an alarm is
	//	triggered.
	//
	{
		for( Int_T j = 0; j < n; j++ )
			if( !ExcludeCols[j] )
			{
				Ptr<Real_T> ac;
				Ptr<Int_T> row;
				Int_T cLen;

				for( LP.GetColumn( j, ac, row, cLen ); cLen;
					--cLen, ++ac, ++row )
					if( !ExcludeRows[ *row ] < 0 )
					{
						Ptr<Real_T> ar;
						Ptr<Int_T> col;
						Int_T rLen;
						Bool_T Found = False;

						for( LP.GetRow( *row, ar, col, rLen ); rLen;
							--rLen, ++ar, ++col )
							if( *col == j )
								if( IsEqual( *ac, *ar ) )
								{
									Found = True;
									break;
								}
#ifndef NDEBUG
								else
									abort();
#endif

						assert( Found );
					}
			}
	}
}
#endif
