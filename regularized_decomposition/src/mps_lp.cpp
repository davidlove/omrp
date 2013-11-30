/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	mps_lp.cpp
CREATED:			1993.09.16
LAST MODIFIED:		1995.12.27

DEPENDENCIES:		stdtype.h, mps_lp.h, myalloc.h, sort_lab.h, work_vec.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	This source file contains constructor and destructor for MPS_LP class
objects.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	MPS_LP::MPS_LP( void )
	MPS_LP::~MPS_LP( void )

PRIVATE FUNCTIONS:
	void MPS_LP::FreeStorage( void )

STATIC FUNCTIONS:
	None.

STATIC DATA:
	None.

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	TOSTD_LOBND,
	TOSTD_UPBND			-	HAVE TO BE DEFINED! Those macros denote
							respectively lower and upper bound on admissible
							non-zero absolute value in linear problem.

------------------------------------------------------------------------------*/

#include <math.h>

#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif
#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif


/*------------------------------------------------------------------------------

	MPS_LP::MPS_LP( void )

PURPOSE:
	This is class MPS_LP constructor. It sets all pointers to NULL and all
problem dimensions to 0. Additionally it sets object status (MPS_LP::OS) of
the newly created object to OS_EMPTY.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

MPS_LP::MPS_LP( void )
	: OS( OS_EMPTY ),
	n( 0 ), m( 0 ), nz( 0 ),
	RowLabels( 50 ), ColLabels( 100 ), ObjRow( -1 ),

	//--------------------------------------------------------------------------
	//  Statistic counters.
	//
	mE( 0 ), mG( 0 ), mL( 0 ), mR( 0 ), mF( 0 ),
	nPL( 0 ), nFX( 0 ), nFR( 0 ), nMI( 0 ), nUP( 0 )
{
	*Name = *RHS_Name = *BoundsName = *RangesName = '\0';
}


/*------------------------------------------------------------------------------

	void MPS_LP::FreeStorage( void )

PURPOSE:
	This is a function that deallocates all storage that may have been
previously allocated. Additionally all problem dimensions are reset to zero
and object staus is set to OS_EMPTY. This effectively enables you to reuse
the object. The function is called when problem reading fails or when object is
deleted.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void MPS_LP::FreeStorage( void )
{
	RowLabels.FreeMemory();
	ColLabels.FreeMemory();

	*Name = *RHS_Name = *BoundsName = *RangesName = '\0';
	n = m = 0;					nz = 0;
	OS = OS_EMPTY;
}


/*------------------------------------------------------------------------------

	void MPS_LP::SetL( Int_T j, Real_T l_j, Bool_T Finite )
	void MPS_LP::SetU( Int_T j, Real_T u_j, Bool_T Finite )
 
PURPOSE:
	Set a new value for any given simple bound on variable.
 
PARAMETERS:
	Bool_T Finite
		Says whether the bound in question is finite or not.

	Real_T u_j / l_j
		The finite value of the bound (ignored if 'Finite' == False).

	Int_T j
		The number of the variable to which the bound refers.
 
RETURN VALUE:
	None.
 
SIDE EFFECTS:
	None.
 
------------------------------------------------------------------------------*/

void MPS_LP::SetL( Int_T j, Real_T l_j, Bool_T Finite )
{
	if( Finite )
	{
		assert( l_j > -INFINITY );
 
		VarType[j] |= VT_LO;
		l[j] = l_j;

		if( !( VarType[j] & VT_FX ) && IsEqual( u[j], l[j] ) )
			VarType[j] |= VT_FX;
	}
	else
	{	
		VarType[j] &= ~((Short_T) VT_LO );  
		l[j] = -INFINITY; 
	}

	if( IsNotEqual( l[j], u[j] ) )
		VarType[j] &= ~((Short_T) VT_FX );
}
 
 
void MPS_LP::SetU( Int_T j, Real_T u_j, Bool_T Finite )
{
	if( Finite ) 
	{
		assert( u_j < +INFINITY );
 
		VarType[j] |= VT_UP;
		u[j] = u_j;

		if( !( VarType[j] & VT_FX ) && IsEqual( u[j], l[j] ) )
			VarType[j] |= VT_FX;
	}  
	else
	{	
		VarType[j] &= ~((Short_T)VT_UP);
		u[j] = +INFINITY; 
	}

	if( IsNotEqual( l[j], u[j] ) )
		VarType[j] &= ~((Short_T) VT_FX );
}


/*------------------------------------------------------------------------------

	Bool_T MPS_LP::CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin,
		Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow )

PURPOSE:
	Creates a linear problem as a fragment of a given problem 'Src' using rows
'RowMin' to 'RowMax' and columns 'ColMin' to 'ColMax' as well as a free row
'CostRow'. All other free rows are omitted.
	If 'CostRow' is somewhere in between 'RowMin' and 'RowMax', then this part
of the matrix is copied without any changes. Otherwise the cost row is placed
as the last row of the '*this' matrix.

PARAMETERS:
	const MPS_LP &Src
		Source MPS_LP object: we take all data from this matrix (and its
		vectors).

	Int_T RowMin, Int_T RowMax
		The range of rows ("RowMin" included, "RowMax" not included)..

	Int_T ColMin, Int_T ColMax
		The range of columns ("ColMin" included, "ColMax" not included)..

	Int_T CostRow
		Objective function row.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T MPS_LP::CreateAsSubmatrix( const MPS_LP &Src, Int_T RowMin, // )
	Int_T RowMax, Int_T ColMin, Int_T ColMax, Int_T CostRow )
{
	OS = OS_EMPTY;

	//--------------------------------------------------------------------------
	//	Check if the arguments make sense.
	//
	assert( RowMin >= 0 && RowMin <= Src.GetM() );
	assert( RowMax >= 0 && RowMax <= Src.GetM() );
	assert( RowMin <= RowMax );

	assert( ColMin >= 0 && ColMin <= Src.GetN() );
	assert( ColMax >= 0 && ColMax <= Src.GetN() );
	assert( ColMin <= ColMax );

	assert( CostRow >= 0 && CostRow <= Src.GetM() );

	Bool_T CostInRng = (CostRow >= RowMin && CostRow < RowMax) ? True : False;

	//--------------------------------------------------------------------------
	//	Reset all the object data (regarding dimensions).
	//
	OS = OS_EMPTY;
	*Name = *RHS_Name = *BoundsName = *RangesName = '\0';
	n = m = 0;
	nz = 0;
	RowLabels.FreeMemory();
	ColLabels.FreeMemory();

	mE = mG = mL = mR = mF = nPL = nFX = nFR = nMI = nUP = 0;

	//--------------------------------------------------------------------------
	//	Calculate the dimensions of the problem by a scan of the 'Src' LP.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;

	n = Int_T( ColMax - ColMin );
	m = Int_T( RowMax - RowMin + ( ( CostInRng ) ? 0 : 1) );
	for( Int_T j = ColMin; j < ColMax; j++ )
	{
		Int_T len;

		for( Src.GetColumn( j, a, row, len ); len; --len, ++row )
			if( ( *row >= RowMin && *row < RowMax &&
				Src.RowType[ *row ] != RT_FR ) || *row == CostRow )
				nz++;
	}

	//--------------------------------------------------------------------------
	//	Resize all arrays to match the dimensions that were just computed.
	//
	RowType.Resize( m );
	ac.Resize( nz );
	Row.Resize( nz );
	ColStart.Resize( n + 1 );
	b.Resize( m );
	VarType.Resize( n );
	l.Resize( n );
	u.Resize( n );
	r.Resize( m );

	//--------------------------------------------------------------------------
	//	Copy the contents of the 'Src' constraint matrix into *this* matrix.
	//
	if( CostInRng )
	{
		Int_T pos, jj;

		for( jj = ColMin, ColStart[0] = 0, pos = 0; jj < ColMax;
			ColStart[ ++jj - ColMin ] = pos )
		{
			Int_T len;

			for( Src.GetColumn( jj, a, row, len ); len; --len, ++row, ++a )
				if( *row >= RowMin && *row < RowMax &&
					( Src.RowType[ *row ] != RT_FR || *row == CostRow ) )
				{
					ac[pos]		= *a;
					Row[pos]	= Int_T( *row - RowMin );
					pos++;
				}
		}

		assert( pos == nz );
	}
	else
	{
		Int_T pos, jj;

		for( jj = ColMin, ColStart[0] = 0, pos = 0; jj < ColMax;
			ColStart[ ++jj - ColMin ] = pos )
		{
			Int_T len;

			for( Src.GetColumn( jj, a, row, len ); len; --len, ++row, ++a )
				if( *row >= RowMin && *row < RowMax &&
					Src.RowType[ *row ] != RT_FR )
				{
					ac[pos]		= *a;
					Row[pos]	= Int_T( *row - RowMin );
					pos++;
				}
				else if( *row == CostRow )
				{
					ac[pos]		= *a;
					Row[pos]	= Int_T( m - 1 );
					pos++;
				}
		}

		assert( pos == nz );
	}

	//--------------------------------------------------------------------------
	//	Copy the contents of all vectors. Mind the fact that copying depends
	//	on whether objective function belongs to the range of rows, or not.
	//
	if( m > 0 )
	{
		if( CostInRng )
		{
			RowType.Copy( 	Src.RowType,	Src.m, m, m,		RowMin );
			b.Copy( 		Src.b,			Src.m, m, m,		RowMin );
			r.Copy( 		Src.r,			Src.m, m, m,		RowMin );
		}
		else
		{
			RowType.Copy( 	Src.RowType,	Src.m, m, m - 1,	RowMin );
			b.Copy( 		Src.b,			Src.m, m, m - 1,	RowMin );
			r.Copy( 		Src.r,			Src.m, m, m - 1,	RowMin );

			RowType[m - 1] = RT_FR;
			b[m - 1] = 0.0;
			r[m - 1] = +INFINITY;
		}
	}

	if( n > 0 )
	{
		l.Copy( Src.l,				Src.n, n, n,	ColMin );
		u.Copy( Src.u,				Src.n, n, n,	ColMin );
		VarType.Copy( Src.VarType,	Src.n, n, n,	ColMin );
	}

	OS = OS_NO_LABELS;

	//--------------------------------------------------------------------------
	//	Copy the labels from the source problem.
	//
	if( n > 0 )
	{
		for( Int_T jj = ColMin; jj < ColMax; jj++ )
			ColLabels.AddLabel( Src.ColLabels.FindLabel( jj ) );
		ColLabels.SortLabels();
	}

	assert( ColLabels.NumberOfLabels() == n );

	if( m > 0 )
	{
		for( Int_T jj = RowMin; jj < RowMax; jj++ )
			RowLabels.AddLabel( Src.RowLabels.FindLabel( jj ) );
		if( !CostInRng )
			RowLabels.AddLabel( Src.RowLabels.FindLabel( CostRow ) );
		RowLabels.SortLabels();
	}

	assert( RowLabels.NumberOfLabels() == m );

	strncpy( Name,			Src.Name,		LAB_LEN+1 );
	strncpy( RHS_Name,		Src.RHS_Name,	LAB_LEN+1 );
	strncpy( BoundsName,	Src.BoundsName,	LAB_LEN+1 );
	strncpy( RangesName,	Src.RangesName,	LAB_LEN+1 );

	Name[LAB_LEN]		= '\0';
	RHS_Name[LAB_LEN]	= '\0';
	BoundsName[LAB_LEN]	= '\0';
	RangesName[LAB_LEN]	= '\0';

	//--------------------------------------------------------------------------
	//	Mark the current object as filled with data.
	//
	OS = OS_FULL;

	return True;
}


Bool_T MPS_LP::CreateAsSubmatrix( const MPS_LP &Src, // )
	const Array<Bool_T> &IncludeRows,
	const Array<Bool_T> &IncludeCols )
{
	OS = OS_EMPTY;

	//--------------------------------------------------------------------------
	//	Reset all the object data (regarding dimensions).
	//
	OS = OS_EMPTY;
	*Name = *RHS_Name = *BoundsName = *RangesName = '\0';
	n = m = 0;
	nz = 0;
	RowLabels.FreeMemory();
	ColLabels.FreeMemory();

	mE = mG = mL = mR = mF = nPL = nFX = nFR = nMI = nUP = 0;

	//--------------------------------------------------------------------------
	//	Calculate the dimensions of the problem by a scan of the 'Src' LP.
	//
	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	Int_T SrcM	= Src.GetM(),
		SrcN	= Src.GetN();

	Int_T i, j;
	for( i = 0; i < SrcM; i++ ) if( IncludeRows[i] ) m++;
	for( j = 0; j < SrcN; j++ ) if( IncludeCols[j] ) n++;

	for( j = 0; j < SrcN; j++ )
		if( IncludeCols[j] )
			for( Src.GetColumn( j, a, row, len ); len; --len, ++row )
				if( IncludeRows[ *row ] )
					nz++;

	//--------------------------------------------------------------------------
	//	Resize all arrays to match the dimensions that were just computed.
	//
	RowType.Resize( m );
	ac.Resize( nz );
	Row.Resize( nz );
	ColStart.Resize( n + 1 );
	b.Resize( m );
	VarType.Resize( n );
	l.Resize( n );
	u.Resize( n );
	r.Resize( m );

	//--------------------------------------------------------------------------
	//	Compute the row permutation vector 'RowNum'.
	//	Copy the contents of the 'Src' constraint matrix into *this* matrix.
	//
	WorkVector<Int_T> RowNum( SrcM );
	Int_T ii;

	RowNum.Fill( -1, SrcM );
	for( i = ii = 0; i < SrcM; i++ )
		if( IncludeRows[i] ) RowNum[i] = ii++;

	assert( ii == m );
		
	Int_T CntCol, CntNZ;

	for( j = 0, ColStart[0] = 0, CntCol = 0, CntNZ = 0; j < SrcN; ++j )
	{
		if( !IncludeCols[j] ) continue;

		for( Src.GetColumn( j, a, row, len ); len; --len, ++row, ++a )
			if( IncludeRows[ *row ] )
			{
				ac[CntNZ]	= *a;
				Row[CntNZ]	= RowNum[ *row ];
				CntNZ++;
			}

		ColStart[ ++CntCol ] = CntNZ;
	}

	assert( CntNZ == nz );

	//--------------------------------------------------------------------------
	//	Copy the contents of all vectors. Mind the fact that copying depends
	//	on whether objective function belongs to the range of rows, or not.
	//
	if( m > 0 )
	{
		Int_T newI = 0;

		for( i = 0; i < SrcM; i++ )
			if( IncludeRows[i] )
			{
				RowType[newI]	= Src.RowType[i];
				b[newI]			= Src.b[i];
				r[newI]			= Src.r[i];
				newI++;
			}
	}

	if( n > 0 )
	{
		Int_T jj = 0;

		for( j = 0; j < SrcN; j++ )
			if( IncludeCols[j] )
			{
				l[jj]		= Src.l[j];
				u[jj]		= Src.u[j];
				VarType[jj]	= Src.VarType[j];
				jj++;
			}
	}

	OS = OS_NO_LABELS;

	//--------------------------------------------------------------------------
	//	Copy the labels from the source problem.
	//
	if( n > 0 )
	{
		for( j = 0; j < SrcN; j++ )
			if( IncludeCols[j] )
				ColLabels.AddLabel( Src.ColLabels.FindLabel( j ) );
		ColLabels.SortLabels();
	}

	assert( ColLabels.NumberOfLabels() == n );

	if( m > 0 )
	{
		for( i = 0; i < SrcM; i++ )
			if( IncludeRows[i] )
				RowLabels.AddLabel( Src.RowLabels.FindLabel( i ) );
		RowLabels.SortLabels();
	}

	assert( RowLabels.NumberOfLabels() == m );

	strncpy( Name,			Src.Name,		LAB_LEN+1 );
	strncpy( RHS_Name,		Src.RHS_Name,	LAB_LEN+1 );
	strncpy( BoundsName,	Src.BoundsName,	LAB_LEN+1 );
	strncpy( RangesName,	Src.RangesName,	LAB_LEN+1 );

	Name[LAB_LEN]		= '\0';
	RHS_Name[LAB_LEN]	= '\0';
	BoundsName[LAB_LEN]	= '\0';
	RangesName[LAB_LEN]	= '\0';

	//--------------------------------------------------------------------------
	//	Mark the current object as filled with data.
	//
	OS = OS_FULL;

	return True;
}


/*------------------------------------------------------------------------------

	void MPS_LP::GetConstraintRange( Int_T row, Real_T &bl, Real_T &bu ) const
	void MPS_LP::SetConstraintRange( Int_T i, Real_T min, Bool_T minFinite,
		Real_T max, Bool_T maxFinite )

PURPOSE:
	Standard MPS file notion of expressing the row activity bounds by the
right hand side and range is equivalent to explicitly giving those bounds.
The two functions named above perform the necessary conversions.

PARAMETERS:
	Int_T row, Real_T &bl, Real_T &bu
		Row (constraint) number and activity bounds. "bl" and "bu" are both
		set by the function.

	Int_T i, Real_T min, Bool_T minFinite, Real_T max, Bool_T maxFinite
		Row (constraint) number and activity values. Flags "min/maxFinite"
		override possible values "min/max" if they are "False".

RETURN VALUE:
	None.

SIDE EFFECTS:
	x

------------------------------------------------------------------------------*/

void MPS_LP::GetConstraintRange( Int_T row, Real_T &bl, Real_T &bu )
	const
{
	assert( row >= 0 && row < m );

	if( RowType[row] & RT_RNG )
	{
		switch( RowType[row] & RT_TYPE )
		{
		case RT_LE:	bl = b[row] - r[row];	bu = b[row];			break;
		case RT_GE:
		case RT_EQ:	bl = b[row];			bu = b[row] + r[row];	break;
		}
	}
	else
	{
		switch( RowType[row] & RT_TYPE )
		{
		case RT_LE: bl = -INFINITY;			bu = b[row];			break;
		case RT_GE: bl = b[row];			bu = +INFINITY;			break;
		case RT_EQ: bl =					bu = b[row];			break;
		case RT_FR:	bl = -INFINITY;			bu = +INFINITY;			break;
		}
	}
}


void MPS_LP::SetConstraintRange( Int_T i, Real_T min, Bool_T minFinite,
	Real_T max, Bool_T maxFinite )
{
	assert( min <= max );

	if( minFinite && maxFinite && IsEqual( min, max ) )
	{
		RowType[i]	= RT_EQ;
		b[i]		= min;
		r[i]		= 0.0;
	}
	else if( minFinite && maxFinite )
	{
		RowType[i]	= RT_LE | RT_RNG;
		b[i]		= max;
		r[i]		= max - min;
	}
	else if( minFinite )
	{
		RowType[i]	= RT_GE;
		b[i]		= min;
		r[i]		= +INFINITY;
	}
	else if( maxFinite )
	{
		RowType[i]	= RT_LE;
		b[i]		= max;
		r[i]		= +INFINITY;
	}
	else
	{
		RowType[i]	= RT_FR;
		b[i]		= 0.0;
		r[i]		= +INFINITY;
	}
}


/*------------------------------------------------------------------------------

	Bool_T MPS_LP::SetMatrixElement( Int_T row, Int_T col, Real_T val )

PURPOSE:
	Sets the value of a matrix element (if it already exists in the matrix). A
new element will not be added. The matrix sparsity structure shall remain
unchanged.

PARAMETERS:
	Int_T row, Int_T col, Real_T val
		Row, column and the value to be set.

RETURN VALUE:
	Success status. The call succeeds only if an entry at the position (row,col)
already exists.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T MPS_LP::SetMatrixElement( Int_T row, Int_T col, Real_T val )
{
	assert( row >= 0 && row < m );
	assert( col >= 0 && col < n );
	assert( val > -INFINITY && val < INFINITY );

	for( Int_T j = ColStart[col], je = ColStart[col+1]; j < je; j++ )
		if( Row[j] == row )
		{
			ac[j] = val;
			return True;
		}

	return False;
}


/*------------------------------------------------------------------------------

	void MPS_LP::CheckLP( void )

PURPOSE:
	After the linear problem is read, it's data is checked for consistency.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void MPS_LP::CheckLP( void )
{
	//--------------------------------------------------------------------------
	//	Check RHS and range.
	//
	for( Int_T i = 0; i < m; i++ )
	{
		Short_T &rt = RowType[i] ;

		if( rt & RT_RNG )
		{
			if( r[i] <= -INFINITY || b[i] > +INFINITY )
			{
				static Bool_T found = False;

				if( !found )
				{
					Error( "RHS/range value too large." );
					found = True;
				}
			}
		}
		else if( rt & ( RT_LE | RT_GE | RT_EQ ) )
		{
			r[i] = ( rt & RT_EQ ) ? 0.0 : +INFINITY;

			if( b[i] >= +INFINITY )
			{
				static Bool_T found = False;

				if( !found )
				{
					Error( "RHS value too large." );
					found = True;
				}
			}
		}
		else
		{
			r[i] = -INFINITY;
			b[i] = +INFINITY;
		}
	}

	//--------------------------------------------------------------------------
	//	Check variable types.
	//
	for( Int_T j = 0; j < n; j++ )
	{
		Short_T &vt = VarType[j];

		if( ( vt & VT_LO ) && l[j] <= -INFINITY )
		{
			static Bool_T found = False;

			vt &= ~(VT_LO|VT_FX);
			l[j] = -INFINITY;

			if( !found )
			{
				Warning( "Infinite lower bound(s). Variable type adjusted." );
				found = True;
			}
		}

		if( ( vt & VT_UP ) && u[j] >= +INFINITY )
		{
			static Bool_T found = False;

			vt &= ~(VT_UP|VT_FX);
			u[j] = +INFINITY;

			if( !found )
			{
				Warning( "Infinite upper bound(s). Variable type adjusted." );
				found = True;
			}
		}
	}

	//--------------------------------------------------------------------------
	//	Check the constraint matrix.
	//
	for( Int_T k = 0; k < nz; k++ )
	{
		if( fabs( ac[k] ) >= +INFINITY )
		{
			static Bool_T found = False;

			if( !found )
			{
				Error( "Infinite constraint matrix coefficient." );
				found = True;
			}
		}
	}
}


/*------------------------------------------------------------------------------

	Array<Int_T> MPS_LP::GetRowLen( void )

PURPOSE:
	Returns an array of row lengths. Needed during writing an LP-DIT output
file.

PARAMETERS:
	None.

RETURN VALUE:
	The array.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Array<Int_T> MPS_LP::GetRowLen( void )
{
	Array<Int_T> RowLen( m, 0.0 );

	Ptr<Real_T> a;
	Ptr<Int_T> row;
	Int_T len;

	for( Int_T j = 0; j < n; j++ )
		for( GetColumn( j, a, row, len ); len; --len, ++row )
			RowLen[*row]++;

	return RowLen;
}
