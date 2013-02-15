/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic data parser and scenario generator.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	changelp.cpp
CREATED:			1995.12.27
LAST MODIFIED:		1996.04.09

DEPENDENCIES:		stdtype.h, smartptr.h, changelp.h
					<assert.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __MAKELAB_H__
#	include "makelab.h"
#endif
#ifndef __CHANGELP_H__
#	include "changelp.h"
#endif


/*------------------------------------------------------------------------------

	virtual Bool_T ChangeCost::Apply( MPS_LP &lp, Bool_T zero ) const
	virtual Bool_T ChangeRHS::Apply( MPS_LP &lp, Bool_T zero ) const
	virtual Bool_T ChangeLoBND::Apply( MPS_LP &lp, Bool_T zero ) const
	virtual Bool_T ChangeUpBND::Apply( MPS_LP &lp, Bool_T zero ) const
	virtual Bool_T ChangeMatrix::Apply( MPS_LP &lp, Bool_T zero ) const

PURPOSE:
	Those virtual member functions apply the changes specified in their
respective objects to the linear problem "lp". If the boolean argument "zero" is
set to "True", then the corresponding value in the object is set to zero.

PARAMETERS:
	MPS_LP &lp
		The linear problem to be modified.

	Bool_T zero
		If "False", the value stored in the object will be used, otherwise a
		zero will be put in the appropriate place. In "changelp.h" a default
		value "False" for this argument is defined.

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ChangeCost::Apply( MPS_LP &lp, Bool_T zero )
const
{
	if( col < 0 )
		return False;
	return lp.SetC( col, zero ? 0.0 : value );
}


Bool_T ChangeRHS::Apply( MPS_LP &lp, Bool_T zero )
const
{
	if( row < 0 )
		return False;
	lp.SetRHS_Elem( row, zero ? 0.0 : value );
	return True;
}


Bool_T ChangeLoBND::Apply( MPS_LP &lp, Bool_T zero )
const
{
	if( col < 0 )
		return False;
	lp.SetL( col, zero ? 0.0 : value );
	return True;
}


Bool_T ChangeUpBND::Apply( MPS_LP &lp, Bool_T zero )
const
{
	if( col < 0 )
		return False;
	lp.SetU( col, zero ? 0.0 : value );
	return True;
}


Bool_T ChangeMatrix::Apply( MPS_LP &lp, Bool_T zero )
const
{
	if( row < 0 || col < 0 )
		return False;

	lp.SetMatrixElement( row, col, zero ? 0.0 : value );

	return True;
}


/*------------------------------------------------------------------------------

	virtual void ChangeCost::Permute( const Array<Int_T> &,
		const Array<Int_T> &ColNum )
	virtual void ChangeRHS::Permute( const Array<Int_T> &RowNum,
		const Array<Int_T> & )
	virtual void ChangeLoBND::Permute( const Array<Int_T> &,
		const Array<Int_T> &ColNum )
	virtual void ChangeUpBND::Permute( const Array<Int_T> &,
		const Array<Int_T> &ColNum )
	virtual void ChangeMatrix::Permute( const Array<Int_T> &RowNum,
		const Array<Int_T> &ColNum )

PURPOSE:
	Initially the row and column numbers stored in the "ChangeLP" objects
correspond to the CORE file form of the linear problem. If however the problem
should be split into periods, the indice need to be changed (not exactly
permuted, but we shall continue to use this term).
	The above listed virtual methods perform the appropriate permutations given
"permutation" tables.

PARAMETERS:
	const Array<Int_T> &RowNum, const Array<Int_T> &ColNum
		Index permutation tables for row and column indice, respectively.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void ChangeCost::Permute( const Array<Int_T> &, const Array<Int_T> &ColNum )
{
	col = ColNum[ col ];
	assert( col >= 0 );
}


void ChangeRHS::Permute( const Array<Int_T> &RowNum, const Array<Int_T> & )
{
	row = RowNum[ row ];
	assert( row >= 0 );
}


void ChangeLoBND::Permute( const Array<Int_T> &, const Array<Int_T> &ColNum )
{
	col = ColNum[ col ];
	assert( col >= 0 );
}


void ChangeUpBND::Permute( const Array<Int_T> &, const Array<Int_T> &ColNum )
{
	col = ColNum[ col ];
	assert( col >= 0 );
}


void ChangeMatrix::Permute( const Array<Int_T> &RowNum, // )
	const Array<Int_T> &ColNum )
{
	col = ColNum[ col ];
	row = RowNum[ row ];
	assert( col >= 0 );
	assert( row >= 0 );
}


/*------------------------------------------------------------------------------

	virtual Bool_T ChangeCost::Write( FILE *fp, Int_T sc, Real_T prob ) const
	virtual Bool_T ChangeRHS::Write( FILE *fp, Int_T sc, Real_T ) const
	virtual Bool_T ChangeLoBND::Write( FILE *fp, Int_T sc, Real_T ) const
	virtual Bool_T ChangeUpBND::Write( FILE *fp, Int_T sc, Real_T ) const
	virtual Bool_T ChangeMatrix::Write( FILE *fp, Int_T sc, Real_T ) const

PURPOSE:
	Write a stochastic data file entry

PARAMETERS:
	FILE *fp, Int_T sc, Real_T prob

RETURN VALUE:
	Boolean success status.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

Bool_T ChangeCost::Write( FILE *fp, Int_T sc, Real_T prob )
const
{
	assert( fp != NULL );

	if( fprintf( fp, "    %-8s  %-8s  %12G\n",
		MakeLabel( col, sc ), "COST", prob * value ) == EOF )
		return False;

	return True;
}


Bool_T ChangeRHS::Write( FILE *fp, Int_T sc, Real_T )
const
{
	assert( fp != NULL );

	if( fprintf( fp, "    %-8s  %-8s  %12G\n",
		"RHS", MakeLabel( row, sc ), value ) == EOF )
		return False;

	return True;
}


Bool_T ChangeLoBND::Write( FILE *, Int_T, Real_T )
const
{
	abort();
	return True;
}


Bool_T ChangeUpBND::Write( FILE *, Int_T, Real_T )
const
{
	abort();
	return True;
}


Bool_T ChangeMatrix::Write( FILE *, Int_T, Real_T )
const
{
	abort();
	return True;
}
