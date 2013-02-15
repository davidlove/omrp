/*------------------------------------------------------------------------------
MODULE TYPE:		Project core code.
PROJECT CODE:		REGULARIXED DECOMPOSITION
PROJECT FULL NAME:	Implementation of the regularized decomposition for two
					stage linear programs.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. Andrzej Ruszczynski

--------------------------------------------------------------------------------

HEADER FILE NAME:	changelp.h
CREATED:			1995.12.27
LAST MODIFIED:		1996.04.09

DEPENDENCIES:		stdtype.h, mps_lp.h
					<stdio.h>, <assert.h>

--------------------------------------------------------------------------------

HEADER CONTENTS:
	x

------------------------------------------------------------------------------*/

#ifndef __CHANGELP_H__
#define __CHANGELP_H__

#include <stdio.h>
#include <assert.h>

#ifndef __MY_DEFS_H__
#	include "my_defs.h"
#endif
#ifndef __MPS_LP_H__
#	include "mps_lp.h"
#endif

//==============================================================================
//	Type definitions for scenario list elements. A two level hierarchy of
//	objects is created with "changeLP" as a common abstract base class.
//

abstract class ChangeLP
{
protected:
	Real_T value;

public:
	enum Type { NONE, RHS, COST, LO_BND, UP_BND, MATRIX };

protected:
	ChangeLP( Real_T v );

public:
	ChangeLP( void );
	virtual ~ChangeLP( void );

	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const pure;
	virtual Type GetType( void ) const pure;
	virtual Int_T GetRow( void ) const pure;
	virtual Int_T GetCol( void ) const pure;
	virtual Real_T GetValue( void ) const;

	virtual void Permute( const Array<Int_T> &RowNum,
		const Array<Int_T> &ColNum ) pure;

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T prob ) const pure;
};


//------------------------------------------------------------------------------

class ChangeCost : public ChangeLP
{
protected:
	Int_T col;

public:
	ChangeCost( void );
	ChangeCost( Int_T j, Real_T v );
	virtual ~ChangeCost( void );
	
	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const;
	virtual Type GetType( void ) const;
	virtual Int_T GetRow( void ) const;
	virtual Int_T GetCol( void ) const;

	virtual void Permute( const Array<Int_T> &, const Array<Int_T> &ColNum );

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T prob ) const;
};


//------------------------------------------------------------------------------

class ChangeRHS : public ChangeLP
{
protected:
	Int_T row;

public:
	ChangeRHS( void );
	ChangeRHS( Int_T i, Real_T v );
	virtual ~ChangeRHS( void );
	
	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const;
	virtual Type GetType( void ) const;
	virtual Int_T GetRow( void ) const;
	virtual Int_T GetCol( void ) const;

	virtual void Permute( const Array<Int_T> &RowNum, const Array<Int_T> & );

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T ) const;
};


//------------------------------------------------------------------------------

class ChangeLoBND : public ChangeLP
{
protected:
	Int_T col;

public:
	ChangeLoBND( void );
	ChangeLoBND( Int_T j, Real_T v );
	virtual ~ChangeLoBND( void );
	
	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const;
	virtual Type GetType( void ) const;
	virtual Int_T GetRow( void ) const;
	virtual Int_T GetCol( void ) const;

	virtual void Permute( const Array<Int_T> &, const Array<Int_T> &ColNum );

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T ) const;
};


//------------------------------------------------------------------------------

class ChangeUpBND : public ChangeLP
{
protected:
	Int_T col;

public:
	ChangeUpBND( void );
	ChangeUpBND( Int_T j, Real_T v );
	virtual ~ChangeUpBND( void );
	
	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const;
	virtual Type GetType( void ) const;
	virtual Int_T GetRow( void ) const;
	virtual Int_T GetCol( void ) const;

	virtual void Permute( const Array<Int_T> &, const Array<Int_T> &ColNum );

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T ) const;
};


//------------------------------------------------------------------------------

class ChangeMatrix : public ChangeLP
{
protected:
	Int_T row, col;

public:
	ChangeMatrix( void );
	ChangeMatrix( Int_T i, Int_T j, Real_T v );
	virtual ~ChangeMatrix( void );
	
	virtual Bool_T Apply( MPS_LP &lp, Bool_T zero = False ) const;
	virtual Type GetType( void ) const;
	virtual Int_T GetRow( void ) const;
	virtual Int_T GetCol( void ) const;

	virtual void Permute( const Array<Int_T> &RowNum,
		const Array<Int_T> &ColNum );

	virtual Bool_T Write( FILE *fp, Int_T sc, Real_T ) const;
};

//
//	End of type definitions for scenario list elements.
//==============================================================================


//==============================================================================
//
//	Inline functions' definitions.
//

/*------------------------------------------------------------------------------

	ChangeLP::ChangeLP( void )
	ChangeLP::ChangeLP( Real_T v )
	virtual ChangeLP::~ChangeLP( void )

PURPOSE:
	Two constructors and a virtual destructor of an abstract base class.

PARAMETERS:
	Real_T v
		A floating point value. Its meaning depends on the type of the object
which inherits the "ChangeLP" object.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeLP::ChangeLP( void )
	: value( 0.0 )
{}


inline
ChangeLP::ChangeLP( Real_T v )
	: value( v )
{}


inline
ChangeLP::~ChangeLP( void )
{}


/*------------------------------------------------------------------------------

	Real_T ChangeLP::GetValue( void ) const

PURPOSE:
	Data access routine. Returns the floating point value stored in the object.

PARAMETERS:
	None.

RETURN VALUE:
	The value.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
Real_T ChangeLP::GetValue( void )
const
{ return value; }



/*------------------------------------------------------------------------------

	ChangeCost::ChangeCost( void )
	ChangeCost::ChangeCost( Int_T j, Real_T v )
	virtual ChangeCost::~ChangeCost( void )

PURPOSE:
	Class "ChangeCost" constructors.

PARAMETERS:
	Int_T j, Real_T v
		Column number and the cost coefficient value.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeCost::ChangeCost( void )
	: ChangeLP(), col( -1 )
{}


inline
ChangeCost::ChangeCost( Int_T j, Real_T v )
	: ChangeLP( v ), col( j )
{ assert( j >= 0 ); }


inline
ChangeCost::~ChangeCost( void )
{}

/*------------------------------------------------------------------------------

	ChangeRHS::ChangeRHS( void )
	ChangeRHS::ChangeRHS( Int_T i, Real_T v )
	virtual ChangeRHS::~ChangeRHS( void )

PURPOSE:
	Class "ChangeRHS" constructors.

PARAMETERS:
	Int_T i, Real_T v
		Row number and the random value.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeRHS::ChangeRHS( void )
	: ChangeLP(), row( -1 )
{}


inline
ChangeRHS::ChangeRHS( Int_T i, Real_T v )
	: ChangeLP( v ), row( i )
{ assert( i >= 0 ); }


inline
ChangeRHS::~ChangeRHS( void )
{}


/*------------------------------------------------------------------------------

	ChangeLoBND::ChangeLoBND( void )
	ChangeLoBND::ChangeLoBND( Int_T j, Real_T v )
	virtual ChangeLoBND::~ChangeLoBND( void )

PURPOSE:
	Class "ChangeLoBND" constructors.

PARAMETERS:
	Int_T j, Real_T v
		Column number and the random value.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeLoBND::ChangeLoBND( void )
	: ChangeLP(), col( -1 )
{}


inline
ChangeLoBND::ChangeLoBND( Int_T j, Real_T v )
	: ChangeLP( v ), col( j )
{ assert( j >= 0 ); }


inline
ChangeLoBND::~ChangeLoBND( void )
{}


/*------------------------------------------------------------------------------

	ChangeUpBND::ChangeUpBND( void )
	ChangeUpBND::ChangeUpBND( Int_T j, Real_T v )
	virtual ChangeUpBND::~ChangeUpBND( void )

PURPOSE:
	Class "ChangeUpBND" constructors.

PARAMETERS:
	Int_T j, Real_T v
		Column number and the random value.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeUpBND::ChangeUpBND( void )
	: ChangeLP(), col( -1 )
{}


inline
ChangeUpBND::ChangeUpBND( Int_T j, Real_T v )
	: ChangeLP( v ), col( j )
{ assert( j >= 0 ); }


inline
ChangeUpBND::~ChangeUpBND( void )
{}


/*------------------------------------------------------------------------------

	ChangeMatrix::ChangeMatrix( void )
	ChangeMatrix::ChangeMatrix( Int_T i, Int_T j, Real_T v )
	virtual ChangeMatrix::~ChangeMatrix( void )

PURPOSE:
	Class "ChangeMatrix" constructors.

PARAMETERS:
	Int_T i, Int_T j, Real_T v
		Row and column number and the random value.

RETURN VALUE:
	Not applicable.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeMatrix::ChangeMatrix( void )
	: ChangeLP(), row( -1 ), col( -1 )
{}


inline
ChangeMatrix::ChangeMatrix( Int_T i, Int_T j, Real_T v )
	: ChangeLP( v ), row( i ), col( j )
{ assert( i >= 0 ); assert( j >= 0 ); }


inline
ChangeMatrix::~ChangeMatrix( void )
{}


/*------------------------------------------------------------------------------

	virtual ChangeLP::Type ChangeCost::GetType( void ) const
	virtual Int_T ChangeCost::GetRow( void ) const
	virtual Int_T ChangeCost::GetCol( void ) const

	virtual ChangeLP::Type ChangeRHS::GetType( void ) const
	virtual Int_T ChangeRHS::GetRow( void ) const
	virtual Int_T ChangeRHS::GetCol( void ) const

	virtual ChangeLP::Type ChangeLoBND::GetType( void ) const
	virtual Int_T ChangeLoBND::GetRow( void ) const
	virtual Int_T ChangeLoBND::GetCol( void ) const

	virtual ChangeLP::Type ChangeUpBND::GetType( void ) const
	virtual Int_T ChangeUpBND::GetRow( void ) const
	virtual Int_T ChangeUpBND::GetCol( void ) const

	virtual ChangeLP::Type ChangeMatrix::GetType( void ) const
	virtual Int_T ChangeMatrix::GetRow( void ) const
	virtual Int_T ChangeMatrix::GetCol( void ) const

PURPOSE:
	Virtual data access routines for the whole family of "ChangeLP" objects.
"GetType" methods identify the type of object. "GetRow" and "GetColumn" methods
return row and column numbers, resp. of the data that is to be changed by the
"Apply" method. If either row or column number is not meaningful for a paticular
type of an "ChangeLP" object, then (-1) is returned.


PARAMETERS:
	None.

RETURN VALUE:
	"GetType":		one of enum "ChangeLP::Type":
					NONE, RHS, COST, LO_BND, UP_BND, MATRIX
	"GetRow":		row number, or (-1).
	"GetColumn":	column number, or (-1).

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

inline
ChangeLP::Type ChangeCost::GetType( void )
const
{ return ChangeLP::COST; }


inline
Int_T ChangeCost::GetRow( void )
const
{ return -1; }


inline
Int_T ChangeCost::GetCol( void )
const
{ return col; }


//------------------------------------------------------------------------------

inline
ChangeLP::Type ChangeRHS::GetType( void )
const
{ return ChangeLP::RHS; }


inline
Int_T ChangeRHS::GetRow( void )
const
{ return row; }


inline
Int_T ChangeRHS::GetCol( void )
const
{ return -1; }


//------------------------------------------------------------------------------

inline
ChangeLP::Type ChangeLoBND::GetType( void )
const
{ return ChangeLP::LO_BND; }


inline
Int_T ChangeLoBND::GetRow( void )
const
{ return -1; }


inline
Int_T ChangeLoBND::GetCol( void )
const
{ return col; }


//------------------------------------------------------------------------------

inline
ChangeLP::Type ChangeUpBND::GetType( void )
const
{ return ChangeLP::UP_BND; }


inline
Int_T ChangeUpBND::GetRow( void )
const
{ return -1; }


inline
Int_T ChangeUpBND::GetCol( void )
const
{ return col; }


//------------------------------------------------------------------------------

inline
ChangeLP::Type ChangeMatrix::GetType( void )
const
{ return ChangeLP::MATRIX; }


inline
Int_T ChangeMatrix::GetRow( void )
const
{ return row; }


inline
Int_T ChangeMatrix::GetCol( void )
const
{ return col; }


//
//	End of inline functions' definitions.
//
//==============================================================================

#endif
