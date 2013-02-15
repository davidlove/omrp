/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming - general purpose.
PROJECT CODE:		---------------
PROJECT FULL NAME:	---------------

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	---------------

--------------------------------------------------------------------------------

SOURCE FILE NAME:	io_dit.cpp
CREATED:			1995.02.01
LAST MODIFIED:		1996.02.03

DEPENDENCIES:		stdtype.h, mps_lp.h, myalloc.h, sort_lab.h, work_vec.h,
					print.h
					<math.h>

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	LP-DIT interface for linear problem input and output is provided here.

--------------------------------------------------------------------------------

PUBLIC INTERFACE:
	x

PRIVATE FUNCTIONS:
	x

STATIC FUNCTIONS:
	x

STATIC DATA:
	x

--------------------------------------------------------------------------------

USED MACROS AND THEIR MEANING:
	x

------------------------------------------------------------------------------*/

#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

#ifdef SUPPORT_LP_DIT

#include <math.h>

#ifndef __SORT_LAB_H__
#	include "sort_lab.h"
#endif
#ifndef __LP_CODES_H__
#	include "lp_codes.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif


#include "lp_dit.h"


static Bool_T ReadLabels	= True,
	Minimize				= True;


Bool_T MPS_LP::ReadLP_DIT( const char *FileName, VerbLevel Verbosity )
{
	if( Verbosity >= V_LOW )
		Print( "Reading an LP-DIT format linear problem from file %s.\n",
			FileName );
	//
	//	Forget all previous data.
	//
	FreeStorage();
	mE = mG = mL = mR = mF = nPL = nFX = nFR = nMI = nUP = 0;

	//
	//	Read in the problem's dimensions.
	//
	lp_init( 0, 1, 0, 0, (char *) FileName, (void *) this );

	//
	//	Allocate memory.
	//
	ac.Resize( nz );
	Row.Resize( nz );
	ColStart.Resize( n + 1 );
	b.Resize( m );			b.Fill( 0.0,					m );
	r.Resize( m );			r.Fill( INFINITY,				m );
	l.Resize( n );			l.Fill( 0.0,					n );
	u.Resize( n );			u.Fill( INFINITY,				n );
	VarType.Resize( n );	VarType.Fill( VTM_UNDEFINED,	n );
	RowType.Resize( m );	RowType.Fill( RT_UNDEFINED,		m );

	ColStart[0] = 0;

	//
	//	Read in the LP.
	//
	lp_def( 0, 0, (void *) this );
	OS = OS_FULL;

	assert( ColStart[n] == nz );

	//
	//	Some more things to be done before the object is ready.
	//
	if( !ReadLabels )
	{
		RowLabels.FreeMemory();
		ColLabels.FreeMemory();

		OS = OS_NO_LABELS;
	}

	if( !Minimize )
	{
		Warning( "LP converted to minimization." );

		for( Int_T k = 0; k < nz; k++ )
			if( Row[k] == ObjRow && IsNonZero( ac[k] ) )
				ac[k] = -ac[k];
	}

	//
	//	The final report on problem statistics.
	//
	if( Verbosity == V_LINE )
	{
		Print(
			/* FORMAT STRING */
			"| %-8s | %10d | %10d | %10ld | %10.4f |", 

			/* DATA */
			Name, (int) m, (int) n, (long int) nz,
			100.0 * double( nz ) / double( m ) / double( n )
			);
	}
	else if( Verbosity >= V_LOW )
	{
		Print(
			"\nLinear problem statistics (as read from %s).\n"
			"\t%-18s%10d (%d E, %d L, %d G, %d N)\n"
			"\t%-18s%10d\n"
			"\t%-18s%10d (%d PL, %d UP, %d FX, %d FR, %d MI)\n"
			"\t%-18s%10ld\n"
			"\t%-18s%10.4f %%\n",

			FileName,
			"No. of rows:",			(int) m, (int) mE, (int) mL, (int) mG,
									(int) mF,
			"No. of ranges:",		(int) mR,
			"No. of variables:",	(int) n, (int) nPL, (int) nUP, (int) nFX,
									(int) nFR, (int) nMI,
			"No. of non-zeros:",	(long int) nz,
			"Density:",				100.0 * (double)nz / (double)m / (double)n
			);
	}
	return True;
}


Bool_T MPS_LP::WriteLP_DIT( const char *FileName, VerbLevel Verbosity )
	const
{
	assert( OS == OS_FULL );

	if( Verbosity >= V_LOW )
		Print( "Writing an LP-DIT file %s....", FileName );

	//--------------------------------------------------------------------------
	//	Declare the problem's dimensions to LP-DIT.
	//
	lp_init( 1, 1, 0, 0, (char *) FileName, (void *) this );

	//--------------------------------------------------------------------------
	//	Store the LP.
	//
	lp_def( 1, 0, (void *) this );

	if( Verbosity >= V_LOW )
		Print( "Finished.\n" );

	return True;
}


//------------------------------------------------------------------------------
//	These are the "C" functions that will be called by LP-DIT functions
//	'lp_init' and 'lp_def'. They in turn will call appropriate MPS_LP member
//	functions. This is not exactly a very safe mechanism - a void pointer is
//	cast to (MPS_LP *). Unfortunately, we see no entirely safe way to obtain the
//	same result.
//

extern "C" {

//
//	Linear problem input functions.
//
void lpo_par( LP_HEAD h, void *lp )
	{ ( (MPS_LP *) lp )->lpo_par( h ); }

void lpo_specs( LP_HEAD h, char **str, void *lp )
	{ ( (MPS_LP *) lp )->lpo_specs( h, str ); }

void lpo_cols( LP_HEAD h, LP_VAR *v, void *lp )
	{ ( (MPS_LP *) lp )->lpo_cols( h, v ); }

void lpo_rows( LP_HEAD h, LP_VAR *v, void *lp )
	{ ( (MPS_LP *) lp )->lpo_rows( h, v ); }

void lpo_vect( LP_HEAD h, LP_MAT *v, void *lp )
	{ ( (MPS_LP *) lp )->lpo_vect( h, v ); }

//
//	Linear problem output functions.
//
void lpi_par( LP_HEAD *h, void *lp )
	{ ( (MPS_LP *) lp )->lpi_par( *h ); }

void lpi_specs( LP_HEAD h, char **str, void *lp )
	{ ( (MPS_LP *) lp )->lpi_specs( h, str ); }

void lpi_cols( LP_HEAD h, LP_VAR *v, void *lp )
	{ ( (MPS_LP *) lp )->lpi_cols( h, v ); }

void lpi_rows( LP_HEAD h, LP_VAR *v, void *lp )
	{ ( (MPS_LP *) lp )->lpi_rows( h, v ); }

void lpi_vect( LP_HEAD h, LP_MAT *v, void *lp )
	{ ( (MPS_LP *) lp )->lpi_vect( h, v ); }

} // end of extern "C"


//------------------------------------------------------------------------------
//
//	THIS IS THE BEGINNING OF THE LP-DIT INTERFACE FUNCTIONS.
//
//------------------------------------------------------------------------------


static Int_T lastCol	= 0,
	lastNZ				= 0;


void MPS_LP::lpo_par( LP_HEAD &h )
{
	ObjRow		= -1;
	Minimize	= True;
	lastCol		= 0;
	lastNZ		= 0;

	if( h.status != LP_INI )
		FatalError( "No initial problem definition from LP-DIT." );

	m	= (Int_T) h.m;
	n	= (Int_T) h.n;
	nz	= (Int_T) h.nz;

	assert( m > 0 && n > 0 && nz > 0 );
	assert( h.min_max == MAXIMIZE || h.min_max == MINIMIZE );

	if( h.min_max == MAXIMIZE )
		Minimize = False;

	ObjRow = (Int_T)h.obj;

	strncpy( Name, h.name, LAB_LEN );
	Name[LAB_LEN] = '\0';
}

void MPS_LP::lpo_specs( LP_HEAD &, char ** )
{
// Do nothing - ignore specifications.
}


void MPS_LP::lpo_cols( LP_HEAD &, LP_VAR *v )
{
	for( Int_T j = 0; j < n; j++ )
	{
		LP_VAR &vv = v[j];

		if( *vv.name == '\0' )
			ReadLabels = False;

		if( ReadLabels )
			ColLabels.AddLabel( vv.name );

		VarType[j] = 0;

		if( vv.is_low_bnd )
		{
			l[j] = vv.low_bnd;
			VarType[j] |= VT_LO;
		}
		else
			l[j] = -INFINITY;

		if( vv.is_upp_bnd )
		{
			u[j] = vv.upp_bnd;
			VarType[j] |= VT_UP;
		}
		else
			u[j] = INFINITY;

		if( IsEqual( l[j], u[j] ) )
			VarType[j] = VT_FIXED;

		if( VarType[j] & VT_FX )						nFX++;
		else if( vv.is_upp_bnd && vv.is_low_bnd )		nUP++;
		else if( vv.is_low_bnd )						nPL++;
		else if( vv.is_upp_bnd )						nMI++;
		else											nFR++;

		assert( l[j] <= u[j] );
	}
}


void MPS_LP::lpo_rows( LP_HEAD &, LP_VAR *v )
{
	for( Int_T i = 0; i < m; i++ )
	{
		LP_VAR &vv = v[i];
		Short_T &rt = RowType[i];

		if( *vv.name == '\0' )
			ReadLabels = False;

		if( ReadLabels )
			RowLabels.AddLabel( vv.name );

		if( vv.is_eq )		{ rt = RT_EQ;	mE++; }
		else if( vv.is_le )	{ rt = RT_LE;	mL++; }
		else if( vv.is_ge )	{ rt = RT_GE;	mG++; }
		else if( vv.is_ne )	{ rt = RT_FR;	mF++; }
		else
			FatalError( "io_dit.cpp: lpo_rows: Invalid row type." );
		
		switch( rt )
		{
		case RT_EQ:
			assert( IsZero( vv.low_bnd - vv.upp_bnd) );

			r[i] = INFINITY;
			b[i] = vv.upp_bnd;
			break;

		case RT_LE:
			assert( vv.upp_bnd > -INFINITY );
			assert( vv.upp_bnd < +INFINITY );

			b[i] = vv.upp_bnd;
			if( vv.is_low_bnd )
			{
				r[i] = fabs( vv.upp_bnd - vv.low_bnd );
				rt |= RT_RNG;
				mR++;
			}
			else
				r[i] = INFINITY;
			break;

		case RT_GE:
			assert( vv.low_bnd > -INFINITY );
			assert( vv.low_bnd < +INFINITY );

			b[i] = vv.low_bnd;
			if( vv.is_upp_bnd )
			{
				r[i] = fabs( vv.upp_bnd - vv.low_bnd );
				rt |= RT_RNG;
				mR++;
			}
			else
				r[i] = INFINITY;
			break;

		case RT_FR:
			b[i] = 0.0;
			r[i] = INFINITY;
			break;
		}
	}

#ifndef NDEBUG
	{ for( Int_T i = 0; i < m; i++ )
		assert( fabs( b[i] ) < INFINITY ); }
#endif
}


void MPS_LP::lpo_vect( LP_HEAD &, LP_MAT *v )
{
	assert( v->index == lastCol );

	lastCol++;
	lastNZ += v->elems;

	assert( lastNZ <= nz );

	Int_T &k	= ColStart[ v->index + 1 ];

	k = ColStart[ v->index ];

	for( Int_T i = 0; i < v->elems; i++ )
	{
		LP_VECT *ve = v->el + i;

		ac[k]	= ve->value;
		Row[k]	= ve->index;

		if( IsNonZero( ac[k] ) )
			k++;
	}
}



void MPS_LP::lpi_par( LP_HEAD &h )
{
	assert( OS == OS_FULL );

	strcpy( h.name, Name );

	h.min_max	= ( !Minimize ) ? MAXIMIZE : MINIMIZE;

	h.m			= m;
	h.obj		= 0;		// To be filled in later in this function...
	h.n			= n;
	h.nint		= 0;
	h.nz		= (long) nz;

	h.specs		= 0;
	h.status	= LP_INI;

	h.feas		= 0.0;
	h.optim		= 0.0;
	h.infty		= +INFINITY;

	//--------------------------------------------------------------------------
	//	Now back to h.obj: we scan the rows for the first free row.
	//
	for( Int_T i = 0; i < m; i++ )
		if( RowType[i] == RT_FR )
		{
			h.obj = i;
			break;
		}
}


void SolvableLP::lpi_par( LP_HEAD &h )
{
	assert( OS == OS_FULL );

	MPS_LP::lpi_par( h );
	h.m++;							// Count the objective row.

	Int_T cnt = 0;
	for( Int_T j = 0; j < n; j++ )
		if( IsNonZero( c[j] ) ) cnt++;

	h.nz += cnt;					// Count the objective row.
	h.obj = m;						// Put the objective at the end.
}


void MPS_LP::lpi_specs( LP_HEAD &, char ** )
{
	assert( OS == OS_FULL );

// Do nothing - ignore specifications.
}


void MPS_LP::lpi_cols( LP_HEAD &, LP_VAR *v )
{
	for( Int_T j = 0; j < n; j++ )
	{
		LP_VAR &vv = v[j];

		strcpy( vv.name, ColLabels.FindLabel( j ) );

		vv.elems		= (LP_IND)( ColStart[j+1] - ColStart[j] );

		vv.is_eq		= 0;
		vv.is_le		= 0;
		vv.is_ge		= 0;
		vv.is_ne		= 0;
		vv.is_low_bnd	= 0;
		vv.is_upp_bnd	= 0;
		vv.mip_type		= 0;
		vv.attr			= 0;

		if( VarType[j] & VT_FX )
		{
			vv.low_bnd		= vv.upp_bnd = u[j];
			vv.is_eq		= 1;
			vv.is_low_bnd	= 1;
			vv.is_upp_bnd	= 1;
		}
		else if( ( VarType[j] & VT_LO ) && ( VarType[j] & VT_UP ) )
		{
			vv.low_bnd		= l[j];
			vv.upp_bnd		= u[j];
			vv.is_low_bnd	= 1;
			vv.is_upp_bnd	= 1;
		}
		else if( VarType[j] & VT_LO )
		{
			vv.low_bnd		= l[j];
			vv.upp_bnd		= +INFINITY;
			vv.is_ge		= 1;
			vv.is_low_bnd	= 1;
		}
		else if( VarType[j] & VT_UP )
		{
			vv.low_bnd		= -INFINITY;
			vv.upp_bnd		= u[j];
			vv.is_le		= 1;
			vv.is_upp_bnd	= 1;
		}
		else
		{
			vv.low_bnd	= -INFINITY;
			vv.upp_bnd	= +INFINITY;
			vv.is_ne	= 1;
		}
	}
}


void SolvableLP::lpi_cols( LP_HEAD &h, LP_VAR *v )
{
	MPS_LP::lpi_cols( h, v );

	for( Int_T j = 0; j < n; j++ )
		if( IsNonZero( c[j] ) )
			v[j].elems++;
}


void MPS_LP::lpi_rows( LP_HEAD &, LP_VAR *v )
{
	for( Int_T i = 0; i < m; i++ )
	{
		LP_VAR &vv = v[i];

		strcpy( vv.name, RowLabels.FindLabel( i ) );

		vv.elems		= 0;
		vv.is_eq		= 0;
		vv.is_le		= 0;
		vv.is_ge		= 0;
		vv.is_ne		= 0;
		vv.mip_type		= 0;
		vv.attr			= 0;

		Short_T rt		= Short_T( RowType[i] & RT_TYPE );

		if( rt & RT_RNG )
			vv.is_ge = 1;
		else if( rt == RT_EQ )
			vv.is_eq = 1;
		else if( rt == RT_LE )
			vv.is_le = 1;
		else if( rt == RT_GE )
			vv.is_ge = 1;
		else
		{
			assert( rt == RT_FR );
			vv.is_ne = 1;
		}

		Real_T ll, uu;

		GetConstraintRange( i, ll, uu );

		vv.low_bnd		= ll;
		vv.upp_bnd		= uu;
		vv.is_low_bnd	= ( ll > -INFINITY ) ? 1 : 0;
		vv.is_upp_bnd	= ( uu < +INFINITY ) ? 1 : 0;

		assert( !vv.is_ge || vv.is_low_bnd );
		assert( !vv.is_eq || vv.is_low_bnd );
		assert( !vv.is_eq || vv.is_upp_bnd );
		assert( !vv.is_le || vv.is_upp_bnd );
		assert( !vv.is_ne || ( !vv.is_low_bnd && !vv.is_upp_bnd ) );
	}

	for( Int_T j = 0; j < nz; j++ )
	{
		assert( Row[j] >= 0 && Row[j] < m );
		v[ Row[j] ].elems++;
	}
}


void SolvableLP::lpi_rows( LP_HEAD &h, LP_VAR *v )
{
	MPS_LP::lpi_rows( h, v );

	LP_VAR &vv	= v[m];

	vv.elems	= 0;
	vv.is_eq	= 0;
	vv.is_le	= 0;
	vv.is_ge	= 0;
	vv.is_ne	= 1;
	vv.mip_type	= 0;
	vv.attr		= 0;

	vv.low_bnd		= -INFINITY;
	vv.upp_bnd		= +INFINITY;
	vv.is_low_bnd	= 0;
	vv.is_upp_bnd	= 0;

	strcpy( vv.name, Obj );

	for( Int_T j = 0; j < n; j++ )
		if( IsNonZero( c[j] ) ) vv.elems++;
}


void MPS_LP::lpi_vect( LP_HEAD &, LP_MAT *v )
{
	assert( v->index < n );
	assert( v->el != NULL );

	Int_T j	= v->index,
		k	= ColStart[j],
		len = Int_T( ColStart[j+1] - k );

	v->elems = (LP_IND)len;

	for( Int_T i = 0; i < len; i++, k++ )
	{
		LP_VECT &ve	= v->el[i];

		ve.index	= (LP_IND)Row[k];
		ve.value	= ac[k];
	}
}


void SolvableLP::lpi_vect( LP_HEAD &h, LP_MAT *v )
{
	MPS_LP::lpi_vect( h, v );

	if( IsNonZero( c[v->index] ) )
	{
		v->el[ v->elems ].index = (LP_IND)m;
		v->el[ v->elems ].value = c[v->index];
		v->elems++;
	}
}


//	End of functions called by LP-dit.
//------------------------------------------------------------------------------

#endif
