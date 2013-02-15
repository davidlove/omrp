/*------------------------------------------------------------------------------
MODULE TYPE:		Stochastic programming core code
PROJECT CODE:		Regularized decomposition
PROJECT FULL NAME:	Implementation of the regularized decomposition of
					Ruszczynski with primal simplex method for subproblem
					solution.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	prof. A. Ruszczynski

--------------------------------------------------------------------------------

SOURCE FILE NAME:	mast_sol.cpp
CREATED:			1994.07.21
LAST MODIFIED:		1996.06.21

DEPENDENCIES:		mast_sol.h, std_math.h
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

#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __MAST_SOL_H__
#	include "mast_sol.h"
#endif
#ifndef __QDX_PUB_H__ 
#	include "qdx_pub.h" 
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif


/*------------------------------------------------------------------------------

	MasterSolver::MasterSolver( const SolvableLP &lp, const Scenarios &Scen )

PURPOSE:
	Allocates memory for the master problem solver of the regularized
decomposition method. Initially checks all problem dimension specified in
the argument list.

PARAMETERS:
	SolvableLP &lp
		First stage constraint matrix as file of packed rows.

	Int_T blocks
		Number of blocks.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

MasterSolver::MasterSolver( const SolvableLP &lp, const Scenarios &Scen )
	: n( lp.GetN() ), l( Int_T( Scen.GetNumberOfScenarios() + 1 ) ),
	mdmat( lp.GetM() ),
	itmax( 1000 ), istop( 0 ), levprt( 3 ), initpen( PENALTY_DEF ),
	dmat( (Real_T *) lp.GetNonZerosByRows() ),
	jcol( (Int_T *) lp.GetColumnNumbers() ),
	iptr( (Int_T *) lp.GetRowStarts() ),
	bmin( NULL ), bmax( NULL ),
	x( NULL ), xmax( NULL ), xmin( NULL ), y( NULL ), yb( NULL ),
	g( n*(n+2*l) ), a( n+2*l ), v( l ), weight( l ),
	q( NULL ), r( NULL ), z( NULL ), w( NULL ),
	pricnb( NULL ), pricba( l ), pi( NULL ), col( NULL ), dpb( l ),
	iblock( n+2*l ), ibasic( l ), icheck( n+2*l ), ieq( n+2*l ), drow( n+2*l), 
	irn( NULL ), istat( NULL ),
	marks( NULL ), status( NULL ),
	Objective( 0.0 ), SubMan( NULL ),
	ExpC( l ), ExpCost (0.0), ExpC2( l ), ExpCost2 (0.0)  
{
	//--------------------------------------------------------------------------
	//	Check the dimensions of the problem.
	//
	assert( n >= 0 );
	assert( mdmat >= 0 );
	assert( l > 1 );

	//--------------------------------------------------------------------------
	//	Recompute the 'jcol' and 'iptr' to match the FORTRAN indexing scheme.
	//
	Int_T i;
	for( i = 0; i < mdmat; i++ )
	{
		Int_T rs	= iptr[i],
			re		= iptr[i+1];

		for( Int_T j = rs; j < re; j++ )
			jcol[j]++;

		iptr[i]++;
	}
	iptr[i]++;
	
	//--------------------------------------------------------------------------
	//	Allocate memory for the master and the coordinator.
	//

	if( mdmat )
	{
		bmin	= new Real_T[mdmat];	Fill( bmin,		mdmat,	0.0 );
		bmax	= new Real_T[mdmat];	Fill( bmax,		mdmat,	0.0 );
		marks	= new Int_T[mdmat];		Fill( marks,	mdmat,	(Int_T) 0 );
		status	= new Int_T[mdmat];		Fill( status,	mdmat,	(Int_T) 0 );

		if( !bmin || !bmax || !marks || !status)
			FatalError( "Not enough memory to initialize the master solver." );
	}

	if( n )
	{
		x		= new Real_T[n];	Fill( x,		n,			0.0 );
		xmin	= new Real_T[n];
		xmax	= new Real_T[n];
		y		= new Real_T[n];	Fill( y,		n,			0.0 );
		yb		= new Real_T[n];	Fill( yb,		n,			0.0 );
				
		q		= new Real_T[n];	Fill( q,		n,			0.0 );
		r		= new Real_T[n*(n+1)/2];
									Fill( r,		n*(n+1)/2,	0.0 );
		z		= new Real_T[n];	Fill( z,		n,			0.0 );
		w		= new Real_T[n];	Fill( w,		n,			0.0 );
		col		= new Real_T[n];	Fill( col,		n,			0.0 );
		irn		= new Int_T[n];		Fill( irn,		n,			(Int_T) 0 );
		istat	= new Int_T[n];		Fill( istat,	n,			(Int_T) 0 );

		if( !x || !xmin || !xmax || !y || !yb || !q || !r || !z || !w ||
			!col || !irn || !istat )
			FatalError( "Not enough memory to initialize the master solver." );
	}
	
	pricnb	= new Real_T[n+1];		Fill( pricnb,	n+1,	0.0 );
	pi		= new Real_T[n+1];		Fill( pi,		n+1,	0.0 );	
	inonba	= new Int_T[n+1];		Fill( inonba,	n+1,	(Int_T) 0 );
	

/*	if( !a || !v || !weight || !pricnb || !pricba || !pi || !dpb ||
		!iblock || !ibasic || !inonba || !icheck ||!ieq ||!drow )
		FatalError( "Not enough memory to initialize the master solver." );
*/


	if( !pricnb || !pi || !inonba )
		FatalError( "Not enough memory to initialize the master solver." );


	//--------------------------------------------------------------------------
	//	Fill the 'bmin', 'bmax', 'xmin' and 'xmax' tables.
	//	Fill 'weight' vector with:
	//		(1)		1.0 for the first stage objective function and
	//		(2)		1.0/(number of scenarios) for each scenario.
	//
	for( i = 0; i < mdmat; i++ )
		lp.GetConstraintRange( i, bmin[i], bmax[i] );

	for( i = 0; i < n; i++ )
	{
		xmin[i] = lp.GetL( i );
		xmax[i] = lp.GetU( i );
	}

	for( i = 0; i < l-1; i++ ) weight[i] = Scen[i].GetProbability();
	weight[ l-1 ] = 1.0;


	//@BEGIN--------------------------------------------------------------------
	//Modifications: The vectors that depend on "l", number of scenarios, 
	//               are declared at the initializer list above
	//				 and here they are filled with 0
	//               -- no need to do this for weight
	
	g.Fill(0.0,		n*(n+2*l)	);
	a.Fill(0.0,		n+2*l		);
	v.Fill(0.0,		l			);
	
	pricba.Fill(0.0,	l		);
	dpb.Fill(0.0,		l		);

	iblock.Fill( (Int_T) 0,	n+2*l );
	ibasic.Fill( (Int_T) 0,	l	  );
	icheck.Fill( (Int_T) 0,	n+2*l );
	ieq.Fill(    (Int_T) 0,	n+2*l );
	drow.Fill(   (Int_T) 0,	n+2*l );


	//Modifications_2:  This is for testing solution quality...	
	//					ExpCost = 0.0 at the initializer list above
	//					and ExpC is set "l" elements and here 
	//					filled with 0

	ExpC.Fill(0.0, l); 
	ExpC2.Fill(0.0, l); 

	//@END----------------------------------------------------------------------
}


void MasterSolver::SetPenalty( Real_T Penalty )
{
	assert( Penalty >= 0 );

	if( Penalty < PENALTY_LO )
	{
		Warning( "Too small initial penalty %g. Assuming %g instead.",
			Penalty, PENALTY_LO );
		initpen = PENALTY_LO;
	}
	else if( Penalty > PENALTY_HI )
	{
		Warning( "Too large initial penalty %g. Assuming %g instead.",
			Penalty, PENALTY_HI );
		initpen = PENALTY_HI;
	}
	else
		initpen = Penalty;
}


/*------------------------------------------------------------------------------

	MasterSolver::~MasterSolver( void )

PURPOSE:
	Deallocates all memory allocated previously.

PARAMETERS:
	None.

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

MasterSolver::~MasterSolver( void )
{
	if( bmin )			delete bmin;
	if( bmax )			delete bmax;
	if( x )				delete x;
	if( xmin )			delete xmin;
	if( xmax )			delete xmax;
	if( y )				delete y;
	if( yb )			delete yb;		
	if( q )				delete q;
	if( r )				delete r;
	if( z )				delete z;
	if( w )				delete w;
	if( pricnb )		delete pricnb;
	if( pi )			delete pi;
	if( col )			delete col;
	if( inonba )		delete inonba;
	if( irn )			delete irn;
	if( istat )			delete istat;
	if( marks )			delete marks;
	
	//@BEGIN----------------------------------------
	//if( weight )		delete weight;
	//if( a )			delete a;
    //if( g )			delete g;
    //if( v )			delete v;
	//if( pricba )		delete pricba;
	//if( dpb )			delete dpb;	
	//if( iblock )		delete iblock;
	//if( ibasic )		delete ibasic;
	//if( icheck )		delete icheck;
	//if( ieq )			delete ieq;

	g.Resize(0); 
	a.Resize(0); 
	v.Resize(0); 
	weight.Resize(0); 
	
	pricba.Resize(0); 
	dpb.Resize(0); 
	
	iblock.Resize(0); 
	ibasic.Resize(0); 
	icheck.Resize(0); 
	ieq.Resize(0); 
	drow.Resize(0); 

	//for testing solution quality part:

	ExpC.Resize(0); 
	ExpC2.Resize(0);

	//@END-------------------------------------------
	   
}


/*------------------------------------------------------------------------------

	void MasterSolver::SetStartingPoint( const Array<Real_T> &x1,
		const Int_T x1n )

PURPOSE:
	Sets the initial starting point. "x1" is used instead of the vector of
zeros. Function will accept ANY point.

PARAMETERS:
	const Array<Real_T> &x1, const Int_T x1n
		Starting point coordinates (a vector and its length).

RETURN VALUE:
	None.

SIDE EFFECTS:
	None.

------------------------------------------------------------------------------*/

void MasterSolver::SetStartingPoint( const Array<Real_T> &x1, const Int_T x1n )
{
	assert( n == x1n );

	for( Int_T i = 0; i < x1n; i++ )
		x[i] = 	y[i] =	 x1[i];
}


/*------------------------------------------------------------------------------

	Int_T MasterSolver::Solve( void )

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

Int_T MasterSolver::Solve( void )
{

	//@BEGIN-------------------------------------------------------------
	//Modifications: The vectors that depend on "l", number of scenarios, 
	//               has been set to "vector_name.start" to get a 
	//				 pointer to the start of that array
	//@END---------------------------------------------------------------
	
	q1cmte_( n, l, x, y, yb, xmin, xmax, v.start, weight.start, mdmat, dmat, jcol, iptr,
		bmin, bmax, marks, status, g.start, a.start, iblock.start, icheck.start, ieq.start, 
		drow.start, ibasic.start, pricba.start, inonba, irn, istat, pricnb, pi, q, r, z, w, 
		col, dpb.start, itmax, &istop, *SubMan, levprt, initpen, ExpC.start, ExpCost, 
		ExpC2.start, ExpCost2);

	return 0;
}


/*------------------------------------------------------------------------------

	StochSolution *MasterSolver::GetSolution( void )

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

StochSolution *MasterSolver::GetSolution( void )
{
	StochSolution *sol = new StochSolution( l, n );

	if( sol == NULL ) FatalError( "Out of memory" );

	Int_T i;

	for( i = 0; i < n; i++ )
		sol->x[i] = ( IsZero( x[i] ) ) ? 0.0 : x[i];

	Objective = 0.0;
	
	
	for( i = 0; i < l; i++ )      
		Objective += 
			( sol->f[i]			= ( IsZero( v[i] ) ) ? 0.0 : v[i] ) *
			( sol->weights[i]	= weight[i] );

	sol->result = Objective;

	return sol;
}



//@BEGIN------------------------------------------------------------------------
//This section contains various functions written for the mastersolver class
//to implement SADSAM and SRP.....
//
//
void MasterSolver::ReInitialize( int InitScen )
{
	//ReInitialize the work vectors that have sizes related to "l"
	//by resizing to correct size 
	//so that when number of scenarios is changed (increased) 
	//they are good to go...


	//takes InitScen, initial number of scenarios as a parameter and assumes, 
	//"l", total number of scenarios, is adjusted to its new value beforehand.

	g.Resize(n*(n+2*l)); 
	a.Resize(n+2*l); 
	v.Resize(l); 
	weight.Resize(l); 

	pricba.Resize(l); 
	dpb.Resize(l);
	
	iblock.Resize(n+2*l);	
    ibasic.Resize(l);
    icheck.Resize(n+2*l);
	ieq.Resize(n+2*l);	
	drow.Resize(n+2*l);

	ExpC.Resize(l); 
	ExpC2.Resize(l); 

	
	//Initialize new elements to 0

	int len, start; 
	len = n+2*l; 
	start = n + 2*(InitScen+1); 

	iblock.Fill(0, len, start); 
	icheck.Fill(0, len, start);
	drow.Fill(  0, len, start);
	ieq.Fill(   0, len, start);

	a.Fill(   0.0, len, start);

	len *= n; 
	start = n*(n+2*(InitScen+1)); 
	g.Fill(0.0, len, start);

	len = l;
	start = InitScen + 1;
	
	v.Fill(      0.0, len, start);
	weight.Fill( 0.0, len, start);
	pricba.Fill(0.0, len, start);
	dpb.Fill(      0, len, start);
	ibasic.Fill(   0, len, start);

	ExpC.Fill(   0.0, len, start);
	ExpC2.Fill(   0.0, len, start);
}


//------------------------------------------------------------------------------

Real_T MasterSolver::CalculateVariance (Real_T& gap)
{
	//This function assumes that for a given xhat, 
	//  1.) ExpC is filled and ExpCost is calculated. 
	//  2.) The zn* problem is solved.
	//
	//Then, the function calcultes the gap and the variance
	//     estimates needed for the single rep. procedure (or,
	//     its variants) for testing solution quality 

	Real_T  var = 0.0;   //stores variance
	Real_T mean1, mean2; //stores means
	int i; 

	
	mean1 = 0.0; 
	mean2 = 0.0; 

	for(i = 0; i<l-1; i++){
		mean1 += ExpC[i];
		mean2 += v[i];
	}

	mean1 /= (Real_T) l-1; 
	mean2 /= (Real_T) l-1; 

	for(i = 0; i<l-1; i++)
		var += pow((ExpC[i] - mean1)  -  (v[i] - mean2), 2); 

	var /=  (Real_T) l-2;		//because l = #of scenarios + 1

	gap = ExpCost - Objective;  //minimization problem
	return var;
}

//------------------------------------------------------------------------------

Real_T MasterSolver::CalculateVariance2 (void)
{
	//Assessing Sol Quality. Calculates different variance

	Real_T  var = 0.0;   //stores variance
	Real_T mean1, mean2; //stores means
	int i; 

	mean1 = 0.0; 
	mean2 = 0.0; 
	for(i = 0; i<l-1; i++){
		mean1 += ExpC[i];
		mean2 += ExpC2[i];
	}
	mean1 /= (Real_T) l-1; 
	mean2 /= (Real_T) l-1; 

	for(i = 0; i<l-1; i++)
		var += pow((ExpC[i] - mean1)  -  (ExpC2[i] - mean2), 2); 

	var /=  (Real_T) l-2;		//because l = #of scenarios + 1

	return var;
}

//------------------------------------------------------------------------------

void MasterSolver::ReInit( void )
{
	//When it starts a news solve, just initialize the two vectors 
	//marks and status to 0, so that it should work fine... 

	Fill( marks,	mdmat,	(Int_T) 0 );
	Fill( status,	mdmat,	(Int_T) 0 );	

}

//------------------------------------------------------------------------------

void MasterSolver::FillSolution( StochSolution *sol )
{
	//This function is like GetSolution but does not create a new
	//StochSolution object, but fills in an already created one, 
	//*sol. If GetSolution is used, since "new" is used, this creates 
	//a memory leak problem.

	Int_T i;

	for( i = 0; i < n; i++ )
		sol->x[i] = ( IsZero( x[i] ) ) ? 0.0 : x[i];

	Objective = 0.0;
	
	
	for( i = 0; i < l; i++ )      
		Objective += 
			( sol->f[i]			= ( IsZero( v[i] ) ) ? 0.0 : v[i] ) *
			( sol->weights[i]	= weight[i] );

	sol->result = Objective;
}
//
//
//@END--------------------------------------------------------------------------



