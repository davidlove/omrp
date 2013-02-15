/* ========================================================================	*/
/*                          FILE Q2											*/
/*																			*/
/*     Q2MSTR - ALGORITHM FOR SOLVING THE REGULARIZED MASTER PROBLEM.		*/
/*     Called by: Q1CMTE													*/
/* ------------------------------------------------------------------------	*/
/*     WRITTEN BY A. RUSZCZYNSKI , INSTITUT FUER OPERATIONS RESEARCH,		*/
/*     UNIVERSITAET ZUERICH, MARCH 1985.									*/
/* ========================================================================	*/
/*      q2.f -- translated by f2c (version 19940705.1).						*/
/*--------------------------------------------------------------------------*/
/*		Modified by Artur Swietanowski										*/
/* ========================================================================	*/


#include <assert.h>
#include <math.h>

#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif

extern Int_T InitScen;


int q2mstr_( Int_T n, Int_T *nfix, Int_T *m, Int_T * mg, Int_T *l, 
	Real_T *g, Real_T *a, Real_T *xmin, Real_T *xmax, Real_T *x, Real_T *y,
	Real_T *yb, Real_T *v, Real_T *weight, Real_T *penlty, Bool_T * newpen,
	Int_T *iblock, Int_T *ibasic, Int_T *inonba, Int_T * icheck, Int_T *ieq,
	Int_T *irn, Int_T *istat, Real_T *q, Real_T *r, Real_T *z, Real_T *w,
	Real_T *pricnb, Real_T *pricba, Real_T *pi, Real_T *col, Real_T *dpb,
	Int_T *index, Real_T *gmax, Int_T *inew, Int_T *jnew,
	Int_T *iter, Real_T *initpen, Real_T &tolcut )
{
	/* System generated locals */
	Int_T g_dim1, g_offset, i_1;
	Real_T d_1;

	/* Local variables */
	static Int_T ibbl,		//
		         idel,		//constraint to be deleted (from the active set).
				 ldel,		//if constraint to be deleted is of a subproblem, gives index to that sub (1..L)
				 nrec = 0;	//?

	static Int_T ires,		//
		         lnew,		//block of the new cut. Takes value 0:if first stage constraint, a number between 1,2,..L:otherwise 
				 i, k, j, 
				 ifdep;		//linear dependence. 0: if linearly independent, 1:linearly dependent
	
	static Real_T sigma,	//these are used to look at QR factorization accuracy.
		         sigsum, 
				 sigmax, 
				 tolold;
	
	static Bool_T nodel;	//False, if there is an active constraint deleted from the active set
	static Int_T itmax;		//
	
	static Bool_T newyb;	//True, if yb is recalculated
	static Int_T ii;
	static Int_T resfac=500,//If QR factorization has not changed for "recfac" times, it is recalculated. 
				 ibl,		//?
				 res = 0,	//Number of times "if reset is desirable check" is done
							//If this number exceeds "resfac", factorization is redone
				 must = 0;	//?

	//Additional notes:     yb=x_B in the paper.... 


	/* Parameter adjustments */
	g_dim1 = n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--a; --xmin; --xmax; --x; --y; --yb; --v; --weight; --iblock; --ibasic;
	--inonba; --icheck; --ieq; --irn; --istat; --q; --r; --z; --w; --pricnb;
	--pricba; --pi; --col; --dpb;


	ifdep = 0;
	if( ! (*newpen) )
		goto L1;

	/*   PENALTY CHANGED BY Q1CMTE (upper level algorithm). RESET Z AND YB. */
	newyb = True;
    idel = 0;
	if ( *inew > 0 )
		if ( iblock[*inew] < 0 ) 
			if ( *gmax < tolcut ) {
				must = 1;
				res = resfac;
				//Print("Feasibility cut in block %d\n",iblock[*inew]);
				tolold = tolcut;
				tolcut = Min( *gmax/2.0, tolcut);
			}
	goto L4;

L1:
	switch( *index )
	{
	case 1: goto L7;    //Cold Start     (Initialize Everything)
	case 2: goto L5;	//Phase 2 Starts (Form original cuts; this contains Cold Start code as well)
	case 3: goto L2;	//Null Step
	case 4: goto L3;	//Serious Step   (Recompute z, yb and the prices)
	case 5: goto L77;	//new scenarios have been added...  GB
	}


	/*     IF NULL STEP, THE WARM START IS EASY.	*/
	/*     ONLY NORMALIZE PRICES FROM TIME TO TIME. */
L2:
	nrec = Int_T( (nrec + 1) % n );
	ires = 0;
	if( nrec == 0 )
		q3resp_( l, m, pricba+1, &pricnb[*nfix + 1], dpb+1,		//normalizes prices, if there are any
			&inonba[*nfix + 1], iblock+1, weight+1, ieq+1 );	//nonbasic cuts that come from a subproblem (objective cut)

	if (*jnew != 0)
	{
		*inew = 0;		//the new cut is a simple bound 
		goto L32;		//Add bound
	}
	else
	{
		assert ( *inew != 0 );
		goto L33;		// new cut is not a simple bound, first stage constraint, Add cut
		
	}

	/*     AFTER A SERIOUS STEP Z MUST BE RECOMPUTED. */


//@BEGIN-------------------------------------------------------------------------------------

L77:	//NEW SCENARIO(S) HAS(HAVE) BEEN ADDED //

	//reset z and yb
	newyb=True;
	q3resz_( &n, nfix, l, m, x+1, yb+1, &g[g_offset], r+1, z+1, pi+1, ibasic+1,
			inonba+1, irn+1, weight+1, penlty, &newyb);

	//calculate pricenb
	q3pric_( &n, nfix, m, y+1, yb+1, g+g_offset, r+1, pi+1, z+1, w+1, penlty,
			inonba+1, istat+1 );

	i_1 = *nfix + *m;
	for (i=1; i<= i_1; i++)
		 pricnb[i] = pi[i];
	
	//normalize prices (priceba)
	q3resp_( l, m, pricba+1, &pricnb[*nfix + 1], dpb+1,	
			&inonba[*nfix + 1], iblock+1, weight+1, ieq+1 );

	//Question: What happens to all of these static variables, nrec, res, ires... etc?
	//ires = 1;
	//nrec = 0;
	//itmax = 10;
	//res = 0;

	//calculate y  :y shouldn't change!! 
/*	q3gety_( &n, nfix, m, y+1, yb+1, pi+1, g+g_offset,
			pricnb+1, inonba+1, irn+1, penlty ); */

	//calculate v: Only calculate the new scenarios' obj func estimates 
    for (j = InitScen + 1; j<*l; j++)
	{
		i = ibasic[j];
		v[j] = a[i] + ddot_( n, &g[i * g_dim1 + 1], y+1);
	}
	//check bounds and cuts..: Bounds and cuts should be OK.  continue the algorithm
	//goto L25;
	++(*iter);
	goto L199; 

//@END---------------------------------------------------------------------------------------

L3:
	*index = 3;
	newyb = False;
L4:
	nodel = False;
	ires = 2;
	res = 0;
	//reset z and yb
	q3resz_( &n, nfix, l, m, x+1, yb+1, &g[g_offset], r+1, z+1, pi+1, ibasic+1,
		inonba+1, irn+1, weight+1, penlty, &newyb);
	//reset prices
	q3resp_( l, m, pricba+1, &pricnb[*nfix + 1], dpb+1, &inonba[*nfix + 1],
		iblock+1, weight+1, ieq+1 );
	
	for (i = 1; i <= *mg; ++i)
		if (icheck[i] == -1)
			icheck[i] = 0;
	goto L85;


L5:		// PHASE 2 STARTS. FORM ORIGINAL CUTS. 
	if( *m >= 1 )
	{
		for( i = 1; i <= *m; ++i )
		{
			ii = inonba[*nfix+i];
			ibl = Int_T( -iblock[ii] );
			if( ibl < 1 ) continue;
			ibl = ibasic[ibl];
			dsuma_( n, &g[ii * g_dim1 + 1], &g[ibl * g_dim1 + 1]);
			a[ii] += a[ibl];
		}
		*nfix = 0;
	}

L7:		// COLD START. EVERYTHING MUST BE INITIALIZED.
	*m = 0;
	ires = 1;
	nrec = 0;
	*penlty = *initpen;
	nodel = True;
	itmax = 10;
	res = 0;
	//resfac = 10;
	for( i = 1; i <= *mg; ++i )
		icheck[i] = 1;

	// FIND BASIC CONSTRAINTS AND RECOVER THE FIRST Y.

	dzero_( n, &yb[1]);
	k = 1;
	for( ii = 1; ii <= *l; ++ii )
	{
		for( i = 1; i <= *mg; ++i )
		{
			if( iblock[i] != k ) continue;   

			icheck[i] = -2;
			ibasic[k] = i;
			pricba[k] = weight[k];
			dstep_( n, &yb[1], &g[i * g_dim1 + 1], -pricba[k] );
			if (++k > *l)					//if not enough, then, basis is not found... 
				goto L12;
		}
	}
	
	// BASIS (Set of Basic Cuts) NOT FOUND. 
	*index = 10;
	goto L200;

L12:
	d_1 = 1. / *penlty;
	dmult_( n, yb+1, &d_1);
	dsuma_( n, yb+1, x+1);
	dcopy_( n, yb+1, y+1);


	
	/*     RECOVER Y (IF DELETION OCCURED) AND V. */

L15:
	if (! nodel)	
		q3gety_( &n, nfix, m, y+1, yb+1, pi+1, g+g_offset,
			pricnb+1, inonba+1, irn+1, penlty );

	for( k = 1; k <= *l; ++k )
	{
		i = ibasic[k];

		v[k] = a[i] + ddot_( n, &g[i * g_dim1 + 1], y+1);
	}

	/*     CHECK WHETHER RESET IS DESIRABLE. */

	++(*iter);
	++res;
	if( res >= resfac )
		goto L205;

	if( ires == 1 || *m < 1 )
		goto L25;

	sigmax = sigsum = 0.0;
	for( i=1; i<=*m; i++ )
	{
		ii = inonba[*nfix+i];
		sigma = fabs( a[ii] + ddot_( n, &g[ii * g_dim1 + 1], y+1) );
		sigsum += sigma;
		if( sigma > sigmax ) sigmax = sigma;
		if( sigma > tolcut ) goto L205;
	}
//	Print( "Factorization accuracy: sum = %g, max = %g\n", sigsum, sigmax );
	if( sigmax < 0.1 * tolcut )
	{
		Real_T NewTolcut = Max( 2.0 * sigmax, 1.0e-8 );

		//if( tolcut > NewTolcut )
			tolcut = NewTolcut;
	}
	goto L25;

	/*     RESET THE FACTORIZATION IF THERE IS A CHANCE OF IMPROVEMENT. */
	/*     OTHERWISE RELAX THE TEST ( THE PROBLEM IS ILL-CONDITIONED ). */

L205:
	if( ires == 2 )
	{
		assert( res <= 1 );
		goto L23;
	}

	Print( "Master factorization reset. Sigma[%d] = %g\n",i,sigma );
	assert( sigma < 1.0 );
	q3resr_( &n, nfix, l, m, &g[g_offset], a+1, x+1, y+1,
		yb+1, weight+1, penlty, q+1, r+1, z+1, 
		w+1, col+1, pi+1, pricba+1, pricnb+1, iblock+1,
		ibasic+1, inonba+1, irn+1);
	q3resp_(l, m, pricba+1, &pricnb[*nfix + 1], dpb+1,
		&inonba[*nfix + 1], iblock+1, weight+1, ieq+1);
	*iter = 1;
	nrec = 0;
	ires = 2;
	nodel = False;
	res = 0;
	goto L85;
L23:
	tolcut = Max( 2.0 * tolcut, sigma );

L25:	/*     CHECK BOUNDS AND CUTS. */

	ires = 0;
	nodel = True;
	*gmax = q3chbd_( Int_T( n - *nfix ), xmin+1, xmax+1, y+1, irn+1, istat+1,
		jnew, tolcut, 1 );
	if ( must)
		if ( *gmax < tolold ) {
			*gmax=0.0;
			*jnew = 0;
		}
	q3chct_( &n, mg, &g[g_offset], a+1, y+1, v+1, iblock+1, icheck+1, inew, gmax,
		tolcut, ieq+1 );
	if ( must ) {
		tolcut = tolold;
		must = 0;
	}

	if( *inew > 0 )
	{
		*jnew = 0;
		goto L33;
	}
	if( *jnew == 0 )
		goto L199;

L32:	/*     ADD THE NEW MEMBER TO THE ACTIVE SET. */
		
	//Add bound
	q4adbd_( &n, nfix, m, xmin+1, xmax+1, &g[g_offset], q+1, r+1, z+1, w+1, y+1,
		yb+1, col+1, pi+1, pricnb+1, inonba+1, irn+1, istat+1, &ifdep, &nodel,
		jnew);
	++(*iter);
	if( ifdep <= 0 )
		goto L85;
	else
		goto L50;

L33://Add cut						
	icheck[*inew] = -2;
	lnew = iblock[*inew];   
	
	//calculate reduced cut
	if( lnew >= 1 )			
	{
		ibbl = ibasic[lnew];
		ddifr_( n, &g[*inew * g_dim1 + 1], &g[ibbl * g_dim1 + 1]);   
		a[*inew] -= a[ibbl];
	}

//	if ( lnew < 0 )
//		Print("Add feasibility cut in block %d with error %G\n",lnew,*gmax);
L45:
	q4adct_( &n, nfix, m, &g[g_offset], a+1, q+1, r+1, z+1, w+1, y+1, yb+1, col+1,
		pi+1, pricnb+1, inonba+1, irn+1, &ifdep, &nodel, inew);
	if( ifdep == 0 )
		goto L85;

	--(*m);

L50:	/*     LINEAR DEPENDENCE. SELECT THE MEMBER TO BE DELETED. */

//	Print( "Linear dependence \n");
	q3ldep_( &n, nfix, m, r+1, pi+1, col+1, g+g_offset, inonba+1, istat+1,
		inew );
	i_1 = Int_T( *nfix + *m + 1 );
	q3corr_( nfix, &i_1, l, &idel, &ldel, &ifdep, pricnb+1, pricba+1, pi+1,
		dpb+1, inonba+1, iblock+1, ieq+1 );
	if( ifdep < 0 )
	{
		/*     THE FEASIBLE SET IS EMPTY. */
		if( ires == 2 ) {
			*index = 20;    //inconsistent constraints, no finite theta (in the paper)
			goto L200;
		}
	}

L55:	/*     DELETE THE SELECTED MEMBER FROM THE ACTIVE SET. */

	nodel = False;
	++(*iter);
	
	//if index to the deleted active constraint is greater than the number of fixed variables,
	//then, delete cut, otherwise delete bound
	if( idel > *nfix )
		q4dtct_( &n, nfix, m, &g[g_offset], a+1, yb+1, weight+1, penlty,
			iblock+1, ibasic+1, inonba+1, icheck+1, r+1, z+1, w+1, pricnb+1,
			pricba+1, col+1, &idel, &ldel, &ifdep);
	else
		q4dtbd_( &n, nfix, m, &g[g_offset], y+1, yb+1, r+1, z+1, w+1, pi+1,
			pricnb+1, inonba+1, irn+1, istat+1, &idel, &ifdep);

	/*     IF DELETION FOLLOWED FROM DEPENDENCE, APPEND AGAIN. */

	if( ifdep == 0 )
		goto L85;

	if( *inew <= 0 )
		goto L32;
	else
		goto L45;


	/*     FIND PRICES FOR THE NEW ACTIVE SET. */

L85:
	q3pric_( &n, nfix, m, y+1, yb+1, g+g_offset, r+1, pi+1, z+1, w+1, penlty,
		inonba+1, istat+1 );
	ifdep = 0;
	ddifr_( Int_T( *nfix + *m ), pi+1, pricnb+1 );  //pi[] = pi[] - pricnb[] 
	i_1 = Int_T( *nfix + *m );
	q3corr_( nfix, &i_1, l, &idel, &ldel, &ifdep, pricnb+1, pricba+1, pi+1,
		dpb+1, inonba+1, iblock+1, ieq+1 );

	/*     IF NEW PRICES ARE NONNEGATIVE, START THE NEXT ITERATION. */
	/*     OTHERWISE RELAX ONE CONSTRAINT. */

	if( idel == 0 )  
		goto L15;	 //start next iteration

	goto L55;		//delete the constraint where non-negativity bound is hit
L199:
	*index = 3;
L200:
	return 0;


}
/* -----END OF Q2MSTR----------------------------------------------------- */

//possible aux function q2newScen()??
