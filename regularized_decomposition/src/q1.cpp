/*==============================================================================
                         FILE Q1
    Q1CMTE - Upper Level Algorithm of the REGULARIZED DECOMPOSITION METHOD

    Q1SLCT - AUXILIARY ROUTINE FOR SELECTING ACTIVE COMMITTE MEMBERS
		and other functions that work with these
--------------------------------------------------------------------------------
    WRITTEN BY A. RUSZCZYNSKI , INSTITUT FUER OPERATIONS RESEARCH,
    UNIVERSITAET ZUERICH, APRIL 1985.
    DATE LAST MODIFIED: NOVEMBER 1985.
--------------------------------------------------------------------------------
    q1.f -- translated by f2c (version 19940705.1). from FORTRAN
--------------------------------------------------------------------------------
    Modified by Artur Swietanowski.
	Then by Guzin Bayraksan.  
	Last modified: Feb 20, 2003 
===============================================================================*/


#include <math.h>
#include <assert.h>
#include <stdlib.h>


#ifndef __QDX_LOC_H__
#	include "qdx_loc.h"
#endif
#ifndef __QDX_PUB_H__
#	include "qdx_pub.h"
#endif
#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __ERROR_H__
#	include "error.h"
#endif
#ifndef __PRINT_H__
#	include "print.h"
#endif
#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif



//@BEGIN-----------
#ifndef __SPTR_DEB_H__
#	include "sptr_deb.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif
extern Real_T PEN;
extern Int_T MGG, NFIX2, M2; 
extern Int_T InitScen;

extern Real_T Opt_Tol; 

extern Real_T whatistolcut_2;
//@END-------------



/*-----------------------------------------------------------------------------------------

	void q1cmte_(...)
	
PURPOSE:
	The upper-level algorithm of the regularized decomposition method. 
	Realizes the phase 1 - phase 2 procedure. 
	At each iteration:
		- verifies bounds on variables and constraints,
		- if constraints satisfied, invokes blocks 1,...,l,
		- updates the committee by deleting inactive members and adding
		  new direct or generated cuts,
		- verifies whether serious or null steps should be made.

PARAMETERS:
	For definition of these parameters, see "mast_sol.h"

RETURN VALUE:
	None.

-----------------------------------------------------------------------------------------*/



void q1cmte_( const Int_T n, Int_T l, Real_T *x, Real_T *y, Real_T *yb,
	Real_T *xmin, Real_T *xmax, Real_T *v, Real_T *weight, Int_T mdmat,
	Real_T *dmat, Int_T *jcol, Int_T *iptr, Real_T *bmin, Real_T *bmax,
	Int_T *marks, Int_T *status, Real_T *g, Real_T *a, Int_T *iblock,
	Int_T *icheck, Int_T *ieq, Int_T *drow,
	Int_T *ibasic, Real_T *pricba, Int_T *inonba, Int_T *irn, 
	Int_T *istat, Real_T *pricnb, Real_T *pi, Real_T *q, Real_T *r, 
	Real_T *z, Real_T *w, Real_T *col, Real_T *dpb, Int_T itmax, 
	Int_T *istop, RD_SubproblemManager &SubMan, Int_T levprt, Real_T initpen, 
	Real_T *ExpC, Real_T& ExpCost, Real_T *ExpC2, Real_T& ExpCost2)
{
	//Print("\nEntered Q1MTE - Upper Level Algorithm\n");
	
	/* System generated locals */
	Int_T g_dim1, g_offset;			//g_dim1 = # of first stage variables, g_offset = g_dim1 + 1
 
	/* Local variables */
	Int_T nsec = 0;					//Number of sections of 1st stage constraints
									//nsec>1 if  no. of 1st st constraints >= 2 * no. of 1st stage variables

	Real_T gmax;					//Maximum difference(infeasibility) in the bounds (or constraints)
	Int_T iter = 0;					//Number of iterations (upper level algorithm)
	Int_T iterq;					//Number of iterations (Solving Regularized Master Algorithm)

	Int_T nfix,						//Number of fixed variables (=at their lower or upper bounds)

		  inew,						//Returns the index of the violated 1st stage constraint, if the violation
									//is larger than the max violation is bounds (CheckContraints in Sections)

		  jnew;						//Returns variable that is out of bounds (CheckBounds),
									//it is -j if lower bound is violated, +j if upper bound is violated
	
	Real_T vsum = 0.0;				//vsum = Sum(i, v[i] * weight[i]), i=1...l,l+1. (l+1 is used for 1st st. obj value)
									//this is the obj func value found by solving the regularized master
									//(it does not include the regularized term!) 
									//vsum <= funold. When two are sufficiently close, the algorithm terminates

	Int_T ifeas;					//0 --> infeasible
	Int_T ifeas1 = 0;				//Infeasibility of the 1st stage
	Int_T ifeas2 = 0;				//Infeasibility of the subproblems (then a feasibility cut is added)
	
	Int_T index;					//Stores different stages/results for solving regularized master problem
	
	Real_T value;					//used in solving subproblems; stores the cut value 

	Real_T f1,						//these are assigned in SolveBlocks and are used in DetermineStepType
		   f2;						//f1 = weighted sum of feas. cut "value"s, f2 = weighted sum of obj cut "value"s
	Real_T fx;						//fx = value of objective at the regularization point

	Real_T gamma = 0.0;				//tolerance to change the regularization point and the penalty
									//0.5 <= gamma <= 1.0  (usually gamma = 0.9)
	
	Int_T i, k, j,
			m;						//m = number of Nonbasic cuts, index to nonbasic cuts starts from nfix+1 in inonba[]
	Int_T j1, j2 = 0;

	Int_T mg;						//size of the committee of cuts
	Int_T iphase = 0;				//phase of problem iphase = 1 --> phase 1, similarly 2
	Real_T funold = 0.0;			//Stores the best upper bound found in the algorithm. 
									//if a better value is found, it can be changed in the DetermineStepType  (funold = f2)
									//funold = Sum(i, a[i]+g[i]x * weight[i]) i=1...l,l+1. (l+1 is used for 1st st. obj value)
									//objective value at the current regularization point
		
	Int_T iser = 0;					//step type -1:Null inf, 0:Null feas, 1:apprx ser, 2:exact serious, 
									//Also, used as: -3: if bounds are violated, -2: if first stage constraints are violated
	Int_T nserap,					//number of approximate serious steps
		  nserex,					//number of exact serious steps
		  nulinf,					//number of null infeasible steps
	       nulf;					//number of null feasible steps

	Real_T penmax = 1.0;			//maximum penalty 
	Bool_T newpen;					//True, if penalty value has been changed. 
	Int_T  nwcmax = 0;				//max. # of violated 1st stage constraints to be added in one iteration = min(n,l)
									//used in CheckAllConstraints 
	
	Real_T funmin = 0.0;			// ?
	Int_T  istart = 0;				// ?  used in DetermineStepType
	Real_T tolcut = 1.0e-7,			//tolcut is used as a tolerace value for numerical differences in calculations
		penlty,						//the value of the penalty (initialized to initpen, later changed)
		wgtsum = 0.0;				//wgtsum = sum of weights. 
	
	
	//@BEGIN-------------------------------------------------------
//	int    ct;						//used for a simple counter
	bool   flag = true;			//used for really exiting q1. 
	Real_T tolcut_2 = 1E-08;   		//to get out with less/more tolerance 
									//original value is tolcut = 1E-08
	Real_T minn;					//minimum of |funold| and |vsum| used for termination check
	bool   lshaped = false;			//true, if want to run l-shaped. 
									//also, to run L-shaped, make flag = true

	whatistolcut_2 = tolcut_2; 
	//@END---------------------------------------------------------

	/* Parameter adjustments */

	g_dim1 = n;
	g_offset = g_dim1 + 1;
	g -= g_offset;
	--x; --y; --yb; --xmin; --xmax; --v; --weight; --a;
	--dmat; --jcol; --iptr; --bmin; --bmax;
	--marks; --status; --iblock; --icheck; --ieq; --drow; --ibasic;
	--pricba; --inonba; --irn; --istat; --pricnb; --pi; --q;
	--r; --z; --w; --col; --dpb;

	//-------------------------------------------------------------------------------------
	//	HOT START.

	if( *istop == 3 ){
		Print ("\nHOT START\n");
		goto L60;
	}


	//@BEGIN-------------------------------------------------------------------------------
	//        PURPOSE : Given xhat, just want obj function estimate
	//
	//if ( *istop == 10 ) goto L60; 
	//@END---------------------------------------------------------------------------------



	//@BEGIN-------------------------------------------------------------------------------
	//        PURPOSE: To Re-Start the algorithm after another scenario 
	//				  has been added.

	if (*istop == 9)
	{
		//Print ("\nRe-Starting the Algorithm...Added another scenario...\n");
		
		//---------------------------------------------------------------------------------
		//Adjust parameters and work vectors
		//---------------------------------------------------------------------------------
		
		//Initialize some variables and Adjust the parameters
		
		iter = 0;
		iterq = 0;
		nserex = 0;
		nserap = 0;
		nulinf = 0;
		nulf = 0;
		inew = 0; 
		jnew = 0; 
		funold = 0.0;
		gamma = .9;
		nsec = Max( mdmat / n, 1 );
		nwcmax = Min( n, l );
			
		mg = MGG;
		nfix = NFIX2;
		m = M2;
	
		iphase = 2;
		ifeas = 1;
		istart = 1;

		penlty = PEN; //initpen;
		newpen = False;	
		Print("\nPENALTY : %10g\n", penlty);
		
		//---------------------------------------------------------------------------------
		// Adjust vectors related to "l", the number of scenarios
		// a[n+2l], g[n*(n+2l)], v[l], weight[l], dpb[l]
		// iblock[n+2l], ibasic[l], priceba[l], icheck[n+2l], ieq[n+2l], drow[n+2l]
		// Note that most of this is done at MasterSolver.ReInitialize 
		//---------------------------------------------------------------------------------
		
		//adjust weight
		for ( i = 1; i < l; i++)
			weight[i] =(double) 1/(l-1);
		weight[l] = 1;

		wgtsum = 1.;
		for( i = 1; i <= l; ++i ) wgtsum += weight[i];

		//adjust priceba
		for (i = 1; i <= l; i++)
			pricba[i] =  weight[i];


		//Below should already be done in MasterSolver.ReInitialize
		//but leaving it anyway... 

		//adjust v and dpb
		v[l] = v[InitScen+1];
		for (i = InitScen+1; i < l; i++)
			v[i] = 0.0;

		for (i = InitScen+2; i<=l; i++)
			dpb[i] = 0.0;
		

		//---------------------------------------------------------------------------------
		// Solve the new scenario(s) at the current point(=optimal point of the prev. solve)
		// and get the new cut, add to the committee of cuts
		//---------------------------------------------------------------------------------
		SubMan.ChangePreviousBlockNumber(0);
		SubMan.ReSizeSolverState(l-1, InitScen);

		for (i = InitScen+1; i < l; i++)
		{
			dzero_( n, &g[(mg+1) * g_dim1 + 1] ); 

			if(!SubMan.SolveSubproblem(i, n, y+1, value, g+((mg+1)*n +1))) {
				//Print("\nworked!!\n");
			
				//Feasibility Cut
					
				//f1 += value * weight[l];
				ifeas2 = 0;
				icheck[mg+1] = 1;
				gmax = value;
				iblock[mg+1] = -i;  //-(l-1);
			}
			else
			{	//Objective Cut
				
				//f2 += value * weight[l];
				//?? if (value <= v[l] + tolcut) continue;
				icheck[mg+1] = 1;        // ?? 
				gmax = value ;
				iblock[mg+1] = i;   //l-1;
			}   

			a[mg+1] = value - ddot_(n, g+((mg+1)*n +1), y+1);
			ieq[mg+1] =0;
			//drow[mg+1] = 0;

			if (inew == 0 ){
				inew = i; //l-1;		
				jnew = 0;
			}

			mg++;
		}  //end of SolveBlocks loop w/ i. 


		//adjust ibasic
		j = ibasic[InitScen+1];    //location of obj cut in G
		ibasic[l] = j;
		iblock [j] = l; 

		for (i=InitScen+1; i<l; i++){
			j = MGG + i - InitScen;
			ibasic[i] = j;
			icheck[j] = -2;
		}

		
	/*	j = ibasic[i];  //l-1];  
		iblock[j]++;
			
		ibasic[l-1] = mg+1; //
		ibasic[l] = j;      //		
		icheck[mg+1] = -2;  //	    
	*/			
	 //	--mg;

		//adjust inonba
		//for (i=1; i<=m; i++) inonba[nfix+i]++;

		//calculate funold 
		for (j = 1; j <= l; j++){  //changed i in for loop to j
			k = ibasic[j];
			funold += (a[k] + ddot_(n, g+(k*g_dim1+1), x+1))*weight[j];
		}
			
		//send to Solve Master	
		index = 5;
		goto L91;   
	}
	//@END---------------------------------------------------------------------------------



	//--------------------------------------------------------------------------
	//	INITIALIZATION
	//
	iter = 0;
	iterq = 0;
	nserex = 0;
	nserap = 0;
	nulinf = 0;
	nulf = 0;
	gamma = .9;
	wgtsum = 1.;
	for( i = 1; i <= l; ++i ) wgtsum += weight[i];

	funold = 0.0;
	vsum = 0.0;


	//
	//	Computing MAX
	//
	nsec = Max( mdmat / n, 1 );
	nwcmax = Min( n, l );
	j2 = 0;
	
	dcopy_( n, x+1, y+1 );
	
	mg = 0;
	m = 0;
	iphase = 1;

	//@BEGIN--------------------------------------------------
	//This is to ensure that when only calculating xhat obj
	//values, v does not change
	if ( *istop != 10 ){
		dzero_( l, v+1 );
	}
	//@END----------------------------------------------------

	nfix = 0;
	for( i = 1; i <= n; ++i )		istat[i]	= 2;
	for( i = 1; i <= n; ++i )		irn[i]		= i;
	for( i = 1; i <= mdmat; ++i )	marks[i]	= 1;

	//--------------------------------------------------------------------------
	//	START THE NEXT ITERATION.
	//
	for(;;)
	{
		ifeas1 = 1;
		ifeas2 = 1;
		newpen = False;
		inew = 0;
		jnew = 0;

		//----------------------------------------------------------------------
		//	CHECK BOUNDS ON VARIABLES.
		//
//		Print("CheckBounds\n");
		gmax = q3chbd_( n-nfix, xmin+1, xmax+1, y+1, irn+1, istat+1, &jnew,
			tolcut, 2);
		//Print("jnew = %d    gmax = %G\n",jnew,gmax);
		if( gmax > tolcut )
		{
			ifeas1 = 0;
			iser = -3;
		}

		//----------------------------------------------------------------------
		//	VERIFY CONSTRAINTS IN SECTIONS.
		//
		if( mdmat < 1 ) goto L60;

//		Print("CheckAllConstraints\n");
		if( CheckAllConstraints( j1, j2, mdmat, n, nsec, mg, iblock, icheck,
			g_dim1, g, drow, a, ieq, marks, iptr, dmat, jcol, y, status, bmin,
			bmax, tolcut, nwcmax, gmax, inew, jnew, ifeas1 ) )
		{
			if( gmax <= tolcut )
				goto L60;

			iser = -2;
		}

		if( iter == 0 ) goto L60;

		ifeas = 0;
		goto L91;

	L60:
		//----------------------------------------------------------------------
		//	SOLVE BLOCKS AT THE CURRENT TRIAL POINT.
		//
		if( ++iter > itmax )
		{
			*istop = 3;
			if( levprt > 0 ) Print( "ITERATION LIMIT.\n" );
			return;
		}

//		Print("Entering SolveBlocks\n");

		SolveBlocks( f1, f2, fx, mg, l, n, y, value, g, weight, ifeas2, iphase,
			iblock, icheck, tolcut, a, ieq, inew, jnew, gmax, g_dim1, v, x,
			SubMan, istop, ExpC, ExpCost, ExpC2, ExpCost2 );

		if( iphase == 1 )
			funold = Max(funold,(Real_T) fx);

		--mg;



		//@BEGIN---------------------------
		if ( (*istop == 10)  || (*istop == 11)){
			 return; 
		}
		
		//Here, we got the "f2" this is our 1/n[sum(i,f(xhat, xi_i))]
		// and we stored the "value"s within the SolveBlocks 
		// as our f(xhat, xi_i). So, we can return to main program.
	
		//@END-----------------------------



		//Print("DetermineStepType\n");
		DetermineStepType( iphase, ifeas1, ifeas2, ifeas, f1, f2, vsum, mg, l,
			n, g_dim1, icheck, iblock, g, a, iter, index, funold, funmin,
			istart, iser, levprt, y, x, pi, yb, istat, irn, penlty, initpen,
			newpen, penmax, nulinf, nserex, nulf, nserap, tolcut, gamma,
			inew, jnew, lshaped );
		
		
	
	//	CountCriticalScenarios( l, m, nfix, iblock, inonba );

		if( levprt >= 3 )
			Print( "%5d  %16.8E  %16.8E  %-12s\n",
				(int) iter, (Real_T) f1, (Real_T) f2,
				( iser == -1 ) ? "NULL INFEAS" :
					( iser == 0 ) ? "NULL FEAS" :
					( iser == 1 ) ? "APPROX SERIOUS" :
					( iser == 2 ) ? "EXACT SERIOUS" : "???" );
		//Print("\nPENALTY : %10g\n", penlty);

		
	L91:
		//----------------------------------------------------------------------
		//	SOLVE THE REGULARIZED MASTER.
		//
		//Print("SolveMaster. Penalty = %G  Tolcut = %G\n",penlty,tolcut);

		//if(flag){
		//	Print("\nv's before entering Solve Master:\n"); 
		//	for(ct = 1; ct<l+1; ct++)  Print("v[%d] = %f\n", ct, v[ct]); 
		//}

//		Print("Entering SolveMaster\n");
//		Print("\nPENALTY : %10g\n", penlty);
		q2mstr_( n, &nfix, &m, &mg, &l, g+g_offset, a+1, xmin+1, xmax+1, x+1,
			y+1, yb+1, v+1, weight+1, &penlty, &newpen, iblock+1, ibasic+1,
			inonba+1, icheck+1, ieq+1, irn+1, istat+1, q+1, r+1, z+1, w+1,
			pricnb+1, pricba+1, pi+1, col+1, dpb+1, &index, &gmax, &inew,
			&jnew, &iterq, &initpen, tolcut );

		//if(flag){
		//	Print("\nv's after exiting from Solve Master:\n"); 
		//	for(ct = 1; ct<l+1; ct++)  Print("v[%d] = %f\n", ct, v[ct]); 
		//} 

		if( index == 20 )
		{
			//	DIRECT CONSTRAINTS ARE INCONSISTENT.
			*istop = 4;
			if( levprt > 0 )
			{
				if( iphase == 1)
					Print( "DIRECT CONSTRAINTS ARE INCONSISTENT\n" );
				else
					Print( "SUBOPTIMAL SOLUTION FOUND\n" );
			}
			return;
		}



		//----------------------------------------------------------------------
		//	COMPRESS THE COMMITTEE, IF THERE ARE MORE THAN N+L MEMBERS.
		//
//		Print("Entering Compress\n");
		CompressCommittee( mg, n, l, m, nfix, istat, inonba, icheck, ieq,
			status, drow, g, a, iblock, ibasic, g_offset );


		//----------------------------------------------------------------------
		//	TEST FOR TERMINATION.
		//

		vsum = ddot_( l, weight+1, v+1 );
	
		//@BEGIN----------------------------------------------------------
		//changing the stopping criterion to suboptimality... 
		
		//find min of |funold| and |vsum|
		//if (fabs(funold) < fabs(vsum)) minn = fabs(funold); 
		//else minn = fabs(vsum); 

		
		//if( ifeas == 0 || iter == 1 ||
		//	funold - vsum >= tolcut_2 * (minn + 1) )  continue; 


		//Below was the original stopping (or, non-stopping) part. 
		if( ifeas == 0 || vsum <= funold - (1.0+wgtsum) * tolcut_2 || iter == 1 )
			continue;

		if (iter == 0) continue;	//I added this part.. GB 
									//I don't remember why exactly but might be to do w/ adding cuts
		//@END------------------------------------------------------------

		if( vsum <= funold + sqrt( (Real_T)l ) * wgtsum * tolcut_2 ||
			vsum <= funold + (fabs(funold) + 1.) * tolcut_2 )
		{
			if( iphase == 2 )
			{
				//--------------------------------------------------------------
				//	SOLUTION FOUND.
				//
				*istop = 2;

				if( levprt > 0 )
					Print( "OPTIMAL POINT FOUND.\n" );
				
				Print("\nPENALTY : %10g", penlty);
				Print("\nTolerace : %10g\n", tolcut_2);
				
				//@BEGIN--------------------------------------------------------
				if(!flag){  //not really exit yet, solve master w/o regularizing term to 
							//get a good lower bound

					penlty = 1.0E-6; //1.0E-15; //original PENALTY_LO = 1.0E-6;  //penlty = (1/sigma) in the paper
					newpen  = True; 
					flag = true; 
					
					while (flag){

						//solve master
//						Print("\nSolve Master\n");
						q2mstr_( n, &nfix, &m, &mg, &l, g+g_offset, a+1, xmin+1, xmax+1, x+1,
								y+1, yb+1, v+1, weight+1, &penlty, &newpen, iblock+1, ibasic+1,
								inonba+1, icheck+1, ieq+1, irn+1, istat+1, q+1, r+1, z+1, w+1,
								pricnb+1, pricba+1, pi+1, col+1, dpb+1, &index, &gmax, &inew,
								&jnew, &iterq, &initpen, tolcut );

						//NOTE: Here, tolcut changes, if there are rounding inaccuracies. 
						//So, it affect the checks after Checkbounds and constraints. 
						
						if( index == 20 ){
							//	DIRECT CONSTRAINTS ARE INCONSISTENT.
							*istop = 4;
							if( levprt > 0 ){
								if( iphase == 1)
									Print( "DIRECT CONSTRAINTS ARE INCONSISTENT\n" );
								else
									Print( "SUBOPTIMAL SOLUTION FOUND\n" );
							}
							return;
						}

						vsum = ddot_( l, weight+1, v+1 );

						//Print("\nRESOLVED MASTER WITH PENALTY %10g", penlty );
						//Print("\nand obtained vsum = %10g\n", vsum );
						//Print ("\nOPTIMAL DECISION VECTOR to Reg Master:\n");
						//for (i = 0; i< g_dim1; i++){
						//	Print ("y*[%1d] = " 
						//		   "%10g\n", 
						//		   i+1, y[i+1] );
						//}
		
						//delete inactive cuts... 
						CompressCommittee( mg, n, l, m, nfix, istat, inonba, icheck, ieq,
											status, drow, g, a, iblock, ibasic, g_offset );
						
						//check if bounds and constraints are satisfied  
						//  (this part is like starting over another iteration). 

						ifeas1 = 1;
						ifeas2 = 1;
						newpen = False;
						inew = 0;
						jnew = 0;

						gmax = q3chbd_( n-nfix, xmin+1, xmax+1, y+1, irn+1, istat+1, &jnew,
							tolcut, 2);
						//Print("jnew = %d    gmax = %G\n",jnew,gmax);
						if( gmax > tolcut ){
							ifeas1 = 0;
							iser = -3;
						}


						if( CheckAllConstraints( j1, j2, mdmat, n, nsec, mg, iblock, icheck,
							g_dim1, g, drow, a, ieq, marks, iptr, dmat, jcol, y, status, bmin,
							bmax, tolcut, nwcmax, gmax, inew, jnew, ifeas1 ) )
						{
							if( gmax <= tolcut )  //then, all constraints and bounds are satisfied
								flag = false; 
							else
								ifeas = 0;
						}

					}  //end of while with flag
	

					Print("\nRESOLVED MASTER WITH PENALTY %10g", penlty );
					Print("\nand obtained vsum = %10g\n", vsum );
					Print ("\nOPTIMAL DECISION VECTOR of Reg Master:\n");
					for (i = 0; i< g_dim1; i++){
						Print ("y*[%1d] = " 
							   "%10g\n", 
							   i+1, y[i+1] );
					}  
					
					//Now, a new twist...  Instead of:
		
					//here, calculate the true optimality gap. . . 
					if (fabs(funold) < fabs(vsum)) minn = fabs(funold); 
					else minn = fabs(vsum); 

					Opt_Tol = (funold - vsum) / minn; 
					
					return;			  //exit q1


					//will have it run L-shaped... 
					//lshaped= true; 
					//flag = true; 
					//goto L60; 

				}  //end of if with (!flag)   */
	 
		

				PEN = penlty; 
				MGG = mg;
				NFIX2 = nfix;
				M2 = m;		
				//@END----------------------------------------------------------

				return;
			}

			if( vsum < wgtsum * tolcut )
				continue;

			//------------------------------------------------------------------
			//	INDUCED AND DIRECT CONSTRAINTS ARE INCONSISTENT.
			//
			*istop = 5;
			if( levprt > 0 )
				Print( "INDUCED AND DIRECT CONSTRAINTS ARE INCONSISTENT\n" );

			return;
		}

		//----------------------------------------------------------------------
		//	IN PHASE 2, VSUM GREATER THAN FUNOLD INDICATES NONCONVEXITY.
		//
		if( iphase == 2 )
		{
			//*istop = 100;
			//Print( "The recourse functions appear non-convex: "
				//"FunOld = %G, vSum = %G.\n", funold, vsum );

			//return;
			continue;
		}

		//----------------------------------------------------------------------
		//	IN PHASE 1 THIS MAY FOLLOW FROM THE CHOICE OF INDUCED CUTS.
		//
		funold = vsum + wgtsum * tolcut;
	}
	
}
/* -----END OF Q1CMTE-----------------------------------------------------	*/



/*--------------------------------------------------------------------------*/
/*			 Q1SLCT															*/
/*  Purpose: TO COMPRESS THE COMMITTEE BY DELETING INACTIVE MEMBERS.		*/
/*																			*/
/*  Called by:	CompressCommittee											*/
/*  Subroutines Called: Dcopy												*/
/*  Return Value:  0														*/
/*--------------------------------------------------------------------------*/


int q1slct_( Int_T *n, Int_T *l, Int_T *m, Int_T *mg,  Real_T *g, // )
	Real_T *a, Int_T *iblock, Int_T *ibasic,  Int_T *inonba, Int_T *icheck,
	Int_T *ieq, Int_T *status, Int_T *drow )
{
	/* System generated locals */
	Int_T g_dim1, g_offset;

	/* Local variables */
	Int_T iold, inew, i;

	/* Parameter adjustments */
	g_dim1 = *n;
	g_offset = Int_T( g_dim1 + 1 );
	g -= g_offset;
	--a;
	--iblock;
	--ibasic;
	--inonba;
	--icheck;
	--ieq;
	--status;
	--drow;

	/* Function Body */
	for (i = 1; i <= *mg; ++i)
		icheck[i] = 0;
	for (i = 1; i <= *l; ++i)
		icheck[ibasic[i]] = Int_T( -i );
	if (*m < 1)
		goto L20;

	for (i = 1; i <= *m; ++i)
		icheck[inonba[i]] = i;
L20:
	inew = 1;
	for (iold = 1; iold <= *mg; ++iold)
	{
		i = icheck[iold];
		if (i == 0) {
			if ( (drow[iold] > 0) && (status[drow[iold]] == 1) ) {
				status[drow[iold]] = 0;
			}
			else if ( (drow[iold] < 0) && (status[-drow[iold]] == -1) ) {
				status[-drow[iold]] = 0;
			}
			drow[iold] = 0;
			continue;
		}
		if (iold != inew)
		{
			dcopy_(*n, &g[iold * g_dim1 + 1], &g[inew * g_dim1 + 1]);
			a[inew] = a[iold];
			iblock[inew] = iblock[iold];
			ieq[inew] = ieq[iold];
			drow[inew] = drow[iold];
			if( i <= 0 )
			{
				i = Int_T( -i );
				ibasic[i] = inew;
			}
			else
				inonba[i] = inew;
		}

		icheck[inew] = -2;
		++inew;
	}
	return 0;
}
/* -----END OF Q1SLCT----------------------------------------------------- */


void CheckConstraintsInSection( Int_T *marks, Int_T j1, Int_T j2, // )
	Int_T *iptr, Real_T *dmat, Int_T *jcol, Real_T *y, Int_T *status,
	Real_T *bmin, Real_T *bmax, Real_T tolcut, Real_T &gsec, Int_T &imax )
{
	for( Int_T i = j1; i <= j2; ++i )
	{
		//------------------------------------------------------------------
		//	CHECK THE I-TH CONSTRAINT, IF PREVIOUSLY VIOLATED.
		//
		if( marks[i] <= 0 )
		{
			marks[i] = 0;
			continue;
		}

		if( status[i] == 2 ) continue;

		Int_T ii = iptr[i];

		assert( ii > 0 );

		Real_T bi = ddots_( Int_T( iptr[i+1] - ii ), dmat+ii, jcol+ii, y+1 );

		if( status[i] != 1 )
		{
			Real_T btem = bi - bmax[i];

			if( btem > tolcut * ( fabs(bi) + 1.0 ) )
			{
				if( btem <= gsec )
					continue;

				gsec = btem;
				imax = i;
				continue;
			}
		}

		if( status[i] != -1 )
		{
			Real_T btem = bmin[i] - bi;

			if( btem <= tolcut * ( fabs(bi) + 1.0 ) || btem <= gsec )
				continue;

			gsec = btem;
			imax = Int_T( -i );
		}
	}
}


Bool_T CheckAllConstraints( Int_T &j1, Int_T &j2, Int_T mdmat, Int_T n, // )
	Int_T nsec, Int_T &mg, Int_T *iblock, Int_T *icheck, Int_T g_dim1,
	Real_T *g, Int_T *drow, Real_T *a, Int_T *ieq, Int_T *marks,
	Int_T *iptr, Real_T *dmat, Int_T *jcol, Real_T *y, Int_T *status,
	Real_T *bmin, Real_T *bmax, Real_T tolcut, Int_T nwcmax,
	Real_T &gmax, Int_T &inew, Int_T &jnew, Int_T &ifeas1 )
{
	const Int_T mgold = mg;

	for( Int_T isec = 1; isec <= nsec; ++isec )
	{
		//----------------------------------------------------------------------
		//	SET THE POINTERS FOR THE CURRENT SECTION.
		//
		j1 = j2 % mdmat + 1;
		j2 = j1 + n - 1;
		if( mdmat - j2 < n ) j2 = mdmat;

		Int_T imax = 0, irun = 0;
		Real_T gsec = 0.0;

		//----------------------------------------------------------------------
		//	GO THROUGH SECTION.
		//
		for(;;)
		{
			CheckConstraintsInSection( marks, j1, j2, iptr, dmat, jcol, y,
				status, bmin, bmax, tolcut, gsec, imax );

			if( gsec > tolcut )	break;
			if( irun > 0 )		break;

			irun = 1;

			//------------------------------------------------------------------
			//	VERIFY CONSTRAINTS SKIPPED IN THE FIRST PASS.
			//
			for( Int_T i = j1; i <= j2; ++i )		//these were the ones that had marks[i]<0
				if( marks[i] == 0 ) marks[i] = 1;
		}
		if( irun > 0 )		continue;

		//----------------------------------------------------------------------
		//	APPEND THE FEASIBILITY CUT TO THE COMMITTEE.
		//
		++mg;
		iblock[mg] = 0;
		icheck[mg] = 1;
		dzero_( n, &g[mg * g_dim1 + 1] );

		Int_T ii;
        //Print("imax = %d    gsec = %G\n",imax,gsec);
		//Print( "Infeasibility Cut\n");
		if( imax > 0 )          //then violation is in the "...>=bmax" side
		{
			status[imax] = 1;   //imax gives the constraint index with the max violation
			drow[mg] = imax;
			a[mg] = -bmax[imax];
			ii = iptr[imax];
			dunpk_( iptr[imax+1]-ii, dmat+ii, jcol+ii, &g[mg*g_dim1+1] );
		}
		else
		{
			imax = Int_T( -imax );
			status[imax] = -1;
			drow[mg] = Int_T( -imax );
			a[mg] = bmin[imax];
			ii = iptr[imax];
			dunne_( iptr[imax+1]-ii, &dmat[ii], &jcol[ii], &g[mg*g_dim1+1] );
		}

		if( IsEqual( bmin[imax], bmax[imax] ) )
		{
			ieq[mg] = 1;
			status[imax] = 2;
		}
		else
			ieq[mg] = 0;

		if( gsec > gmax )
		{
			gmax = gsec;
			inew = mg;
			jnew = 0;
			ifeas1 = 0;
		}

		//----------------------------------------------------------------------
		//	IF SUFFICIENTLY MANY CUTS APPENDED, EXIT FROM THE LOOP.
		//
		if( mg >= mgold + nwcmax ) break;
	}

	return ( mg < mgold + nwcmax ) ? True : False;
}


void SolveBlocks( Real_T &f1, Real_T &f2, Real_T &fx, Int_T &mg, Int_T l, // )
	Int_T n, Real_T *y, Real_T &value, Real_T *g, Real_T *weight, Int_T &ifeas2,
	Int_T &iphase, Int_T *iblock, Int_T *icheck, Real_T tolcut, Real_T *a,
	Int_T *ieq, Int_T &inew, Int_T &jnew, Real_T &gmax, Int_T g_dim1, Real_T *v,
	Real_T *x, RD_SubproblemManager &SubMan, Int_T *istop, Real_T *ExpC, Real_T& ExpCost, 
	Real_T *ExpC2, Real_T& ExpCost2 )
{
	f1 = f2 = fx = 0.0;

	++mg;

	for( Int_T i = 1; i <= l; ++i )
	{
		Real_T gi;

		if( !SubMan.SolveSubproblem( i, n, y+1, value, g+(mg*g_dim1+1) ) )
		{
			//	INDUCED CONSTRAINT.
			//
			f1 += value * weight[i];
			ifeas2 = 0;

			if( iphase == 1 )
			{
				icheck[mg]	= 1;
				gi			= value - v[i];
				iblock[mg]	= i;
			}
			else
			{
				icheck[mg]	= 1;
				gi			= value;
				iblock[mg]	= -i;
			}
		}
		else
		{
			//	OBJECTIVE CUT.
			//
			//Print ( "Objective Cut\n" ); 
			f2 += value * weight[i];

			if( iphase == 2 )
			{
				if (value <= v[i] + tolcut) continue;

				icheck[mg]	= 1;
				gi			= value - v[i];
				iblock[mg]	= i;
			}
			else
			{
				icheck[mg]	= -1;
				gi			= -v[i];
				iblock[mg]	= i;
			}
		}

		a[mg]	= value - ddot_( n, g+(mg*g_dim1+1), y+1 );
		ieq[mg]	= 0;
		fx		+= weight[i] * ( a[mg] + ddot_( n, g+(mg*g_dim1+1), x+1 ) );
		++mg;


		//@BEGIN--------------------------------------------------------------
		//		 This is to store the "value"s and "f2" for the 
		//			the Single Replication Procedure for testing solution
		//			quality... given xhat...=
		//			(assuming we only have objective cuts...)
		
		if (*istop == 10){
			ExpC[i-1] = value; 
			if( i == l)	ExpCost = f2; 
		}
		
		if (*istop == 11){
			ExpC2[i-1] = value; 
			if( i == l)	ExpCost2 = f2; 
		}
		//@END----------------------------------------------------------------



		if( gi <= gmax && inew != 0 ) continue;

		inew = mg - 1;
		jnew = 0;
		gmax = gi;
	}
}


void DetermineStepType( Int_T &iphase, Int_T ifeas1, Int_T ifeas2, // )
	Int_T &ifeas, Real_T f1, Real_T f2, Real_T vsum,
	Int_T mg, Int_T l, Int_T n, Int_T g_dim1, Int_T *icheck, Int_T *iblock,
	Real_T *g, Real_T *a, Int_T iter, Int_T &index, Real_T &funold,
	Real_T &funmin, Int_T &istart, Int_T &iser, Int_T levprt, Real_T *y,
	Real_T *x, Real_T *pi, Real_T *yb, Int_T *istat, Int_T *irn,
	Real_T &penlty, Real_T initpen, Bool_T &newpen, 
	Real_T penmax, Int_T &nulinf, Int_T &nserex, Int_T &nulf,
	Int_T &nserap, Real_T tolcut, Real_T gamma, Int_T inew, Int_T jnew, bool lshaped)
{
	//const  ::omega1 and 2 used to be constants, but I am changing that! 
	Real_T omega1 = 0.5, omega2 = 1.5; 
	//Used to be 0.7 and 1.3 resp.  -->this is a comment from before. 

	if(lshaped){
		omega1 = 1; 
		omega2 = 1; 
	//To run L-shaped fron scratch, set Decompoptions:DecompOptions InitPen (1e-6) in main.h
	}


	Real_T fun;

	/*    CONTINUE PHASE 1. */
	if( iphase != 2 && ifeas1 <= 0 )
	{
		ifeas = ifeas1;

		/*    FOR MORE FREEDOM IN GENERATING AND CHOOSING INDUCED CONSTRAINTS */
		/*    FORCE F1 TO LOOK LIKE A CONVEX FUNCTION. */

		for( Int_T i = mg - l + 1; i <= mg; ++i )
		{
			if( icheck[i] > 0 ) continue;

			dzero_( n, &g[i * g_dim1 + 1] );
			a[i] = 0.0;
		}

		fun = Max( f1, vsum );
		if( iter > 1 )
			goto L81;

		index = 1;
		funold = funmin = fun;
		istart = ifeas1;
		iser = 1;
		return;
	}

	/*     INITIALIZE PHASE 2. */
	if( iphase != 2 && ifeas2 > 0 )
	{
		if( levprt > 0 )
			Print( "FEASIBLE POINT FOUND\n" );

		iphase = 2;
		funold = f2;
		ifeas = 1;
		istart = 1;
		dcopy_( n, y+1, x+1 );
		index = 2;

		Int_T i;
		for( i = 1; i <= n; ++i ) istat[i] = 2;
		for( i = 1; i <= n; ++i ) irn[i] = i;

		Int_T ii = mg - l;
		if( ii < 1 )
			return;

		for( i = 1; i <= ii; ++i ) iblock[i] = -iblock[i];
		for( i = 1; i <= ii; ++i ) icheck[i] = 1;
		penlty = initpen;
		newpen = False;
		iser = 1;
		return;
	}

	//	CONTINUE PHASE 2.
	ifeas = ifeas1 * ifeas2;
	//fun = f2;
		fun = Max( f2, vsum );

	//	MAKE A NULL OR SERIOUS STEP AND CORRECT THE PENALTY.

L81:
	if( ifeas <= 0 )
	{
		//	INFEASIBLE POINT.

		iser = -1;


		penlty *= omega2;
		newpen = True;
		++nulinf;
		return;
	}

	if( istart <= 0 )
	{
		//	THE FIRST POINT SATISFYING DIRECT CONSTRAINTS.

		istart = 1;
		iser = 1;
		++nserex;
		MakeStep( n, x, y, yb, pi, funold, fun, index );
		return;
	}

	if( iphase != 2 )
	{
		iser = 0;
		if( vsum > tolcut )
		{
			penlty *= omega1;
			//penlty = Max(penlty,1.0e-7);
			newpen = True;
		}

		if( fun > funmin )
			return;

		funmin = fun;
		iser = 1;
		dcopy_( n, y+1, x+1 );
		funold = fun;
		index = 4;

		return;
	}

	if( fun - vsum <= (1.0 - gamma) * ( funold - vsum ) )
	{
		newpen = True;
		penlty *= omega1;
		//penlty = Max(penlty,1.0e-7);

		//	EXACT PREDICTION.

		++nserex;
		iser = 2;
		newpen = True;
		dcopy_( n, y+1, x+1 );
		funold = fun;
		index = 4;
		return;
	}

	if ( fun < funold || ( inew == 0 && jnew == 0 ) )
		goto L87;

	/*       LONG STEPS. INCREASE PENALTY. */
	if( penlty < penmax )
	{
		penlty *= omega2;
		newpen = True;
	}

	++nulf;
	iser = 0;
	return;

L87:
	if( fun - vsum > gamma * (funold - vsum) && iphase == 2 )
		if ( inew != 0 || jnew != 0 )
		{
			++nulf;
			iser = 0;
			return;
		}

	/*     SERIOUS STEP. */
	++nserap;
	iser = 1;
	MakeStep( n, x, y, yb, pi, funold, fun, index );
}



//CountCriticalScenarios is not used anywhere! 

void CountCriticalScenarios( Int_T l, Int_T m, Int_T nfix, Int_T *iblock, // )
	Int_T *inonba )
{
	//--------------------------------------------------------------------------
	//	count critical scenarios:
	//	nctot is the total number of critical scenarios
	//	ncsum is the total amount of criticality 
	//
	Array<Int_T> ncrit( l+1, 0 );
	Int_T ncsum	= 0,
		nctot	= 0;

	for( Int_T i = 1; i <= m; i++)
	{
		Int_T ii = iblock[inonba[nfix+i]];
		if( ii == 0 ) continue;

		if( ii < 0 ) ii = -ii;
		if( ncrit[ii] == 0 )
		{
			nctot++;
			ncsum++;
		}
		ncrit[ii]++;
	}
	Print( "Critical scenarios: %6d, Degree of criticality: %6d\n",
		nctot, ncsum );
}


void MakeStep( Int_T n, Real_T *x, Real_T *y, Real_T *yb, Real_T *pi, // )
	Real_T &funold, Real_T fun, Int_T &index )
{
	dcopy_( n, y+1, pi+1 );
	ddifr_( n, pi+1, x+1 );
	dsuma_( n, yb+1, pi+1 );
	dcopy_( n, y+1, x+1 );
	funold = fun;
	index = 4;
}


/*--------------------------------------------------------------------------*/
/*			 CompressCommittee												*/
/*  Purpose: TO COMPRESS THE COMMITTEE BY DELETING INACTIVE MEMBERS.		*/
/*																			*/
/*  Called by:	Q1CMTE														*/
/*  Subroutines Called: Q1SLCT												*/
/*--------------------------------------------------------------------------*/


void CompressCommittee( Int_T &mg, Int_T n, Int_T l, Int_T m, Int_T nfix, // )
	Int_T *istat, Int_T *inonba, Int_T *icheck, Int_T *ieq, Int_T *status,
	Int_T *drow, Real_T *g, Real_T *a, Int_T *iblock, Int_T *ibasic,
	Int_T g_offset )
{
	Int_T nmem = mg, i;
	for( i = 1; i <= n; ++i )
		if( istat[i] < 2 )
			++nmem;

	if( nmem <= n + l )
		return;

	q1slct_( &n, &l, &m, &mg, g+g_offset, a+1, iblock+1, ibasic+1,
		inonba+nfix+1, icheck+1, ieq+1, status+1, drow+1 );
	for( i = 1; i <= n; ++i )
		if( istat[i] == 1 )
			istat[i] = 2;
	Print("Compress Committee...\n"); 
	mg = m + l;
}
