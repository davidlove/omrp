/*------------------------------------------------------------------------------
MODULE TYPE:		Linear programming, linear problem presolving.
PROJECT CODE:		PRESOLVER
PROJECT FULL NAME:	Linear problem presolver implementation.

MODULE AUTHOR:		Artur Swietanowski.

PROJECT SUPERVISOR:	dr Jacek Gondzio, prof. Andrzej P. Wierzbicki

--------------------------------------------------------------------------------

SOURCE FILE NAME:	gnu_fix.cpp
CREATED:			1994.01.12
LAST MODIFIED:		1996.10.07

DEPENDENCIES:		stdtype.h, std_tmpl.h, smartptr.h, error.h, memblock.h,
					smartdcl.h, sptr_deb.h, sptr_ndb.h, work_vec.h, vec_pool.h,
					myalloc.h

--------------------------------------------------------------------------------

SOURCE FILE CONTENTS:
	Some code necessary for the GNU C++ compiler v. 2.6.x to compile and link
the project properly. All templates used in the project are "declared" here.

------------------------------------------------------------------------------*/

#if defined( explicit_templates )

#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif

#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SMARTPTR_H__
#	include "smartptr.h"
#endif
#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif
#ifndef __MYALLOC_H__
#	include "myalloc.h"
#endif


//------------------------------------------------------------------------------
//	File "smartptr.h" templates' instatiations.
//
template class SmartPointerBase<Int_T>;
template class Array<Int_T>;
template class Ptr<Int_T>;

template class SmartPointerBase<Real_T>;
template class Array<Real_T>;
template class Ptr<Real_T>;

template class SmartPointerBase<Bool_T>;
template class Array<Bool_T>;
template class Ptr<Bool_T>;

template class SmartPointerBase<Short_T>;
template class Array<Short_T>;
template class Ptr<Short_T>;

template class SmartPointerBase<unsigned int>;
template class Array<unsigned int>;
template class Ptr<unsigned int>;

//------------------------------------------------------------------------------
//	File "std_tmpl.h" templates' instatiations.
//
template Real_T Max( Real_T, Real_T );
template Int_T Max( Int_T, Int_T );
template unsigned int Max( unsigned int, unsigned int );

template Real_T Min( Real_T, Real_T );
template Int_T Min( Int_T, Int_T );

//------------------------------------------------------------------------------
//	File "myalloc.h" templates instatiations.
//
template Short_T *MALLOC( Short_T *& Table, size_t len );
template Short_T *REALLOC( Short_T *& Table, size_t len );
template void FREE( Short_T *& Table );

template Real_T *MALLOC( Real_T *& Table, size_t len );
template Real_T *REALLOC( Real_T *& Table, size_t len );
template void FREE( Real_T *& Table );

template Int_T *MALLOC( Int_T *& Table, size_t len );
template Int_T *REALLOC( Int_T *& Table, size_t len );
template void FREE( Int_T *& Table );

template Bool_T *MALLOC( Bool_T *& Table, size_t len );
template Bool_T *REALLOC( Bool_T *& Table, size_t len );
template void FREE( Bool_T *& Table );

template unsigned int *MALLOC( unsigned int *& Table, size_t len );
template unsigned int *REALLOC( unsigned int *& Table, size_t len );
template void FREE( unsigned int *& Table );

template void FREE( char *& Table );

//------------------------------------------------------------------------------
//	File "work_vec.h" templates' instatiations.
//
template class WorkVector<Bool_T>;
template class WorkVector<Int_T>;
template class WorkVector<Real_T>;

#endif

#if defined( explicit_templates )

#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif

int WorkVector<Bool_T>::VecType = WorkVectorPool::NO_TYPE;
int WorkVector<Int_T>::VecType = WorkVectorPool::NO_TYPE;
int WorkVector<Real_T>::VecType = WorkVectorPool::NO_TYPE;

#endif
