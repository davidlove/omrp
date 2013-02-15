//------------------------------------------------------------------------------
//	Explicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )

#include <stdlib.h>

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

template class Array<long>;
template class SmartPointerBase<long>;
template class Ptr<long>;

//------------------------------------------------------------------------------
//	File "std_tmpl.h" templates' instatiations.
//
template Real_T Max( Real_T, Real_T );
template Int_T Max( Int_T, Int_T );
template Long_T Max( Long_T, Long_T );
template unsigned int Max( unsigned int, unsigned int );
template size_t Max( size_t, size_t );

template Real_T Min( Real_T, Real_T );
template Int_T Min( Int_T, Int_T );

template Real_T Sqr( const Real_T );

template void Swap( Real_T &, Real_T & );
template void Swap( Int_T &, Int_T & );

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

template long *MALLOC( long *& Table, size_t len );
template long *REALLOC( long *& Table, size_t len );
template void FREE( long *& Table );

template void FREE( char *& Table );

//------------------------------------------------------------------------------
//	File "work_vec.h" templates' instatiations.
//
template class WorkVector<Bool_T>;
template class WorkVector<Int_T>;
template class WorkVector<Short_T>;
template class WorkVector<Real_T>;

#endif

#if defined( explicit_templates )

#ifndef __WORK_VEC_H__
#	include "work_vec.h"
#endif

int WorkVector<Bool_T>::VecType = WorkVectorPool::NO_TYPE;
int WorkVector<Int_T>::VecType = WorkVectorPool::NO_TYPE;
int WorkVector<Short_T>::VecType = WorkVectorPool::NO_TYPE;
int WorkVector<Real_T>::VecType = WorkVectorPool::NO_TYPE;


#ifndef __SC_TREE_H__
#	include "sc_tree.h"
#endif
#ifndef __SOLV_LP_H__
#	include "solv_lp.h"
#endif
#ifndef __PERIODS_H__
#	include "periods.h"
#endif
#ifndef __CHANGELP_H__
#	include "changelp.h"
#endif
	template class ExtendArray<ChangeLP *>;
	template class Array<ChangeLP *>;
	template class SmartPointerBase<ChangeLP *>;
	template class Ptr<ChangeLP *>;

	template class ExtendArray<ExtendArray<ChangeLP *> *>;
	template class Array<ExtendArray<ChangeLP *> *>;
	template class SmartPointerBase<ExtendArray<ChangeLP *> *>;
	template class Ptr<ExtendArray<ChangeLP *> *>;

	template class ExtendArray<ScenarioTree::ScStruct>;
	template class Array<ScenarioTree::ScStruct>;
	template class SmartPointerBase<ScenarioTree::ScStruct>;
	template class Ptr<ScenarioTree::ScStruct>;

	template class ExtendArray<ExtendArray<ScenarioTree::ScStruct> *>;
	template class Array<ExtendArray<ScenarioTree::ScStruct> *>;
	template class SmartPointerBase<ExtendArray<ScenarioTree::ScStruct> *>;
	template class Ptr<ExtendArray<ScenarioTree::ScStruct> *>;

	template ChangeLP **MALLOC( ChangeLP **& Table, size_t len );
	template ChangeLP **REALLOC( ChangeLP **& Table, size_t len );
	template void FREE( ChangeLP **& Table );

	template ExtendArray<ChangeLP *> **MALLOC( ExtendArray<ChangeLP *> **& Table, size_t len );
	template ExtendArray<ChangeLP *> **REALLOC( ExtendArray<ChangeLP *> **& Table, size_t len );
	template void FREE( ExtendArray<ChangeLP *> **& Table );

	template ScenarioTree::ScStruct *MALLOC( ScenarioTree::ScStruct *& Table, size_t len );
	template ScenarioTree::ScStruct *REALLOC( ScenarioTree::ScStruct *& Table, size_t len );
	template void FREE( ScenarioTree::ScStruct *& Table );

	template ExtendArray<ScenarioTree::ScStruct> ** MALLOC( ExtendArray<ScenarioTree::ScStruct> ** & Table, size_t len );
	template ExtendArray<ScenarioTree::ScStruct> ** REALLOC( ExtendArray<ScenarioTree::ScStruct> ** & Table, size_t len );
	template void FREE( ExtendArray<ScenarioTree::ScStruct> ** & Table );

#endif
