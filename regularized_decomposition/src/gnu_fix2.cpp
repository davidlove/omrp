#ifndef __STD_TMPL_H__
#	include "std_tmpl.h"
#endif
#ifndef __SCENARIO_H__
#	include "scenario.h"
#endif

//------------------------------------------------------------------------------
//	Implicit template instantiation (for GNU C++  ver. 2.6.2 or later only).
//
#if defined( explicit_templates )

//------------------------------------------------------------------------------
//	Implementation of some function templates from the "std_tmpl.h".
//

template void Fill( Real_T *array, Int_T len, const Real_T &val );
template void Fill( Int_T *array, Int_T len, const Int_T &val );

//
//	Array<Delta *>
//
template class SmartPointerBase<Delta *>;
template class Array<Delta *>;
template class Ptr<Delta *>;

template Delta **MALLOC( Delta **& Table, size_t len );
template Delta **REALLOC( Delta **& Table, size_t len );
template void FREE( Delta **& Table );

//
//	Array<Distribution *>
//
template class SmartPointerBase<Distribution *>;
template class Array<Distribution *>;
template class Ptr<Distribution *>;

template Distribution **MALLOC( Distribution **& Table, size_t len );
template Distribution **REALLOC( Distribution **& Table, size_t len );
template void FREE( Distribution **& Table );

//
//	Array<StochasticDataBlock *>
//
template class SmartPointerBase<StochasticDataBlock *>;
template class Array<StochasticDataBlock *>;
template class Ptr<StochasticDataBlock *>;

template StochasticDataBlock **MALLOC( StochasticDataBlock **& Table,
	size_t len );
template StochasticDataBlock **REALLOC( StochasticDataBlock **& Table,
	size_t len );
template void FREE( StochasticDataBlock **& Table );

//
//	Array<Scenario *>
//
template class SmartPointerBase<Scenario *>;
template class Array<Scenario *>;
template class Ptr<Scenario *>;

template Scenario **MALLOC( Scenario **& Table, size_t len );
template Scenario **REALLOC( Scenario **& Table, size_t len );
template void FREE( Scenario **& Table );

//
//	Array<Array<Real_T> *>
//
template class SmartPointerBase<Array<Real_T> *>;
template class Array<Array<Real_T> *>;
template class Ptr<Array<Real_T> *>;

template Array<Real_T> **MALLOC( Array<Real_T> **& Table, size_t len );
template Array<Real_T> **REALLOC( Array<Real_T> **& Table, size_t len );
template void FREE( Array<Real_T> **& Table );
//
//	End of explicit template instantiation.
//------------------------------------------------------------------------------

#endif
