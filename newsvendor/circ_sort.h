#ifndef __CIRCSORT_H__
#define __CIRCSORT_H__

// #include <stdio>
// #include <stdlib>
// #include <math>
#include <iostream>
#include <fstream>
#include <algorithm>
// #include "main.h"
//
using namespace std;

class compare_circ_class
{
   public:
      bool operator()( double* a, double* b )
      {
         return ( *a < *b );
      }
};

class circ_sort
{
   public:
      circ_sort( void );
      circ_sort( int arraySize );
      ~circ_sort( void );
      bool Initialize( int arraySize );
      bool AddEntry( double newEntry );
      // Return the index'th smallest value, i.e., from the sorted array
      double operator[]( size_t index );
      // Return unsorted values (i.e., according to the order they were added)
      double operator()( size_t index );
      // int Compare( const void* arg1, const void* arg2 )
      // static bool Compare( double* a, double* b )
      // {
         // return ( *a < *b );
      // }
      void Sort( void );
      void PrintQueue( void );
      void PrintQueue( ofstream& debugFile );
      void PrintSorted( void );
      void PrintSorted( ofstream& debugFile );

   private:
      int                  size,
                           nextEntryLocation;
      double            *  circBuffer,
                        ** sortedBuffer;
      compare_circ_class   Compare;
};

#endif
