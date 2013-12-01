#include "circ_sort.h"

circ_sort::circ_sort( void )
{
   size = 0;
   nextEntryLocation = 0;
   circBuffer = NULL;
   sortedBuffer = NULL;

   return;
}

circ_sort::circ_sort( int arraySize )
{
   size = 0;
   nextEntryLocation = 0;
   circBuffer = NULL;
   sortedBuffer = NULL;
   Initialize( arraySize );

   return;
}

bool circ_sort::Initialize( int arraySize )
{
   if( circBuffer )
   {
      return false;
   }

   size = arraySize;
   nextEntryLocation = 0;
   circBuffer = new double[size];
   sortedBuffer = new double*[size];

   for( int ii = 0; ii < size; ii++ )
   {
      sortedBuffer[ii] = circBuffer + ii;
   }

   return true;
}

circ_sort::~circ_sort( void )
{
   if( circBuffer )
      delete [] circBuffer;
   if( sortedBuffer )
      delete [] sortedBuffer;
   size = 0;

   return;
}

bool circ_sort::AddEntry( double newEntry )
{
   circBuffer[nextEntryLocation++] = newEntry;
   if( nextEntryLocation >= size )
      nextEntryLocation = 0;

   return true;
}

double circ_sort::operator[]( size_t index )
{
   return *(sortedBuffer[index]);
}

double circ_sort::operator()( size_t index )
{
   switch( nextEntryLocation )
   {
      case 0:
         return circBuffer[nextEntryLocation + index];
      default:
         return circBuffer[ (nextEntryLocation + index) % size ];
   }

   return 0.0;
}

void circ_sort::Sort( void )
{
   // qsort( sortedBuffer, (size_t) size, sizeof( double* ), Compare );
   sort( sortedBuffer, sortedBuffer + size, Compare );

   return;
}

void circ_sort::PrintQueue( void )
{
   for( int ii = 0; ii < size; ii++ )
      cout << "circBuffer[" << ii << "] = " << circBuffer[ii] << endl;
   
   return;
}

void circ_sort::PrintQueue( ofstream& debugFile )
{
   for( int ii = 0; ii < size; ii++ )
      debugFile << "circBuffer[" << ii << "] = " << circBuffer[ii] << endl;
   
   return;
}

void circ_sort::PrintSorted( void )
{
   for( int ii = 0; ii < size; ii++ )
      cout << "sortedBuffer[" << ii << "] = " << *(sortedBuffer[ii]) << endl;
   
   return;
}

void circ_sort::PrintSorted( ofstream& debugFile )
{
   for( int ii = 0; ii < size; ii++ )
      debugFile << "sortedBuffer[" << ii << "] = " << *(sortedBuffer[ii]) << endl;
   
   return;
}
