#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <vector>

#define LINE_LENGTH 800

using namespace std;

int main( int argc, char** argv )
{
   int            batchSize,
                  numBatches,
                  degreesFreedom,
                  gamma,
                  time;
   double         coverage, 
                  // bias,
                  avgVar,
                  varVar,
                  widthCI,
                  ciDenom,
                  varWidthCI,
                  sizeCI,
                  varSizeCI;
                  // probOverlap;
   char           c,
                  inputLine[LINE_LENGTH + 1];
   bool           fileError = false,
                  appendTag = false,
                  firstLoop = true;
   string         outputFileName,
                  firstWordTest;
   vector<string> inputFileName;
   ifstream       currentFile; 
   ofstream       outputFile;

   while ((c = getopt (argc, argv, "ao:")) != -1) 
   {
      switch (c) 
      { 
         case 'a':
            appendTag = true;
            break;
         case 'o':
            outputFileName = optarg;
            break;
      }
   }
   for( int ii = optind; ii < argc; ii++ )
   {
      inputFileName.push_back( argv[ii] );
   }
   if( outputFileName.empty() )
   {
      cerr << "Must specify an output file name" << endl;
      abort();
   }
   if( appendTag )
      outputFile.open( outputFileName.c_str(), ofstream::app );
   else
      outputFile.open( outputFileName.c_str() );
   if( !outputFile ) {
      cerr << "Could not open output file " << outputFileName << endl;
      abort();
   }

   for( vector<string>::size_type fileIndex = 0; fileIndex < inputFileName.size(); fileIndex++ )
   {
      currentFile.open( inputFileName[fileIndex].c_str() );
      if( !currentFile )
      {
         cerr << "Could not open file " << argv[fileIndex] << endl;
         fileError = true;
         break;
      }
      cout << "File: " << inputFileName[fileIndex] << ' ';
      do
      {
         if( currentFile.eof() )
         {
            cout << "End of file " << inputFileName[fileIndex] << "." << endl;
            return 0;
         }
         currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      }
      while( strcmp( inputLine, "BEGIN MRP" ) && !currentFile.eof() );
      if( currentFile.eof() )
      {
         cout << "End of file " << inputFileName[fileIndex] << "." << endl;
         currentFile.close();
         currentFile.clear();
         continue;
      }
      printf( "%s\n", inputLine );
      // Skip the first 3 lines of input
      // Get the batch size
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Batch Size = %d", &batchSize );
      // Get the batch nonoverlap
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "gamma = %d", &gamma );
      // Get the CI width
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "CI Width = %lf", &widthCI );
      // Get the bias
      //currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      //sscanf( inputLine, "Bias                      =  %lf", &bias );
      // Get the coverage probability -- one blank line before getting the good stuff
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Coverage, MRP: = %lf", &coverage ); 
      // Get the average variance
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Average Variance =  %lf", &avgVar );
      // Get the variance of the variance
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Variance of Variance =  %lf", &varVar );
      
      // NEW STUFF -- Enhanced look at CI widths
      // Get the number of batches
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Number Batches = %d\n", &numBatches);
      // Get degrees of freedom
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Degrees Freedom = %d\n", &degreesFreedom);
      // Get denominator of CI (i.e., n/m)
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "CI Denom = %lf\n", &ciDenom);
      // Get variance of the CI Width
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Var CI Width = %lf\n", &varWidthCI);
      // Get the size of the CI (i.e., z_m* + widthCI)
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "CI Sze = %lf\n", &sizeCI);
      // Get variance of CI size
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Var CI Size = %lf\n", &varSizeCI);
      
      // Get computation time
      currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      sscanf( inputLine, "Time = %d\n", &time );

      // Get the probability of solution overlap
      //currentFile.getline( inputLine, LINE_LENGTH, '\n' );
      //sscanf( inputLine, "Prob of Solution Overlap  =  %lf", &probOverlap );
      // Close the file
      currentFile.close();
      currentFile.clear();

      if( firstLoop && !appendTag )
      {
         outputFile << "% Newsvendor Problem with m = " << batchSize << endl;
         outputFile << "\%Gamma" 
            //<< '\t' << "Prob Soln Overlap"
            << "\t" << "Cover Prob"
            << "\t" << "Width of CI" 
            //<< '\t' << "Bias"
            << "\t" << "Average of Var"
            << '\t' << "Variance of Var"
            << '\t' << "Num Batches"
            << '\t' << "Degrees Freedom"
            << '\t' << "CI Denominator"
            << '\t' << "Var of CI Width"
            << '\t' << "Size of CI"
            << '\t' << "Var of CI Size"
            << "\t" << "Computation Time"
            << endl;
      }
      firstLoop = false;
      outputFile << gamma
         //<< '\t' << showpoint << setprecision(7) << probOverlap
         << "\t" << showpoint << setprecision(4) << coverage
         << "\t" << setprecision(4) << widthCI
         //<< "\t" << setprecision(2) << bias
         << "\t" << setprecision(7) << avgVar
         << "\t" << varVar
         << "\t" << numBatches
         << "\t" << degreesFreedom
         << "\t" << ciDenom
         << "\t" << varWidthCI
         << "\t" << sizeCI
         << "\t" << varSizeCI
         << "\t" << time
         << endl;
   }
   
   if( fileError )
   {
      exit(1);
   }

   outputFile.close();
   outputFile.clear();

   return 0;
}
