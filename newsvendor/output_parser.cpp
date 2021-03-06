#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <vector>

using namespace std;

int main( int argc, char** argv )
{
   int            batchSize,
                  gamma,
                  numBatches;
   double         bias,
                  coverage,
                  avgVar,
                  varVar,
                  widthCI,
                  probOverlap,
                  runTime;
   char           c,
                  inputLine[81];
   bool           fileError = false,
                  appendTag = false;
   string         outputFileName;
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
      // Skip the first 3 lines of input
      for( int ii = 0; ii < 3; ii++ )
      {
         currentFile.getline( inputLine, 80, '\n' );
      }
      // Get the batch size
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Batch Size        :    m  = %d", &batchSize );
      // Get the batch nonoverlap
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Batch Non-overlap : gamma = %d", &gamma );
      // Get the number of batches
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Number Batches    :       =  %d", &numBatches );
      // Skip the next three lines of input
      for( int ii = 0; ii < 3; ii++ )
      {
         currentFile.getline( inputLine, 80, '\n' );
      }
      // Get the CI width
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Width of CI               =  %lf", &widthCI );
      // Get the bias
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Bias                      =  %lf", &bias );
      // Get the average variance
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Average Variance          =  %lf", &avgVar );
      // Get the variance of the variance
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Variance of VG            =  %lf", &varVar );
      // Get the probability of solution overlap
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Prob of Solution Overlap  =  %lf", &probOverlap );
      // Get the coverage probability -- one blank line before getting the good stuff
      currentFile.getline( inputLine, 80, '\n' );
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "%lf", &coverage ); 
      // Get the solution time -- one blank line before
      currentFile.getline( inputLine, 80, '\n' );
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Time to run the program is :%lf seconds.", &runTime ); 
      // Close the file
      currentFile.close();
      currentFile.clear();

      if( !fileIndex && !appendTag )
      {
         outputFile << "% Newsvendor Problem with m = " << batchSize << endl;
         outputFile << "\%Gamma" 
            << '\t' << "Prob Soln Overlap"
            << '\t' << "Cover Prob"
            << '\t' << "Width of CI" 
            << '\t' << "Bias"
            << '\t' << "Average of Var"
            << '\t' << "Variance of Var"
            << '\t' << "Run Time"
            << '\t' << "Run Time/Batch"
            << endl;
      }
      outputFile << gamma
         << '\t' << showpoint << setprecision(7) << probOverlap
         << '\t' << setprecision(4) << coverage
         << '\t' << setprecision(4) << widthCI
         << '\t' << setprecision(2) << bias
         << '\t' << setprecision(7) << avgVar
         << "\t" << varVar
         << "\t" << runTime
         << "\t" << runTime / numBatches
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
