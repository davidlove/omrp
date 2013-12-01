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
                  gamma;
   double         bias,
                  coverage,
                  avgVar,
                  varVar,
                  widthCI,
                  probOverlap;
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
      // Skip the first 4 lines of input
      for( int ii = 0; ii < 4; ii++ )
      {
         currentFile.getline( inputLine, 80, '\n' );
      }
      // Get the batch size
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Batch Size        :    m  = %d", &batchSize );
      cout << inputLine << endl;
      // Get the batch nonoverlap
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "Batch Non-overlap : gamma = %d", &gamma );
      // Skip the next four lines of input
      for( int ii = 0; ii < 4; ii++ )
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
      // currentFile.getline( inputLine, 80, '\n' );
      // sscanf( inputLine, "Prob of Solution Overlap  =  %lf", &probOverlap );
      // Get the coverage probability -- one blank line before getting the good stuff
      currentFile.getline( inputLine, 80, '\n' );
      currentFile.getline( inputLine, 80, '\n' );
      sscanf( inputLine, "%lf", &coverage ); 
      // Close the file
      currentFile.close();
      currentFile.clear();

      if( !fileIndex && !appendTag )
      {
         outputFile << "% Newsvendor Problem with m = " << batchSize << endl;
         outputFile << "\%Gamma" 
            // << '\t' << "Prob Soln Overlap"
            << '\t' << "Cover Prob"
            << '\t' << "Width of CI" 
            << '\t' << "Bias"
            << "\t\t" << "Average of Var"
            << '\t' << "Variance of Var"
            << endl;
      }
      outputFile << gamma
         // << '\t' << showpoint << setprecision(7) << probOverlap
         << "\t\t" << setprecision(4) << coverage
         << "\t\t" << setprecision(4) << widthCI
         << "\t\t" << setprecision(2) << bias
         << "\t\t" << setprecision(7) << avgVar
         << "\t" << varVar
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
