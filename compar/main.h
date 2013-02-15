#ifndef __MAIN_H__
#define __MAIN_H__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <unistd.h>

#include "rand01_twister.h"

using namespace std;

/*
#ifndef __STDTYPE_H__
#	include "stdtype.h"
#endif
#ifndef __SIMPLEX_H__
#	include "simplex.h"
#endif
#ifndef __SUB_MAN_H__
#	include "sub_man.h"
#endif
#ifndef __PARSESPC_H__
#	include "parsespc.h"
#endif
*/

//==============================================================================
//
//	Function prototypes.
//

void run( int argc, char *argv[] );
void ReGenerateScenarios( int scennum, int gamma );
double CalculateGap( int batchSize );
void ReadScenarios( string filename );

//
//	End of function prototypes.
//
//==============================================================================

#endif
