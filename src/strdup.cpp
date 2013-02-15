#include <string.h>
#include <stdlib.h>

#ifndef __STRDUP_H__
#	include "strdup.h"
#endif

char *DuplicateString( const char *s )
{
	char *c = (char *)malloc( strlen(s) + 1 );

	if( c ) strcpy( c, s );
	return c;
}
