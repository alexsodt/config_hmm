// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int readNInts( char *buffer, int *vals, int nvalues );
int readNDoubles( char *buffer, double *vals, int nvalues );
void getLine( FILE *theFile, char *theBuffer );
void print5( int val, char *str );
const char *advance_string( const char *t, int nadv );
int decodeString( char *buf, char **out, int nmax );
