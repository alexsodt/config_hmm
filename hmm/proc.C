#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include <string.h>

int main( int argc, char **argv )
{
	FILE * theFile = fopen(argv[1],"r");

	char *buffer = (char *)malloc( sizeof(char) * 40000 );

	int step = 0;
	while( !feof(theFile) )
	{

		getLine( theFile, buffer );
		if( feof(theFile) ) break;
	
		for( int s = 0; s < strlen(buffer); s++ )
			printf("%d,%d,%c\n", 1+s, step, (buffer[s] == '1' ? '1' : '0' ) ); 	
		step++;
	}
}
