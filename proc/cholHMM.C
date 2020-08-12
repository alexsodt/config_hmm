#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

int main( int argc, char **argv )
{
	if( argc != 2 )
	{
		printf("Syntax: cholHMM.exe file\n");
		exit(1);
	}
	
	int max_chol = 1000;
	double vals[max_chol];

	FILE *theFile = fopen(argv[1],"r");
	
	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	getLine( theFile, buffer );

	int nr = readNDoubles( buffer, vals, max_chol );

	if( nr == max_chol )
	{
		printf("WARNING: the program is hard coded for maximum 1000 cholesterols. The number in the file either met or exceeded that level.\n");
	}

	int nchl = nr;

	int nframes=0;
	int nframesSpace = 10;
	double *zvals = (double *)malloc( sizeof(double) * nchl * nframesSpace );
	
	memcpy( zvals, vals, sizeof(double) * nchl );

	while( !feof(theFile) )
	{
		getLine( theFile, buffer );

		if( feof(theFile) ) break;
		int nr = readNDoubles( buffer, vals, max_chol );
	
		if( nframes == nframesSpace )
		{
			nframesSpace *= 2;

			zvals = (double *)realloc( zvals, sizeof(double) * nframesSpace * nchl );	
		}

		memcpy( zvals + nframes*nchl, vals, sizeof(double) * nchl );

		nframes++;
	} 

	double cut_val = 8.0;

	printf("COMMAND TRAIN\n");
	for( int c = 0; c < nchl; c++ )
	{
		for( int t = 0; t < nframes; t++ )
		{
			double z = zvals[t*nchl+c];
			if( z < -cut_val )
				printf("A");
			else if( z< 0 )
				printf("B");
			else if( z< cut_val )
				printf("C");
			else
				printf("D");
		}
		printf("\n");
	}
	printf("STOP\n");
}
