#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main( int argc, char **argv )
{
	int nvals = 0;
	int nvalSpace = 10;
	int *vals = (int *)malloc( sizeof(int) * nvalSpace );

	FILE *theFile = fopen(argv[1],"r");
	
	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1] );
		return -1;
	}

	while( !feof(theFile) )
	{
		int val;
		int nr = fscanf(theFile, "%d\n", &val );

		if( nr == 1 )
		{
			if( nvals == nvalSpace )
			{
				nvalSpace *= 2;
				vals = (int *)realloc( vals, sizeof(int) * nvalSpace );
			}

			vals[nvals] = val;
			nvals++;
		}
	}

	int done = 0;
		// bubble sort. if you care: congratulations. do not email me.
	while(!done)
	{
		done = 1;

		for( int t = 0; t < nvals-1; t++ )
		{
			if( vals[t] > vals[t+1] )
			{
				int temp = vals[t];
				vals[t] = vals[t+1];
				vals[t+1] = temp;
				done = 0;
			}
		}
	}

	for( int t = 0; t < nvals; t++ )
	{
		if( t == 0 )
			printf("%d\n", vals[t] );
		else if( vals[t] != vals[t-1] )
			printf("%d\n", vals[t] );
	}
}
