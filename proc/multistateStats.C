//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "proc_definition.h"

int main( int argc, char **argv )
{
	char *buffer = (char *)malloc( sizeof(char) * 1000000 );
	char *buffer2 = (char *)malloc( sizeof(char) * 1000000 );

	if( argc != 2 )
	{
		printf("Syntax: multistateStats.exe multistate_hmm_assignment\n");
		exit(1);
	}		
	int max_states = 11;
	int max_simul  = 5;

	double *raw_cnt = (double *)malloc( sizeof( double) * max_states * max_simul );
	memset( raw_cnt, 0, sizeof(double) * max_states * max_simul );
	
	double *pair_cnt = (double *)malloc( sizeof( double) * max_states * max_states );
	memset( pair_cnt, 0, sizeof(double) * max_states * max_states );

	FILE *stateFile = fopen(argv[1],"r");

	while( !feof(stateFile) )
	{
		getLine(stateFile, buffer );

		char *tread = buffer;

		int done = 0;
		
		while( *tread )
		{
			int nsimul = *tread - '0'; 

			if( nsimul == 0 )
			{
				//  zero state, zero simul.
				raw_cnt[0] += 1;		
				tread += 1;
			}
			else
			{
				tread += 1;

				int state_list[nsimul];
				for( int s = 0; s < nsimul; s++ )
				{
					int state = 1 + *tread - '0';
					if( state >= max_states )
					{
						printf("ERROR. State '%c' out of range.\n", *tread );
						exit(1);
					}
					
					state_list[s] = state;
					tread += 1;

					if( nsimul < max_simul )
						raw_cnt[state*max_simul+nsimul] += 1;
				}
				if( nsimul > 1 )
				{
					for( int s = 0; s < nsimul; s++ )
					for( int s2 = 0; s2 < nsimul; s2++ )
					{
						if( s == s2 ) continue;

						pair_cnt[state_list[s]*max_states+state_list[s2]] += 1;
					}
				}
			}
		}
	}
	
	printf("Total counts.\n");
	for( int simul = 0; simul < max_simul; simul++ )
	{
		for( int s = 0; s < max_states; s++ )
			printf("State %c simultaneous %d cnt %lf\n", (s == 0 ? 'X' : '0' + s - 1), simul, raw_cnt[s*max_simul+simul] );
	}
	printf("Pair correlation:\n");
	for( int s = 0; s < max_states; s++ )
	{
		double sum = 0;

		for( int s2 = 0; s2 < max_states; s2++ )
			sum += pair_cnt[s*max_states+s2];
		printf(" For state '%c':", (s == 0 ? 'X' : '0'+s-1) );
		for( int s2 = 0; s2 < max_states; s2++ )		
			printf(" %lf", pair_cnt[s*max_states+s2]/sum );
		printf("\n");
	}

}




