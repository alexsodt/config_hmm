#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include <math.h>

int main( int argc, char ** argv )
{
	if( argc < 4 )
	{
		printf("Syntax: bestFitC <outputFromFourierExtract_mode_curvature> qlow qhigh [state_select]\n");
		return 0;
	}
	FILE *theFile = fopen(argv[1],"r");

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1]);
		return 0;
	}

	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	
	getLine( theFile, buffer );

	int state_select = -1;
	int nmodes;
	int nstates;
	int ebar_mode = 0;
	double q_low = atof(argv[2]);
	double q_high = atof(argv[3]);

	if( argc > 4 )
	{
		state_select = atoi(argv[4]);
		if( argv[4][0] >= 'A' && argv[4][0] <= 'Z' )
			state_select = argv[4][0] - 'A';
		if( argv[4][0] >= 'a' && argv[4][0] <= 'z' )
			state_select = 26  + argv[4][0] - 'a';
	}

	int nr = sscanf( buffer, "nmodes %d nstates %d", &nmodes, &nstates );

	if( nr != 2 )
	{
		nr = sscanf( buffer, "nmodes %d nconfigs %d", &nmodes, &nstates );

		if( nr != 2 )
		{
			printf("Expected, but failed, to read the number of states and modes from the first line of the file.\n");
			return 0;
		}
	}

	double *state_pop = (double *)malloc( sizeof(double) * nstates );

	char *read_states = buffer + strlen("nmodes ");
	while( *read_states && *read_states != ' ' && *read_states != '\t' ) read_states += 1; // skip through the number
	while( *read_states && (*read_states == ' ' || *read_states == '\t') ) read_states += 1; // skip to next
	while( *read_states && *read_states != ' ' && *read_states != '\t' ) read_states += 1; // skip through the word "nstates"
	while( *read_states && (*read_states == ' ' || *read_states == '\t') ) read_states += 1; // skip to next
	while( *read_states && *read_states != ' ' && *read_states != '\t' ) read_states += 1; // skip through the number of states
	read_states += strlen(" num_per_state");
	nr = readNDoubles( read_states, state_pop, nstates );

	if( nr != nstates )
	{
		printf("Failed to read nstates %d populations from line '%s'.\n", nstates, read_states );
		return 0;
	}

	double *qvals = (double *)malloc( sizeof(double) * nmodes );

	getLine( theFile, buffer );

	nr = readNDoubles( buffer + strlen("mode q vals "), qvals, nmodes );
	
	if( nr != nmodes )
	{
		printf("Failed to read %d qvals from mode line.\n", nmodes );
		return 0;
	}

	double *state_curvature = (double *)malloc( sizeof(double) * nstates * nmodes );
	double *state_ebars = (double *)malloc( sizeof(double) * nstates * nmodes );

	for( int s = 0; s < nstates; s++ )
	{
		getLine( theFile, buffer );
		char *t = buffer + strlen("state ");
		while( *t && *t != ' ' ) t += 1;
		t += strlen(" curv ");
				
		nr = readNDoubles( t, state_curvature+s*nmodes, nmodes );
		if( nr != nmodes )
		{
			printf("Read failure in curvature.\n");
			exit(1);
		}
	}

	int *mode_sorter = (int *)malloc( sizeof(int) * nmodes );

	for( int i = 0; i < nmodes; i++ )
		mode_sorter[i]=i;
	int done = 0;
	while(!done)
	{
		done = 1;

		for( int t = 0; t < nmodes-1; t++ )
		{
			if( qvals[mode_sorter[t]] > qvals[mode_sorter[t+1]] )
			{
				int x = mode_sorter[t];
				mode_sorter[t] = mode_sorter[t+1];
				mode_sorter[t+1] = x;
				done = 0;
			}
		}	
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		getLine( theFile, buffer );
		char *t = buffer + strlen("state ");
		while( *t && *t != ' ' && *t != '\t' ) t += 1;
		t += strlen(" ebar ");
				
		nr =readNDoubles( t, state_ebars+s*nmodes, nmodes );
		if( nr != nmodes )
		{
			printf("Read failure in error bars.\n");
			exit(1);
		}
	}

	int mode_count[nstates];
	memset( mode_count, 0, sizeof(int)* nstates );
	double c[nstates];
	double ce[nstates];
	memset( c, 0, sizeof(double) * nstates );
	memset( ce, 0, sizeof(double) * nstates );

	for( int mode_iq = 1; mode_iq < nmodes; mode_iq++ )
	{	
		if( qvals[mode_iq] <= 1e-30 ||  qvals[mode_iq] < q_low || qvals[mode_iq] > q_high )
			continue;

		for( int s = 0; s < nstates; s++ )
		{
			double alt_err = state_ebars[s*nmodes+mode_iq];
	//		double use_err = alt_err_per_pops / sqrt( state_pop[s]);
			//printf("alt_err: %le use_err: %le\n", alt_err, use_err );
			double use_err = 1;
			double w = 1.0 / pow( use_err, 2.0 );
			c[s] += state_curvature[s*nmodes+mode_iq];
			ce[s] += state_curvature[s*nmodes+mode_iq]*state_curvature[s*nmodes+mode_iq]; 
			mode_count[s] += 1;
		}
	}

	if( state_select >= 0 )
	{
		if( state_select >= nstates )
		{
			printf("State selection greater than the number of states.\n");
			exit(1);
		}

		printf("%le\n", c[state_select] / mode_count[state_select] );
	}
	else
	{
		for( int s = 0; s < nstates; s++ )
			printf(" %le", c[s] );
		for( int s = 0; s < nstates; s++ )
			printf(" %d", mode_count[s] );
		printf(" err");
		for( int s = 0; s < nstates; s++ )
		{
			double x2 = ce[s] / mode_count[s];
			double xm = c[s] / mode_count[s];
			
			printf(" %le", sqrt(x2-xm*xm)/sqrt(mode_count[s]) );
		}
		printf("\n");
	}
	return 0;
}
