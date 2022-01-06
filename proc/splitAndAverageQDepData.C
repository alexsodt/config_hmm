#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

double q_bin_width = 0.030;

int main( int argc, char **argv )
{
	if( argc < 2 )
	{
		printf("Syntax: splitAndAverageQDepData [qwidth=#] file1 ... fileN\n");
		return 0;
	}

	int arg_off = 0;

	if( !strncasecmp(argv[1],"qwidth=", 7 ) )
	{
		sscanf( argv[1], "qwidth=%lf", &q_bin_width );
		if( q_bin_width < 1e-14 )
		{
			printf("ERROR: q_bin_width %lf\n", q_bin_width );
			exit(1);
		}
		arg_off += 1;
	}

	double q_max = 0.3;
	int t = q_max / q_bin_width + 1;
	q_max = t * q_bin_width;
	int n_q_bins = lround(q_max / q_bin_width);


	double *sum_c0_q = ( double *)malloc( sizeof(double) * n_q_bins );
	double *sum_c02_q = ( double *)malloc( sizeof(double) * n_q_bins );
	double *nav_c0_q = (double *)malloc( sizeof(double) * n_q_bins );

	memset( sum_c0_q, 0, sizeof(double) * n_q_bins );
	memset( sum_c02_q, 0, sizeof(double) * n_q_bins );
	memset( nav_c0_q, 0, sizeof(double) * n_q_bins );

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	for( int c = 1+arg_off; c < argc; c++ )
	{
		FILE *theFile = fopen( argv[c], "r" );
		if( !theFile )
		{
			printf("Couldn't open file '%s'.\n", argv[c] );
			exit(1);
		}

		int nm = 0;

		getLine( theFile, buffer );	

		while( !feof(theFile) && buffer[0] == '#' )
			getLine( theFile, buffer );	

		int nr = sscanf( buffer, "nmodes %d nstates", & nm );

		double *qvals = (double*)malloc( sizeof(double) * nm );
		memset( qvals, 0, sizeof(double) * nm );

		int nrq, nrm;
		getLine( theFile, buffer );
		if( strlen(buffer) > strlen("mode q vals ") )
			nrq = readNDoubles( buffer + strlen("mode q vals "), qvals, nm );	
		else
		{
			printf("File read error '%s'\n", argv[c]);
			exit(1);
		}
		
		double *curvs = (double*)malloc( sizeof(double) * nm );
		memset( curvs, 0, sizeof(double) * nm );
		
		getLine( theFile, buffer );
		if( strlen(buffer) > strlen("state 0 curv ") )
			nrm = readNDoubles( buffer + strlen("state 0 curv "), curvs, nm );	
		else
		{
			printf("File read error '%s'\n", argv[c]);
			exit(1);
		}
	
		for( int iq = 0; iq < nrq && iq < nrm; iq++ )
		{
			if( qvals[iq] < 1e-10 ) continue;
			int qbin = qvals[iq] / q_bin_width; 
			if( qbin >= 0 && qbin < n_q_bins )
			{
				sum_c0_q[qbin] += curvs[iq];
				sum_c02_q[qbin] += curvs[iq]*curvs[iq];
				nav_c0_q[qbin] += 1;
			}
			printf("scatter %lf %le\n", qvals[iq], curvs[iq] );	
		}

		fclose(theFile);		
	}

	for( int iq = 0; iq < n_q_bins; iq++ )
	{	
		sum_c0_q[iq] /= nav_c0_q[iq];
		sum_c02_q[iq] /= nav_c0_q[iq];
	}


	for( int iq = 0; iq < n_q_bins; iq++ )
	{
		double symmetry_factor = 0.5; // the pos/neg modes are pre-averaged and duplicated so appear twice --- makes the error appear smaller. correct here
		double err = sqrt( sum_c02_q[iq] - sum_c0_q[iq] * sum_c0_q[iq] ) / sqrt( nav_c0_q[iq] * symmetry_factor);
		if( nav_c0_q[iq] > 0 )
			printf("hist %le %le %le %lf\n", (iq+0.5)*q_bin_width, sum_c0_q[iq], err, nav_c0_q[iq] );
	}
}
