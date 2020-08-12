#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

int main( int argc, char **argv )
{
	double kc = 9.55;
	if( argc < 3 )
	{
		printf("Syntax: estimatedFdR.exe modelC.out states.out [nPOPC] [POPC F']\n");
		printf("Takes the c_0 of a lipid dimer model (modelC.out) and a state assignment (states.out) and produces an estimate of dFdC\n");
		printf("nPOPC default assumes there are 200 total lipids, with the rest being POPC.\n");
		printf("POPC F' set to be 0.038.\n");
		printf("Monolayer kc = 9.55 kcal/mol.\n");
		return 0;
	}
	FILE *theFile = fopen(argv[1],"r");

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	#define MAX_STATES 48
	
	double stateC[MAX_STATES];
	double stateE[MAX_STATES];
	int nstates = 0;

	double nPOPC = -1;

	if( argc > 3 ) 
		nPOPC = atof(argv[3]);

	while( !feof(theFile) && nstates < MAX_STATES )
	{
		getLine(theFile, buffer );
		if( feof(theFile) ) break;

		int tstate;
		double c, std;
		int nr = sscanf( buffer, "State %d <c> %le std %le", &tstate, &c, &std );

		if( nr == 3 )
		{
			stateC[nstates] = c;
			stateE[nstates] = std;
			nstates++;
		}
		else
			printf("Failed to read three fields from line '%s'.\n", buffer );
	}
	fclose(theFile);

	theFile = fopen(argv[2],"r");
	if( !theFile )
	{
		printf("Couldn't open state file '%s'.\n", argv[2] );
		return -1;
	}

	double av_Fp = 0;
	double av_Fp_err2 = 0;
	double nav = 0;

	double Fp_POPC = 0.038;
	double Fp_POPC_err = 0.003;

	double av_n1 = 0;

	while( !feof(theFile) )
	{
		getLine(theFile, buffer );
		if( feof(theFile) ) break;

		double nPOPC_inst = nPOPC;
		if( nPOPC_inst < 0 )
			nPOPC_inst = 200 - strlen(buffer);

		double nLipids = nPOPC_inst + strlen(buffer);

		double inst_Fp = nPOPC_inst * Fp_POPC;
		double inst_Fp_err2 = pow( nPOPC_inst * Fp_POPC_err, 2.0);
	
		for( int t = 0; t < strlen(buffer); t++ )
		{
			int state = buffer[t] - '0';

			if( state == 2 )
				av_n1 += 1;
		
			inst_Fp -= stateC[state] * kc;
			inst_Fp_err2 += pow( stateE[state] * kc, 2.0 );	
		}

		inst_Fp /= nLipids;
		inst_Fp_err2 /= nLipids * nLipids;
		
		av_Fp += inst_Fp;
		av_Fp_err2 += inst_Fp_err2;

		nav += 1;
	}

	av_n1 /= nav;

	av_Fp /= nav;
	av_Fp_err2 /= nav*nav;

	printf("Expected: %le err: %le n_1: %lf\n", av_Fp, sqrt(av_Fp_err2), av_n1);
}
