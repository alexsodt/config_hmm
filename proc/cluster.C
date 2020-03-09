// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"
#include "ctype.h"
#include "alignSet.h"
#include <sys/time.h>
#include "proc_definition.h"

#define N_SITE_PER_CHAIN 6

#include "comparison.h"
#define SUB_PRESS
#define LEN 100

#define OUTPUT_PS
#define REMOVE_COM

#define MAX_LIPID 500
#define MAX_LATOMS 200

int N_TALLY = 6; 

#define N_GRID 50

#define MIN_X -30
#define MAX_X  30

#define MIN_Y -30
#define MAX_Y  30

#define MIN_Z -30
#define MAX_Z  30


int getAtCode( char * res )
{
	if( strlen(res) == 4 )
	{
		if( !strncasecmp( res, "DNPC", 4 ) ) 
			return 2; 
		if( !strncasecmp( res, "DEPC", 4 ) ) 
			return 2; 

		if( res[0] == 'D' && res[2] == 'P' )
			return 1;
	}	

	return 0;
}
void getShellEnhancement( double * cur_pos_top, double * cur_pos_bottom, 
			int ntop, int nbot, const char *unique, 
			double Lx, double Ly, int * at_code_top, int * at_code_bottom, 
			int *shell_top, int *shell_bottom,
			double *redArea, double *purpleArea, double *blueArea, // these are broken down by shell. 
			double *redTop, double *redBottom, 
			double *blueTop, double *blueBottom, int *res1_top, int *res2_bottom );

void printMap( double * cur_pos_top, double * cur_pos_bottom, int ntop, int nbot, const char *unique, double Lx, double Ly, int * at_code_top, int * at_code_bottom, 
char *fileName, double *redArea, double *purpleArea, double *blueArea, double *redTop, double *redBottom, double *blueTop, double *blueBottom );

char code( int cnts[3] )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];

	if( done == 0 )
	{
		array = (int *)malloc( sizeof(int) * (N_TALLY+1)*(N_TALLY+1) );

		char cur = 'A';

		for( int x = 0; x <= N_TALLY; x++ )
		for( int y = 0; y <= N_TALLY; y++ )
			array[x*(N_TALLY+1)+y] = 'a';

		for( int Neither = 0; Neither <= N_TALLY; Neither++ )
		for( int y = 0; y <= Neither; y++ )
		{
			int x = Neither - y;

			array[y*(N_TALLY+1)+x] = cur;

			if( cur == 'Z' )
				cur = 'a';
			else
				cur += 1;
		}

		done = 1;
	}	

	return array[cnts[0]*(N_TALLY+1)+cnts[1]];
}

struct elem
{
	int i, j,k;
	double PXX,PYY,PZZ;
};



int main( int argc, char **argv )
{
        struct timeval tp;
 
        gettimeofday( &tp, NULL );
 
	char buffer[4096];

	if( argc < 3)
	{
		printf("Syntax: cluster NMED definition.inp PDB [PDB2 ...]\n");	
		return 0;
	}


	int kmedoids = atoi( argv[1] );
	
	int arg_offset = 0;
	
	load_pair_definition( argv[2] );

	binary_penalty = pair_binary_penalty();
	binary_benefit = pair_binary_benefit();
	hbond0_penalty = pair_hbond0_penalty();
	hbond0_benefit = -hbond0_penalty;

	printf("benefit: %lf penalty: %lf hbond0: %lf\n",
		binary_benefit, binary_penalty, hbond0_penalty );

	int nat = 0;
	{
		FILE *thePDB = fopen(argv[3+arg_offset],"r");
		
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) ) nat++;
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}

		fclose(thePDB);
	}


	int nstructsSpace = 100;
	aStructure *structs = (aStructure *)malloc( sizeof(aStructure) * nstructsSpace );

	for( int i = 0; i < nstructsSpace; i++ )
		structs[i].atoms = (double *)malloc( sizeof(double) * 3 * nat );

//	double *structs = (double *)malloc( sizeof(double) *3 * nat * nstructsSpace  );
	char **remarks = (char **)malloc( sizeof(char *) * nstructsSpace );
	char **remarksc = (char **)malloc( sizeof(char *) * nstructsSpace );
	int nstructs = 0;
	int cat = 0;

	for( int x = 0; x < nstructsSpace; x++ )
	{
		remarks[x] = NULL;	
		remarksc[x] = NULL;	
		structs[x].binary_data = 0;
	}
	
	
	for( int c = 3+arg_offset; c < argc; c++ )
	{
		FILE *thePDB = fopen(argv[c],"r");
	
		int nstructs_local = 0;
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) {
				 break;	
			}
	
			if( !strncasecmp( buffer, "REMARK CODE", 11 ) )
			{
				remarksc[nstructs] = (char *)malloc( sizeof(char) * (1+strlen(buffer) ) );
				strcpy( remarksc[nstructs], buffer );

				unsigned long code;

				sscanf(buffer, "REMARK CODE %lu", &code );

				structs[nstructs].binary_data = code;
			}
			else if( !strncasecmp( buffer, "REMARK", 6 ) )
			{
				remarks[nstructs] = (char *)malloc( sizeof(char) * (1+strlen(buffer) ) );
				strcpy( remarks[nstructs], buffer );
			}
	
			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
				if( nat == cat )
				{
					printf("ATOM # inconsistency in file %s, structure %d.\n", argv[c], nstructs_local );
					exit(1);
				}
				struct atom_rec tat;
	
				readATOM( buffer, &tat );

				structs[nstructs].atoms[3*cat+0] = tat.x;
				structs[nstructs].atoms[3*cat+1] = tat.y;
				structs[nstructs].atoms[3*cat+2] = tat.z;

				cat++;
			}

//			printf("buffer: %s\n", buffer );
			if( !strncasecmp( buffer, "END", 3) || !strncasecmp( buffer, "TER", 3) )
			{
				nstructs_local++;
				nstructs++;
//				printf("new struct ! %d\n", nstructs );

				if( nstructs == nstructsSpace )
				{
					nstructsSpace *= 2;
	
					structs = (aStructure *)realloc( structs, sizeof(aStructure) * nstructsSpace );
					remarks = (char **)realloc( remarks, sizeof(char *) * nstructsSpace );
					remarksc = (char **)realloc( remarksc, sizeof(char *) * nstructsSpace );
					for( int x = nstructs; x < nstructsSpace; x++ )
					{
						structs[x].binary_data = 0;
						remarks[x] = NULL;	
						remarksc[x] = NULL;	
						structs[x].atoms = (double *)malloc( sizeof(double) * 3 * nat );
					}
				}	
				cat = 0;
			}
		}
	}
	printf("Read %d structures.\n", nstructs);
	
	int medoids[kmedoids];
	
	doKMedoidClusteringLR( structs, nstructs, nat, kmedoids, medoids, 0.8, 0.5, 500); 
	//doKMedoidClustering( lipid_configs, nconfigs, lipid_length, kmedoids, medoids ); 

	FILE *cenFile = fopen("centers.pdb","w");

	for(int x =0; x < kmedoids; x++ )
	{
		FILE *thePDB = fopen(argv[3+arg_offset],"r");

		fprintf(cenFile, "%s\n", remarks[medoids[x]] );
		fprintf(cenFile, "%s\n", remarksc[medoids[x]] );
	
		int tp = 0;	
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
				struct atom_rec atp;

				readATOM( buffer, &atp );

				atp.x = structs[medoids[x]].atoms[tp*3+0];
				atp.y = structs[medoids[x]].atoms[tp*3+1];
				atp.z = structs[medoids[x]].atoms[tp*3+2];

				printATOM( cenFile, atp.bead, atp.res, &atp );

				tp++;
			}
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}
		fprintf(cenFile, "END\n");

		fclose(thePDB);

		if( remarks[medoids[x]] )
			printf("medoid %d struct %d %s %s\n", x, medoids[x], remarks[medoids[x]], remarksc[medoids[x]] );
		else
			printf("medoid %d struct %d\n", x, medoids[x]  );
	}	

	fclose(cenFile);	
}

