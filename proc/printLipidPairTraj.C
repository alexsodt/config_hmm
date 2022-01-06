#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "proc_definition.h"

const double cutoff = 20.0;


int main( int argc, char **argv )
{
	if( argc < 5 )
	{
		printf("Syntax: printLipidPairTraj.exe psf definition.inp LINE# hmm.inp/out\n");
		return 0;
	}

	FILE *psfFile = fopen(argv[1], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		return 0;
	}

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) )
		loadPSFfromPDB( psfFile );	
	else
		loadPSF( psfFile );

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

	load_pair_definition( argv[2] );

	int init_done = 0;
	int nat = curNAtoms();

	int line_number = atoi(argv[3]);
	
	FILE *decodeFile = fopen(argv[4],"r");
	char *buffer = (char *)malloc( sizeof(char) * 1000000 );

	for( int l = 0; l < line_number && !feof(decodeFile); l++ )
		getLine( decodeFile, buffer );

	if( feof(decodeFile) )
	{
		printf("Failed to read to line %d.\n", line_number );
		exit(1);
	}

	const char *t = buffer;
	
	t = advance_string(t, 2 );

	char *segid1 = (char *)malloc( sizeof(char) * ( strlen(t)+1) );
	char *segid2 = (char *)malloc( sizeof(char) * ( strlen(t)+1) );
	char *resName1 = (char *)malloc( sizeof(char) * (strlen(t)+1) );
	char *resName2 = (char *)malloc( sizeof(char) * (strlen(t)+1) );
	
	int res1,res2;	
	int nextra;
	int nr;
	if( is_seg(0) )
	{
		nr = sscanf( t, "# %s", segid1 );
		t = advance_string(t, 2 );
	}
	else
	{
		nr = sscanf( t, "# %s %d", resName1, &res1 );
		t = advance_string(t, 3 );
	}
	
	if( is_seg(1) )
	{
		nr = sscanf( t, "%s", segid2 );
		t = advance_string(t, 1 );
	}
	else
	{
		nr = sscanf( t, "%s %d", resName2, &res2 );
		t = advance_string(t, 2 );
	}

	nr = sscanf( t, "nextra %d", &nextra );
	t = advance_string(t,2);	
	printf("Reading %s %d / %s %d.\n", resName1, res1, resName2, res2 );

/*
	if( nr != 5 )
	{
		printf("Failed to read residues and residue #s from string '%s'.\n", t );
		exit(1);
	}
*/

	char *extra_seg[nextra];
	char *extra_resName[nextra];
	int extra_resNum[nextra];

	for( int e = 0; e < nextra; e++ )
	{
		extra_seg[e] = (char *)malloc( sizeof(char) * (1+strlen(t)) );
		extra_resName[e] = (char *)malloc( sizeof(char) * (1+strlen(t)) );

		int nr = sscanf(t, "%s %s %d", extra_seg[e], extra_resName[e], extra_resNum+e );

		printf("Extra: %s %s %d\n", extra_seg[e], extra_resName[e], extra_resNum[e] );

		t = advance_string(t, 3);
	}

	int start_frame;
	int stop_frame;

	t = advance_string( t, 3 );

	char *fileName_start = (char *)malloc( sizeof(char) * (1+strlen(t)) ); 
	char *fileName_stop = (char *)malloc( sizeof(char) * (1+strlen(t)) ); 

	nr = sscanf(t, "%s frame %d", fileName_start, &start_frame );
	t = advance_string(t, 6 );
	nr = sscanf(t, "%s frame %d", fileName_stop, &stop_frame );

	if( strcasecmp( fileName_start, fileName_stop ) )
	{
		printf("As currently implemented, the start and stop frames must be in the same dcd (here: %s and %s).\n", fileName_start, fileName_stop);
		exit(1);
	}

	printf("Start: %s %d\n", fileName_start, start_frame );
	printf("Stop: %s %d\n", fileName_stop, stop_frame );

	FILE *dcdFile = fopen(fileName_start,"r" );
	if( !dcdFile )
	{
		printf("Couldn't open file '%s'. It must be in the current directory.\n");
		exit(1);
	}

	readDCDHeader(dcdFile);
	setAligned();
	int nframes = curNFrames();

	off_t fp = ftello(dcdFile);
	loadFrame( dcdFile, at );

	for( int a= 0; a < curNAtoms(); a++ )
		at[a].zap();

	off_t fp2 = ftello(dcdFile);

	off_t del_per = fp2-fp;
	off_t offset = fp + del_per * (long)start_frame;

	fseeko( dcdFile, offset, SEEK_SET );	
	int do_break = 0;

	int nframes_traj = stop_frame-start_frame;

	for( int f = 0; f < nframes_traj; f++ )
	{
		double La, Lb, Lc;
		double alpha,beta,gamma;
		
		loadFrame( dcdFile, at );
		
		PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

		int init = 1;
		double last[3] = { 0,0,0};

		for( int unit = 0; unit < 3; unit++ )
		for( int a = 0; a < curNAtoms(); a++ )
		{
			int match = 0;
			if( is_seg(0) && !strcasecmp( at[a].segid, segid1) && unit == 0)
				match = 1;
			if( is_seg(1) && !strcasecmp( at[a].segid, segid2) && unit == 1)
				match = 1;
			if( !is_seg(0) && !strcasecmp(at[a].resname, resName1) && at[a].res == res1 && unit == 0)
				match = 1;
			if( !is_seg(1) && !strcasecmp(at[a].resname, resName2) && at[a].res == res2 && unit == 1 )
				match = 1;

			if( nextra > 0 && unit == 2)
			{
				for( int e = 0; e < nextra; e++ )
				{
					if( strcasecmp( at[a].segid, extra_seg[e] ) )
						continue;
					if( strcasecmp( at[a].resname, extra_resName[e] ) )
						continue;
					if( at[a].res != extra_resNum[e] )
						continue;

					match = 1;
				}
			}
			

			if( !match ) continue;

			if( !init )
			{
				while( at[a].x - last[0] < -La/2 ) at[a].x += La;
				while( at[a].x - last[0] > La/2 ) at[a].x -= La;
				while( at[a].y - last[1] < -Lb/2 ) at[a].y += Lb;
				while( at[a].y - last[1] > Lb/2 ) at[a].y -= Lb;
				while( at[a].z - last[2] < -Lc/2 ) at[a].z += Lc;
				while( at[a].z - last[2] > Lc/2 ) at[a].z -= Lc;
			}
			else
			{
				last[0] = at[a].x;
				last[1] = at[a].y;
				last[2] = at[a].z;
				init = 0;
			}
		
			printATOM( stdout, at[a].bead, at[a].res, at+a );
		}
		printf("END\n");

		for( int a= 0; a < curNAtoms(); a++ )
			at[a].zap();
	}

	fclose(dcdFile);
}





















