#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "proc_definition.h"

static double cutoff = 14.0;


int main( int argc, char **argv )
{
	if( argc < 4 )
	{
		printf("Syntax: printFullDimers psf definition.inp centers.pdb\n");
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

	cutoff = getCutoff();

	int init_done = 0;
	int nat = curNAtoms();

	FILE *thePDB = fopen(argv[3],"r");

	char *buffer = (char *)malloc( sizeof(char) * 1000000 );

	while( !feof(thePDB) )
	{
		getLine(thePDB, buffer );
		if( feof(thePDB) ) break;
	
		if( !strncasecmp( buffer, "REMARK CODE", 11 ) )
		{
		}
		else if( !strncasecmp( buffer, "REMARK", 6 ) )
		{
			char segid1[256];
			char segid2[256];
			char resName1[256];
			char resName2[256];
			int res1,res2;
			char fileName[256];
			int frame;

			int nextra = 0;
			const char *extra_string = NULL;

			const char *t = buffer;

			sscanf( buffer, "REMARK %s", fileName );
	
			int err = 0;
			int arg_advance = 0;

			// advance past REMARK and filename.
			t = advance_string(t, 2 );

			if( is_seg(0) )
			{
				int nr = sscanf(t, "%s ", segid1 );
				if( nr != 1 ) err = 1;
				arg_advance = 1;
			}
			else
			{
				int nr = sscanf(t, "%s %d ", resName1, &res1 );
				if( nr != 2 ) err = 1;
				arg_advance = 2;
			}
	
			t = advance_string( t, arg_advance );
			
			if( is_seg(1) )
			{
				int nr = sscanf(t, "%s ", segid2 );
				if( nr != 1 ) err = 1;
				arg_advance = 1;
			}
			else
			{
				int nr = sscanf(t, "%s %d ", resName2, &res2 );
				if( nr != 2 ) err = 1;
				arg_advance = 2;
			}
		
			t = advance_string( t, arg_advance );	

			int nr = sscanf(t, "frame %d", &frame );

			t = advance_string( t, 2 );

			if( !strncasecmp( t, "nextra", 6 ) )
			{
				sscanf( t, "nextra %d", &nextra );
				t = advance_string(t, 2 );

				extra_string = t;
			}

			FILE *dcdFile = fopen(fileName,"r" );
			if( !dcdFile )
			{
				printf("Couldn't open file '%s'. It must be in the current directory.\n");
				exit(1);
			}
		
			readDCDHeader(dcdFile);
			setAligned();
			int nframes = curNFrames();

#define SEEK_IT

#ifdef SEEK_IT
			off_t fp = ftello(dcdFile);
			loadFrame( dcdFile, at );

			for( int a= 0; a < curNAtoms(); a++ )
				at[a].zap();

			off_t fp2 = ftello(dcdFile);
		
			off_t del_per = fp2-fp;
			off_t offset = fp + del_per * (long)frame;

			fseeko( dcdFile, offset, SEEK_SET );	
#endif
			int do_break = 0;

			int f = frame;

#ifndef SEEK_IT
			for( int f = 0; f < nframes; f++ )
#endif
			{
				double La, Lb, Lc;
				double alpha,beta,gamma;
				
				loadFrame( dcdFile, at );
				
				PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

				if( f == frame )
				{
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
							const char * t = extra_string;

							for( int e = 0; e < nextra && t; e++ )
							{
								if( strncasecmp( at[a].segid, t, strlen(at[a].segid) ) )
									continue;
								while( *t && *t != '_' ) t+= 1;
								if( !*t ) continue;
								t+=1;
								if( strncasecmp( at[a].resname, t, strlen(at[a].resname) ) )
									continue;
								while( *t && *t != '_' ) t+= 1;
								if( !*t ) continue;
								t+=1;
								if( at[a].res != atoi(t) )
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
#ifndef SEEK_IT	
					do_break = 1;
#endif
				}

				for( int a= 0; a < curNAtoms(); a++ )
					at[a].zap();

			if( do_break ) break;

			}
			fclose(dcdFile);
		}
	}

}





















