#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "proc_definition.h"
#include "comparison.h"
	


const double cutoff = 14.0;
double info_hbond_cutoff = 2.5;

struct foundABond
{
	int resNum[2];
	int atIndex[2];
	char atName1[256];
	char atName2[256];

	int example_abs[2];
	int example_f;
	char example_dcd[256];
	int count;
	
	struct foundABond *next;
};

int main( int argc, char **argv )
{
	if( argc < 5 )
	{
		printf("Syntax: extractDimers psf/pdb stride definition.inp dcd1 [dcd2 ...]\n");
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

	int init_done = 0;
	int nat = curNAtoms();

	double n_frame_add = 0;

	int *atom_list[2];
		atom_list[0]  = (int *)malloc( sizeof(int) * curNAtoms() );
		atom_list[1]  = (int *)malloc( sizeof(int) * curNAtoms() );
	int *first_atom[2];
		first_atom[0]  = (int *)malloc( sizeof(int) * curNAtoms() );
		first_atom[1]  = (int *)malloc( sizeof(int) * curNAtoms() );
	int *num_atom[2]; 
		num_atom[0]  = (int *)malloc( sizeof(int) * curNAtoms() );
		num_atom[1]  = (int *)malloc( sizeof(int) * curNAtoms() );
	int *unit_ptr[2];
		unit_ptr[0]  = (int *)malloc( sizeof(int) * curNAtoms() );
		unit_ptr[1]  = (int *)malloc( sizeof(int) * curNAtoms() );
	int *unit_size[2];
		unit_size[0]  = (int *)malloc( sizeof(int) * curNAtoms() );
		unit_size[1]  = (int *)malloc( sizeof(int) * curNAtoms() );
	int nunits[2] = {0,0};
	int natoms[2] = {0,0};

	double *unit_r[2] = { NULL, NULL };

	int stride = atoi(argv[2]);

	int cntr = 0;

	struct foundABond *allBonds = NULL;

	int arg_offset = 0;

	binary_penalty = 0;
	binary_benefit = 0;
	hbond0_penalty = 0;

	load_pair_definition( argv[3] );

	for( int c = 4; c < argc; c++ )
	{
		FILE *dcdFile = fopen( argv[c], "r");
	
		readDCDHeader(dcdFile);
		setAligned();

		int nframes = curNFrames();
	
		fprintf(stderr, "Processing %s.\n", argv[c] );
		for( int f = 0; f < nframes; f++,  n_frame_add++, cntr++)
		{
			double La, Lb, Lc;
			double alpha,beta,gamma;
			
			loadFrame( dcdFile, at );
			
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

			if( cntr % stride != 0 )
			{
				for( int a = 0; a < curNAtoms(); a++ )
					at[a].zap();
				continue;
			}

			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}
		
			if( !init_done )
			{
				int pres = -1;	
				int cur_hex = 0;

				char psegid[256];

				for( int pass = 0; pass < 2; pass++ )
				{
					pres = -1;
					psegid[0] = '\0';

					int cur_unit = -1;
					int op = 0;
		
					for( int a = 0; a < curNAtoms(); a++ )
					{
						int is_new = 0;
						if( is_seg(pass) )
						{
							if( strcasecmp( at[a].segid, psegid) )
							{
								op = 0;
								is_new = 1;
							}
						}	
						else
						{
							if( at[a].res != pres )
							{
								op = 0;
								is_new = 1;
							}
						}
			
						if( !use_atom( USE_ANYWHERE, pass, at[a].res, at[a].atname, at[a].resname, at[a].segid ) )
							continue;
	
						if( is_new )
						{
							cur_unit++;
							op=1;
							
							first_atom[pass][cur_unit] = a;
							num_atom[pass][cur_unit] = 1;
							unit_ptr[pass][cur_unit] = natoms[pass];
							unit_size[pass][cur_unit] = 0;
							nunits[pass]++;
					
						}
						else if( op )
						{
							num_atom[pass][cur_unit]++;
						}
						
						if( use_atom( USE_CHI2, pass, at[a].res, at[a].atname, at[a].resname, at[a].segid ) )
						{
							atom_list[pass][natoms[pass]] = a;
							unit_size[pass][cur_unit]++;
							natoms[pass]++;	
						}
					
						pres = at[a].res;	
						strcpy( psegid, at[a].segid );
					}		
				}

				unit_r[0] = (double *)malloc( sizeof(double) * 3 * nunits[0] );
				unit_r[1] = (double *)malloc( sizeof(double) * 3 * nunits[1] );

				init_done = 1;
			} 
	
			for( int pass = 0; pass < 2; pass++ )
			for( int u = 0; u < nunits[pass]; u++ )
			{
				unit_r[pass][3*u+0] = 0;
				unit_r[pass][3*u+1] = 0;
				unit_r[pass][3*u+2] = 0;
				for( int ax = 0; ax < unit_size[pass][u]; ax++ )
				{
					int a = atom_list[pass][unit_ptr[pass][u]+ax];

					unit_r[pass][3*u+0] += at[a].x;
					unit_r[pass][3*u+1] += at[a].y;
					unit_r[pass][3*u+2] += at[a].z;
				}
	
				unit_r[pass][3*u+0] /= unit_size[pass][u];
				unit_r[pass][3*u+1] /= unit_size[pass][u];
				unit_r[pass][3*u+2] /= unit_size[pass][u];
			}

			for( int u = 0; u < nunits[0]; u++ )
			{
				for( int u2 = u+1; u2 < nunits[1]; u2++ )
				{
					if( is_sym() && u2 <= u ) continue;

					double dr[3] = { 
						unit_r[1][3*u2+0] - unit_r[0][3*u+0],
						unit_r[1][3*u2+1] - unit_r[0][3*u+1],
						unit_r[1][3*u2+2] - unit_r[0][3*u+2] };

					double shift[3] = { 0,0,0};
					while( dr[0]+shift[0] < -La/2 ) shift[0] += La;
					while( dr[0]+shift[0] > La/2 ) shift[0] -= La;
					while( dr[1]+shift[1] < -Lb/2 ) shift[1] += Lb;
					while( dr[1]+shift[1] > Lb/2 ) shift[1] -= Lb;
					while( dr[2]+shift[2] < -Lc/2 ) shift[2] += Lc;
					while( dr[2]+shift[2] > Lc/2 ) shift[2] -= Lc;
					
					dr[0] += shift[0];
					dr[1] += shift[1];
					dr[2] += shift[2];

					double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

					if( r < cutoff )
					{
						long binary_hbond = 0;

						for( int a1 = first_atom[0][u]; a1 < first_atom[0][u] + num_atom[0][u]; a1++ ) 
						for( int a2 = first_atom[1][u2]; a2 < first_atom[1][u2] + num_atom[1][u2]; a2++ ) 
						{
							if( !(at[a1].atname[0] == 'O' && at[a2].atname[0] == 'H' ) &&
							    !(at[a2].atname[0] == 'O' && at[a1].atname[0] == 'H' ) )
								continue;
							
							double adr[3] = { at[a1].x - at[a2].x, at[a1].y - at[a2].y, at[a1].z - at[a2].z };	

							while( adr[0] < -La/2 )adr[0] += La;
							while( adr[0] > La/2 ) adr[0]-= La;
							while( adr[1] < -Lb/2 )adr[1] += Lb;
							while( adr[1] > Lb/2 ) adr[1]-= Lb;
							while( adr[2] < -Lc/2 )adr[2] += Lc;
							while( adr[2] > Lc/2 ) adr[2]-= Lc;

							double r = sqrt(adr[0]*adr[0]+adr[1]*adr[1]+adr[2]*adr[2]);

							if( r < info_hbond_cutoff )
							{
								int bond_index = binaryHBonder( at+a1, at+a2 );
								if( bond_index >= 0 && !(binary_hbond & (1<<bond_index)) )
									binary_hbond += (1<<bond_index); 
							}
						}

						if( info_mode() )
						{
							for( int a1 = first_atom[0][u]; a1 < first_atom[0][u] + num_atom[0][u]; a1++ ) 
							for( int a2 = first_atom[1][u2]; a2 < first_atom[1][u2] + num_atom[1][u2]; a2++ ) 
							{
								if( !(at[a1].atname[0] == 'O' && at[a2].atname[0] == 'H' ) &&
								    !(at[a1].atname[0] == 'H' && at[a2].atname[0] == 'O' ) ) 
									continue;
	
								double adr[3] = { at[a1].x - at[a2].x, at[a1].y - at[a2].y, at[a1].z - at[a2].z };	
	
								while( adr[0] < -La/2 )adr[0] += La;
								while( adr[0] > La/2 ) adr[0]-= La;
								while( adr[1] < -Lb/2 )adr[1] += Lb;
								while( adr[1] > Lb/2 ) adr[1]-= Lb;
								while( adr[2] < -Lc/2 )adr[2] += Lc;
								while( adr[2] > Lc/2 ) adr[2]-= Lc;
	
								double r = sqrt(adr[0]*adr[0]+adr[1]*adr[1]+adr[2]*adr[2]);
	
								if( r < 1.5 )
									continue; // OH bond
	
								if( r < info_hbond_cutoff )
								{
									int gotit = 0;
									int rel_index1 = a1 - first_atom[0][u];
									int rel_index2 = a2 - first_atom[1][u2];
									for( struct foundABond *bond = allBonds; bond; bond = bond->next )
									{
										if( (bond->atIndex[0] == rel_index1 && bond->atIndex[1] == rel_index2 ) ||
										    (bond->atIndex[1] == rel_index1 && bond->atIndex[0] == rel_index2 ) )
										{
											gotit =1;
											bond->count++;
										}	
									}
	
									if( !gotit )
									{
										foundABond *bond = (foundABond *)malloc( sizeof(foundABond) );
										bond->atIndex[0] = rel_index1;
										bond->atIndex[1] = rel_index2;
										bond->resNum[0] = at[a1].res;
										bond->resNum[1] = at[a2].res;
										bond->example_abs[0] = a1;
										bond->example_abs[1] = a2;
										bond->example_f = f;
										strcpy( bond->example_dcd, argv[c] );
										strcpy( bond->atName1, at[a1].atname );
										strcpy( bond->atName2, at[a2].atname );
										bond->count=1;
										bond->next = allBonds;
										allBonds = bond;
									}
								}
							}
						}
						printf("REMARK %s", argv[c] );

						if( is_seg(0) )
							printf(" %s", at[atom_list[0][unit_ptr[0][u]]].segid );
						else
							printf(" %s %d", at[atom_list[0][unit_ptr[0][u]]].resname, at[atom_list[0][unit_ptr[0][u]]].res );
						
						if( is_seg(1) )
							printf(" %s", at[atom_list[1][unit_ptr[1][u]]].segid );
						else
							printf(" %s %d", at[atom_list[1][unit_ptr[1][u]]].resname, at[atom_list[1][unit_ptr[1][u]]].res );
							
						printf(" frame %d\n", f );

						printf("REMARK CODE %lu\n", binary_hbond );
						for( int ax = 0; ax < unit_size[0][u]; ax++ )
							printATOM( stdout,
							        at[atom_list[0][unit_ptr[0][u]+ax]].bead, 
							        at[atom_list[0][unit_ptr[0][u]+ax]].res, 
								at+atom_list[0][unit_ptr[0][u]+ax] );
	
			

						for( int ax = 0; ax < unit_size[1][u2]; ax++ )
						{
							int a = atom_list[1][unit_ptr[1][u2]+ax];

							at[a].x += shift[0];
							at[a].y += shift[1];
							at[a].z += shift[2];
							printATOM( stdout,
							        at[a].bead, 
							        at[a].res, 
								at+a );
							at[a].x -= shift[0];
							at[a].y -= shift[1];
							at[a].z -= shift[2];
						}
						printf("END\n");
					} 
				}
			}



			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}

		fclose(dcdFile);
	}

#ifdef INFO_MODE
	for( foundABond *bond = allBonds; bond; bond = bond->next )
	{
		printf("hbond %s %s %d %d res %d %d count: %d dcd %s frame %d a1/2 %d %d\n",
			bond->atName1, bond->atName2, bond->atIndex[0], bond->atIndex[1], bond->resNum[0], bond->resNum[1], bond->count,	
				bond->example_dcd, bond->example_f, bond->example_abs[0], bond->example_abs[1] );	
	}
#endif
}





















