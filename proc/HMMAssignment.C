#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "proc_definition.h"

int micelle_nthresh = 6;
int thresh_ncontacts_interior = 30;
int thresh_ncontacts_connection = 3;

static double cutoff = 14.0;

double getChi2Symm( double *, double *, int ); 

struct tracked_pair
{
	char unit1[256];
	char unit2[256];
	char *state_seq;
	char *obs_seq;
	struct tracked_pair *next;

	int active;
	int inactive;
	int res1;	
	int res2;

	int len;
	int cntr;

	char *start_file;
	int start_frame;

	char *end_file;
	int end_frame;
};


int main( int argc, char **argv )
{
	if( argc < 4 )
	{
		printf("Syntax: HMMAssignment psf hmm.decoded stride definition.inp  [mode=config] [mode=multistate] [mode=only_donor] dcd1 [dcd2 ...]\n");
		printf("(Only donor only outputs a state if you are a hydrogen bond donor).\n");
		return 0;
	}

	char *buffer = (char *)malloc( sizeof(char)* 1000000 );
	FILE *psfFile = fopen(argv[1], "r" );
	int only_donor = 0;
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

	struct tracked_pair *pairs = NULL;

	FILE *hmm = fopen(argv[2],"r");
		
	int max_len = 1000000;
	char *hmm_seq_read = (char *)malloc( sizeof(char) * max_len );
	char *class_seq_read = (char *)malloc( sizeof(char) * max_len );


	int prioritize_by_kinetics = 0;

	rewind(hmm);
	
	getLine(hmm,buffer);	
	getLine(hmm,buffer);	
	int nstates;
	int nr = sscanf( buffer, "Setting up system with %d states", &nstates );

	double *kinetics=(double*)malloc( sizeof(double) * nstates );

	const char *search_string = "Final sorted transition probabilities";
	while( !feof(hmm) && strncasecmp( buffer, search_string, strlen(search_string) ) )
		getLine(hmm,buffer);
	for( int s = 0; s < nstates; s++ )
	{
		getLine( hmm, buffer );
		double ss[nstates];
		int nr = readNDoubles( buffer, ss, nstates );
		if( nr != nstates )
		{
			printf("Failed to read transition probabilities from line '%s'.\n", buffer );
			exit(1);
		}
		
		kinetics[s] = ss[s];
	}
	rewind(hmm);

	int sorter[nstates];

	for( int s = 0; s < nstates; s++ )
		sorter[s] =s;

	int done_sorting = 0;

	while( 0 && !done_sorting )
	{
		done_sorting = 1;

		for( int s = 1; s < nstates-1; s++ )
		{
			if( kinetics[sorter[s+1]] > kinetics[sorter[s]] )
			{
				int t = sorter[s];
				sorter[s] = sorter[s+1];
				sorter[s+1] = t;
				done_sorting = 0;
			}
		}
	}
	while( !feof(hmm) )
	{
		getLine( hmm, buffer );

		if( feof(hmm) ) break;

		int found_pound = 0;

		const char * t = buffer;

		while( *t && *t != '#' ) t+=1;

		if( *t != '#' )
			continue;
	
		char * h = buffer;

		int hmm_len = 0;
		int class_len = 0;

		while( *h && *h != ' ' && *h != '\t' && hmm_len < max_len )
		{
			hmm_seq_read[hmm_len] = *h;
			hmm_seq_read[hmm_len+1] = '\0';
			h+=1;
			hmm_len += 1;
		}

		while( *h == ' ' || *h == '\t' ) h+= 1;

		while( *h && *h != ' ' && *h != '\t' && hmm_len < max_len )
		{
			class_seq_read[class_len] = *h;
			class_seq_read[class_len+1] = '\0';
			h+=1;
			class_len += 1;
		}

		if( class_len != hmm_len )
		{
			printf("at line '%s' class_len %d != hmm_len %d.\n", buffer, class_len, hmm_len );
			exit(1);
		}

		char segid1[256];
		char segid2[256];
		segid1[0] = 'X'; segid1[1] = '\0';
		segid2[0] = 'X'; segid2[1] = '\0';
		char dcdNameS[256];
		char dcdNameE[256];
		char resname1[256], resname2[256];
		int frameS, frameE;
		int res1=-1, res2=-1;

		int err = 0;
		int arg_advance = 0;
		if( is_seg(0) )
		{
			int nr = sscanf(t, "# %s ", segid1 );
			if( nr != 1 ) err = 1;
			arg_advance = 2;
		}
		else
		{
			int nr = sscanf(t, "# %s %d ", resname1, &res1 );
			if( nr != 2 ) err = 1;
			arg_advance = 3;
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
			int nr = sscanf(t, "%s %d ", resname2, &res2 );
			if( nr != 2 ) err = 1;
			arg_advance = 2;
		}
		
		t = advance_string( t, arg_advance );	

		if( !strncasecmp( t, "nextra", 6 ) )
		{
			int nex;
			int nr = sscanf( t, "nextra %d", &nex );
			t = advance_string( t, 2 );	
	
			for( int e = 0; e < nex; e++ )
				t = advance_string( t, 3 );	
		}

		int nr = sscanf(t, "starting in dcd %s frame %d, ending in dcd %s frame %d.",
			dcdNameS, &frameS, dcdNameE, &frameE );			

		if( nr != 4 )
			err = 1;

		if( err )
		{
			printf("Format error in line '%s'.\n", buffer );
			exit(1);
		}

		struct tracked_pair *apair = (struct tracked_pair *)malloc( sizeof(struct tracked_pair) );

		apair->next = pairs;
		pairs = apair;

		apair->res1 = res1;
		apair->res2 = res2;
		apair->active = 0;
		apair->inactive = 0;
		apair->len = class_len;
		apair->cntr = 0;

		if( is_seg(0) ) 
			strcpy( apair->unit1, segid1 );
		else
			strcpy( apair->unit1, resname1 );

		if( is_seg(1) ) 
			strcpy( apair->unit2, segid1 );
		else
			strcpy( apair->unit2, resname2 );

		apair->obs_seq = (char *)malloc( sizeof(char)*(1+strlen(class_seq_read) ) );	
		apair->state_seq = (char *)malloc( sizeof(char)*(1+strlen(class_seq_read) ) );	
		strcpy( apair->state_seq, class_seq_read );
		strcpy( apair->obs_seq, hmm_seq_read );

		apair->start_file = (char *)malloc( sizeof(char)*(1+strlen(dcdNameS) ) );
		strcpy( apair->start_file, dcdNameS);
		apair->end_file = (char *)malloc( sizeof(char)*(1+strlen(dcdNameE) ) );
		strcpy( apair->end_file, dcdNameE);

		apair->start_frame = frameS;
		apair->end_frame = frameE;

	}
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
	int *unit_map = (int *)malloc( sizeof(int) * curNAtoms() ); 

	for( int u = 0; u < curNAtoms(); u++ )
		unit_map[u] = -1;

	int stride = atoi(argv[3]);
	
	int arg_offset = 0;
	int mode = 0; // standard mode, the HMM state is written

	load_pair_definition( argv[4] );
	
	cutoff = getCutoff();

	if( argc > 5 && !strncasecmp( argv[5], "mode=config", 11 ) )
	{	
		// config mode: a special format is used:
		// first a number indicating the number of simultaneous complexes
		// then, the letter codes, ie.
		// 001M2AA1B00 ...
	
		mode = 1;		
		arg_offset += 1;
	}
	
	if( argc> 5 && !strncasecmp( argv[5], "mode=multistate", strlen("mode=multistate") ) )
	{	
		// config mode: a special format is used:
		// first a number indicating the number of simultaneous complexes
		// then, the number codes. confusing to look at because it must be interpreted sequentially
		// 0011201321012...
		// is 0 0 1:1 2:01 3:210 1:2	

		mode = 2;		
		arg_offset += 1;
	}
	
	if( argc> 5 && !strncasecmp( argv[5], "mode=only_donor", strlen("mode=only_donor") ) )
	{	
		only_donor = 1;
		arg_offset += 1;
	}

	int tot_frame = 0;

	// for asymmetric pairs we don't want to label a pair twice, for example, in A <-> A|B|C. This marks the pair as being used.
	int *marked = (int *)malloc( sizeof(int) * nat );
	for( int u = 0; u < nat; u++ )
		marked[u] = -1;

	int real_units_2 = 0;

	for( int c = 5+arg_offset; c < argc; c++ )
	{
		FILE *dcdFile = fopen( argv[c], "r");
	
		readDCDHeader(dcdFile);
		setAligned();

		int nframes = curNFrames();
		fprintf(stderr, "Processing %s.\n", argv[c] );
		for( int f = 0; f < nframes; f++, tot_frame++)
		{
			double La, Lb, Lc;
			double alpha,beta,gamma;
			
			loadFrame( dcdFile, at );
			
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

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
							//printf("Adding unit %s.\n", at[a].resname );	

							if( marked[a] >= 0 ) // which element of unit[0] this is.
								unit_map[cur_unit] = marked[a];
							else
							{
								if( pass == 1 )
									real_units_2 += 1; 
								marked[a] = cur_unit;
							}
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


				init_done = 1;
			} 

			int nunits_tot = nunits[0] + real_units_2; 
		
			if( is_sym() ) nunits_tot = nunits[0];

			int max_simul = 9;
			char *cur_state = (char *)malloc( sizeof(char) * (nunits_tot+1) );
			char *hmm_state = (char *)malloc( sizeof(char) * (nunits_tot+1) * max_simul );
			char *multi_state = (char *)malloc( sizeof(char) * (nunits_tot+1) * max_simul );
			int *nhmm = (int *)malloc( sizeof(int) * nunits_tot );
			for( int u = 0; u < nunits_tot; u++ )
				nhmm[u] = 0;
			for( int u = 0; u < nunits_tot; u++ )
				cur_state[u] = 'X';
			cur_state[nunits_tot] = '\0';

			// activate/deactivate
			for( struct tracked_pair *apair = pairs; apair; apair = apair->next )
			{
				if( apair->inactive ) continue;

				if( !apair->active && !strcasecmp( argv[c], apair->start_file) )
				{
					if( f >= apair->start_frame )
					{
						apair->active = 1;
					}
				}				
				
				if( !apair->inactive && !strcasecmp( argv[c], apair->end_file) )
				{
					if( f >= apair->end_frame )
						apair->active = 0;
				}

			}

			for( struct tracked_pair *apair = pairs; apair; apair = apair->next )
			{
				if( !apair->active || apair->inactive ) continue;

				int npasses = 2;
		//		if( is_sym() ) npasses = 1;

				int u_cntr = 0;

				for( int pass = 0; pass < npasses; pass++ )
				{
					if( is_sym() )
						u_cntr = 0;

				for( int u = 0; u < nunits[pass]; u++ )
				{
					int uoff = u_cntr;
//					if( pass == 1 ) uoff += nunits[0];


					if( is_sym() )
					{
						u_cntr++;
					}
					else
					{
						if( pass == 1 )
						{
							if( unit_map[u] >= 0 )
								uoff = unit_map[u];
							else
								u_cntr++;
						}
						else 
							u_cntr++;
					}

					int match = 0;

					if( is_seg(pass) )
					{
						if( pass ==0 && !strcasecmp( apair->unit1, at[atom_list[pass][unit_ptr[pass][u]]].segid ) ) match = 1; 	
						if( pass ==1 && !strcasecmp( apair->unit2, at[atom_list[pass][unit_ptr[pass][u]]].segid ) ) match = 1; 	
					}
					else
					{
						if( pass ==0 && !strcasecmp( apair->unit1, at[atom_list[pass][unit_ptr[pass][u]]].resname ) &&	
							apair->res1 == at[atom_list[pass][unit_ptr[pass][u]]].res ) match = 1; 	
						if( pass ==1 && !strcasecmp( apair->unit2, at[atom_list[pass][unit_ptr[pass][u]]].resname ) && 
							apair->res2 == at[atom_list[pass][unit_ptr[pass][u]]].res ) match = 1; 	
					}
	
					if( !match ) continue;

					if( only_donor )
					{
#if 0
						for( int a1 = first_atom[0][u]; a1 < first_atom[0][u] + num_atom[0][u]; a1++ ) 
						for( int a2 = first_atom[1][u2]; a2 < first_atom[1][u2] + num_atom[1][u2]; a2++ ) 
						{
							if( !(at[a1].atname[0] == 'O' && at[a2].atname[0] == 'H' ) &&
							    !(at[a2].atname[0] == 'O' && at[a1].atname[0] == 'H' ) )
								continue;
						
							// only donors:
							if( pass == 0 && at[a1].atname[0] != 'H' ) continue;
							if( pass == 1 && at[a2].atname[0] != 'H' ) continue;
	
							double adr[3] = { at[a1].x - at[a2].x, at[a1].y - at[a2].y, at[a1].z - at[a2].z };	

							while( adr[0] < -La/2 )adr[0] += La;
							while( adr[0] > La/2 ) adr[0]-= La;
							while( adr[1] < -Lb/2 )adr[1] += Lb;
							while( adr[1] > Lb/2 ) adr[1]-= Lb;
							while( adr[2] < -Lc/2 )adr[2] += Lc;
							while( adr[2] > Lc/2 ) adr[2]-= Lc;

							double r = sqrt(adr[0]*adr[0]+adr[1]*adr[1]+adr[2]*adr[2]);

							if( r < hbond_cut )
							{
								int bond_index = binaryHBonder( at+a1, at+a2 );
								if( bond_index >= 0 && !(binary_hbond & (1<<bond_index)) )
									binary_hbond += (1<<bond_index); 
							}
						}
				
						if( binary_hbond == 0 ) continue;
#endif
					}

					if( apair->cntr < apair->len )
					{
						if( nhmm[uoff] < max_simul )
						{
							hmm_state[uoff*max_simul+nhmm[uoff]]  = apair->obs_seq[apair->cntr];
							multi_state[uoff*max_simul+nhmm[uoff]] = apair->state_seq[apair->cntr];

							nhmm[uoff]++;	
							
						}
/*						else
						{
							printf("Exceeded max_simul.\n");
							exit(1);
						}
*/
						switch( apair->state_seq[apair->cntr] )
						{	
							case '0':
								if( cur_state[uoff] == 'X' )
									cur_state[uoff] = apair->state_seq[apair->cntr];
								break;
							default:
								int state_off = cur_state[uoff] - '0';	
								int new_state_off = apair->state_seq[apair->cntr] - '0';	

								if( cur_state[uoff] == 'X' || cur_state[uoff] == '0' || (sorter[new_state_off] < sorter[state_off] && new_state_off != 0) )
									cur_state[uoff] = apair->state_seq[apair->cntr];
								break;
						}
					}
				}
				}
			}

			for( struct tracked_pair *apair = pairs; apair; apair = apair->next )
			{
				if( apair->active && !apair->inactive )
					apair->cntr += 1;
			}

			if( mode == 0 )
			{
			if( tot_frame % stride == 0 )
				printf("%s\n",cur_state );
			}
			else if( mode == 1 ) 
			{
				for( int u = 0; u < nunits_tot; u++ )
				{
					printf("%d", nhmm[u] );
					for( int t = 0; t < nhmm[u]; t++ )
						printf("%c", hmm_state[u*max_simul+t]); 
				}
				printf("\n");
			}
			else if( mode == 2 ) 
			{
				for( int u = 0; u < nunits_tot; u++ )
				{
					if( nhmm[u] >= 10 )
					{
						printf("OH NO!\n");
						exit(1);
					}
					printf("%d", nhmm[u] );
					for( int t = 0; t < nhmm[u]; t++ )
						printf("%c", multi_state[u*max_simul+t]); 
				}
				printf("\n");
			}
			free(cur_state);
			free(hmm_state);
			free(nhmm);
			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}

		fclose(dcdFile);
	}
				
}



















