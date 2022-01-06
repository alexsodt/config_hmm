#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dcd.h"
#include "util.h"
#include "comparison.h"
#include "proc_definition.h"

extern int debug_trigger;
extern int current_letter;

double normalize( double *dr )
{
	double lr = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

	if( fabs(lr) > 0  )
	{
		dr[0] /= lr;
		dr[1] /= lr;
		dr[2] /= lr;
	}

	return lr;
}
double cross( double *dr1, double *dr2, double *cp )
{
	cp[0] = (dr1[1] * dr2[2] - dr1[2] * dr2[1]);
	cp[1] =-(dr1[0] * dr2[2] - dr1[2] * dr2[0]);
	cp[2] = (dr1[0] * dr2[1] - dr1[1] * dr2[0]);

	double l = sqrt( cp[0]*cp[0]+cp[1]*cp[1]+cp[2]*cp[2]);

	return l;
}


static double cutoff = 14.0;

double getChi2Symm( double *, double *, int ); 

const char *require_amide_donor = "HNF";
const char *require_amide_acceptor = "OF";
const char *require_hydroxyl_donor = "HO3";
const char *require_hydroxyl_acceptor = "O3";

double hbond_cut = 2.5;

#define MAX_EXTRA	10

struct tracked_pair
{
	int unit1;
	int unit2;
	int res1;
	int res2;
	char *hmm_seq;
	int hmm_len;
	int hmm_space;
	int nextra;
	int extra[MAX_EXTRA];
	struct tracked_pair *next;
	int flag;
	char *start_file;
	int start_frame;
};


int main( int argc, char **argv )
{
	if( argc < 5 )
	{
		printf("Syntax: extractDimerHMM psf centers.pdb definition.inp dcd1 [dcd2 ...]\n");
		return 0;
	}

	char *buffer = (char *)malloc( sizeof(char)* 1000000 );
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

	struct atom_rec *dummy = (struct atom_rec *)malloc( sizeof( struct atom_rec ) * curNAtoms() );
	loadDummy( dummy );

	aStructure *structs = NULL;
	int nstructs = 0;

	char **struct_atom_names = NULL;
	
	int hmm_nat = 0;
	{
		FILE *thePDB = fopen(argv[2],"r");
		
		if( !thePDB )
		{
			printf("Couldn't open centers file '%s'.\n", argv[2] );
			return -1;
		}

		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;
	
			if( !strncasecmp( buffer, "ATOM", 4 ) ) hmm_nat++;
	
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}

		struct_atom_names = (char **)malloc( sizeof(char*) * hmm_nat );
		rewind(thePDB);	
		int ta =0;
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
	
			if( feof(thePDB) ) break;

			if( !strncasecmp( buffer, "REMARK CODE", 11) )
			{
				
			}
	
			if( !strncasecmp( buffer, "ATOM", 4 ) ){
				struct atom_rec tat;
				readATOM( buffer, &tat);
				struct_atom_names[ta] = (char *)malloc( sizeof(char) * (1 + strlen(tat.atname) ) );
				strcpy( struct_atom_names[ta], tat.atname );
			 ta++;
			}
			if( !strncasecmp( buffer, "TER", 3) || 
			    !strncasecmp( buffer, "END", 3) )
				break;
		}

		fclose(thePDB);
	}

	use_atoms = hmm_nat;

	{
		int nstructsSpace = 100;
		structs = (aStructure *)malloc( sizeof(aStructure) * nstructsSpace  );
		char **remarks = (char **)malloc( sizeof(char *) * nstructsSpace );
		nstructs = 0;
		int cat = 0;
	
		for( int x = 0; x < nstructsSpace; x++ )
		{
			remarks[x] = NULL;	
			structs[x].atoms = ( double*)malloc( sizeof(double) * 3 * hmm_nat );
		}
	
		FILE *thePDB = fopen(argv[2],"r");
		
		while( !feof(thePDB) )
		{
			getLine( thePDB, buffer );
		
			if( feof(thePDB) ) {
				 break;	
			}
			
			if( !strncasecmp( buffer, "REMARK CODE", 11 ) )
			{
				long the_code;
				sscanf( buffer, "REMARK CODE %lu", &the_code );

				structs[nstructs].binary_data = the_code;
			}	
			else if( !strncasecmp( buffer, "REMARK", 6 ) )
			{
				remarks[nstructs] = (char *)malloc( sizeof(char) * (1+strlen(buffer) ) );
				strcpy( remarks[nstructs], buffer );
			}
		
			if( !strncasecmp( buffer, "ATOM", 4 ) )
			{
	
				if( hmm_nat == cat )
				{
					printf("ATOM # inconsistency in file %s, structure %d.\n", argv[2], nstructs );
					exit(1);
				}
				struct atom_rec tat;
		
				readATOM( buffer, &tat );
		
				structs[nstructs].atoms[3*cat+0] = tat.x;	
				structs[nstructs].atoms[3*cat+1] = tat.y;	
				structs[nstructs].atoms[3*cat+2] = tat.z;	
				cat++;
			}
	
			if( !strncasecmp( buffer, "END", 3) || !strncasecmp( buffer, "TER", 3) )
			{
				nstructs++;
	
				if( nstructs == nstructsSpace )
				{
					nstructsSpace *= 2;
		
					structs = (aStructure *)realloc( structs, sizeof(aStructure) * nstructsSpace );
					remarks = (char **)realloc( remarks, sizeof(char *) * nstructsSpace );
					for( int x = nstructs; x < nstructsSpace; x++ )
					{
						structs[x].atoms = (double *)malloc( sizeof(double) * 3 * hmm_nat );
						remarks[x] = NULL;	
					}
				}	
				cat = 0;
			}
		}
	}

	int *has_amide_hbond = (int *)malloc( sizeof(int) * nstructs );
	memset( has_amide_hbond, 0, sizeof(int) * nstructs );
	int *has_mixed_hbond = (int *)malloc( sizeof(int) * nstructs );
	memset( has_mixed_hbond, 0, sizeof(int) * nstructs );
	int *has_hydroxyl_hbond = (int *)malloc( sizeof(int) * nstructs );
	memset( has_hydroxyl_hbond, 0, sizeof(int) * nstructs );

	for( int s = 0; s < nstructs; s++ )
	{
		for( int x1 = 0; x1 < hmm_nat; x1++ )
		for( int x2 = x1+1; x2 < hmm_nat; x2++ ) 
		{
			if( 
				(!strcasecmp( struct_atom_names[x1], require_amide_donor) && !strcasecmp( struct_atom_names[x2], require_amide_acceptor)) || 
				(!strcasecmp( struct_atom_names[x2], require_amide_donor) && !strcasecmp( struct_atom_names[x1], require_amide_acceptor)) )
			{ 
				double dr[3] = { 
					structs[s].atoms[x1*3+0] - structs[s].atoms[x2*3+0], 
					structs[s].atoms[x1*3+1] - structs[s].atoms[x2*3+1], 
					structs[s].atoms[x1*3+2] - structs[s].atoms[x2*3+2] };
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				if( r < hbond_cut )
					has_amide_hbond[s] = 1; 
			}
			
			if( 
				(!strcasecmp( struct_atom_names[x1], require_hydroxyl_donor) && !strcasecmp( struct_atom_names[x2], require_hydroxyl_acceptor)) || 
				(!strcasecmp( struct_atom_names[x2], require_hydroxyl_donor) && !strcasecmp( struct_atom_names[x1], require_hydroxyl_acceptor)) )
			{ 
				double dr[3] = { 
					structs[s].atoms[x1*3+0] - structs[s].atoms[x2*3+0], 
					structs[s].atoms[x1*3+1] - structs[s].atoms[x2*3+1], 
					structs[s].atoms[x1*3+2] - structs[s].atoms[x2*3+2] };
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				if( r < hbond_cut && r > 1.3 )
					has_hydroxyl_hbond[s] = 1; 
			}
			
			if( 
				(!strcasecmp( struct_atom_names[x1], require_amide_donor) && !strcasecmp( struct_atom_names[x2], require_hydroxyl_acceptor)) || 
				(!strcasecmp( struct_atom_names[x1], require_hydroxyl_donor) && !strcasecmp( struct_atom_names[x2], require_amide_acceptor)) || 
				(!strcasecmp( struct_atom_names[x2], require_amide_donor) && !strcasecmp( struct_atom_names[x1], require_hydroxyl_acceptor)) || 
				(!strcasecmp( struct_atom_names[x2], require_hydroxyl_donor) && !strcasecmp( struct_atom_names[x1], require_amide_acceptor)) )
			{ 
				double dr[3] = { 
					structs[s].atoms[x1*3+0] - structs[s].atoms[x2*3+0], 
					structs[s].atoms[x1*3+1] - structs[s].atoms[x2*3+1], 
					structs[s].atoms[x1*3+2] - structs[s].atoms[x2*3+2] };
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				if( r < hbond_cut )
					has_mixed_hbond[s] = 1; 
			}
		}
	}


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
	
	// qualifies as an extra atom for some unit:
	int *qual_extra = (int *)malloc( sizeof(int) * curNAtoms() );
	memset( qual_extra, 0, sizeof(int) * curNAtoms() );
	int nqual = 0; // this way we can loop over this set of atoms only, looking for nearby atoms.


	int cntr = 0;

	struct tracked_pair *pairs = NULL;

	FILE *hmmFile = fopen("hmm.inp","w");


	fprintf(hmmFile, "COMMAND TRAIN\n");

	int last_nframes = 0;
	if( argc < 4 )
	{
		char amide_class[nstructs+1];
		char hydroxyl_class[nstructs+1];
		char mixed_class[nstructs+1];
	
		int namide = 0;
		int nhydroxyl = 0;
		int nmixed = 0;

		amide_class[namide] = '\0';
		hydroxyl_class[namide] = '\0';
		mixed_class[namide] = '\0';

		for( int s = 0; s < nstructs; s++ )
		{
			printf("State %c amide? %c hydroxyl? %c mixed? %c.\n", (s < 26 ? 'A'+s : 'a'+s-26 ), 
				has_amide_hbond[s] ? 'Y' : 'N',
				has_hydroxyl_hbond[s] ? 'Y' : 'N',
				has_mixed_hbond[s] ? 'Y' : 'N'
			 );
		
			if( has_amide_hbond[s] )
			{
				amide_class[namide] = (s < 26 ? 'A'+s : 'a'+s-26 );
				namide++;
				amide_class[namide] = '\0';
			}
			if( has_hydroxyl_hbond[s] )
			{
				hydroxyl_class[nhydroxyl] = (s < 26 ? 'A'+s : 'a'+s-26 );
				nhydroxyl++;
				hydroxyl_class[nhydroxyl] = '\0';
			}
			if( has_mixed_hbond[s] )
			{
				mixed_class[nmixed] = (s < 26 ? 'A'+s : 'a'+s-26 );
				nmixed++;
				mixed_class[nmixed] = '\0';
			}
		}
		printf("Amide:    %s\n", amide_class );
		printf("Hydroxyl: %s\n", hydroxyl_class );
		printf("Mixed:    %s\n", mixed_class );
		return -1;
	}
	
	int arg_offset = 0;
	
	load_pair_definition( argv[3] );

	cutoff = getCutoff();

	binary_penalty = pair_binary_penalty();
	binary_benefit = pair_binary_benefit();
	hbond0_penalty = pair_hbond0_penalty();
	hbond0_benefit = -hbond0_penalty;

	/**** DOPE DOPE hack: remove for distribution ****/
	double *av_angle = (double *)malloc( sizeof(double) * 52  );
	double *nav_angle = (double *)malloc( sizeof(double) *  52 );
	memset( av_angle, 0, sizeof(double) * 52 );
	memset( nav_angle, 0, sizeof(double) * 52 );
	int activate_hack = 0;
	/**** end hack */

	for( int c = 4; c < argc; c++ )
	{
		FILE *dcdFile = fopen( argv[c], "r");
	
		if( !dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return -1;
		}

		readDCDHeader(dcdFile);
		setAligned();

		int nframes = curNFrames();
		last_nframes = nframes;
		fprintf(stderr, "Processing %s.\n", argv[c] );
		for( int f = 0; f < nframes; f++,  n_frame_add++, cntr++)
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
				
				for( int a = 0; a < curNAtoms(); a++ )
				{
					if( used_as_ion( at[a].res, at[a].resname, at[a].segid, NULL ) >= 0 ) 
						qual_extra[nqual++] = a;	
				}

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

				for( int pass = 0; pass < 2; pass ++ )
				for( int u = 0; u < nunits[pass]; u++ )
				{
					double n_pos = 0;
					double av_c_pos = 0;
					double nav_c_pos = 0;

					for( int a = first_atom[pass][u]; a < first_atom[pass][u] + num_atom[pass][u]; a++ )
					{
						if( at[a].atname[0] == 'C' && at[a].atname[1] != '1' && strlen(at[a].atname) > 2 )
						{
							av_c_pos += at[a].z;	
							nav_c_pos += 1;
						}
						if( at[a].atname[0] == 'N' )
							n_pos = at[a].z;
						
					}

					av_c_pos /= nav_c_pos;

				}


				if( nextra() > MAX_EXTRA) { printf("Increase MAX_EXTRA and recompile.\n"); exit(1); }

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

			for( struct tracked_pair * pair = pairs; pair; pair = pair->next )
				pair->flag = 0;

			for( int u = 0; u < nunits[0]; u++ )
			{
				for( int u2 = 0; u2 < nunits[1]; u2++ )
				{
					if( is_sym() && u2 <= u ) continue;
					
					if( first_atom[0][u] == first_atom[1][u2] ) continue;								

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
					
					int ions_covered = 1;
					int num = nextra();
					int extra_to_print[1+num];
					double cenp[3] = { unit_r[0][3*u+0] + dr[0]/2, 
							   unit_r[0][3*u+1] + dr[1]/2,		
							   unit_r[0][3*u+2] + dr[2]/2 };		
	

					if( r < cutoff && num > 0 )
					{
						for( int i = 0; i < num; i++ )
							extra_to_print[i] = -1;


						for( int qx = 0; qx < nqual; qx++ )
						{
							int q = qual_extra[qx];

							double max_r;
							double qmr = -1;
							int ex = used_as_ion( at[q].res, at[q].resname, at[q].segid, &max_r );

							double dr[3] = { cenp[0] - at[q].x, cenp[1] - at[q].y, cenp[2] - at[q].z };
	
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
		
							if( qmr < 0 || r < qmr )
							{
								qmr = r;
 								if( qmr < max_r )
									extra_to_print[ex] = q;
							}
						}

						for( int t = 0; t < num; t++ )
						{
							if( extra_to_print[t] < 0 )
								ions_covered = 0;
						}

					}

					if( r < cutoff && (ions_covered || allowNoIons())  )
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

							if( r < hbond_cut )
							{
								int bond_index = binaryHBonder( at+a1, at+a2 );
								if( bond_index >= 0 && !(binary_hbond & (1<<bond_index)) )
									binary_hbond += (1<<bond_index); 
							}
						}

						if( force_hbonds() && binary_hbond == 0 )
							continue;

						double *tcen = (double *) malloc( sizeof(double) * 3 * hmm_nat );
	
						int tp = 0;
						for( int ax = 0; ax < unit_size[0][u]; ax++ )
						{
							tcen[3*tp+0] = at[atom_list[0][unit_ptr[0][u]+ax]].x;
							tcen[3*tp+1] = at[atom_list[0][unit_ptr[0][u]+ax]].y;
							tcen[3*tp+2] = at[atom_list[0][unit_ptr[0][u]+ax]].z;
							tp++;
						}
						
						for( int ex = 0; ex < num; ex++ ) 
						{
							int q = extra_to_print[ex];

							if( q >= 0 )
							{
							double dr[3] = { 
								at[q].x - cenp[0],	
								at[q].y - cenp[1] ,	
								at[q].z - cenp[2] };	
								
							while( dr[0] < -La/2 )dr[0] += La;
							while( dr[0] > La/2 ) dr[0]-= La;
							while( dr[1] < -Lb/2 )dr[1] += Lb;
							while( dr[1] > Lb/2 ) dr[1]-= Lb;
							while( dr[2] < -Lc/2 )dr[2] += Lc;
							while( dr[2] > Lc/2 ) dr[2]-= Lc;

							double t[3] = { at[q].x, at[q].y, at[q].z };

							tcen[3*tp+0] = cenp[0] + dr[0];
							tcen[3*tp+1] = cenp[1] + dr[1];
							tcen[3*tp+2] = cenp[2] + dr[2];
							}
							else
							{
								tcen[3*tp+0] = 9999.99;
								tcen[3*tp+1] = 9999.99;
								tcen[3*tp+2] = 9999.99;
							}
							tp++;
						} 

						for( int ax = 0; ax < unit_size[1][u2]; ax++ )
						{
							tcen[3*tp+0] = at[atom_list[1][unit_ptr[1][u2]+ax]].x + shift[0];
							tcen[3*tp+1] = at[atom_list[1][unit_ptr[1][u2]+ax]].y + shift[1];
							tcen[3*tp+2] = at[atom_list[1][unit_ptr[1][u2]+ax]].z + shift[2];
							tp++;
						}


						aStructure tcen_struct;

						tcen_struct.binary_data = binary_hbond;
						tcen_struct.atoms = tcen;
						tcen_struct.cluster = 0;

						double best_chi2 = 1e10;
						int best_match = -1;

						int r1 = at[first_atom[0][u]].res;
						int r2 = at[first_atom[1][u2]].res;

						if( f == 132 && (r1 == 1108 && r2 == 1225) || (r1 == 1225 && r2 == 1108 ) )
						{
//							printf("Debug this.\n");
//							debug_trigger = 1;
						} 

						for( int pass = 0; pass < 2; pass++ )
						{
							for( int sc = 0; sc < nstructs; sc++ )
							{
								current_letter = sc;

//								double lchi2 = getChi2Symm( tcen, structs + sc * 3 * hmm_nat, hmm_nat );	
								double lchi2 = value_comparison( &tcen_struct, structs+sc, is_sym()  ); 

								if( lchi2 < best_chi2 )
								{
									best_chi2 = lchi2;
									best_match = sc;
								}
							}	

							if( best_match >= 0 ) break;
						}
						struct tracked_pair *gotit = NULL;

						for( struct tracked_pair *pair = pairs; pair; pair = pair->next )
						{
							if( is_sym() )
							{
								if( 
									(pair->unit1 == u && pair->unit2 == u2) ||
									(pair->unit1 == u2 && pair->unit2 == u) )
									gotit = pair;
							}
							else
							{
								if( 
									pair->unit1 == u && pair->unit2 == u2 )
									gotit = pair;
							}
						}				

						if( !gotit )
						{
							gotit = (struct tracked_pair *)malloc( sizeof(struct tracked_pair) );
							gotit->next = pairs;
							pairs = gotit;
							gotit->res1 = at[atom_list[0][unit_ptr[0][u]]].res;
							gotit->res2 = at[atom_list[1][unit_ptr[1][u2]]].res;
							gotit->unit1 = u;
							gotit->unit2 = u2;
							gotit->hmm_space = 10;
							gotit->hmm_len = 0;
							gotit->hmm_seq = (char *)malloc( sizeof(char) * (1 + gotit->hmm_space) );
							gotit->start_file = (char *)malloc( sizeof(char) * (1 + strlen(argv[c]) ) );
							strcpy( gotit->start_file, argv[c] );
							gotit->start_frame = f;

							gotit->nextra = num;
							for( int e = 0; e < num; e++ )
								gotit->extra[e] = extra_to_print[e];
			
						}
							
						gotit->flag = 1;	
							
						// Find which cluster center we are closest to.

						if( gotit->hmm_len == gotit->hmm_space )
						{
							gotit->hmm_space *= 2;
							gotit->hmm_seq = (char *)realloc( gotit->hmm_seq, sizeof(char) * ( 1 + gotit->hmm_space ) );
						}

						if( best_match >= 26 )
							gotit->hmm_seq[gotit->hmm_len] = 'a' + (best_match - 26);
						else
							gotit->hmm_seq[gotit->hmm_len] = 'A' + best_match;
			
						if( gotit->hmm_seq[gotit->hmm_len] == '@' )
						{
						
						}

						gotit->hmm_seq[gotit->hmm_len+1] = '\0';
						gotit->hmm_len += 1;

						free(tcen);
					} 
				}
			}

			struct tracked_pair *prev = NULL;
			struct tracked_pair *next = NULL;
			for( struct tracked_pair *pair = pairs; pair; pair = next )
			{
				next = pair->next;

				if( pair->flag == 0 )
				{
					if( prev ) 
						prev->next = next;
					else
						pairs = next;

					if( strlen(pair->hmm_seq) > 20 )
					{
						char *state = (char *)malloc( sizeof(char) * (1 + pair->hmm_len) ) ;

						for( int t = 0; t < pair->hmm_len; t++ )
							state[t] = '.';
						state[pair->hmm_len-1] = 'x';
						state[pair->hmm_len] = '\0';

						fprintf(hmmFile, "%s %s #", pair->hmm_seq, state );

						if( is_seg(0) == 0 )
							fprintf(hmmFile, " %s %d", 
								at[atom_list[0][unit_ptr[0][pair->unit1]]].resname,
								at[atom_list[0][unit_ptr[0][pair->unit1]]].res );
						else
							fprintf(hmmFile, " %s", 
								at[atom_list[0][unit_ptr[0][pair->unit1]]].segid );

						if( is_seg(1) == 0 )
							fprintf(hmmFile, " %s %d", 
								at[atom_list[1][unit_ptr[1][pair->unit2]]].resname,
								at[atom_list[1][unit_ptr[1][pair->unit2]]].res );
						else
							fprintf(hmmFile, " %s", 
								at[atom_list[1][unit_ptr[1][pair->unit2]]].segid );

						fprintf(hmmFile, " nextra %d", nextra() ); 
						
						for( int i = 0; i <  nextra(); i++ )
						{
							if( pair->extra[i] >= 0 )
								fprintf(hmmFile, " %s %s %d", at[pair->extra[i]].segid, at[pair->extra[i]].resname, at[pair->extra[i]].res );						
							else
								fprintf(hmmFile, " NULL NULL 0");
	
						}
						fprintf(hmmFile, " starting in dcd %s frame %d, ending in dcd %s frame %d.\n",
							pair->start_file,
							pair->start_frame,
							argv[c],
							f );
						free(state);
					}
					free( pair->start_file );
					free( pair );
				}
				else
					prev = pair;
			}


			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		}

		fclose(dcdFile);
	}
				
	for( struct tracked_pair *pair = pairs; pair; pair = pair->next)
	{
		char *segname1;
		char *segname2;

		char *resname1;
		char *resname2;

		getSegName( atom_list[0][unit_ptr[0][pair->unit1]], &segname1 );
		getSegName( atom_list[1][unit_ptr[1][pair->unit2]], &segname2 );
		
		getResName( atom_list[0][unit_ptr[0][pair->unit1]], &resname1 );
		getResName( atom_list[1][unit_ptr[1][pair->unit2]], &resname2 );

		int res1 = pair->res1; 
		int res2 = pair->res2; 

		if( pair->hmm_len > 20 )
		{
			char *state = (char *)malloc( sizeof(char) * (1 + pair->hmm_len) ) ;

			for( int t = 0; t < pair->hmm_len; t++ )
				state[t] = '.';
			state[pair->hmm_len] = '\0';

			fprintf(hmmFile, "%s %s #", pair->hmm_seq, state );

			if( is_seg(0) == 0 )
				fprintf(hmmFile, " %s %d", resname1, res1 );
			else
				fprintf(hmmFile, " %s", 
					segname1 );

			if( is_seg(1) == 0 )
				fprintf(hmmFile, " %s %d", resname2, res2 ); 
			else
				fprintf(hmmFile, " %s",  segname2 );
						
			fprintf(hmmFile, " nextra %d", nextra() ); 
						
			for( int i = 0; i <  nextra(); i++ )
			{
				if( pair->extra[i] >= 0 )
					fprintf(hmmFile, " %s %s %d", dummy[pair->extra[i]].segid, dummy[pair->extra[i]].resname, dummy[pair->extra[i]].res );						
				else
					fprintf(hmmFile, " NULL NULL 0");
			}	
			fprintf(hmmFile, " starting in dcd %s frame %d, ending in dcd %s frame %d (unterminated).\n",
				pair->start_file,
				pair->start_frame,
				argv[argc-1],
				last_nframes );
			free(state);	
		}
		free(segname1);
		free(segname2);	
	}
	fprintf(hmmFile, "STOP\n");			
	
	fclose(hmmFile);

	if( activate_hack  )
	{
		FILE *dopeData = fopen("dope_data.txt","w");

		for( int  i = 0;  i < 52; i++ )
			fprintf(dopeData,"%c %lf %lf\n", (i < 26 ? 'A' + i : 'a' + (i-26) ), av_angle[i] / nav_angle[i], nav_angle[i] );
		fclose(dopeData);
	}
}





















