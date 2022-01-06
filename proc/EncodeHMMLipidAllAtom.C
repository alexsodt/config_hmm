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
#define SUB_PRESS
#define LEN 100

#define LEAFLET_BOTH	0
#define LEAFLET_UPPER	1
#define LEAFLET_LOWER	2

//#define FORCE_LOLD 

#define REMOVE_COM
int N_TALLY = 6; 
#define TYPE_PROT 3
static char prot_state;

char code( int cnts[4] )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];

	if( done == 0 )
	{
		// extra state for any protein.
		array = (int *)malloc( sizeof(int) * ((N_TALLY+1)*(N_TALLY+1)+1) );

		char cur = 'A';

		for( int x = 0; x <= N_TALLY; x++ )
		for( int y = 0; y <= N_TALLY; y++ )
			array[x*(N_TALLY+1)+y] = 'a';

		// x y other, fast slow implicit, 1 0 2
		// 0 0 6
		// 1 0 5
		// 0 1 5

		int max_state = 0;

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

		prot_state = cur;
		done = 1;
	}	
	
	if( cnts[TYPE_PROT] > 0 )
		return prot_state;

	return array[cnts[0]*(N_TALLY+1)+cnts[1]];
}

char ext_code( int *cnts, int nres_types )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];

	int total = (N_TALLY+1);

	for( int x = 1; x < nres_types; x++ )
		total *= (N_TALLY+1);

	if( done == 0 )
	{
		// extra state for any protein.
		array = (int *)malloc( sizeof(int) * ((N_TALLY+1)*(N_TALLY+1)+1) );

		char cur = 'A';

		for( int x = 0; x <= N_TALLY; x++ )
		for( int y = 0; y <= N_TALLY; y++ )
			array[x*(N_TALLY+1)+y] = 'a';

		// x y other, fast slow implicit, 1 0 2
		// 0 0 6
		// 1 0 5
		// 0 1 5

		int max_state = 0;

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

		prot_state = cur;
		done = 1;
	}	
	
	if( cnts[TYPE_PROT] > 0 )
		return prot_state;

	return array[cnts[0]*(N_TALLY+1)+cnts[1]];
}

struct elem
{
	int i, j,k;
	double PXX,PYY,PZZ;
};

int main( int argc, char **argv )
{
	char buffer[4096];

	if( argc < 5 )
	{
		printf("Syntax: EncodeHMMLipid LipidsPerConcentrationPoint [upper|lower|both] psf dcd\n");	
		return 0;
	}

	N_TALLY = atoi(argv[1]);
	FILE *psfFile = fopen(argv[3], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[3] );
		return 0;
	}

	int leaflet_code = LEAFLET_BOTH;

	if( !strcasecmp( argv[2], "upper") )
		leaflet_code = LEAFLET_UPPER;
	else if( !strcasecmp( argv[2], "lower") )
		leaflet_code = LEAFLET_LOWER;

	if( !strcasecmp( argv[3] + strlen(argv[3])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );
	fclose(psfFile);	
	

	char *res_select = NULL;


	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );


	int np_link[100];
	int n_np = 0;

	int init_done = 0;

	int nat = curNAtoms();
	int type[nat];

	int grain = 20;

//	nframes = 1000;

	int nx = grain;
	int ny = grain;
	int nz = 1;

struct binnit
{
	double *pos;
	int *id;
	int npos;
	int nposSpace;
};
	binnit *bins = (binnit *)malloc( sizeof(binnit) * nx * ny * nz );
	for( int x = 0; x < nx*ny*nz; x++ )
	{
		bins[x].nposSpace = 100;
		bins[x].pos = (double*)malloc( sizeof(double) * bins[x].nposSpace * 3);
		bins[x].id = (int*)malloc( sizeof(int) * bins[x].nposSpace );
		bins[x].npos = 0; 
	}

	double *last_pos = (double *)malloc( sizeof(double) * 3 * curNAtoms() );
	memset( last_pos, 0, sizeof(double) * curNAtoms() * 3 );
	
	int did_init = 0;

	int *bin_lookup_x = (int *)malloc( sizeof(int) * curNAtoms() );
	int *bin_lookup_y = (int *)malloc( sizeof(int) * curNAtoms() );


	int tf = 0;

#define RESMAX 16
	char resNames[RESMAX][256];	
	int nres = 0;

#ifdef FORCE_LOLD
	nres = 3;
	strcpy( resNames[0], "CHL1");
	strcpy( resNames[1], "PSM");
	strcpy( resNames[2], "PLPC");
//	strcpy( resNames[3], "PROT");
#endif

	int *use_atoms = (int *)malloc( sizeof(int) * curNAtoms() );
	int nuse = 0;

	int state_size = curNAtoms();

	char *states = (char *)malloc( sizeof(char)* curNAtoms() * state_size  );

	int used_twice = 0;
	
	int nframes_use = 0;

	int *leaflet = (int *)malloc( sizeof(int) * curNAtoms() );

	for( int c = 4; c < argc; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");
	
		if( ! dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return 0;
		}
		readDCDHeader(dcdFile);
		setFractional();
		int nframes = curNFrames();


		for( int f = 0; f < nframes; f++, nframes_use++  )
		{
			if( nframes_use * nuse >= state_size )
			{	
				state_size = nframes_use * nuse * 2;
				states = (char *)realloc( states, sizeof(char) * state_size );
			} 

			double cur_com[2] = {0,0};
	
			double La, Lb, Lc;
			double alpha,beta,gamma;
			
			loadFrame( dcdFile, at );
	
	
			if( !DCDsuccess() )
			{
				nframes = f;
				break;
			}
	
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
	
			if( !init_done )
			{
				strcpy( resNames[TYPE_PROT], "PROT" );

#define N_BINS_MOLDIST 100
			
			        double best_chi2 = 1e10;
				double wrapto = 0;
				double moldist[N_BINS_MOLDIST];
				memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );

				for( int ax = 0; ax < curNAtoms(); ax++ )
				{
					if( at[ax].atname[0] != 'C' ) continue;

					int p = ax;

					double tz = at[p].z*Lc;
					
					while( tz < 0 ) tz += Lc;
					while( tz >= Lc ) tz -= Lc;
					
					int zb = N_BINS_MOLDIST * tz / Lc; // this is right
					if( zb < 0 ) zb = 0;
					if( zb >= N_BINS_MOLDIST ) zb = N_BINS_MOLDIST-1;
					moldist[zb] += 1;
				}
	
				for( int zb = 0; zb < N_BINS_MOLDIST; zb++ )
				{
				        double zv = Lc * (zb+0.5) / (double)N_BINS_MOLDIST;
				
				         int zlow  = zb- N_BINS_MOLDIST/2;
				         int zhigh = zlow + N_BINS_MOLDIST;
				
				         double Lchi2 = 0;
				         for( int iz = zlow; iz < zhigh; iz++ )
				         {
				                 double dz = Lc * (iz+0.5) / N_BINS_MOLDIST - zv;
				
				                 int iiz = iz;
				                 while( iiz < 0 ) iiz += N_BINS_MOLDIST;
				                 while( iiz >= N_BINS_MOLDIST ) iiz -= N_BINS_MOLDIST;
				
				                 Lchi2 += moldist[iiz] * (dz) * (dz);
				         }
				
				         if( Lchi2 < best_chi2 )
				         {
				                 best_chi2 = Lchi2;
				                 wrapto = zv/Lc;
				         }
				}
	
				// wrap around z periodic dimension

				for( int p = 0; p < curNAtoms(); p++ )	
				{
					while( at[p].z - wrapto < -1./2 ) at[p].z += 1;
					while( at[p].z - wrapto > 1./2 ) at[p].z -= 1;
	
					if( at[p].z > wrapto )
						leaflet[p] = 1;
					else
						leaflet[p] = -1;
				} 

	
				for( int a = 0; a < curNAtoms(); a++ )
				{

					if( !strcasecmp( at[a].resname, "TIP3" ) ) continue;
					if( !strcasecmp( at[a].resname, "ZMA" ) ) continue; // ligand for A2A
					if( !strcasecmp( at[a].resname, "POT" ) ) continue;
					if( !strcasecmp( at[a].resname, "CLA" ) ) continue;
					if( !strcasecmp( at[a].resname, "SOD" ) ) continue;
	
					if( !strcasecmp( at[a].resname, "BGLC") ) continue;
					if( !strcasecmp( at[a].resname, "BGAL") ) continue;
					if( !strncasecmp( at[a].resname, "ANE5", 4) ) continue;

					if( !strcasecmp(at[a].resname, "CHL1" ) )
					{
						if( strcasecmp( at[a].atname, "O3") )
							continue;
					}
					else if( !strncasecmp(at[a].segid, "PRO", 3 ) )
					{
						if( strcasecmp(at[a].atname, "CA" ) )
							continue;			
					}
					else if( strcasecmp( at[a].atname, "C2" ) && strcasecmp( at[a].atname, "C1F" ))
						continue;
				
					if( leaflet[a] == -1 && leaflet_code == LEAFLET_UPPER ) continue;	
					if( leaflet[a] == 1 && leaflet_code == LEAFLET_LOWER ) continue;	

					use_atoms[nuse] = a;
					nuse++;
	
					int isProt = 0;
					if( !strncasecmp( at[a].segid, "PRO", 3) )
						isProt = 1;	
					int gotIt = -1;
	
					// Prot is automatically type TYPE_PROT.
	
					for( int b = 0; b < nres; b++ )
					{
						if( isProt && !strcasecmp(resNames[b],"PROT") )
							gotIt = b;
	
						if( !strcasecmp( at[a].resname, resNames[b] ) )
						{
							gotIt = b;
						}
					}
	
					if( isProt )
						gotIt = TYPE_PROT;
			
					if( gotIt == -1 )
					{
						if( nres == RESMAX )
						{
							printf("Too many different residues for the HMM: was expecting %d or fewer.\n", RESMAX );		
							exit(1);
						}
	
						gotIt = nres;
	
						if( gotIt == TYPE_PROT )
						{
							nres++;
							gotIt = nres;
						}
	
						if( isProt )
							strcpy( resNames[nres], "PROT" );
						else
							strcpy( resNames[nres], at[a].resname );
						nres++;		
					}
	
					type[a] = gotIt;
				}
	
				if( nres <= TYPE_PROT )
				{
					for( int x = nres; x < TYPE_PROT; x++ )
						strcpy(resNames[x], "UNK");
	
					nres = TYPE_PROT+1;
				}
				init_done = 1;
		
		printf("COMMENT");
		for( int r = 0; r < nres; r++ )
		{
			int pa = -1;
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
	
				if( type[a] == r )
				{
					pa = a;
					break;
				}
			}
			if( pa == -1 ) 
			{	//printf(" UNK"); 
			}
			else if( r == TYPE_PROT ) printf(" PROT"); else printf(" %s", at[pa].resname );
		}
		printf("\n");
		fflush(stdout);
			}
		
			if( did_init )
			{
				//for( int a = 0; a < curNAtoms(); a++ )
				for( int ax = 0; ax < nuse; ax++ )
				{
					int a = use_atoms[ax];
	
					double dr[3] = { at[a].x - last_pos[3*a+0], at[a].y - last_pos[3*a+1], at[a].z - last_pos[3*a+2] };
	
					while( dr[0] < -0.5 ) dr[0] += 1;
					while( dr[1] < -0.5) dr[1] += 1;
					while( dr[2] < -0.5 ) dr[2] += 1;
	
					while( dr[0] > 0.5 ) dr[0] -= 1;
					while( dr[1] > 0.5 ) dr[1] -= 1;
					while( dr[2] > 0.5 ) dr[2] -= 1;
	
					at[a].x = last_pos[3*a+0] + dr[0];
					at[a].y = last_pos[3*a+1] + dr[1];
					at[a].z = last_pos[3*a+2] + dr[2];
				}	
			}
	//		printf("%d %lf %lf %lf\n", f, at[0].x, at[0].y, at[0].z );
	
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
				last_pos[3*a+0] = at[a].x;
				last_pos[3*a+1] = at[a].y;
				last_pos[3*a+2] = at[a].z;
			}
	
			did_init = 1;
	
	
			int was_used[nat];
			memset( was_used, 0, sizeof(int) * nat );
		
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
				cur_com[0] += at[a].x;
				cur_com[1] += at[a].y;
			}
	
			cur_com[0] /= nuse;
			cur_com[1] /= nuse;
	
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
				at[a].x -= cur_com[0];
				at[a].y -= cur_com[1];
			}
			
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
	
				while( at[a].x < 0 ) at[a].x += 1.0;
				while( at[a].x > 1.0 ) at[a].x -= 1.0;
				while( at[a].y < 0 ) at[a].y += 1.0;
				while( at[a].y > 1.0 ) at[a].y -= 1.0;
	
				int bx = nx * (at[a].x );
	 			while( bx < 0 ) bx += nx;
	 			while( bx >= nx ) bx -= nx;
				int by = ny * (at[a].y );
	 			while( by < 0 ) by += ny;
	 			while( by >= ny ) by -= ny;
				int bz = nz * (at[a].z );
	 			while( bz < 0 ) bz += nz;
	 			while( bz >= nz ) bz -= nz;
	
				int tbin = bz+nz*(by+ny*bx);
	
				if( bins[tbin].nposSpace == bins[tbin].npos )
				{
					bins[tbin].nposSpace *= 2;
					bins[tbin].pos = (double *)realloc( bins[tbin].pos, sizeof(double) * 3 * bins[tbin].nposSpace );
				}
	
				bins[tbin].pos[3*bins[tbin].npos+0] = at[a].x;
				bins[tbin].pos[3*bins[tbin].npos+1] = at[a].y;
				bins[tbin].pos[3*bins[tbin].npos+2] = at[a].z;
	
				bins[tbin].id[bins[tbin].npos] = a;
				bin_lookup_x[a] = bx;
				bin_lookup_y[a] = by;
	
				bins[tbin].npos += 1;
			}
	
			int a_use = 0;
	
			for( int ax = 0; ax < nuse; ax++ )
			{
				int a = use_atoms[ax];
				if( type[a] == TYPE_PROT) continue;
				int px = bin_lookup_x[a];
				int py = bin_lookup_y[a];
				double pt[2] = { La * at[a].x, Lb * at[a].y };
				int done = 0;
				double dist[nat];
				int n_at_use = 0;
				int sorter[nat];
	
				int grain_expand = 4;
	
				while( !done )
				{
					double L_for_sure;
					double La_for_sure = La * (grain_expand+0.5) / grain;
					double Lb_for_sure = Lb * (grain_expand+0.5) / grain;
	
					if( La_for_sure < Lb_for_sure )
						L_for_sure = La_for_sure;
					else
						L_for_sure = Lb_for_sure;
	
					n_at_use = 0;
	
					for( int dx = -grain_expand; dx <= grain_expand; dx++ )
					for( int dy = -grain_expand; dy <= grain_expand; dy++ )
					{		
						int tbin_x = px + dx;
						int tbin_y = py + dy;
						int tbin_z = 0;
		
						if( tbin_x < 0 ) tbin_x += nx;
						if( tbin_y < 0 ) tbin_y += ny;
		
						if( tbin_x >= nx ) tbin_x -= nx;
						if( tbin_y >= ny ) tbin_y -= ny;
	
						int bin = tbin_x * ny * nz + tbin_y * nz + tbin_z;
						for( int ax = 0; ax < bins[bin].npos; ax++ )
						{ 
	//						if( bins[bin].id[ax] == a ) continue;
	
							double dr[3] = { La * at[bins[bin].id[ax]].x - pt[0], Lb * at[bins[bin].id[ax]].y - pt[1], 0 };	
		
							if( at[bins[bin].id[ax]].z < 0 && at[a].z > 0 ) continue;
							if( at[bins[bin].id[ax]].z > 0 && at[a].z < 0 ) continue;
	
						
	
							while( dr[0] < -La / 2 ) dr[0] += La;
							while( dr[1] < -Lb / 2 ) dr[1] += Lb;
							while( dr[0] > La / 2 ) dr[0] -= La;
							while( dr[1] > Lb / 2 ) dr[1] -= Lb;
		
							double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]);
	
							if( r < L_for_sure )
							{
								sorter[n_at_use] = bins[bin].id[ax];
								dist[bins[bin].id[ax]] = r;
								n_at_use++;	
							}
						}
					}
	
					if( n_at_use >= N_TALLY )
						done = 1;
	
					grain_expand += 1;	
				}
				done = 0;
	
				while( !done )
				{
					done = 1;
	
					for( int x = 0; x < n_at_use-1; x++ )
					{
						if( dist[sorter[x+1]] < dist[sorter[x]] )
						{
							done = 0;
							int t = sorter[x];
							sorter[x] = sorter[x+1];
							sorter[x+1] = t;
						}
					}
				}
	
				int nto[4] = {0,0,0,0};
	
				for( int x = 0; x < N_TALLY; x++ )
				{
					was_used[sorter[x]] += 1;
					nto[type[sorter[x]]] += 1;
				}
				states[a_use+nuse*nframes_use] = code(nto);	
				a_use++;
	
				//if( px == 0 && py == 0 )
	
	
	//			printf("State %c %d contains: %d of POT %d of SOD %d of CLA.\n", 
	//				states[(px*grain+py)*nframes+f], (int)(states[(px*grain+py)*nframes+f]-'A'), nto[0], nto[1], nto[2] ); 
				
	//			printf("%c", states[(px*grain+py)*nframes+f] );
	//			if( f % 80 == 0 && f != 0 ) printf("\n");
	
	//				$iprintf("%c", states[f], nto[0], nto[1], nto[2] );
			}
		
			int viable = 0;
			int used = 0;
			int used_twice = 0;
			for( int a = 0; a < curNAtoms(); a++ )
			{
				viable++;
				if( was_used[a] )
					used++;
				if( was_used[a] > 1 )
					used_twice++;
			}
	
	//		printf("%lf%% were used, %lf%% were used twice.\n",
	//			100.0*(double)used / (double)viable, 100.0 * (double)used_twice/(double)viable );
			for( int x = 0; x < nx*ny*nz; x++ )
				bins[x].npos = 0; 
	
			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
			tf++;
		}
		fclose(dcdFile); 
	}

	
	printf("COMMAND TRAIN\n");
	for( int a = 0; a < nuse; a++ )
	{
		for( int f = 0; f < nframes_use; f++ )
			printf("%c", states[f*nuse+a] );
		printf("\n");
	}
	printf("STOP\n");
	
}






