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


struct selection
{
	char *resName;
	struct selection *next;
};

struct elem
{
	int i, j,k;
	double PXX,PYY,PZZ;
};


struct an_atom
{
	char atname[256];
	double sum_r;
	double sum_r_real;
	double n_r;
	int index;
	int tag;
	struct an_atom *next;
};

/*
const char *sorted_list[] =
{
	"NC3",
	"NH3",
	"PO4",
	"GL1",
	"GL2",
	"C1A",
	"C2A",
	"C3A",
	"D3A",
	"C4A",
	"C5A",
	"C1B",
	"C2B",
	"C3B",
	"D3B",
	"C4B",
	"C5B"
};

*/
const char *sorted_list[] =
{
	"C1F",
	"C2F",
	"C3F",
	"C4F",
	"C5F",
	"C6F",

	"C1S", 
	"C2S",
	"C3S",

	"C15",
	"C14",
	"C13",
	"N",
	"C12",
	"C11",
	"O12",
	"P",
	"O13",
	"O14",
	"O11",
	"C1",
	"C2",
	"C3",
	"O21",
	"O31",
	"C21",
	"C31",
	"O32",
	"O22",
	"C22",
	"C32",
	"C23",
	"C33",
	"C24",
	"C34",
	"C25",
	"C35",
	"C26",
	"C36",
	"C27",
	"C37",
	"C28",
	"C37",
	"C29",
	"C39",
	"C210",
	"C310",
	"C211",
	"C311",
	"C212",
	"C312",
	"C213",
	"C313",
	"C214",
	"C314",
	"C215",
	"C315",
	"C216",
	"C316",
	"C217",
	"C317",
	"C218",
	"C318"
};


/*
const char *sorted_list[] =
{
	"N",
	"HN1",
	"HN2",
	"HN3",
	"C12",
	"C11",
	"O12",
	"P",
	"O13",
	"O14",
	"O11",
	"C1",
	"C2",
	"C3",
	"O21",
	"O31",
	"C21",
	"C31",
	"O32",
	"O22",
	"C22",
	"C32",
	"C23",
	"C33",
	"C24",
	"C34",
	"C25",
	"C35",
	"C26",
	"C36",
	"C27",
	"C37",
	"C28",
	"C38",
	"C29",
	"C39",
	"C210",
	"C310",
	"C211",
	"C311",
	"C212",
	"C312",
	"C213",
	"C313",
	"C214",
	"C314",
	"C215",
	"C315",
	"C216",
	"C316",
	"C217",
	"C317",
	"C218",
	"C318"
};
*/


void MinImage( double *f_target, double *f_pt )
{
	double t_t[3] = { f_target[0], f_target[1], f_target[2] };
	double t_f[3] = { f_pt[0], f_pt[1], f_pt[2] };

	TransformFractional(t_t);

	double rmin = 1e10;
	double mt[3] = {0,0,0};
	for( int dx = -1; dx <= 1; dx ++ )
	for( int dy = -1; dy <= 1; dy ++ )
	{
		double t_f[3] = { f_pt[0] + dx, f_pt[1] + dy, f_pt[2] };
		TransformFractional(t_f);
		double dr[3] = { t_f[0] - t_t[0], t_f[1] - t_t[1], t_f[2] - t_t[2] };
		double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];

		if( r2 < rmin )
		{
			rmin = r2;
			mt[0] = f_pt[0] + dx;
			mt[1] = f_pt[1] + dy;
			mt[2] = f_pt[2];
		}
	}

	f_pt[0] = mt[0];
	f_pt[1] = mt[1];
	f_pt[2] = mt[2];
}

void MinImageR( double *r_target, double *f_pt )
{
	double t_t[3] = { r_target[0], r_target[1], r_target[2] };
	double t_f[3] = { f_pt[0], f_pt[1], f_pt[2] };

	double rmin = 1e10;
	double mt[3] = {0,0,0};
	for( int dx = -1; dx <= 1; dx ++ )
	for( int dy = -1; dy <= 1; dy ++ )
	{
		double t_f[3] = { f_pt[0] + dx, f_pt[1] + dy, f_pt[2] };
		TransformFractional(t_f);
		double dr[3] = { t_f[0] - t_t[0], t_f[1] - t_t[1], t_f[2] - t_t[2] };
		double r2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];

		if( r2 < rmin )
		{
			rmin = r2;
			mt[0] = f_pt[0] + dx;
			mt[1] = f_pt[1] + dy;
			mt[2] = f_pt[2] ;
		}
	}

	f_pt[0] = mt[0];
	f_pt[1] = mt[1];
	f_pt[2] = mt[2];
}


int main( int argc, char **argv )
{
	char buffer[4096];

	if( argc < 5 )
	{
		printf("Syntax: AtomHeight psf refpdb RESIDUE dcd [dcd2 ...]\n");	
		return 0;
	}


	FILE *psfFile = fopen(argv[1], "r" );
	if( ! psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[2] );
		return 0;
	}

	struct an_atom *allAtoms = NULL;

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) )
		loadPSFfromPDB( psfFile );	
	else
		loadPSF( psfFile );
	fclose(psfFile);	
	
	struct atom_rec *ref_pdb = (atom_rec *)malloc( sizeof(atom_rec) * curNAtoms() );
	FILE *refFile = fopen(argv[2],"r");
	loadPDB( refFile, ref_pdb );
	fclose(refFile);
	int isPureWater = 1;
	int init_done = 0;

	int nat = curNAtoms();

	int noxSpace = 100;
	int nox = 0;
	int *oxlist = (int *)malloc( sizeof(int) * noxSpace );
	
	struct an_atom **atom_link = (struct an_atom **)malloc( sizeof(struct an_atom *) * curNAtoms() );

	double *lipid_sign = (double *)malloc( sizeof(double) * curNAtoms() );
	for( int x = 0; x < curNAtoms(); x++ )
		lipid_sign[x] = 1.0;
	
	int *lipid_wrap = (int *)malloc( sizeof(double) * curNAtoms() );
	for( int x = 0; x < curNAtoms(); x++ )
		lipid_wrap[x] = -1;

	double prev_z = 0;

	const char *res_select = argv[3];

	selection *theSelection = NULL;

	int done = 0;
	char *get = argv[3];
	char *tbuf = (char *)malloc( sizeof(char)*(1+strlen(argv[3])));

	while( !done )
	{
		int s = 0;
		while( *get && *get != '+' )
		{
			tbuf[s] = *get;
			tbuf[s+1] = '\0';
			get += 1;
			s++;
		}

		if( *get == '+' ) get += 1;

		if( s > 0 )
		{
			selection *aSel = (selection *)malloc( sizeof(selection) );
			aSel->resName = (char *)malloc( sizeof(char)* (1+strlen(tbuf)) );
			strcpy( aSel->resName, tbuf );
			aSel->next = theSelection;
			theSelection = aSel;
			printf("Selecting: %s\n", aSel->resName );
		}
		else done = 1;
	
	}

	double frames_tot = 0;

	double dZ = 0.5;
	double Z_MIN = -50.0;
	double Z_MAX = 50.0;
	int nbins_hist = (Z_MAX-Z_MIN)/dZ;
	double *histogram[5] ={NULL,NULL,NULL,NULL,NULL};
	for( int t = 0; t < 5; t++ )
	{
		histogram[t] = (double *)malloc( sizeof(double) * nbins_hist );
		memset( histogram[t], 0, sizeof(double) * nbins_hist );	
	}

	double av_Lx = 0;
	double  av_Ly = 0;

	for( int c = 4; c < argc; c++ )
	{
		FILE *dcdFile = fopen(argv[c], "r");
	
		if( ! dcdFile )
		{
			printf("Couldn't open dcd file '%s'.\n", argv[c] );
			return 0;
		}
	
		int isPDB = 0;

		if( !strcasecmp( argv[c]+strlen(argv[c])-3,"pdb") )
			isPDB=1;

		int nframes = 1;

		if( !isPDB )
		{
			readDCDHeader(dcdFile);
			setSymmetric();
			nframes = curNFrames();
		}
		struct atom_rec *at = (struct atom_rec *)malloc( sizeof(struct atom_rec) * curNAtoms() );

		for( int f = 0; f < nframes; f++, frames_tot += 1)
		{
			if( isPDB )
				loadPDB( dcdFile, at );
			else
				loadFrame( dcdFile, at );
	
			double La,Lb,Lc,alpha,beta,gamma;
			PBCD(&La,&Lb,&Lc,&alpha,&beta,&gamma);
	
			av_Lx += La;
			av_Ly += Lb;

			if( !isPDB && !DCDsuccess() )
			{
				nframes = f;
				break;
			}
	
			int pres = -1;
			int cres = 0;

			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( at[a].res != pres )
				{
					pres = at[a].res;
					cres++;
				}

				at[a].res = cres;
			}

			if( !init_done )
			{
				for( int a = 0; a < curNAtoms(); a++ )
				{
					int gotit=0;
					for( selection *aSel = theSelection; aSel; aSel = aSel->next )
					{
						if( !strcasecmp( aSel->resName, "PEPC") && (!strcasecmp(at[a].resname, "DOPE") || !strcasecmp(at[a].resname, "DOPC") ) )
							gotit = 1;
						if( !strcasecmp( at[a].resname, aSel->resName) ) 
							gotit=1;
					} 	
					if( !gotit ) continue;
//					if( strcasecmp( at[a].resname, res_select) ) continue;

					atom_link[a] = NULL;	
	
					char *tat = at[a].atname;
					char *trt = at[a].resname;		
			
					while( *tat == ' ' ) tat += 1;
					while( *trt == ' ' ) trt += 1;
	
					int np_mol = at[a].res-1;
					
					if( (!strcasecmp( trt, "TIP3" ) && tat[0] == 'O') || (tat[0] == 'W'))
					{
						if( nox == noxSpace )
						{
							noxSpace *= 2;
							oxlist = (int *)realloc( oxlist, sizeof(int) * noxSpace );
						}
						
						oxlist[nox] = a;
						nox++;
					}
					else if( strcasecmp( trt, "TIP3" ) )	
						isPureWater = 0;
	
					if( !strcasecmp( trt, "TIP3" ) || !strcasecmp( trt, "TETD") || !strcasecmp(trt, "W") )
						continue;
/*
					if( strcasecmp( trt, "SDPC" ) &&
					    strcasecmp( trt, "SDPE" ) &&
					    strcasecmp( trt, "DAPC" ) && 
					    strcasecmp( trt, "DOPC" ) && 
					    strcasecmp( trt, "POPC" ) && 
					    strcasecmp( trt, "POPS" ) && 
					    strcasecmp( trt, "DOPE" ) && 
					    strcasecmp( trt, "DOPS" ) && 
					    strncasecmp( trt, "CER1", 4 ) && 
					    strncasecmp( trt, "PHPC", 4 ) && 
					    strncasecmp( trt, "EHPC", 4 ) && 
					    strncasecmp( trt, "PSM", 3 ) && 
					    strcasecmp( trt, "DPPC" ) && 
					    strcasecmp( trt, "DPPE" ) && 
					    strcasecmp( trt, "SOPC" ) 
							) continue;
*/
					struct an_atom *the_link = NULL;

					for( struct an_atom *link = allAtoms; link; link = link->next )
					{
						if( !strcasecmp( link->atname, tat ) )
							{ the_link = link; break; }
					}
			
					if( !the_link )
					{
						the_link = (struct an_atom *)malloc( sizeof(struct an_atom) );
						strcpy( the_link->atname, tat );
						the_link->sum_r = 0;
						the_link->sum_r_real = 0;
						the_link->n_r = 0;
						the_link->next = allAtoms;
						allAtoms = the_link;
					}
					
					atom_link[a] = the_link;
				}
			
				for( int a = 0; a < nat; a++ )
				{
					if( !atom_link[a] ) continue;

					if( lipid_wrap[at[a].res-1] == -1 && (!strcasecmp( at[a].atname, "C216" )||!strcasecmp(at[a].atname, "C16S")) ) lipid_wrap[at[a].res-1] = a;
				}	


			}

#define N_BINS_MOLDIST 100
	                 double best_chi2 = 1e10;
			double wrapto = 0;
	                 int nbins = N_BINS_MOLDIST;
	                 double moldist[N_BINS_MOLDIST];
	                 memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );
	
//			for( int p = nuse; p < nuse+nuset; p++ )
			for( int a = 0; a < nat; a++ )
			{
		//		if( !atom_link[a] ) continue;
	//			if( at[a].atname[0] != 'C' ) continue;
				if( strcasecmp(at[a].atname, "C14" ) && strcasecmp( at[a].atname, "C14S") && strcasecmp( at[a].atname, "C214") && strcasecmp( at[a].atname, "C314") && strcasecmp( at[a].atname, "C4B") && strcasecmp( at[a].atname, "D4B") ) continue;

	                              double tz = at[a].z/Lc;

	                              while( tz < 0 ) tz += 1;
	                              while( tz >= 1 ) tz -= 1;
	
	                              int zb = N_BINS_MOLDIST * tz; // this is right
	                              if( zb < 0 ) zb = 0;
	                              if( zb >= N_BINS_MOLDIST ) zb = N_BINS_MOLDIST-1;
	                              moldist[zb] += 1;
			}

                 for( int zb = 0; zb < nbins; zb++ )
                 {
                         double zv = (zb+0.5) / (double)N_BINS_MOLDIST;
//			printf("%lf %lf\n", zv, moldist[zb] );
                          int zlow  = zb - nbins/2;
                          int zhigh = zlow + nbins;

                          double lchi2 = 0;
                          for( int iz = zlow; iz < zhigh; iz++ )
                          {
                                  double dz =  (iz+0.5) / nbins - zv;

                                  int iiz = iz;
                                  while( iiz < 0 ) iiz += nbins;
                                  while( iiz >= nbins ) iiz -= nbins;

                                  lchi2 += moldist[iiz] * (dz) * (dz);
                          }

                          if( lchi2 < best_chi2 )
                          {
                                  best_chi2 = lchi2;
                                  wrapto = zv;
                          }
                 }


		while( wrapto < -0.5 ) wrapto += 1;
		while( wrapto >  0.5 ) wrapto -= 1;
	
		wrapto *= Lc;
	
		printf("Wrapto: %lf\n", wrapto );	
	
		for( int x = 0; x < nat; x++ )
		{
			while( at[x].z -wrapto < -Lc/2 ) at[x].z += Lc;
			while( at[x].z -wrapto >  Lc/2 ) at[x].z -= Lc;

			at[x].z -= wrapto;
		}

			double avz = 0;
			double navz = 0;
	
			for( int a = 0; a < nat; a++ )
			{
				if( strcasecmp(at[a].atname, "C14" ) && strcasecmp( at[a].atname, "C14S") && strcasecmp( at[a].atname, "C214") && strcasecmp( at[a].atname, "C314") && strcasecmp( at[a].atname, "C4B") && strcasecmp( at[a].atname, "D4B") ) continue;

//				if( atom_link[a] )
				{
//					while( at[a].z - prev_z < -Lc/2 ) at[a].z += Lc;
//					while( at[a].z - prev_z >  Lc/2 ) at[a].z -= Lc;
					avz += at[a].z;	
					navz += 1;
				}
			}

			avz /= navz;
			prev_z = avz;

//			printf("wrapto: %le avz: %le\n", wrapto, avz );

/*
			for( int a = 0; a < nat; a++ )
			{
				if( atom_link[a] )
				{
					avz += at[a].z;	
					navz += 1;
				}
			}

			avz /= navz;
*/

			if( init_done == 0 )
			{
				int nup = 0, ndown = 0;
				for( int a = 0; a < nat; a++ )
				{
					if( atom_link[a]  )
					{
#ifdef OLD
					if( !strncasecmp( at[a].segid, "TPE", 3 ) ) 
						lipid_sign[at[a].res-1] = 1;
					else if( !strncasecmp( at[a].segid, "BPE", 3 ) )
						lipid_sign[at[a].res-1] = -1;
					else if( at[a].atname[0] == 'O' && at[a].z > avz )
						lipid_sign[at[a].res-1] = 1;
					else if( at[a].atname[0] == 'O' )
						lipid_sign[at[a].res-1] = -1;
					else if( at[a].atname[0] == 'P' && at[a].z > avz )
					{	lipid_sign[at[a].res-1] = 1; nup++; 
//printf("at up %s z: %lf\n", at[a].atname, at[a].z); 
					}
					else if( at[a].atname[0] == 'P' )
					{	lipid_sign[at[a].res-1] = -1; ndown++; 
//printf("at down %s z: %lf\n", at[a].atname, at[a].z); 
					}
					else  if( (!strcasecmp(at[a].atname,"NH3") || !strcasecmp(at[a].atname,"NC3")) && at[a].z > avz  )
						lipid_sign[at[a].res-1] = 1;
					else if( (!strcasecmp(at[a].atname,"NH3") || !strcasecmp(at[a].atname,"NC3")) && at[a].z < avz  )
						lipid_sign[at[a].res-1] = -1;
#else
					if( !strncasecmp( at[a].segid, "TPE", 3 ) ) 
						lipid_sign[at[a].res-1] = 1;
					else if( !strncasecmp( at[a].segid, "BPE", 3 ) )
						lipid_sign[at[a].res-1] = -1;
					else if( at[a].atname[0] == 'O' && ref_pdb[a].z > 0 )
						lipid_sign[at[a].res-1] = 1;
					else if( at[a].atname[0] == 'O' )
						lipid_sign[at[a].res-1] = -1;
					else if( !strcasecmp( at[a].resname, "CPC") && at[a].atname[0] == 'N' && ref_pdb[a].z > 0 )
					{	lipid_sign[at[a].res-1] = 1; }
					else if( !strcasecmp( at[a].resname, "CPC") && at[a].atname[0] == 'N' && ref_pdb[a].z <= 0 )
					{	lipid_sign[at[a].res-1] = -1; }
					else if( at[a].atname[0] == 'P' && ref_pdb[a].z > 0 )
					{	lipid_sign[at[a].res-1] = 1; nup++; 
						//printf("at up %s z: %lf\n", at[a].atname, at[a].z); 
					}
					else if( at[a].atname[0] == 'P' )
					{	lipid_sign[at[a].res-1] = -1; ndown++; 
//						printf("at down %s z: %lf\n", at[a].atname, at[a].z); 
					}
					else  if( (!strcasecmp(at[a].atname,"NH3") || !strcasecmp(at[a].atname,"NC3")) && ref_pdb[a].z > 0  )
						lipid_sign[at[a].res-1] = 1;
					else if( (!strcasecmp(at[a].atname,"NH3") || !strcasecmp(at[a].atname,"NC3")) && ref_pdb[a].z < 0  )
						lipid_sign[at[a].res-1] = -1;
					}
#endif
				}
//				printf("nup: %d ndown: %d\n", nup, ndown );
			}
				
			init_done = 1;

			double av_n = 0;

			for( int a = 0; a < curNAtoms(); a++ )
			{
				if( atom_link[a] )
				{
					double com_use[3] = {0,0,0};

					atom_link[a]->sum_r += lipid_sign[at[a].res-1] * (at[a].z - avz);
					atom_link[a]->n_r += 1;
				
					int zbin = nbins_hist  * (at[a].z-avz-Z_MIN)/(Z_MAX-Z_MIN);
					if( zbin < 0 ) zbin = 0;
					if( zbin >= nbins_hist ) zbin = nbins_hist-1;
					if( lipid_sign[at[a].res-1] < 0 )
					{	
						if( at[a].atname[0] != 'C' && at[a].atname[0] != 'D' ) 
							histogram[2][zbin] += 1;
						else if( at[a].atname[0] == 'C' || at[a].atname[0] == 'D' )
							histogram[3][zbin] += 1;
					}
					else
					{	
						if( at[a].atname[0] != 'C' && at[a].atname[0] != 'D' ) 
							histogram[0][zbin] += 1;
						else if( at[a].atname[0] == 'C' || at[a].atname[0] == 'D' )
							histogram[1][zbin] += 1;
					}
				}
				else
				{	
					int zbin = nbins_hist  * (at[a].z-avz-Z_MIN)/(Z_MAX-Z_MIN);
					if( zbin < 0 ) zbin = 0;
					if( zbin >= nbins_hist ) zbin = nbins_hist-1;
							
					histogram[4][zbin] += 1;
				}				
			}
	
//			printf("%d %lf\n", f, av_n ); 
		

			for( int a = 0; a < curNAtoms(); a++ )
				at[a].zap();
		} 

		free(at);
		fclose(dcdFile);
	}

	int nentry = 0;
	for( struct an_atom *at = allAtoms; at; at = at->next )
		nentry++;
	struct an_atom **alist = (struct an_atom **)malloc( sizeof(struct an_atom *) * nentry );
	nentry = 0;
	for( struct an_atom *at = allAtoms; at; at = at->next )
		{ alist[nentry] = at; nentry++; }

	int sorter[nentry];
	for( int i = 0; i < nentry; i++ )
		sorter[i] = i;
	done = 0;

	while( !done )
	{
		done = 1;

		for( int i = 0; i < nentry-1; i++ )
		{
			if( alist[sorter[i]]->sum_r < alist[sorter[i+1]]->sum_r )
			{
				int t = sorter[i];
			
				sorter[i] = sorter[i+1];
				sorter[i+1] = t;

				done = 0;
			}
		}
	}

//	for( int i = 0; i < nentry; i++ )
//	{
//		printf("%s %lf %lf\n", alist[sorter[i]]->atname, alist[sorter[i]]->sum_r / alist[sorter[i]]->n_r, alist[sorter[i]]->n_r );
//	}

//	if( 1 ||( !strncasecmp( res_select, "CER1", 4 ) || !strncasecmp( res_select, "PSM", 3 ) ) )
	{
		for( int a = 0; a < nentry; a++ )
		{
			int ta = sorter[a];
			if( alist[ta]->atname[0] == 'H' ) continue;
				printf("%s %lf %lf\n", alist[ta]->atname, alist[ta]->sum_r / alist[ta]->n_r, alist[ta]->sum_r_real / alist[ta]->n_r );
		}
	}
/*	else
	{ 
		for( int x = 0; x < sizeof(sorted_list)/sizeof(char*); x++ )
		{
			int ta = -1;
	
			for( int a = 0; a < nentry; a++ )
			{
				if( !strcasecmp( alist[a]->atname, sorted_list[x] ) )
					ta = a;
			}
	
			if( ta == -1 )
			{
				continue;
				printf("Error. Couldn't find atom '%s'.\n", sorted_list[x] );
				return 0;
			}
	
			printf("%s %lf %lf\n", alist[ta]->atname, alist[ta]->sum_r / alist[ta]->n_r, alist[ta]->sum_r_real / alist[ta]->n_r );
		}
	}*/

//	for( struct an_atom *at = allAtoms; at; at = at->next )
//	{
//		printf("%s %lf %lf\n", at->atname, at->sum_r / at->n_r, at->n_r );
//	}

	av_Lx /= frames_tot;
	av_Ly /= frames_tot;
#if 0	
	for( int b = 0; b < nbins_hist; b++ )
	{
		double z1 = Z_MIN +     b*dZ;
		double z2 = Z_MIN + (b+1)*dZ;
		double vol = av_Lx * av_Ly * (z2-z1);
		double zav = (z1+z2)/2;
		printf("%lf %lf %lf %lf %lf %lf\n", zav, 
			histogram[0][b] / (frames_tot*vol), 
			histogram[1][b] / (frames_tot*vol), 
			histogram[2][b] / (frames_tot*vol), 
			histogram[3][b] / (frames_tot*vol),  
			histogram[4][b] / (frames_tot*vol)  
			);
	}

#endif
}






