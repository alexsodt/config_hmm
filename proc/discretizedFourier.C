//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "pdb.h"
#include "dcd.h"

#define BILAYER_MODE
//#define COMPUTE_D
//#define DENSITY_CURVATURE_CORR

//#define RANDOMIZE
#define FIXED_SEED
#define FORCE_PASS 	0
#define TOLERANT_PASS 	1
#define RESTRICTED_PASS 2
#define RANDOM_PASS	3
//#define SKIP_OPT

//#define APPLY_EXTERNAL


int main( int argc, char **argv )
{
	double bin_width = 3.0;

#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];
	
	if( argc < 6 )
	{
		printf("Syntax: fourier.exe psf pdb_sets_leaflet dcd SA1[+SA2][+...] TA1[+TA2][+...] [name]\n");
		printf("Put a lower case 'x' before the atom name to stop it from being included in the Fourier transform.\n");
		return -1;
	}

	       
	FILE *psfFile = fopen(argv[1],"r");
	
	if( !psfFile )
	{
		printf("Couldn't open PSF file '%s'.\n", argv[1] );
		exit(1);
	}

	if( !strcasecmp( argv[1] + strlen(argv[1])-3, "pdb" ) ) 
		loadPSFfromPDB( psfFile );    
        else
		loadPSF( psfFile );

	int nat = curNAtoms();
	
	int *atom_list = (int *)malloc( sizeof(int) * curNAtoms() );
	int *x_atom_list = (int *)malloc( sizeof(int) * curNAtoms() );
	memset( x_atom_list, 0, sizeof(int) * curNAtoms() );
	int *tatom_list = (int *)malloc( sizeof(int) * curNAtoms() );
	int nuse = 0;
	int nuset = 0;
	FILE *dcdFile = fopen( argv[3], "r");
	
	if( !dcdFile )
	{
		printf("Couldn't open dcd file '%s'.\n", argv[3] );
		exit(1);
	}

        readDCDHeader(dcdFile);
        setFractional();

	int nframes = curNFrames();

	struct atom_rec *at = (struct atom_rec*)malloc( sizeof(struct atom_rec) * nat );
	struct atom_rec *at_leaflet = (struct atom_rec*)malloc( sizeof(struct atom_rec) * nat );
	int *assign_leaflet = (int *)malloc( sizeof(int) * nat );
	memset( assign_leaflet, 0, sizeof(int) * nat );
	FILE *pdbFile = fopen(argv[2],"r");

	loadPDB( pdbFile, at_leaflet ); 

	int nox=0,noy=0,noz=0;
	int doQII = 0;
	int doHex = 0;

	loadFrame( dcdFile, at );
	if( !DCDsuccess() )
	{
		printf("Could not load intended frames\n");
		exit(1);
	}
	double last_PBC[3];
	double cur_PBC[3];
	double PBC_vecs[3][3];
	double La,Lb,Lc,alpha,beta,gamma;	
	PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );

	printf("La: %lf Lb: %lf Lc: %lf alpha: %lf beta: %lf gamma: %lf\n",
		La, Lb, Lc, alpha, beta, gamma );

	cur_PBC[0] = La;
	cur_PBC[1] = Lb;
	cur_PBC[2] = Lc;
	
	double MAX_SPACING = 15.0;

	int N_QX = ceil(La / MAX_SPACING);
	int N_QY = ceil(Lb / MAX_SPACING);

	int dqx = N_QX/2;
	int dqy = N_QY/2;

	N_QX = 2*dqx+1;
	N_QY = 2*dqy+1;

	double *hq = (double *)malloc( sizeof(double) * N_QX * N_QY * 2 * 2 );
	memset( hq, 0, sizeof(double) * N_QX * N_QY * 2 * 2 );
	

	int nx = 1 + La / bin_width;
	int ny = 1 + Lb / bin_width;
	int nz = 1 + Lc / bin_width;
	memcpy( last_PBC, cur_PBC, sizeof(double) * 3 );	
	
	// load starting lattice.
	
	int max_s=64;
	char *surfaceAtoms[max_s];
	char *tailAtoms[max_s];
	char *res_surfaceAtoms[max_s];
	char *res_tailAtoms[max_s];


	int n_atoms_check = decodeString(argv[4],surfaceAtoms,res_surfaceAtoms,max_s);
	int n_tatoms_check = decodeString(argv[5],tailAtoms,res_tailAtoms,max_s);

	int *leaflet = (int *)malloc( sizeof(int) * nat );
	
	char fileName[256];
	
	if( argc > 6 )
		sprintf(fileName, "%s_ldata.txt", argv[6] );
	else
		sprintf(fileName, "ldata.txt");
	
	FILE *ldata = (FILE *)fopen(fileName, "w");

	fprintf(ldata, "# %s", argv[3] );

	// get leaflet from P or O atoms, not acyl chains.

	int pres = at_leaflet[0].res;
	int res_atom_start = 0;

	// center of mass from input pdb file. need not be zero centered.
	double leaflet_av_z = 0;
	double leaflet_n = 0;
	for( int ta = 0; ta < nat; ta++ )
	{
		for( int xa = 0; xa < n_atoms_check; xa++ )
		{
			if( !strcasecmp( at[ta].atname, surfaceAtoms[xa] ) && (!strcasecmp( "any", res_surfaceAtoms[xa]) || !strcasecmp( at[ta].resname, res_surfaceAtoms[xa])) )
			{
				leaflet_av_z += at_leaflet[ta].z;
				leaflet_n += 1;
			}
		}
	}
	leaflet_av_z /= leaflet_n;

	for( int a = 0; a <= nat; a++ )
	{
		if( a == nat || at[a].res != pres )
		{
			int cur_leaflet = 0, alt_cur_leaflet = 0;

			for( int ta = res_atom_start; ta < a; ta++ )
			{
				// get the leaflet from the phosphorous atom if it's there.
				if( at_leaflet[ta].atname[0] == 'P' ) 
				{
					if( at_leaflet[ta].z < leaflet_av_z )
						cur_leaflet = -1;
					else
						cur_leaflet = 1;
				}	
				// alternatively, get the leaflet from the surfaceAtom since we know it's there.
				for( int xa = 0; xa < n_atoms_check; xa++ )
				{
					if( !strcasecmp( at[ta].atname, surfaceAtoms[xa] ) && (!strcasecmp( "any", res_surfaceAtoms[xa]) || !strcasecmp( at[ta].resname, res_surfaceAtoms[xa])) )
					{
						if( at_leaflet[ta].z < leaflet_av_z ) 
							alt_cur_leaflet = -1;
						else
							alt_cur_leaflet = 1;
					}
				}
			}	


			if( cur_leaflet == 0 ) cur_leaflet = alt_cur_leaflet;
			
			for( int ta = res_atom_start; ta < a; ta++ )
				assign_leaflet[ta] = cur_leaflet;

			if( a < nat ) res_atom_start = a;
		}

		if( a == nat ) break;

		pres = at[a].res;
	}

	for( int a = 0; a < nat; a++ )
	{
		for( int xa = 0; xa < n_atoms_check; xa++ )
		{
			char *puse = surfaceAtoms[xa];
			if( *puse == 'x' )
				puse = puse + 1;
			int use = 0;

			if( !strncasecmp( at[a].segid, "GLPA", 4 ) && strncasecmp( at[a].resname, "CER", 3) ) 
				continue;

			if( !strcasecmp( at[a].atname, puse ) && (!strcasecmp( "any", res_surfaceAtoms[xa]) || !strcasecmp( at[a].resname, res_surfaceAtoms[xa])))
			{
				if( assign_leaflet[a] < 0 ) 
					leaflet[nuse] = -1;
				else
					leaflet[nuse] = 1;

				fprintf(ldata, " %s,%s,%d,%d,%s", at[a].resname, at[a].atname,leaflet[nuse], at[a].res, at[a].segid );

				if( surfaceAtoms[xa][0] == 'x' )
					x_atom_list[nuse] = 1;

				atom_list[nuse] = a;
				nuse++;
				break;
			}
		}
	}
	fprintf(ldata, "\n");
	
	for( int xa = 0; xa < n_tatoms_check; xa++ )
	{
		for( int a = 0; a < nat; a++ )
		{
			int use = 0;
			
			if( !strncasecmp( at[a].segid, "GLPA", 4 ) && strncasecmp( at[a].resname, "CER", 3) ) 
				continue;

			if( !strcasecmp( at[a].atname, tailAtoms[xa] ) && (!strcasecmp( "any", res_tailAtoms[xa]) || !strcasecmp( at[a].resname, res_tailAtoms[xa])) )
			{
				tatom_list[nuset] = a;
				nuset++;
			}
		}
	}
	
	double fr[3]={0,0,0};

	double *last_align = (double *)malloc( sizeof(double) * (nuse+nuset) *3 );
	double *cur_set = (double *)malloc( sizeof(double) * (nuse+nuset) *3 );
	int init_done = 0;
	

	if( argc > 6 )
		sprintf(fileName, "%s_hq.txt", argv[6] );
	else
		sprintf(fileName, "hq.txt");

	FILE *fq    = (FILE *)fopen(fileName, "w");

#ifdef BILAYER_MODE
	fprintf(fq,"mode=bilayer\n");
#endif

	fprintf(fq, "# %d %d %s", N_QX, N_QY, argv[3] );

	for( int iqx = 0; iqx < N_QX; iqx++ )
	for( int iqy = 0; iqy < N_QY; iqy++ )
	{
		double qx = 2 * M_PI * (-dqx+iqx) / La;
		double qy = 2 * M_PI * (-dqy+iqy) / Lb;
	
		fprintf(fq, " %lf,%lf", qx, qy );
	}	
	fprintf(fq, "\n");

	double *hsum = (double *)malloc( sizeof(double) * N_QX * N_QY*3 );
	double *nsum = (double *)malloc( sizeof(double) * N_QX * N_QY*3 );

	double pwrapto = -1e10;
	int pwrapto_set = 0;

	for( int f = 0; f < nframes; f++ )
	{
		memset( hq, 0, sizeof(double) * 2 * N_QX * N_QY );

		memset( hsum, 0, sizeof(double) * N_QX*N_QY*3);
		memset( nsum, 0, sizeof(double) * N_QX*N_QY*3);
		
		int pres = at[0].res;
		int wrap_at = 0;
		for( int a = 0; a < nat; a++ )
		{
			if( at[a].res == pres )
			{
				while( at[a].x - at[wrap_at].x < -0.5 ) at[a].x += 1.0;
				while( at[a].y - at[wrap_at].y < -0.5 ) at[a].y += 1.0;
				while( at[a].z - at[wrap_at].z < -0.5 ) at[a].z += 1.0;
				while( at[a].x - at[wrap_at].x >  0.5 ) at[a].x -= 1.0;
				while( at[a].y - at[wrap_at].y >  0.5 ) at[a].y -= 1.0;
				while( at[a].z - at[wrap_at].z > 0.5 ) at[a].z -= 1.0;
			}
			else
			{
				pres = at[a].res;	
				wrap_at = a;
			}
		}

		for( int ax = 0; ax < nuse; ax++ )
		{
			cur_set[3*ax+0] = at[atom_list[ax]].x;	
			cur_set[3*ax+1] = at[atom_list[ax]].y;	
			cur_set[3*ax+2] = at[atom_list[ax]].z;	
		}
		
		for( int ax = nuse; ax < nuse+nuset; ax++ )
		{
			cur_set[3*ax+0] = at[tatom_list[ax-nuse]].x;	
			cur_set[3*ax+1] = at[tatom_list[ax-nuse]].y;	
			cur_set[3*ax+2] = at[tatom_list[ax-nuse]].z;	
		}


		if( init_done )
		{
			for( int ax = 0; ax < nuse+nuset; ax++ )
			{
				while( cur_set[3*ax+0] - last_align[3*ax+0] < -0.5 ) cur_set[3*ax+0] += 1.0;
				while( cur_set[3*ax+0] - last_align[3*ax+0] >  0.5 ) cur_set[3*ax+0] -= 1.0;
					
				while( cur_set[3*ax+1] - last_align[3*ax+1] < -0.5 ) cur_set[3*ax+1] += 1.0;
				while( cur_set[3*ax+1] - last_align[3*ax+1] >  0.5 ) cur_set[3*ax+1] -= 1.0;
				
				while( cur_set[3*ax+2] - last_align[3*ax+2] < -0.5 ) cur_set[3*ax+2] += 1.0;
				while( cur_set[3*ax+2] - last_align[3*ax+2] >  0.5 ) cur_set[3*ax+2] -= 1.0;
			}		
		}
	
		memcpy( last_align, cur_set, sizeof(double) * 3 * (nuse+nuset) );

#define N_BINS_MOLDIST 100
                 double best_chi2 = 1e10;
		double wrapto = 0;
                 int nbins = N_BINS_MOLDIST;
                 double moldist[N_BINS_MOLDIST];
                 memset( moldist, 0, sizeof(double) * N_BINS_MOLDIST );

		for( int p = nuse; p < nuse+nuset; p++ )
		{
                              double tz = cur_set[3*p+2];

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

		double dwrap = wrapto - pwrapto;
		while( dwrap < -0.5 ) dwrap += 1.0;
		while( dwrap > 0.5 ) dwrap -= 1.0;

		if( pwrapto_set && fabs(dwrap) > 0.2  )
		{
			printf("wrapto: %le pwrapto: %le dwrap: %le\n", wrapto, pwrapto, dwrap );
			printf("Too big a step in wrapto. resetting to %le\n", pwrapto );
			wrapto = pwrapto;
		}

		pwrapto = wrapto;
		pwrapto_set = 1;

	
		for( int x = 0; x < nuse; x++ )
		{
			while( cur_set[3*x+2] -wrapto < -0.5 ) cur_set[3*x+2] += 1.0;
			while( cur_set[3*x+2] -wrapto >  0.5 ) cur_set[3*x+2] -= 1.0;
		}

		double com[3] = { 0,0,0};

		for( int x = 0; x < nuse; x++ )
		{
			com[0] += cur_set[3*x+0]; 
			com[1] += cur_set[3*x+1]; 
			com[2] += cur_set[3*x+2]; 
		}

		com[0] /= nuse;
		com[1] /= nuse;
		com[2] /= nuse;

		init_done = 1;

		double leaflet_n[2] = { 0,0};

		double mean[2] = { 0,0};
		
		for( int ax = 0; ax < nuse; ax++ )
		{
			if( x_atom_list[ax] ) continue;

#ifdef DISABLE_SPHINGO
			if( !strncasecmp(at[atom_list[ax]].segid, "GLPA",4) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "CER",3) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "PSM",3) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "SSM",3) )
				continue;
#endif

			double x = cur_set[3*ax+0] * La;
			double y = cur_set[3*ax+1] * Lb;
#ifdef RANDOMIZE
			double tz = com[2]*Lc + 3.0 * 2*(-0.5+rand()/(double)RAND_MAX);
			cur_set[3*ax+2] = tz/Lc+com[2];
#endif			
			double z = (cur_set[3*ax+2]-com[2]) * Lc;
			if( leaflet[ax] > 0 )
				mean[0] += z;
			else
				mean[1] += z;
			
			if( leaflet[ax] > 0 ) 
				leaflet_n[0] += 1;
			else
				leaflet_n[1] += 1;
		}
	
		mean[0] /= leaflet_n[0];
		mean[1] /= leaflet_n[1];
		
		for( int ax = 0; ax < nuse; ax++ )
		{
			double x = cur_set[3*ax+0];
			double y = cur_set[3*ax+1];
			double z = (cur_set[3*ax+2]-com[2]) * Lc;

			if( x_atom_list[ax] ) continue;

#ifdef DISABLE_SPHINGO
			if( !strncasecmp(at[atom_list[ax]].segid, "GLPA",4) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "CER",3) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "PSM",3) )
				continue;
			if( !strncasecmp(at[atom_list[ax]].resname, "SSM",3) )
				continue;
#endif

			while( x < 0 ) x += 1.0;
			while( x >= 1.0 ) x -= 1.0; 
			
			while( y < 0 ) y += 1.0;
			while( y >= 1.0 ) y -= 1.0; 

			int ix = x * N_QX;
			int iy = y * N_QY;

			while( ix < 0 ) ix += N_QX;
			while( ix >= N_QX ) ix -= N_QX;
			
			while( iy < 0 ) iy += N_QY;
			while( iy >= N_QY ) iy -= N_QY;

			hsum[(ix*N_QY+iy)*3+0] += z;

			if( leaflet[ax] > 0 )
				hsum[(ix*N_QY+iy)*3+1] += z;
			else
				hsum[(ix*N_QY+iy)*3+2] += z;

			nsum[(ix*N_QY+iy)*3+0] += 1;

			if( leaflet[ax] > 0 )
				nsum[(ix*N_QY+iy)*3+1] += 1;
			else
				nsum[(ix*N_QY+iy)*3+2] += 1;
		}


		for( int ax = 0; ax < nuse; ax++ )
		{
			double x = cur_set[3*ax+0] * La;
			double y = cur_set[3*ax+1] * Lb;
			double z = (cur_set[3*ax+2]-com[2]) * Lc;
			fprintf(ldata, " %lf,%lf,%lf", x,y,
#ifdef BILAYER_MODE
z
#else
(leaflet[ax] > 0 ? z-mean[0] : z - mean[1] )
#endif
);
		}

		fprintf(ldata,"\n");
			

		for( int iqx = 0; iqx < N_QX; iqx++ )
		for( int iqy = 0; iqy < N_QY; iqy++ )
		{
			double qx = 2 * M_PI * (-dqx+iqx) / La;
			double qy = 2 * M_PI * (-dqy+iqy) / Lb;

			for( int ix = 0; ix < N_QX; ix++ )
			for( int iy = 0; iy < N_QY; iy++ )
			{
				double x = La * (ix+0.5) /N_QX;
				double y = Lb * (iy+0.5) /N_QY;
				double z[3] = {0,0,0};
				
				for( int type = 0; type < 3; type++ )
				{
					if( nsum[(ix*N_QY+iy)*3+type] == 0 )
					{
						if( iqx == 0 && iqy == 0 )
						{
//							printf("Interpolating.\n");
						}
						double w = 0;
						
						for( int dx = -1; dx <= 1; dx++ )
						for( int dy = -1; dy <= 1; dy++ )
						{
							int bx = ix + dx;
							int by = iy + dy;
							if( bx < 0 ) bx += N_QX; if( bx >= N_QX ) bx -= N_QX;	
							if( by < 0 ) by += N_QY; if( by >= N_QY ) by -= N_QY;	

//							if( nsum[(bx*N_QY+by)*3+type] > 0 && iqx == 0 && iqy == 0 ) 
//								printf("using %le (%lf)\n", hsum[(bx*N_QY+by)*3+type]/nsum[(bx*N_QY+by)*3+type], nsum[(bx*N_QY+by)*3+type] );

#ifdef INTERP_METHOD_1
							z[type] += hsum[(bx*N_QY+by)*3+type];
							w += nsum[(bx*N_QY+by)*3+type];
#else
							if( nsum[(bx*N_QY+by)*3+type] > 0 )
							{
								z[type] += hsum[(bx*N_QY+by)*3+type]/nsum[(bx*N_QY+by)*3+type];
								w += 1;
							}
#endif
						}			

						if( w < 1 )
						{
							printf("The surface was poorly covered at a spot. This has always indicated that particular lipids were left out. Stopping.\n");
							exit(1);
						}

						z[type] /= w;
	
//						if( iqx == 0 && iqy == 0 ) 
//						printf("interpolating to %le\n", (type == 1 ? 1 : -1 ) * z[type] );	
					}
					else
						z[type] =  (hsum[(ix*N_QY+iy)*3+type]/nsum[(ix*N_QY+iy)*3+type]);
				}

		//		if( iqx == 2 && iqy == 0 )
		//		printf("%lf %lf %lf %lf\n", x, y, z[1], z[2] );
#ifdef COMPUTE_D
				hq[(iqx*N_QY+iqy)*2+0] +=  (z[1]-z[2]) * cos( qx * x + qy * y ); 
				hq[(iqx*N_QY+iqy)*2+1] -=  (z[1]-z[2]) * sin( qx * x + qy * y ); 
#else
				hq[(iqx*N_QY+iqy)*2+0] +=  (z[1]+z[2])/2 * cos( qx * x + qy * y ); 
				hq[(iqx*N_QY+iqy)*2+1] -=  (z[1]+z[2])/2 * sin( qx * x + qy * y ); 
//				hq[(iqx*N_QY+iqy)*2+0] +=  (z[0]) * cos( qx * x + qy * y ); 
//				hq[(iqx*N_QY+iqy)*2+1] -=  (z[0]) * sin( qx * x + qy * y ); 
#endif
			}
		}

		for( int iq = 0; iq < N_QX *N_QY; iq++ )
		{
			hq[2*iq+0] /= (N_QX*N_QY/2);
			hq[2*iq+1] /= (N_QX*N_QY/2);
		}
	
		double Upper[10];
		double NUpper[10];
		double Lower[10];
		double NLower[10];
	
		memset( Upper, 0, sizeof(double) * 10 );	
		memset( NUpper, 0, sizeof(double) * 10 );	
		memset( Lower, 0, sizeof(double) * 10 );	
		memset( NLower, 0, sizeof(double) * 10 );	

#ifdef DENSITY_CURVATURE_CORR
		for( int ul = 0; ul < 2; ul++ )
		for( int ix = 0; ix < N_QX; ix++ )
		for( int iy = 0; iy < N_QY; iy++ )
		{
			double rho = nsum[(ix*N_QY+iy)*3+1+ul];
			double x = ix * La / N_QX;
			double y = iy * Lb / N_QY;

			double ctot = 0;
			for( int iqx = 0; iqx < N_QX; iqx++ )
			for( int iqy = 0; iqy < N_QY; iqy++ )
			{

				double qx = 2 * M_PI * (-dqx+iqx) / La;
				double qy = 2 * M_PI * (-dqy+iqy) / Lb;
				double q = sqrt(qx*qx+qy*qy);

				if( q > 0.15 ) continue;

				double use_hq[2] = { hq[2*(iqx*N_QY+iqy)+0], hq[2*(iqx*N_QY+iqy)+1] };

				double cr_c   =  qx*qx * cos( qx * x + qy * y ) * use_hq[0]; 
				cr_c  +=  qy*qy * cos( qx * x + qy * y ) * use_hq[0]; 
				double cr_s   =- qx*qx * sin( qx * x + qy * y ) * use_hq[1];
				cr_s  +=- qy*qy * sin( qx * x + qy * y ) * use_hq[1];

				ctot += cr_c + cr_s;
			}
		
			int irho = (int)lround(rho);
		
			if( ul == 0 )
			{
				Upper[irho] += ctot;
				NUpper[irho] += 1;
			}
			else 
			{
				Lower[irho] += ctot;
				NLower[irho] += 1;
			}

//			printf("ul %d rho %lf c %lf\n", ul, rho, ctot );
		}

		printf("UPPER %le %le %le %le %le\n", Upper[0] / NUpper[0], Upper[1] / NUpper[1], Upper[2] / NUpper[2], Upper[3] / NUpper[3], Upper[4] / NUpper[4] );
		printf("LOWER %le %le %le %le %le\n", Lower[0] / NLower[0], Lower[1] / NLower[1], Lower[2] / NLower[2], Lower[3] / NLower[3], Lower[4] / NLower[4] );
#endif

	
		fprintf(fq, "# %lf %lf %lf", La, Lb, Lc ); 

		for( int iqx = 0; iqx < N_QX; iqx++ )
		for( int iqy = 0; iqy < N_QY; iqy++ )
		{
			double qx = 2 * M_PI * (-dqx+iqx) / La;
			double qy = 2 * M_PI * (-dqy+iqy) / Lb;

#ifdef BILAYER_MODE
			fprintf(fq, " %lf,%lf", hq[2*(iqx*N_QY+iqy)+0], hq[2*(iqx*N_QY+iqy)+1] );
#else		
			fprintf(fq, " %lf,%lf,%lf,%lf", hq[4*(iqx*N_QY+iqy)+0], hq[4*(iqx*N_QY+iqy)+1], hq[4*(iqx*N_QY+iqy)+2], hq[4*(iqx*N_QY+iqy)+3] );
#endif
		}	
		fprintf(fq, "\n");

		for( int x = 0; x < curNAtoms(); x++ )
			at[x].zap();

		if( f < nframes-1)
		{
			loadFrame( dcdFile, at );

			if( !DCDsuccess() )
			{
				printf("Could not load intended frames\n");
				exit(1);
			}
		
			PBCD( &La, &Lb, &Lc, &alpha, &beta, &gamma );
		}
	}


}
