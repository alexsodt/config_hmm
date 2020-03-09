#include <stdio.h>
#include <stdlib.h>
#include "pdb.h"
#include "dcd.h"
#include <math.h>
#include <sys/time.h>
#define FORM_1

//#define ONE_STEP_TEST

//#define DO_ONE_MODE
//#define JUST_COS

double lipid_area = 65;

int main( int argc, char **argv )
{
	if( argc < 5	)
	{
		printf("Syntax: make_dbg_data ngridx ngridy kc nsteps c0 [indep_leaflet=0]\n");
		exit(1);
	}

	timeval theTime;
	gettimeofday( &theTime, NULL);
	srand(theTime.tv_usec );	
	int ngridx = atoi(argv[1]);
	int ngridy = atoi(argv[2]);
	double kc = atof(argv[3]);
	double kcm = kc/2; // monolayer kc
	int nsteps = atoi(argv[4]);
	int indep = 0;
	double c0 = 0;

	if( argc > 5 )
		c0 = atof(argv[5]);
	double one_step_mag = 1.0;
#ifdef ONE_STEP_TEST
	one_step_mag = atof(argv[5]);
	c0 = 0;
#endif

	int iat = 1;
	int ires = 1;

	char atname1[256];
	char atname2[256];
	char atname3[256];
	char resname[256];
	char aresname[256];
	char segname[256];

	sprintf(atname1, "C21");
	sprintf(atname2, "C22");
	sprintf(atname3, "C23");
	sprintf(resname, "DBG");
#ifdef ONE_STEP_TEST
	sprintf(aresname, "DBG");
#else
	sprintf(aresname, "TRG");
#endif
	sprintf(segname, "DBG");

	FILE *dbgFile = fopen("dbg.pdb","w");

	struct atom_rec *at = (struct atom_rec *)malloc( sizeof(atom_rec) * 2 * ngridx * ngridy *3 );

#ifdef ONE_STEP_TEST
	double cc_u =  one_step_mag;
	double sc_u =  0;
	double cc_l =  -one_step_mag;
	double sc_l =  0;
#else
	double cc_u = 0;
	double sc_u = 0;
	double cc_l = 0;
	double sc_l = 0;
#endif
	int icntr=0;

	double dx = 7.5;
	double dy = 7.5;

	double *xy_sites = (double *)malloc( sizeof(double) * 2 * ngridx * ngridy );
	

	for( int ix = 0; ix < ngridx; ix++ )
	for( int iy = 0; iy < ngridy; iy++ )
	{
		xy_sites[(ix*ngridy+iy)*2+0] = (ix+0.5) * dx;
		xy_sites[(ix*ngridy+iy)*2+1] = (iy+0.5) * dy;
	}

	for( int leaflet = 0; leaflet < 2; leaflet++ )
	for( int ix = 0; ix < ngridx; ix++ )
	for( int iy = 0; iy < ngridy; iy++ )
	{
		for( int t = 0; t < 3; t++ )
		{

			at[icntr].bead = iat;
			at[icntr].res = ires;
		
			if( t == 0 ) at[icntr].atname = atname1;
			if( t == 1 ) at[icntr].atname = atname2;
			if( t == 2 ) at[icntr].atname = atname3;

			at[icntr].segid = segname;

			if( ix == 0 && iy == 0 && leaflet == 1)
				at[icntr].resname = aresname;
			else
				at[icntr].resname = resname;
	
			at[icntr].altloc = ' ';
			at[icntr].chain  = ' ';

			at[icntr].x = ix * dx;
			at[icntr].y = iy * dy;
			at[icntr].z = 0;


			icntr++;
			iat++;
		}
		ires++;
	}

#ifdef ONE_STEP_TEST
	nsteps=1;
	int ninner = 1;
	int nequil = 0;
#else
	int ninner = 100;
	int nequil = 100;
#endif
	double mag = 1;

#ifdef ONE_STEP_TEST
	mag = 0;
#endif

	double cure = 0;

	double area = ngridx*ngridy*dx*dy;

	int m_x = 2;
	int m_y = 0;
	
	double Lx = dx*ngridx;
	double Ly = dy*ngridy;
	double Lz = 1000;
	double qx = 2 * M_PI * m_x / Lx;
	double qy = 2 * M_PI * m_y / Ly;
	double q = sqrt(qx*qx+qy*qy);

	double beta = 1.0/0.616;
	
	FILE *dcdFile = fopen("dbg.dcd", "wb");

	writeDCDHeader( dcdFile, 3*ngridx*ngridy*2, nsteps );

	int init_pdb = 1;

	int target_site[2] = {0,0};

	double expec_hq = 0;
	double nexpec_hq = 0;
	double fac = 2;
			
	double px = xy_sites[(target_site[0]*ngridy+target_site[1])*2+0];
	double py = xy_sites[(target_site[0]*ngridy+target_site[1])*2+1];
	double c_at_site = -fac*qx*qx*cc_u*cos(qx*px+qy*py) - fac*qx*qx*sc_u*sin(qx*px+qy*py);
	double pe = 0.5*(kcm)*lipid_area*((c_at_site-c0)*(c_at_site-c0) - (c_at_site)*(c_at_site));

	cure += pe;

	for( int s = 0; s < nequil + nsteps; s++ )
	{
		for( int inner = 0; inner < 100; inner++ )
		{
			double save[4] = { cc_u,sc_u,cc_l,sc_l };
	
			double dc_u = mag * 2*(rand()/(double)RAND_MAX-0.5);	
			double ds_u = mag * 2*(rand()/(double)RAND_MAX-0.5);	
		
			double dc_l = mag * 2*(rand()/(double)RAND_MAX-0.5);	
			double ds_l = mag * 2*(rand()/(double)RAND_MAX-0.5);	

			int delx = rand() % 3 -1; 
			int dely = rand() % 3 -1; 

#ifdef ONE_STEP_TEST
			delx = 0;
			dely = 0;
#endif

			if( s == 0 && inner == 0 )
			{
				// must start in 0, 0 to be consistent with pdb.
				delx = 0;
				dely = 0;
			}

			target_site[0] += delx;
			target_site[1] += dely;
			
			if( target_site[0] < 0 ) target_site[0] += ngridx;
			if( target_site[0] >= ngridx ) target_site[0] -= ngridx;
			if( target_site[1] < 0 ) target_site[1] += ngridy;
			if( target_site[1] >= ngridy ) target_site[1] -= ngridy;

			if( !indep )	
			{	
				dc_l = -dc_u;
				ds_l = -ds_u;
			}

#ifdef DO_ONE_MODE
			ds_u = dc_u;
			ds_l = dc_l;
#endif

#ifdef JUST_COS
			ds_u = 0;
			ds_l = 0;
#endif

			cc_u += dc_u;
			cc_l += dc_l;
			
			sc_u += ds_u;			
			sc_l += ds_l;
	
			double px = xy_sites[(target_site[0]*ngridy+target_site[1])*2+0];
			double py = xy_sites[(target_site[0]*ngridy+target_site[1])*2+1];
			double c_at_site = -fac*qx*qx*cc_u*cos(qx*px+qy*py) - fac*qx*qx*sc_u*sin(qx*px+qy*py);
			double pe = 0.5*(kcm)*lipid_area*((c_at_site-c0)*(c_at_site-c0) - (c_at_site)*(c_at_site));
			
			double cc = (cc_u-cc_l)/2;
			double sc = (sc_u-sc_l)/2;

			double energy = 
				fac*fac*0.25 * (kc) * q*q*q*q * sc * sc * area + 	
				fac*fac*0.25 * (kc) * q*q*q*q * cc * cc * area;

			energy += pe;

			double dE = energy - cure;
	
			double pr = exp(-beta * dE );
	
			double rn = (double)rand()/(double)RAND_MAX;
	
			if( rn < pr )
			{		
				cure = energy;
			}
			else
			{
				cc_u = save[0];
				sc_u = save[1];
				cc_l = save[2];
				sc_l = save[3];

				target_site[0] -= delx; 
				target_site[1] -= dely; 

				if( target_site[0] < 0 ) target_site[0] += ngridx;
				if( target_site[0] >= ngridx ) target_site[0] -= ngridx;
				if( target_site[1] < 0 ) target_site[1] += ngridy;
				if( target_site[1] >= ngridy ) target_site[1] -= ngridy;
			}
		
			cc = (cc_u-cc_l)/2;
			sc = (sc_u-sc_l)/2;
		
			if( s >= nequil )
			{
				expec_hq += cc*cc + sc*sc;
				nexpec_hq += 1;
			}
		}
		
		int icntr = 0;
		for( int leaflet = 0; leaflet < 2; leaflet++ )
		for( int ix = 0; ix < ngridx; ix++ )
		for( int iy = 0; iy < ngridy; iy++ )
		{
			double x = xy_sites[2*(ix*ngridy+iy)+0];
			double y = xy_sites[2*(ix*ngridy+iy)+1];
			
			double dz[3] = { -1, 0, 1  };

			double cc = cc_u;
			double sc = sc_u;

			if( leaflet == 0)
			{
				cc = cc_l;
				sc = sc_l;
			}

			double h = fac*cc * cos( qx * x + qy * y ) + fac*sc * sin(qx * x + qy * y ); 
			double nrm[3] = { -fac*sc * qx * cos( qx * x + qy * y) + fac*cc * qx * sin( qx * x + qy * y), 
					  -fac*sc * qy * cos( qx * x + qy * y) + fac*cc * qy * sin( qx * x + qy * y),
					1 - fac *pow( sc * qx * cos( qx *x ) - cc * qx * sin(qx*x +qy*y), 2.0) };

			if( leaflet == 0 )
			{
				// sign of x/y coordinates reversed by sign in sc/cc
				nrm[2] *= -1;
				h *= -1;
			}
 
			double hoff[2] = { -10, 10 };
			for( int t = 0; t < 3; t++ )
			{

				double r[3] = { x + (dz[t] * nrm[0]), y + (dz[t] * nrm[1]), hoff[leaflet] + ( h + dz[t] * nrm[2]) };
			
				at[icntr].x = r[0]; 
				at[icntr].y = r[1]; 
				at[icntr].z = r[2];
	
		
				if( init_pdb )
					printATOM( dbgFile, at[icntr].bead, at[icntr].res, at+icntr );		
				
				icntr++;
				iat++;		
			}
			ires++;
		}
		init_pdb = 0;

//		printf("%le %le\n", cc_u, sc_u );

		if( s >= nequil )
		{

			int base_site = 3* ngridx*ngridy;

			for( int t = 0; t < 3; t++ )
			{
				double tr[3] = { at[base_site+t].x, at[base_site+t].y, at[base_site+t].z };
				at[base_site+t].x = at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].x;
				at[base_site+t].y = at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].y;
				at[base_site+t].z = at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].z;
				at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].x=tr[0];
				at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].y=tr[1];
				at[base_site+(target_site[0]*ngridy+target_site[1])*3+t].z=tr[2];
			}
			writeFrame( dcdFile, at, Lx, Ly, Lz );
		}
	} 
	fclose(dbgFile);

	printf("expec_hq^2: %le\n", expec_hq / nexpec_hq );
	double expec_u2 = expec_hq/nexpec_hq;
#ifdef DO_ONE_MODE
	double expec_kc = 2 /(beta*area*q*q*q*q*fac*fac*expec_u2);
#else
	double expec_kc = 4 /(beta*area*q*q*q*q*fac*fac*expec_u2);
#endif
	printf("kc obs %le\n", expec_kc );
}
