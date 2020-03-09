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
#include "comparison.h"
#include "proc_definition.h"

#define IJ_ACCELERATOR
			

static double *dmat = NULL;

int qsort_compar( const void * frame1, const void *frame2 )
{
	int i1 = *(int *)frame1;
	int i2 = *(int *)frame2;

	if( dmat[i1] < dmat[i2] )
		return -1;
	else if( dmat[i1] > dmat[i2] )
		return 1;
	return 0;
}

void doKMedoidClusteringLR( 
	aStructure *xyz, int nxyz, int nat, int KMED, int *outputMedoids, double r, double theta, int Kmax)
{
	use_atoms = nat;

//	FILE *tempDump = fopen("temp.dump","w");
//	printf("Dumping xyzs to a temp file for debugging.\n");
//	fwrite( &nat, sizeof(int), 1, tempDump );
//	fwrite( &nxyz, sizeof(int), 1, tempDump );
//	fwrite( xyz, sizeof(double), nxyz*nat*3, tempDump );
//	nxyz = 1000;
//	fclose(tempDump);
//	printf("Done.\n");
//	fflush(stdout);
	double bestLowerBound = -1e10;
	double bestUpperBound = 1e10;
	printf("nxyz: %d\n", nxyz );
	double alpha = 0.1;
	int do_swap_config = is_sym();

//	int maxK = 1000;
	// the distance between medoids.
	double *dij = (double *)malloc( sizeof(double) * nxyz * nxyz );
	char *xij = (char *)malloc( sizeof(char) * nxyz * nxyz );
	memset( xij, 0, sizeof(char) * nxyz * nxyz );

	double *gi_cur = (double *)malloc( sizeof(double) * nxyz );
	double *gi_prev = (double *)malloc( sizeof(double) * nxyz );
	memset( gi_prev, 0, sizeof(double) * nxyz ); 
	int *curMedoids = (int *)malloc( sizeof(int) * KMED ); 
	int *bestMedoids = (int *)malloc( sizeof(int) * KMED );

	double *muik = (double *)malloc( sizeof(double) * nxyz );
	memset( muik, 0, sizeof(double) * nxyz );
	double *muik_prev = (double *)malloc( sizeof(double) * nxyz );
	memset( muik_prev, 0, sizeof(double) * nxyz );

	int *nz_ij = (int *)malloc( sizeof(int) * nxyz * nxyz );
	int *nnz_ij = (int *)malloc( sizeof(int) * nxyz );

	for( int i = 0; i < nxyz; i++ )
	{
		printf("Processing structure %d.\n", i );
		fflush(stdout);
		for( int j = i; j < nxyz; j++ )
		{
			if( i == j ) 
				dij[i*nxyz+j] = value_comparison( xyz+i, xyz+j, do_swap_config );
			else
				dij[i*nxyz+j] = value_comparison( xyz+i, xyz+j, do_swap_config );  
	
//			if( dij[i*nxyz+j] > 1000 )
//			{
//				printf("eh?\n");
//				getChi2Symm( xyz+i*3*nat, xyz+j*3*nat, nat ); 
//				exit(1);
//			}
	
			dij[j*nxyz+i] = dij[i*nxyz+j];
		}	
	}	

#ifdef IJ_ACCELERATOR	
	int *ij_sorter = (int *)malloc( sizeof(int) * nxyz * nxyz );
	for( int i = 0; i < nxyz; i++ )
	{

		int *t_sorter = ij_sorter + i * nxyz;	
		for( int j = 0; j < nxyz; j++ )
			t_sorter[j]= j;
		printf("Sorting %d.\n", i );
#ifdef BUBBLE_SORT
		int done = 0;


		while( !done )
		{
			done = 1;

			for( int j = 0; j < nxyz-1; j++ )
			{
				if( dij[i*nxyz+t_sorter[j]] > dij[i*nxyz+t_sorter[j+1]] )
				{
					int t = t_sorter[j];
					t_sorter[j] = t_sorter[j+1];
					t_sorter[j+1] = t; 
					done = 0;
				}
			}
		}
#else
		dmat = dij + i * nxyz;
	
		qsort( t_sorter, nxyz, sizeof(int), qsort_compar ); 

		if( dmat[t_sorter[0]] > dmat[t_sorter[1]] )
		{
			printf("CRAP!\n");
			exit(1);
		}
#endif
	}
#endif

//	for( int j = 0; j < nxyz; j++ )
//		muik_prev[j] = 0.25;

	int iter_tot = 500;
	//double r = 0.95;

	int maxK = Kmax;

	double nit = iter_tot / ( (1-pow(r,maxK))/(1-r) );;
	
	printf("maxK: %d\n", maxK);

	double lambda = 0.01;
	double gmag = 1.0;
				
	double *jvals = (double *)malloc( sizeof(double) * nxyz );

	int n_iters_real = 0;
	for( int k = 0; k < maxK; k++ )
	{
		double alpha = lambda;
	
		if( k > 0 && fabs(gmag) > 1e-7 )
		{
			alpha = lambda * (bestUpperBound-bestLowerBound) / gmag;
			lambda *= r;
			nit *= r;

			if( nit < 1.01 ) nit = 1.01;
		}

		for( int sub_iter = 0; sub_iter < lround(nit) && fabs(bestUpperBound-bestLowerBound) > 0.01; sub_iter++, n_iters_real++ )
	//	while( 1 )
		{
			// Maximize Lu
		
			memset( nnz_ij, 0, sizeof(int) * nxyz );
		
			int *curBest = (int *)malloc( sizeof(int) * KMED );
			double *valBest = (double *)malloc( sizeof(double) * KMED );
			int nbGot = 0;
	
			memset( xij, 0, sizeof(char) * nxyz * nxyz );
		
			if( k == 0 && 0 )
			{
				for( int j = 0; j < KMED; j++ )
					curBest[j] = j * (nxyz / KMED);
	
				for( int i = 0; i < nxyz; i++ )
				{
					double smallestD = 1e10;
					int bestm = 0;
					for( int m = 0; m < KMED; m++ )
					{
						if( dij[i*nxyz+curBest[m]] < smallestD )
						{
							bestm = m;
							smallestD = dij[i*nxyz+curBest[m]];
						}
					}

					nz_ij[i*nxyz+nnz_ij[i]] = curBest[bestm];
					nnz_ij[i] += 1; 
					//xij[i*nxyz+curBest[bestm]] = 1;
		
				}
			}
			else
			{
#ifdef IJ_ACCELERATOR
				memset( jvals, 0, sizeof(double) * nxyz );

				for( int i = 0; i < nxyz; i++ )
				{
					// find loop closure by bisection.

					int low = 0;
					int high = nxyz-1;

					int done = 0;
					double *t_dij = dij+i*nxyz;
					int *t_sorter = ij_sorter+i*nxyz;
					while( high-low > 1 )
					{
						int test = (low+high)/2;
		
						if( muik[i] > t_dij[t_sorter[test]] )
							low = test;
						else
							high = test;
					}

					for( int jt = 0; jt < high; jt++ )
						jvals[t_sorter[jt]] += t_dij[t_sorter[jt]] - muik[i];

/*					int jt = 0;
					for( jt = 0; jt < nxyz; jt++ )
					{
						int j = ij_sorter[i*nxyz+jt];

						if( muik[i] < dij[j*nxyz+i] )
							break;

						jvals[j] += dij[j*nxyz+i] - muik[i];
					}

					printf("terminated at %d would have terminated at %d.\n", jt, high );
*/
				}
			
#endif
				// assign medoid centers.
				for( int j = 0; j < nxyz; j++ )
				{
#ifdef IJ_ACCELERATOR
					double valj = jvals[j];
#else
					double valj = 0;
					for( int i = 0; i < nxyz; i++ )
					{
						if( dij[j*nxyz+i] - muik[i] < 0 )
							valj += dij[j*nxyz+i] - muik[i];
					}
#endif
		//			printf("struct %d valj: %le\n", j, valj );
		
					if( nbGot < KMED )
					{
						curBest[nbGot] = j;
						valBest[nbGot] = valj;
						nbGot++;
					}
					else
					{
						double cur_largest = -1e10;
						int tpos = 0;
	
						for( int l = 0; l < nbGot; l++ )
						{
							if( valBest[l] > cur_largest )
							{
								cur_largest = valBest[l];
								tpos = l;
							}
						}
	
						if( valj < cur_largest )
						{
							valBest[tpos] = valj;
							curBest[tpos] = j;
						}
	
		/*				printf("cur set:");
						for( int z = 0; z < nbGot; z++ )
							printf(" %lf", valBest[z] );
						printf(" :");
						for( int z = 0; z < nbGot; z++ )
							printf(" %d", curBest[z] );
						printf("\n");
		*/	
					}
				}
			//	exit(1);
			}
	
			memcpy( curMedoids, curBest, sizeof(int) * KMED );
	//		if( k != 0  && 0 )
			{
				for( int i = 0; i < nxyz; i++ )
				{
					double *t_dij = dij+i*nxyz;
					for( int tj = 0; tj < KMED; tj++ )
					{
						int j = curMedoids[tj];
		
						// xij only non-zero if j is a medoid.
		
						if( i == j )
						{
						//	xij[i*nxyz+j] = 1;
							nz_ij[i*nxyz+nnz_ij[i]] = j;
							nnz_ij[i] += 1; 
						}
						else	
						{ 
							if( t_dij[j] - muik[i] < 0 )
							{
							//	xij[i*nxyz+j] = 1;						
								nz_ij[i*nxyz+nnz_ij[i]] = j;
								nnz_ij[i] += 1; 
							}
						}
					}				
				}
			}
			
			double L = 0;
	
			for( int i = 0; i < nxyz; i++ )
			{
				L += muik[i];
	
				for( int tj = 0; tj < nnz_ij[i]; tj++ )
				{
					int j = nz_ij[i*nxyz+tj];

					L += dij[i*nxyz+j] - muik[i];
				}
/*
				for( int j = 0; j < nxyz; j++ )
				{
					if( xij[i*nxyz+j] != 0 )
						L += dij[i*nxyz+j] - muik[i]; 
				}*/
			}	
			if( L > bestLowerBound )
				bestLowerBound = L;
	
			double curL = 0;
	
			for( int i = 0; i < nxyz; i++ )
			{
				double *t_dij = dij+i*nxyz;
				double smallestD = 1e10;
				int bestm = 0;
				for( int m = 0; m < KMED; m++ )
				{
					if( t_dij[curMedoids[m]] < smallestD )
					{
						bestm = m;
						smallestD = dij[i*nxyz+curMedoids[m]];
					}
				}
	
				curL += dij[i*nxyz+curMedoids[bestm]];
			}
	
			if( curL < bestUpperBound )
			{
				bestUpperBound = curL;
				memcpy( bestMedoids, curMedoids, sizeof(int) * KMED ); 
			}
	
			gmag = 0;
	
			for( int i = 0; i < nxyz; i++ )
			{
				gi_cur[i] = theta + (1-theta) * gi_prev[i] - theta*nnz_ij[i];
				

/*				for( int j = 0; j < nxyz; j++ )
				{
					if( xij[i*nxyz+j] )
						gi_cur[i] -= 1;
				}
*/	
				gmag += gi_cur[i] * gi_cur[i];
			}
	
	
	
			memcpy( gi_prev, gi_cur, sizeof(double) * nxyz );
	
			double max_mu = -1e10;

			for( int i = 0; i < nxyz; i++ )
			{
				muik[i] = muik_prev[i] + alpha * gi_cur[i];

				if( muik[i] > max_mu )
					max_mu = muik[i];
			}
			
/*			printf("mui: ");
			for( int z = 0; z < nxyz; z++ )
				printf("%lf ", muik[z] );
			printf("\n");
*/	
			memcpy( muik_prev, muik, sizeof(double) * nxyz );
	
//			for( int x = 0; x < KMED; x++ )
//				printf("%d ", curMedoids[x] );
		}

		printf("Upper bound: %lf Lower bound: %lf\n", bestUpperBound, bestLowerBound );
		fflush(stdout);
	}	
	
	printf("Niters real: %d\n", n_iters_real ); 

	int sorter[KMED];
	double pop[KMED];
	double sumd[KMED];
	double maxd[KMED];

	for( int i = 0; i < KMED; i++ )
	{
		sumd[i] = 0;
		maxd[i] = -1e10;
		sorter[i] = i;
		pop[i] = 0;
	}	
	for( int i = 0; i < nxyz; i++ )
	{
		double smallestD = 1e10;
		int bestm = 0;
		for( int m = 0; m < KMED; m++ )
		{
			if( dij[i*nxyz+curMedoids[m]] < smallestD )
			{
				bestm = m;
				smallestD = dij[i*nxyz+curMedoids[m]];
			}
		}
	
		pop[bestm] += 1;
		sumd[bestm] += dij[i*nxyz+curMedoids[bestm]];
		if( dij[i*nxyz+curMedoids[bestm]] > maxd[bestm] )
			maxd[bestm] = dij[i*nxyz+curMedoids[bestm]];
	}

	for( int i = 0; i < KMED; i++ )
		sumd[i] /= pop[i];

	int done = 0;

	while( !done )
	{
		done = 1;

		for( int i = 0; i < KMED-1; i++ )
		{
			if( pop[sorter[i]] < pop[sorter[i+1]] )
			{
				int t = sorter[i];
				sorter[i] = sorter[i+1];
				sorter[i+1] = t;
				done = 0;
			}
		}
	}

	printf("mat:\n");
	for( int i = 0; i < KMED; i++ )
	{
		for( int j = 0; j < KMED; j++ )
		{
			printf("%lf ", dij[curMedoids[sorter[i]]*nxyz+curMedoids[sorter[j]]] ); 
		}
		printf("\n");
	}

	for( int i = 0; i < KMED; i++ )
	{
		outputMedoids[i] = curMedoids[sorter[i]];
	
		printf("Cluster centered on %d has pop %lf avd %lf maxd %lf\n",
			outputMedoids[i], pop[sorter[i]], sumd[sorter[i]], maxd[sorter[i]] );
	}


//	memcpy( outputMedoids, bestMedoids, sizeof(int) * KMED );
}



















