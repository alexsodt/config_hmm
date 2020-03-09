#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <sys/time.h>

#define ALPHA	0.05

#include "HMMDirective.h"
#include "HMM.h"
#include "Fuzzy.h"
#include "nr_lin_min.h"

double pseudo_add = 0.0;

#define DO_PROT
int howManyObservables( int NLipids )
{
	return (NLipids+2)*(NLipids+1)/2
#ifdef DO_PROT
	+1
#endif
;
}


int globalN_OBSERVABLES = 0;//howManyObservables(6);
extern int deactivate_N;

extern double global_f;
void setupModel( HMMSystem &theSystem );
void makeSchemeTrainer( HMMSystem &theSystem, int scheme );
int my_isnan( double a);

#define MAX_TALLY 10

void codeToCounts( int obs, int cnts[3], int ndo )
{
	static int done = 0;
	static int *array=  NULL;//[(N_TALLY+1)*(N_TALLY+1)];
	static int maxdone = 0;

	if( done == 0 )
	{
		array = (int *)malloc( sizeof(int) * (MAX_TALLY+1)*(MAX_TALLY+1) * 3 );

		int cur = 0;

		for( int Neither = 0; Neither <= ndo; Neither++ )
		for( int y = 0; y <= Neither; y++ )
		{
			int x = Neither - y;

			int n_cmp_3 = ndo - Neither;
			int n_cmp_2 = x;
			int n_cmp_1 = y; 

			array[cur*3+0] = n_cmp_1; // MOL1 (chol)
			array[cur*3+1] = n_cmp_2; // MOL2 (DPPC)
			array[cur*3+2] = n_cmp_3; // MOL3 (DOPC)
			maxdone = cur;
			cur++;
		}

		array[cur*3+0] = 0;
		array[cur*3+1] = 0;
		array[cur*3+2] = 0;
		cur++;

		done = 1;
	}	

	if( obs < 0 || obs > maxdone )		
	{
		printf("exceeded state code!\n");
	}

	cnts[0] = array[obs*3+0];
	cnts[1] = array[obs*3+1];
	cnts[2] = array[obs*3+2];
}

int main( int argc, char **argv )
{
	int ndo = 6;
	FILE *read_n_do = fopen("read.ndo","r");
	if( read_n_do )
	{
		fscanf(read_n_do, "%d", &ndo );
		fclose(read_n_do);
	}


        struct timeval tp;

       int ierr;
#ifdef PARALLEL
          if ( (ierr = MPI_Init(&argc, &argv)) ) {
            exit(1);
          }
#endif
        
	int nprocs=1, taskid=0;
#ifdef PARALLEL
        /* Change the random number if this is part of a parallel run */
        if ( (ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs)) ) {
          exit(1);
        }
        if ( (ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid)) ) {
          exit(1);
        }
#endif
        
	gettimeofday( &tp, NULL );

	int seed = tp.tv_usec / nprocs;
	int seed_out = 0;	

	//seed = 5;

#ifdef PARALLEL
	MPI_Allreduce( &seed, &seed_out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
#else
	seed_out = seed;
#endif
	seed_out /= nprocs;

	printf("Using seed %d.\n", seed_out );
	
	srand(seed_out);
	//srand(0);

	if( argc < 3 )
	{
		printf("Syntax: LipidShellHMM input_file network_file\n");
		return 0;
	}

	int getNObservables( const char *fileName );
	globalN_OBSERVABLES = getNObservables( argv[2] );
	
	FILE *theFile = fopen( argv[1], "r" );

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1] );
		return 0;
	}

	HMMCommands *theCommands = new HMMCommands( theFile );	
	HMMSystem *theSystem = new HMMSystem();

	if( (taskid != 0) )	
		freopen( "/dev/null", "w", stdout ); 


	int activate_copy = 0;
	int clear_aqueous = 0;
	int copy_from = 0;
	int did_assess = 0;

	double cur_sc_weighting = 1.0;
	double cur_g_weighting = 1.0;
	double cur_s_weighting = 1.0;
	
	int MakeGeneralNetwork( HMMSystem &theSystem, int ignore, const char *fileName );

	MakeGeneralNetwork( *theSystem, theCommands->allDirectives[0]->network, argv[2] );
			
	int nstates = theSystem->getNStates();
	double *component_values = (double *)malloc(sizeof(double) * 3 * nstates );
	memset( component_values, 0, sizeof(double) *3 * nstates );

	double last_val = 0;
	const char *tempSave = "tmp.model.save";

	int start_at = 0;
	int stop_at  = 1000000000;
	int do_randomize = 0;

	for( int c = 0; c < theCommands->ndirectives; c++ )
	{
		if( theCommands->allDirectives[c]->command == COMMAND_TRAIN )
		{
			if( do_randomize )
				theSystem->randomizeTimeSequence( theCommands->allDirectives[c] );
			
	
			check_surface_consistency( *theSystem, *(theCommands->allDirectives[c]) ); 
			surface_grad_opt( *theSystem, *(theCommands->allDirectives[c]) ); 

/*
			double A = 0.0;
			double prev_free = 0;
	
			int doBaumWelch = 1;
				
			theSystem->setPseudocounts( theCommands->allDirectives[c] );
			
				for( int s = 0; s < nstates; s++ )
				{
					for( int o = 0; o < globalN_OBSERVABLES; o++ )
						printf("%lf ", theSystem->p_obs[s*globalN_OBSERVABLES+o] );
					printf("\n");
				}
				for( int s1 = 0; s1 < nstates; s1++ )
				{
					for( int s2 = 0; s2 < nstates; s2++ )
					{	
						printf(" %lf", theSystem->t_ab[s1*nstates+s2] );
					}
					printf("\n");
				}

			int nPreIter = 10;
			int best_init = 0;
			double best_val = -1e20;	
			int switch_to_best = 1;
			int lenPreIter = 25;

			for( int iter = 0; iter < (nPreIter*lenPreIter+50); iter++, A *= 0.95 ) //( iter < 500 ? 0.99 : (iter < 750 ? 0.95 : 0.0) ) )
			{
				theSystem->zeroNextPMat();
	
				if( iter != 0 && (iter % lenPreIter == 0) && (iter/lenPreIter) < nPreIter )
				{
					if( !best_init || (last_val > best_val) )
					{
						printf("best_init: %d last_val: %lf best_val: %lf, saving.\n", best_init, last_val, best_val );
						best_val = last_val;
						theSystem->SaveModel(tempSave);
						best_init = 1;
					}
					
					theSystem->randomizeObs();
				}
				else if( switch_to_best && ((iter/lenPreIter) >= nPreIter) )
				{
					printf("Loading model.\n");
					last_val = best_val;
					theSystem->ReadModel(tempSave);
					switch_to_best = 0;
				}

				double log_sum = 0;
				double log_sum_clamped = 0;
				double log_sum_free = 0;


	
				for( int d = taskid; d < theCommands->allDirectives[c]->n_data; d += nprocs )
				{
					int s_off = start_at;

					int len = theCommands->allDirectives[c]->translated_len[d];
					if( len < s_off ) continue;


					len -= s_off;					
					
					if( len > stop_at - start_at ) len = stop_at - start_at;	

					char *wild_string = (char *)malloc( sizeof(char) * (len+1) );
					for( int s = 0; s < len; s++ )
						wild_string[s] = '.';
					wild_string[len] = '\0';
					double *alpha_constrained = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha = (double *)malloc( sizeof(double) * len * nstates );
					double *beta  = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_norm  = (double *)malloc( sizeof(double) * len * nstates );
	

#if 0 // not sure what I was doing here.	
	//				calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d], theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
	//				calculateBetaForStringConstrained( beta_constrained, beta_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d],  theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
#endif					

#if 0 //unconstrained
					calculateAlphaForStringConstrained( alpha, alpha_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
					calculateBetaForStringConstrained( beta, beta_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
#endif
					calculateAlphaForStringConstrained( alpha, alpha_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d]+s_off,  len, *theSystem );
					calculateBetaForStringConstrained( beta, beta_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d]+s_off,  len, *theSystem );
	//				printf("%le %le\n", alpha_constrained_norm[len-1], alpha_norm[len-1] );	
					if( theCommands->allDirectives[c]->command != COMMAND_CHECK )
					{
						theSystem->incrementNextPMat_Free( alpha, beta, theCommands->allDirectives[c]->translated_observable_string[d]+s_off,  len, theCommands->allDirectives[c]->weights[d] );				
		
						double pfree = 0;
						double pclamped = 0;
		
						for( int s = 0; s < nstates; s++ )
						{
							pfree += alpha[(len-1)*nstates+s] * theCommands->allDirectives[c]->weights[d];
							pclamped += alpha_constrained[(len-1)*nstates+s] * theCommands->allDirectives[c]->weights[d];
						}
	
						log_sum += alpha_constrained_norm[len-1] * theCommands->allDirectives[c]->weights[d] - alpha_norm[len-1] * theCommands->allDirectives[c]->weights[d];			
		
						log_sum_clamped += alpha_constrained_norm[len-1] * theCommands->allDirectives[c]->weights[d];
						log_sum_free += alpha_norm[len-1] * theCommands->allDirectives[c]->weights[d];
		
						double logFac = alpha_constrained_norm[len-1]* theCommands->allDirectives[c]->weights[d] - alpha_norm[len-1] * theCommands->allDirectives[c]->weights[d];
			
					}
					free(wild_string);
					free(alpha);
					free(beta);
					free(alpha_constrained);
					free(beta_constrained);
					free(alpha_norm);
					free(beta_norm);
					free(alpha_constrained_norm);
					free(beta_constrained_norm);
				}
				
				theSystem->addNoise(A);
	
			
				// Now combine em all 
				double f_sum = 0;
				MPI_Allreduce( &log_sum_clamped, &f_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				double F_sum = 0;
				MPI_Allreduce( &log_sum, &F_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				double l_s_c_t = 0;
				MPI_Allreduce( &log_sum_clamped, &l_s_c_t, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				double l_s_f_t = 0;
				MPI_Allreduce( &log_sum_free, &l_s_f_t, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	#ifdef FULL_COMM
				double *recv_t_ab_next = (double *)malloc( sizeof(double) * nstates * nstates ); 
				double *recv_p_obs_next = (double *)malloc( sizeof(double) * nstates * globalN_OBSERVABLES ); 
				memset( recv_t_ab_next, 0, sizeof(double) * nstates * nstates );
				memset( recv_p_obs_next, 0, sizeof(double) * nstates * globalN_OBSERVABLES );
	
				int ar_err;
	
				MPI_Allreduce( theSystem->t_ab_next_free, recv_t_ab_next, nstates*nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( theSystem->p_obs_next_free, recv_p_obs_next, nstates*globalN_OBSERVABLES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
			
				memcpy( theSystem->t_ab_next_free, recv_t_ab_next, nstates*nstates*sizeof(double) );
				memcpy( theSystem->p_obs_next_free, recv_p_obs_next, globalN_OBSERVABLES*nstates*sizeof(double) );
				
				memset( recv_t_ab_next, 0, sizeof(double) * nstates * nstates );
				memset( recv_p_obs_next, 0, sizeof(double) * nstates * globalN_OBSERVABLES );
				
				MPI_Allreduce( theSystem->t_ab_next_clamped, recv_t_ab_next, nstates*nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( theSystem->p_obs_next_clamped, recv_p_obs_next, nstates*globalN_OBSERVABLES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				
				memcpy( theSystem->t_ab_next_clamped, recv_t_ab_next, nstates*nstates*sizeof(double) );
				memcpy( theSystem->p_obs_next_clamped, recv_p_obs_next, globalN_OBSERVABLES*nstates*sizeof(double) );
	
				free( recv_t_ab_next);
				free( recv_p_obs_next);
	#else
				theSystem->CommunicatePMats();
	#endif
	
				theSystem->addPseudocounts( pseudo_add );
				theSystem->finishNextPMat();
				
				if( theCommands->allDirectives[c]->command == COMMAND_CHECK )
					break;
	
	//			printf("logSum: %.12le [clamped: %.12le free: %.12le]\n", log_sum, log_sum_clamped, log_sum_free );
	
	
				memcpy( theSystem->t_ab, theSystem->t_ab_next_free, nstates*nstates*sizeof(double) );
				memcpy( theSystem->p_obs, theSystem->p_obs_next_free, globalN_OBSERVABLES*nstates*sizeof(double) );
	
//#ifdef PRINT_INTERMEDIATE
				printf("Iteration emission probabilities:\n");
				for( int s = 0; s < nstates; s++ )
				{
					for( int o = 0; o < globalN_OBSERVABLES; o++ )
						printf("%le ", theSystem->p_obs[s+nstates*o] );
					printf("\n");
				}
				printf("Iteration transition probabilities:\n");
				for( int s = 0; s < nstates; s++ )
				{	
					for( int q = 0; q < nstates; q++ )
						printf("%lf ", theSystem->t_ab[s*nstates+q]);
					printf("\n");
				}
//#endif
				printf("iteration log_probability: %lf A: %le\n", l_s_f_t, A);
	
				last_val = l_s_f_t;
				MPI_Barrier(MPI_COMM_WORLD);
	
			}
				int maxl = 0;
	
				for( int d = taskid; d < theCommands->allDirectives[c]->n_data; d += nprocs )
				{
					int len = theCommands->allDirectives[c]->translated_len[d] - start_at;
					if( len > stop_at - start_at ) len = stop_at - start_at;	
					if( len > maxl )
						maxl = len;
				}		
	
				char *tstate = (char *)malloc( sizeof(char) * maxl * theCommands->allDirectives[c]->n_data );
				memset( tstate, 0, sizeof(char) * maxl * theCommands->allDirectives[c]->n_data );
	
	
				for( int d = 0; d < theCommands->allDirectives[c]->n_data; d ++ )
				{
					int s_off = start_at;

					int len = theCommands->allDirectives[c]->translated_len[d];
					if( len < s_off ) continue;


					len -= s_off;					
					
					if( len > stop_at - start_at ) len = stop_at - start_at;	
		
					char *wild_string = (char *)malloc( sizeof(char) * (len+1) );
					for( int s = 0; s < len; s++ )
						wild_string[s] = '.';
					wild_string[len] = '\0';
					double *alpha_constrained = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha = (double *)malloc( sizeof(double) * len * nstates );
					double *beta  = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_norm  = (double *)malloc( sizeof(double) * len * nstates );
	
	
					calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
					calculateBetaForStringConstrained( beta_constrained, beta_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off,  theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
					
//					calculateAlphaForStringConstrained( alpha, alpha_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
//					calculateBetaForStringConstrained( beta, beta_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
	
					int *predicted_state = (int *)malloc( sizeof(int) * len );
					
//					theSystem->estimateClassViterbi( alpha_constrained, beta_constrained, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string, len, predicted_state  );
					theSystem->estimateClassViterbi( alpha_constrained, beta_constrained, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d]+s_off, len, predicted_state  );
	
					for( int t = 0; t < len; t++ )
					{
						int o = theCommands->allDirectives[c]->translated_observable_string[d][s_off+t];
						int nok[3];
						codeToCounts(o,nok,ndo);

						int nmol1 =nok[0];
						int nmol2 =nok[1];
						int nmol3 =nok[2];
	
						component_values[predicted_state[t]*3+0] += nmol1;
						component_values[predicted_state[t]*3+1] += nmol2;
						component_values[predicted_state[t]*3+2] += nmol3;
					
//						if( d == 0 )
//							printf("stategrep %d %d %d %d\n", predicted_state[t], ndppc, ndopc, ncholesterol ); 
					}
	
	
					for( int t = 0; t < len; t++ )
						tstate[t*(theCommands->allDirectives[c]->n_data)+d] = (char)predicted_state[t];
	
					free(wild_string);
					free(alpha);
					free(beta);
					free(alpha_constrained);
					free(beta_constrained);
					free(alpha_norm);
					free(beta_norm);
					free(alpha_constrained_norm);
					free(beta_constrained_norm);
				}

			int do_switch = 0;

//			if( component_values[1]  < component_values[4] )
//				do_switch = 1;

			printf("State trajectory:\n");
			for( int t = 0; t < maxl; t++ )
			{
				printf("%d ", t);
				for( int s = 0; s < theCommands->allDirectives[c]->n_data; s++ )
					printf("%d", (do_switch ? 1 - tstate[t*theCommands->allDirectives[c]->n_data+s] : tstate[t*theCommands->allDirectives[c]->n_data+s] ) );
				printf("\n");
			} 
	
			theSystem->ReportEmissionGroups( theCommands->allDirectives[0]->interpreter );
			*/
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_LOAD )
		{
			printf("Loading model parameters.\n");
			int d= 0;
			char trimmed_file_name[256];
			strcpy(trimmed_file_name, theCommands->allDirectives[c]->raw_observable_string[d] );
			
			char *tptr = trimmed_file_name;

			while( *tptr == ' ' || *tptr == '\t' ) tptr += 1;

			while( tptr[strlen(tptr)-1] == ' ' || tptr[strlen(tptr)-1] == '\t' )
				tptr[strlen(tptr)-1] = '\0';

			theSystem->ReadModel(tptr);
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_SAVE )
		{
			printf("Saving model parameters.\n");
			int d= 0;
			char trimmed_file_name[256];
			strcpy(trimmed_file_name, theCommands->allDirectives[c]->raw_observable_string[d] );
			
			char *tptr = trimmed_file_name;

			while( *tptr == ' ' || *tptr == '\t' ) tptr += 1;

			while( tptr[strlen(tptr)-1] == ' ' || tptr[strlen(tptr)-1] == '\t' )
				tptr[strlen(tptr)-1] = '\0';

			theSystem->SaveModel(tptr);
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_PSEUDO ) 
		{
			pseudo_add = atof( theCommands->allDirectives[c]->raw_observable_string[0] );
			printf("Setting pseudo add to %lf.\n", pseudo_add );
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_START ) 
		{
			start_at = atoi( theCommands->allDirectives[c]->raw_observable_string[0] );
			printf("Starting at a minimum of %d.\n", start_at );
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_STOP ) 
		{
			stop_at = atoi( theCommands->allDirectives[c]->raw_observable_string[0] );
			printf("Stopping at a maximum of %d.\n", stop_at  );
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_RANDOMT ) 
		{
			printf("Randomizing.\n");
			do_randomize = 1;
		}
		else if( theCommands->allDirectives[c]->command == COMMAND_DECODE )
		{
			if( do_randomize )
				theSystem->randomizeTimeSequence( theCommands->allDirectives[c] );
	
				int maxl = 0;
	
				for( int d = taskid; d < theCommands->allDirectives[c]->n_data; d += nprocs )
				{
					int len = theCommands->allDirectives[c]->translated_len[d] - start_at;
					if( len > stop_at - start_at ) len = stop_at - start_at;	
					if( len > maxl )
						maxl = len;
				}		
	
				char *tstate = (char *)malloc( sizeof(char) * maxl * theCommands->allDirectives[c]->n_data );
				memset( tstate, 0, sizeof(char) * maxl * theCommands->allDirectives[c]->n_data );
	
	
				for( int d = 0; d < theCommands->allDirectives[c]->n_data; d ++ )
				{
					int s_off = start_at;

					int len = theCommands->allDirectives[c]->translated_len[d];
					if( len < s_off ) continue;


					len -= s_off;					
					
					if( len > stop_at - start_at ) len = stop_at - start_at;	
		
					char *wild_string = (char *)malloc( sizeof(char) * (len+1) );
					for( int s = 0; s < len; s++ )
						wild_string[s] = '.';
					wild_string[len] = '\0';
					double *alpha_constrained = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha = (double *)malloc( sizeof(double) * len * nstates );
					double *beta  = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
					
					double *alpha_norm = (double *)malloc( sizeof(double) * len * nstates );
					double *beta_norm  = (double *)malloc( sizeof(double) * len * nstates );
	
	
					calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
					calculateBetaForStringConstrained( beta_constrained, beta_constrained_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off,  theCommands->allDirectives[c]->raw_class_string[d],  len, *theSystem );
					
//					calculateAlphaForStringConstrained( alpha, alpha_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
//					calculateBetaForStringConstrained( beta, beta_norm, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, wild_string,  len, *theSystem );
	
					int *predicted_state = (int *)malloc( sizeof(int) * len );
					
					theSystem->estimateClassViterbi( alpha_constrained, beta_constrained, theCommands->allDirectives[c]->translated_observable_string[d]+s_off, theCommands->allDirectives[c]->raw_class_string[d], len, predicted_state  );
	
					for( int t = 0; t < len; t++ )
					{
						int o = theCommands->allDirectives[c]->translated_observable_string[d][s_off+t];
						int nok[3];
						codeToCounts(o,nok,ndo);
						int nmol1 =nok[0];
						int nmol2 =nok[1];
						int nmol3 =nok[2];
	
						component_values[predicted_state[t]*3+0] += nmol1;
						component_values[predicted_state[t]*3+1] += nmol2;
						component_values[predicted_state[t]*3+2] += nmol3;
					}
	
	
					for( int t = 0; t < len; t++ )
						tstate[t*(theCommands->allDirectives[c]->n_data)+d] = (char)predicted_state[t];
	
					free(wild_string);
					free(alpha);
					free(beta);
					free(alpha_constrained);
					free(beta_constrained);
					free(alpha_norm);
					free(beta_norm);
					free(alpha_constrained_norm);
					free(beta_constrained_norm);
				}

			int do_switch = 0;

//			if( component_values[1]  < component_values[4] )
//				do_switch = 1;

			for( int t = 0; t < maxl; t++ )
			{
				printf("trajgrep ");
				for( int s = 0; s < theCommands->allDirectives[c]->n_data; s++ )
					printf("%d", (do_switch ? 1 - tstate[t*theCommands->allDirectives[c]->n_data+s] : tstate[t*theCommands->allDirectives[c]->n_data+s] ) );
				printf("\n");
			} 
	
			theSystem->ReportEmissionGroups( theCommands->allDirectives[0]->interpreter );
		}
	}
	printf("last_val: %lf ", last_val );
	for( int s = 0; s < nstates; s++ )
	{
		double tsum = component_values[s*3+0] + component_values[s*3+1] + component_values[s*3+2];
		printf("State %d has %lf%% MOL1 %lf%% MOL2 %lf%% MOL3 ",
			s, 100*component_values[s*3+0]/tsum, 100*component_values[s*3+1]/tsum, 100*component_values[s*3+2]/tsum );
	}
	printf("\n");

#ifdef PARALLEL	
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	fflush(stdout);

#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	return 0;
}

