#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "HMM.h"
#define CUT_PRINT 4700
extern int globalN_OBSERVABLES;
//#define NANCHECK
//#define CHECK_IMPROBABLE
int my_isnan( double val );

int deactivate_N = 1;

char recode( int i , int j )
{
	return 'A' + i;
}

char class_match( char src, char target )
{
	if( src == ' ' || target == ' ' || src == '.' || target == '.' ) return 1;
	return src == target;
}

void calculateAlphaForStringConstrained( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	int do_SC = 1;

	if( len != 40 ) do_SC = 0;

	int   nstates = theSystem.getNStates();
	double *p_obs = theSystem.getObservableP();
	double *t_ab  = theSystem.getTransitionP();
	double *init_p = theSystem.getInitP();

	// list of observables.	

	// alpha[0] 
/*
	if( do_SC )
	{
		double n_SC = 0;

		for( int s = 0; s < nstates; s++ )
		{
		if( !strcasecmp( theSystem.states[s]->label + strlen(theSystem.states[s]->label)-3, "_SC") )
//		if( !strcasecmp( theSystem.states[s]->label, "InnerStartState_SC") && class_match( classes[0], theSystem.states[s]->theClass) )
			n_SC++;
		}
		for( int s = 0; s < nstates; s++ )
		{
//		if( !strcasecmp( theSystem.states[s]->label, "InnerStartState_SC") && class_match( classes[0], theSystem.states[s]->theClass) )
		if( !strcasecmp( theSystem.states[s]->label + strlen(theSystem.states[s]->label)-3, "_SC" ) )
			alpha[s] = p_obs[observables[0] * nstates + s] * 1.0 / n_SC;
		else
			alpha[s] = 0;
		}
	}
	else
	{
		for( int s = 0; s < nstates; s++ )
		{
			if( !class_match( classes[0], theSystem.states[s]->theClass)  ) // classes[0] != '.' && theSystem.states[s]->theClass != classes[0] )
				alpha[s] = 0;
			else
			{
				alpha[s] = init_p[s] * p_obs[ observables[0] * nstates + s];
			}
		}
		
	}
*/
	int do_err = 0;
	for( int s = 0; s < nstates; s++ )
	{
		/*if( strcasecmp( theSystem.states[s]->label + strlen(theSystem.states[s]->label)-3, "_SC" ) )
		{
			alpha[s] = 0;
			continue;
		}*/
			

		if( !class_match( classes[0], theSystem.states[s]->theClass)  ) // classes[0] != '.' && theSystem.states[s]->theClass != classes[0] )
			alpha[s] = 0;
		else
		{
			alpha[s] = init_p[s] * p_obs[ observables[0] * nstates + s];
		}
	}

	double ptot= 0;

	for( int x = 0; x < nstates; x++ )
		ptot += alpha[x];
	for( int x = 0; x < nstates; x++ )
	{
		alpha[x] /= ptot;	
	//	if( fabs(alpha[x]) > 0.00001 )
	//		printf("State %s is occupied at the outset.\n", theSystem.states[x]->label );
	}
	alpha_norm[0] = log(ptot);

	FILE *err_file = NULL;
	char err_file_name[256];
	sprintf(err_file_name, "err.%d.file", taskid );

	for( int t = 1; t < len; t++ )
	{
		double prev_sum = 0;
/*
		for( int s = 0; s < nstates; s++ )
		{
			prev_sum += alpha[(t-1)*nstates+s];

			if( t > CUT_PRINT && fabs(alpha[(t-1)*nstates+s]) > 1e-20 )
			{
				printf("allowed state %s at t=%d %le\n", theSystem.states[s]->label, t-1, alpha[(t-1)*nstates+s] );
				for( int c = 0; c < globalN_OBSERVABLES; c++ )
				{
					printf("emissionp %d %le\n", c, theSystem.p_obs[c*nstates+s] );
				}	
			}
		}*/

		for( int s = 0; s < nstates; s++ )
		{
			int state_list[MAX_TRANSITIONS];
			int ntrans;

			alpha[t*nstates+s] = 0;

			//theSystem.getTransitionsTo( s, state_list, &ntrans, classes[t-1] );
			theSystem.getTransitionsTo( s, state_list, &ntrans, '.' );
/*		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];
				if( fabs(t_ab[s2*nstates+s]) < 1e-10 )
				{
//					printf("Allowed transition with zero probability %s to %s.\n", theSystem.states[s2]->label, theSystem.states[s]->label );
				}
				else if( t > CUT_PRINT && alpha[(t-1)*nstates+s2] > 1e-10 && classes[t] != theSystem.states[s]->theClass && classes[t] != '.' )
				{
					printf("Allowed transition at time t=%d with probability %lf %s to %s and source probability %lf, but denied because required class '%c' does not match '%c'.\n", t, t_ab[s2*nstates+s], theSystem.states[s2]->label, theSystem.states[s]->label,
				 			alpha[(t-1)*nstates+s2], classes[t], theSystem.states[s]->theClass );
				}
	
			}
*/
			if( !class_match( classes[t], theSystem.states[s]->theClass) ) // classes[t] != '.' && classes[t] != theSystem.states[s]->theClass )
				continue;
			
	
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];
#ifdef CHECK_STUFF
				if( fabs(t_ab[s2*nstates+s]) < 1e-10 )
				{
			//		printf("Allowed transition with zero probability %s to %s.\n", theSystem.states[s2]->label, theSystem.states[s]->label );
				}
				else if( t > CUT_PRINT && alpha[(t-1)*nstates+s2] > 1e-10 )
				{
//					printf("Allowed transition with probability %le %s to %s and source probability %le!.\n", t_ab[s2*nstates+s], theSystem.states[s2]->label, theSystem.states[s]->label,
//				 			alpha[(t-1)*nstates+s2]);
				}
				else 
				{
			//		printf("Allowed transition with probability %lf %s to %s and no source probability %le.\n", t_ab[s2*nstates+s], theSystem.states[s2]->label, theSystem.states[s]->label,
			//	 			alpha[(t-1)*nstates+s2]);
				}
				
#endif
				alpha[t*nstates+s] += alpha[(t-1)*nstates+s2] * t_ab[s2*nstates+s] * p_obs[ observables[t] * nstates + s];
#ifdef NANCHECK
				if( my_isnan( alpha[(t-1)*nstates+s2] ) || my_isnan( t_ab[s2*nstates+s]) || my_isnan( p_obs[ observables[t] * nstates + s] ) )
				{
					err_file = fopen( err_file_name, "w");
					fprintf(err_file,"NAN error!!!!\n");
					fprintf(err_file,"classes: %s\n", classes );
					do_err = 1;
//					exit(1);
				}
#endif

			}
		}
		

			ptot =0;
			for( int x = 0; x < nstates; x++ )
			ptot += alpha[t*nstates+x];

#ifdef CHECK_IMPROBABLE
		if( fabs(ptot) < 1e-300 || do_err )
		{
			if( !err_file) err_file = fopen( err_file_name, "w");
			fprintf(err_file,"Fatal error, improbable alpha step [%le] at t=%d\n", ptot, t );
			int tmark = t;

			fprintf(err_file, "initp: ");
			for( int x = 0; x < nstates; x++ )
				fprintf(err_file, "%le ", theSystem.init_p[x] );
			fprintf(err_file,"\n");

			for( int t = tmark - 30; t <= tmark; t++ )
			{
				if( t < 1 ) continue;
	
			for( int s = 0; s < nstates; s++ )
			{
				if(  fabs(alpha[(t-1)*nstates+s]) > 1e-20 )
				{
					fprintf(err_file,"allowed state %s at t=%d %le\n", theSystem.states[s]->label, t-1, alpha[(t-1)*nstates+s] );
					for( int c = 0; c < globalN_OBSERVABLES; c++ )
					{
			//			printf("emissionp %d %le\n", c, theSystem.p_obs[c*nstates+s] );
					}	
				}
			}
		
			for( int s = 0; s < nstates; s++ )
			{
				int state_list[MAX_TRANSITIONS];
				int ntrans;
	
				alpha[t*nstates+s] = 0;
	
				//theSystem.getTransitionsTo( s, state_list, &ntrans, classes[t-1] );
				theSystem.getTransitionsTo( s, state_list, &ntrans, '.' );
			
				for( int xs = 0; xs < ntrans; xs++ )
				{
					int s2 = state_list[xs];
					if( fabs(t_ab[s2*nstates+s]) < 1e-10 )
					{
	//					printf("Allowed transition with zero probability %s to %s.\n", theSystem.states[s2]->label, theSystem.states[s]->label );
					}
					else if( alpha[(t-1)*nstates+s2] > 1e-10 && !class_match( classes[t], theSystem.states[s]->theClass ) ) //classes[t] != theSystem.states[s]->theClass && classes[t] != '.' )
					{
						fprintf(err_file, "Allowed transition at time t=%d with probability %lf %s to %s and source probability %lf, but denied because required class '%c' does not match '%c'.\n", t, t_ab[s2*nstates+s], theSystem.states[s2]->label, theSystem.states[s]->label,
					 			alpha[(t-1)*nstates+s2], classes[t], theSystem.states[s]->theClass );
					}
					else if( class_match( classes[t], theSystem.states[s]->theClass ) )//classes[t] == '.' || theSystem.states[s]->theClass == classes[t] ) //if( alpha[(t-1)*nstates+s2] > 1e-15 ) 
					{
						fprintf(err_file,"Allowed transition at time t=%d with probability %le [%d to %d] %s to %s and source probability %le, allowed because required class '%c' matches '%c'.\n", t, t_ab[s2*nstates+s], s2, s, theSystem.states[s2]->label, theSystem.states[s]->label,
					 			alpha[(t-1)*nstates+s2], classes[t], theSystem.states[s]->theClass );
					}
		
				}
	
				if( !class_match( classes[t], theSystem.states[s]->theClass ) ) // classes[t] != '.' && classes[t] != theSystem.states[s]->theClass )
					continue;
				
		
				for( int xs = 0; xs < ntrans; xs++ )
				{
					int s2 = state_list[xs];
	
					if( fabs(t_ab[s2*nstates+s]) < 1e-10 )
					{
				//		printf("Allowed transition with zero probability %s to %s.\n", theSystem.states[s2]->label, theSystem.states[s]->label );
					}
					else if( t > CUT_PRINT && alpha[(t-1)*nstates+s2] > 1e-10 )
					{
						fprintf(err_file,"Allowed transition with probability %le [%d to %d] %s to %s and source probability %le!.\n", t_ab[s2*nstates+s], s2, s, theSystem.states[s2]->label, theSystem.states[s]->label,
					 			alpha[(t-1)*nstates+s2]);
					}
					else 
					{
				//		printf("Allowed transition with probability %lf %s to %s and no source probability %le.\n", t_ab[s2*nstates+s], theSystem.states[s2]->label, theSystem.states[s]->label,
				//	 			alpha[(t-1)*nstates+s2]);
					}
					
	
					alpha[t*nstates+s] += alpha[(t-1)*nstates+s2] * 
						t_ab[s2*nstates+s] * p_obs[ observables[t] * nstates + s];
				}
			}

			ptot = 0;
			for( int x = 0; x < nstates; x++ )
				ptot += alpha[t*nstates+x];

		for( int x = 0; x < nstates; x++ )
			alpha[t*nstates+x] /= ptot;	

			}
			fflush(err_file);
			system("sleep 120");

			MPI_Finalize();
			exit(1);
		}
#endif
		for( int x = 0; x < nstates; x++ )
			alpha[t*nstates+x] /= ptot;	
		alpha_norm[t] = log(ptot) + alpha_norm[t-1];
	}

	// reweight alphas based on 
}


// beta(T), the probability of observing the 
	
void calculateBetaForStringConstrained( double *beta, double *beta_norm, int *observables, char *classes, int len, HMMSystem &theSystem )
{
	int   nstates = theSystem.getNStates();
	double *p_obs = theSystem.getObservableP();
	double *t_ab  = theSystem.getTransitionP();
	double *init_p = theSystem.getInitP();

	// list of observables.	

	double n_matching = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( !class_match( classes[len-1], theSystem.states[s]->theClass ) ) //classes[len-1] != '.' && theSystem.states[s]->theClass != classes[len-1] )
			{}
		else
			n_matching += 1;
	}
/*
	for( int s = 0; s < nstates; s++ )
	{
		beta[(len-1)*nstates+s] = 1.0 / nstates;
	}
*/

	for( int s = 0; s < nstates; s++ )
	{
		if( !class_match( classes[len-1], theSystem.states[s]->theClass ) ) //classes[len-1] != '.' && theSystem.states[s]->theClass != classes[len-1] )
			beta[(len-1)*nstates+s] = 0.0;
		else
			beta[(len-1)*nstates+s] = 1.0 / n_matching;
	}

	double ptot =0;

	for( int x = 0; x < nstates; x++ )
		ptot += beta[(len-1)*nstates+x];

	for( int x = 0; x < nstates; x++ )
		beta[(len-1)*nstates+x] /= ptot;	

	beta_norm[(len-1)] = log(ptot);
	
	for( int t = len-2; t >= 0; t-- )
	{
		for( int s = 0; s < nstates; s++ )
		{

			int state_list[MAX_TRANSITIONS];
			int ntrans;

			beta[t*nstates+s] = 0;

			theSystem.getTransitionsFrom( s, state_list, &ntrans, '.');
#ifdef CHECK_STUFF			
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];
				if( fabs(t_ab[s*nstates+s2]) < 1e-10 )
				{
//					printf("Allowed transition with zero probability %s to %s.\n", theSystem.states[s]->label, theSystem.states[s2]->label );
				}
				else if( beta[(t+1)*nstates+s2] > 1e-10 && classes[t] != theSystem.states[s]->theClass )
				{
//					printf("Allowed transition with probability %lf %s to %s and source probability %lf, but denied because required class '%c' does not match '%c'.\n", t_ab[s*nstates+s2], theSystem.states[s]->label, theSystem.states[s2]->label,
//				 			beta[(t+1)*nstates+s2], classes[t], theSystem.states[s]->theClass );
				}

			}
#endif	
			if( !class_match( classes[t], theSystem.states[s]->theClass) ) //classes[t] != '.' && classes[t] != theSystem.states[s]->theClass )
				continue;
	
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];

				beta[t*nstates+s] += beta[(t+1)*nstates+s2] * t_ab[s*nstates+s2] * p_obs[ observables[t+1] * nstates + s2];
#ifdef NANCHECK
				if( my_isnan( beta[(t+1)*nstates+s2] ) || my_isnan( t_ab[s*nstates+s2]) || my_isnan( p_obs[ observables[t+1] * nstates + s2] ) || beta[t*nstates+s] < 0)
				{
					printf("NAN error.\n");
					exit(1);
				}
#endif
			}
		}
		
		ptot =0;
	
		for( int x = 0; x < nstates; x++ )
			ptot += beta[t*nstates+x];

#ifdef CHECK_IMPROBABLE
		if( fabs(ptot) < 1e-100 )
		{
			printf("Improbable backward transition: %le at t=%d\n", ptot, t );
			exit(1);
		}
#endif
	
		for( int x = 0; x < nstates; x++ )
			beta[t*nstates+x] /= ptot;	
	
		beta_norm[t] = log(ptot) + beta_norm[t+1];
	}

}


void estimateClass( double *alphas, double *betas, int *observables, char *theClass, int len, HMMSystem &theSystem )
{
/*
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	int nstates = theSystem.getNStates();
	
	double *p_obs = theSystem.getObservableP();
	double class_p[N_CLASSES*len];
	memset( class_p, 0, sizeof(double) * N_CLASSES * len );


	if( !theClass )
	{
		printf("SEQ ");
		for( int c = 0; c < N_CLASSES; c++ )
			printf("%c ", TM_CLASSES[c] );
	}

	for( int t = 0; t < len; t++ )
	{
		for( int s = 0; s < nstates; s++ )
		{
			int state_class = interpretTMClass( theSystem.states[s]->theClass );

			double p = alphas[t*nstates+s] * betas[t*nstates+s] * p_obs[observables[t]*nstates+s];

			if( t == len-1)
				p = alphas[t*nstates+s] * p_obs[observables[t]*nstates+s];

			class_p[t*N_CLASSES+state_class] += p;	 
		}

		double ptot = 0;
		for( int c = 0; c < N_CLASSES; c++ )
			ptot += class_p[t*N_CLASSES+c];

		double large_p = 0;
		char likely_class = ' ';

		for( int c = 0; c < N_CLASSES; c++ )
		{
			class_p[t*N_CLASSES+c] /= ptot;
			if( class_p[t*N_CLASSES+c] > large_p )
			{
				large_p = class_p[t*N_CLASSES+c];
				likely_class = TM_CLASSES[c];
			}
		}

		if( theClass && theClass[t] == '.' )
			theClass[t] = likely_class;

		if( !theClass )
		{
			printf("%c ", all_res[observables[t]] );
			for( int c = 0; c < N_CLASSES; c++ )
				printf("%lf ", class_p[t*N_CLASSES+c] );
			printf("\n");
		}
		//%lf %lf %lf %lf\n", all_res[observables[t]], class_p[t*N_CLASSES+0], class_p[t*N_CLASSES+1], class_p[t*N_CLASSES+2], class_p[t*N_CLASSES+3] );
	}
	
	if( !theClass )
 		return;

	// post-process impossibilities.
	int done = 0;
*/
}


double HMMSystem::estimateClass1Best( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
/*
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	int * one_best_hypothesis = (int *)malloc( sizeof(int) * len * len );
	double *cur_p = (double *)malloc( sizeof(double) * len );
	
	double class_p[N_CLASSES*len];
	memset( class_p, 0, sizeof(double) * N_CLASSES * len );

	for( int x = 0; x < len; x++ )
		cur_p[x] = 1.0;


	for( int t2 = 0; t2 < len; t2++ )
	{
		int best_hypothesis = 0;
		double best_local_p = 0;
		int best_local_state = 0;

		if( t2 == 0 )
		{
			// we have not generated any hypotheses yet...
			
			double best_p = 0;
			int best_state;

			for( int s = 0; s < nstates; s++ )
			{
				double p = init_p[s] * p_obs[observables[0]*nstates+s] * betas[s];// * p_obs[observables[0]*nstates+s];
				
				if( p > best_p )
				{
					best_state = s;
			
					one_best_hypothesis[0] = s;

					cur_p[0] = log(init_p[s] * p_obs[observables[0]*nstates+s]);

					best_p = p;
				}
			}
		}


		for( int hyp = 0; hyp < t2; hyp++ )
		{
			int trans_list[MAX_TRANSITIONS];
			int ntrans;
			int s = one_best_hypothesis[hyp*len+(t2-1)];

			getTransitionsFrom( one_best_hypothesis[hyp*len+(t2-1)], trans_list, &ntrans, '.' );
		
			double best_p = -1e300;		
			double best_check = 1.0;
			int best_state = 0;
			int s_1, s_2;
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = trans_list[xs];
				double p = cur_p[hyp] + log( betas[t2*nstates+s2] * p_obs[observables[t2]*nstates+s2] * t_ab[s*nstates+s2] );
				double check_it_p = betas[t2*nstates+s2] * p_obs[observables[t2]*nstates+s2] * t_ab[s*nstates+s2];
				double local_p = alphas[t2*nstates+s2] * betas[t2*nstates+s2] * p_obs[observables[t2]*nstates+s2];

				if( hyp == 0 && xs == 0 )
				{
					best_hypothesis = hyp;
					memcpy( one_best_hypothesis + t2 * len, one_best_hypothesis + best_hypothesis * len, sizeof(int) * len );
					one_best_hypothesis[t2*len+t2] = s2;
//					cur_p[t2] = p;
					cur_p[t2] = cur_p[hyp] + log( p_obs[observables[t2]*nstates+s2] * t_ab[s*nstates+s2] );
				}


				
				if( local_p > best_local_p )
				{
					// this is the hypothesis for this state.
					best_hypothesis = hyp;
					memcpy( one_best_hypothesis + t2 * len, one_best_hypothesis + best_hypothesis * len, sizeof(int) * len );
					one_best_hypothesis[t2*len+t2] = s2;
					cur_p[t2] = cur_p[hyp] + log( p_obs[observables[t2]*nstates+s2] * t_ab[s*nstates+s2] );
					best_local_p = local_p;
				}

				if( p > best_p )
				{	
					s_1 = s;
					s_2 = s2;
					// this is the best state for this hypothesis.
					one_best_hypothesis[hyp*len+t2] = s2;
					best_p = p;
					best_check = check_it_p;
				}
			}
					
			cur_p[hyp] = cur_p[hyp] + log( p_obs[observables[t2]*nstates+s_2] * t_ab[s_1*nstates+s_2] ); // best_p;

		}
	}

	double best_prob = -1e300;
	int best_hyp = 0;

	for( int t = 0; t < len-2; t++ )
	{
		if( cur_p[t] > best_prob )
		{
			best_prob = cur_p[t];
			best_hyp = t;			
		}
	}
	
	if( !theClass && taskid == 0 )
	{
		printf("SEQ ");
		for( int c = 0; c < N_CLASSES; c++ )
			printf("%c ", TM_CLASSES[c] );
	}

	for( int t = 0; t < len; t++ )
	{
//		printf("%d %s\n", t, states[one_best_hypothesis[best_hyp*len+t]]->label );
		if( ML_state )
			ML_state[t] = one_best_hypothesis[best_hyp*len+t];

		if( !theClass )
		{
			if( taskid == 0 )
			printf("%c %c\n", all_res[observables[t]], states[one_best_hypothesis[best_hyp*len+t]]->theClass );
		}
		else
			theClass[t] = states[one_best_hypothesis[best_hyp*len+t]]->theClass;

//		states[one_best_hypothesis[best_hyp*len+t]]->obs_population[observables[t]] += 1;
	}	


	char ctm3 = theClass[len-3];

	if( ctm3 != 'M' )
	{
		for( int t = len-2; t < len; t++ )
			theClass[t] = ctm3;		
	}

	free( one_best_hypothesis );
	free( cur_p );
	

	if( !theClass )
 		return best_prob;


	return best_prob;

*/
	return 1.0;
}

double grn( void )
{
	double v = ((double)rand() ) / ((double)RAND_MAX) - 1e-20;
	
	return v;
}

void HMMSystem::generateSequence( int *observables, char *classes, int len )
{
	// list of observables.	

	// alpha[0] 

	int *ostates = (int *)malloc( sizeof(int) * len );

	double r = grn();

	double cum = 0;
	
	ostates[0] = nstates-1;

	for( int s = 0; s < nstates; s++ )
	{
		cum += init_p[s];
		
		if( r < cum )
		{
			ostates[0] = s;
			break;
		}
	}

	for( int t = 1; t < len; t++ )
	{
		double prev_sum = 0;

		r = grn();
		cum = 0;

		int state_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( ostates[t-1], state_list, &ntrans, '.' );

		ostates[t] = state_list[ntrans-1];

		for( int xs = 0; xs < ntrans; xs++ )
		{
			int s2 = state_list[xs];

			cum += t_ab[ostates[t-1]*nstates+s2];

			if( r < cum )
			{
				ostates[t] = s2;
				break;
			}
		}
	}

	for( int t = 0; t < len; t++ )
	{
		cum = 0;
		r = grn();

		observables[t] = globalN_OBSERVABLES-1;
		classes[t] = states[ostates[t]]->theClass;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			cum += p_obs[o*nstates+ostates[t]];

			if( r < cum )
			{
				observables[t] = o;
				break;
			}	
		}	
	}
}

/*
int HMMSystem::estimateClassViterbi( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.
	
	int do_SC = (len == 40 );

		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[0], states[s]->theClass) )
			{
					cur_prob[s] = 0;
			}
			else
				cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
		}	


	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	

		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
	}
	for( int t = 0; t < len; t++ )
	{
		if( ML_state )
			ML_state[t] = the_path[t];

		if( !theClass )
		{
			if( taskid == 0 )
				printf("%c %c\n", 'A' + observables[t], states[the_path[t]]->theClass );
				//printf("%c %c\n", all_res[observables[t]], states[the_path[t]]->theClass );
		}
		else
		{
			if( !class_match( theClass[t], states[the_path[t]]->theClass ) )
			{
				printf("Fatal error in class match!!!!\n");
			//	exit(1);
				free(back_ptr);
				return 0;
			}
			theClass[t] = states[the_path[t]]->theClass;
		}
		states[the_path[t]]->obs_population[observables[t]] += 1;
	}	

	free(back_ptr);

	return 1;
}
*/

int HMMSystem::estimateClassViterbi( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	double *T1 = (double *)malloc( sizeof(double) * nstates * len ); 
	int    *T2 = (int *)malloc( sizeof(int) * nstates * len ); 

	for( int s = 0; s < nstates; s++ )
	{
		T1[s] = init_p[s] * p_obs[observables[0]*nstates+s];
		T2[s] = 0;
	}

	for( int t = 1; t < len; t++ )
	{
		for( int s = 0; s < nstates; s++ )
		{
			double maxp = 0;

			for( int k = 0; k < nstates; k++ )
			{
				double valp = T1[k+nstates*(t-1)] * t_ab[k*nstates+s] * p_obs[observables[t]*nstates+s];

				if( valp > maxp )
				{
					maxp = valp;
					T1[s+nstates*t] = valp;
					T2[s+nstates*t] = k;
				}
			}			
		}
		double norm = 0;
		for( int s = 0; s < nstates; s++ )
			norm += T1[s+nstates*t];
		for( int s = 0; s < nstates; s++ )
			T1[s+nstates*t] /= norm;
	}
		
	double zt = 0;
	int xt = 0;
	for( int k = 0; k < nstates; k++ )
	{
		if( T1[k+nstates*(len-1)] > zt )
		{
			zt = T1[k+nstates*(len-1)];
			xt = k;
//			xt = T2[k+nstates*(len-1)];
		}
	}
		
	ML_state[len-1] = xt;

	for( int t = len-2; t >= 0; t-- )
	{
		ML_state[t] = T2[xt+nstates*(t+1)];
		xt = ML_state[t];	
	}
/*
	double zt = 0;
	int xt = 0;

	for( int k = 0; k < nstates; k++ )
	{
		if( T1[k+nstates*(len-1)] > zt )
		{
			zt = T1[k+nstates*(len-1)];
			xt = k;
		}
	}

	printf("T1 end: %lf and %lf\n",
		T1[0+nstates*(len-1)], T1[1+nstates*(len-1)] );
	for( int t = len-1; t >= 0; t-- )
	{
		double zt = 0;
		for( int k = 0; k < nstates; k++ )
		{
			if( T1[k+nstates*(t)] > zt )
			{
				zt = T1[k+nstates*(t)];
				xt = k;
			}
		}		
		ML_state[t] = xt;
	}
*/

	free(T1);
	free(T2);

	return 1;
}


int HMMSystem::estimateClassPlayground( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.

	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
		if( t_ab[the_path[t]*nstates+the_path[t+1]] < 1e-20 )
		{
			printf("Fatal error, improbable transition chosen.\n");
//			free(back_ptr);
			return 0;
		} 
	}


/* REDO!!!*/
/*
	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		char recode( int, int );
		char o = recode( observables[t], 0 );

		printf("At t=%d, with observable %c, the three most likely states are:\n", t, o);
		for( int x = 0; x < 3; x++ )
		{
			int s1 = back_ptr[t*nstates+sorter[x]];
			int s2 = sorter[x];

			double p_o = p_obs[observables[t]*nstates+s2];
			double p_t = t_ab[s1*nstates+s2];

			
			printf("State '%s' with probability %le, most likely from state '%s' with p_obs %le p_trans %le.\n",
				states[sorter[x]]->label, new_probs[sorter[x]], states[back_ptr[t*nstates+sorter[x]]]->label, p_o, p_t );
		}
		
		{
			int s1 = the_path[t-1];
			int s2 = the_path[t];

			double p_o = p_obs[observables[t]*nstates+s2];
			double p_t = t_ab[s1*nstates+s2];
			printf("The state we're actually going to use: '%s', with probability %le, most likely from state '%s' with p_obs %le p_trans %le.\n",
				states[s2]->label, new_probs[s2], states[s1]->label, p_o, p_t );
		}

		printf("greppath %d %le %le\n", new_probs[sorter[0]], new_probs[the_path[t]] );

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}
*/

	// Compare the most likely membrane path to the most likely globular path, find the divergence point.

	int globular_path[len];
	int membrane_path[len];

	int g_state = -1;
	int m_state = -1;

	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->theClass == 'O' )
		{
			g_state = s;	
			printf("State selected as globular state: '%s'.\n", states[g_state]->label );
			break;
		}
	}
	
	double bp = 0;
	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > bp && s != g_state )
		{
			m_state = s;
			bp = cur_prob[s];
		}
	}
	
	globular_path[len-1] = g_state;
	membrane_path[len-1] = m_state;

	for( int t = len-2; t >= 0; t-- )
	{
		globular_path[t] = back_ptr[(t+1)*nstates+globular_path[t+1]];
		membrane_path[t] = back_ptr[(t+1)*nstates+membrane_path[t+1]];
	}

	double best_increment = 1;
	double best_decrement = 1;
	int p_dec = 0;
	int p_inc = 0;

	for( int t = 1; t < len; t++ )
	{
		double p_m_o = p_obs[observables[t]*nstates+membrane_path[t]];
		double p_m_t = t_ab[membrane_path[t-1]*nstates+membrane_path[t]];
		
		double p_g_o = p_obs[observables[t]*nstates+globular_path[t]];
		double p_g_t = t_ab[globular_path[t-1]*nstates+globular_path[t]];

		double fac = (p_g_o * p_g_t) / (p_m_o * p_m_t);

		if( fac > best_increment )
		{
			best_increment = fac;
			p_inc = t;
		}
		if( fac < best_decrement )
		{
			best_decrement = fac;
			p_dec = t;
		}
	}
		
	char recode( int, int );

	char o1 = recode( observables[p_dec-1], 0 ); 
	char o2 = recode( observables[p_dec], 0 ); 

	printf("Best de-enhancement: %le\n", best_decrement);
	printf("OBS1: %c OBS2 %c\n", o1, o2);
	printf("M: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le.\n", states[membrane_path[p_dec-1]]->label, states[membrane_path[p_dec]]->label,
		p_obs[observables[p_dec-1]*nstates+membrane_path[p_dec-1]],
		p_obs[observables[p_dec]*nstates+membrane_path[p_dec]],
		t_ab[membrane_path[p_dec-1]*nstates+membrane_path[p_dec]]
		);
	printf("G: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le\n", states[globular_path[p_dec-1]]->label, states[globular_path[p_dec]]->label,
		p_obs[observables[p_dec-1]*nstates+globular_path[p_dec-1]],
		p_obs[observables[p_dec]*nstates+globular_path[p_dec]],
		t_ab[globular_path[p_dec-1]*nstates+globular_path[p_dec]] );
	o1 = recode( observables[p_inc-1], 0 ); 
	o2 = recode( observables[p_inc], 0 ); 
	
	printf("Best enhancement: %le\n", best_increment);
	printf("OBS1: %c OBS2 %c\n", o1, o2);
	printf("M: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le.\n", states[membrane_path[p_inc-1]]->label, states[membrane_path[p_inc]]->label, 
		p_obs[observables[p_inc-1]*nstates+membrane_path[p_inc-1]],
		p_obs[observables[p_inc]*nstates+membrane_path[p_inc]],
		t_ab[membrane_path[p_inc-1]*nstates+membrane_path[p_inc]]
		);
	printf("G: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le\n", states[globular_path[p_inc-1]]->label, states[globular_path[p_inc]]->label,
		p_obs[observables[p_inc-1]*nstates+globular_path[p_inc-1]],
		p_obs[observables[p_inc]*nstates+globular_path[p_inc]],
		t_ab[globular_path[p_inc-1]*nstates+globular_path[p_inc]]
		);
	if( !strncasecmp( states[membrane_path[p_inc]]->label, "CoreState", 9 ) )
	{
		int spot;
		
		sscanf( states[membrane_path[p_inc]]->label, "CoreState_%d", &spot ); 

		printf("grepspot %d %c\n", spot, o2 );
	}



















	for( int t = 0; t < len; t++ )
	{
		if( ML_state )
			ML_state[t] = the_path[t];

		if( !theClass )
		{
			if( taskid == 0 )
				printf("%c %c\n", 'A' + observables[t], states[the_path[t]]->theClass );
//				printf("%c %c\n", all_res[observables[t]], states[the_path[t]]->theClass );
		}
		else
		{
			if( !class_match( theClass[t], states[the_path[t]]->theClass ) )
			{
				printf("Fatal error in class match!!!!\n");
			//	exit(1);
				free(back_ptr);
				return 0;
			}
			theClass[t] = states[the_path[t]]->theClass;
		}
		states[the_path[t]]->obs_population[observables[t]] += 1;
	}	

	printf("The following path was used:\n");
	for( int t = 1; t < len; t++ )
	{
		int s1 = the_path[t-1];
		int s2 = the_path[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];

		printf("t=%d state: '%s' p_transition: %le\n", t, states[the_path[t]]->label, p );
			
	}
	
//	printf("Done.\n");
//	printf("Type return to continue.\n");
//	char tc;
//	scanf("%c", &tc);

	free(back_ptr);

	return 1;
}

int HMMSystem::estimateClassPlayground2( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.

	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
		if( t_ab[the_path[t]*nstates+the_path[t+1]] < 1e-20 )
		{
			printf("Fatal error, improbable transition chosen.\n");
//			free(back_ptr);
			return 0;
		} 
	}


/* REDO!!!*/
/*
	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		char recode( int, int );
		char o = recode( observables[t], 0 );

		printf("At t=%d, with observable %c, the three most likely states are:\n", t, o);
		for( int x = 0; x < 3; x++ )
		{
			int s1 = back_ptr[t*nstates+sorter[x]];
			int s2 = sorter[x];

			double p_o = p_obs[observables[t]*nstates+s2];
			double p_t = t_ab[s1*nstates+s2];

			
			printf("State '%s' with probability %le, most likely from state '%s' with p_obs %le p_trans %le.\n",
				states[sorter[x]]->label, new_probs[sorter[x]], states[back_ptr[t*nstates+sorter[x]]]->label, p_o, p_t );
		}
		
		{
			int s1 = the_path[t-1];
			int s2 = the_path[t];

			double p_o = p_obs[observables[t]*nstates+s2];
			double p_t = t_ab[s1*nstates+s2];
			printf("The state we're actually going to use: '%s', with probability %le, most likely from state '%s' with p_obs %le p_trans %le.\n",
				states[s2]->label, new_probs[s2], states[s1]->label, p_o, p_t );
		}

		printf("greppath %d %le %le\n", new_probs[sorter[0]], new_probs[the_path[t]] );

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}
*/

	// Compare the most likely membrane path to the most likely globular path, find the divergence point.


	int m1_state = -1;
	double m1_prob = 0;
	int m2_state = -1;
	double m2_prob = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( !strncasecmp( states[s]->label+strlen(states[s]->label)-2, "_0", 2 ) )
		{
			printf("Considering state '%s' [%le].\n", states[s]->label, cur_prob[s] );
			if( cur_prob[s] > m1_prob )
			{
				m1_prob = cur_prob[s];
				m1_state = s;
			}
		}
		if( !strncasecmp( states[s]->label+strlen(states[s]->label)-2, "_1", 2 ) )
		{
			if( cur_prob[s] > m2_prob )
			{
				m2_prob = cur_prob[s];
				m2_state = s;
			}
		}
	}
	
	int membrane_path2[len];
	int membrane_path1[len];
	
	membrane_path1[len-1] = m1_state;
	membrane_path2[len-1] = m2_state;

	for( int t = len-2; t >= 0; t-- )
	{
		membrane_path1[t] = back_ptr[(t+1)*nstates+membrane_path1[t+1]];
		membrane_path2[t] = back_ptr[(t+1)*nstates+membrane_path2[t+1]];
	}

	double best_increment = 1;
	double best_decrement = 1;
	int p_dec = 0;
	int p_inc = 0;

	for( int t = 1; t < len; t++ )
	{
		double p_m_o = p_obs[observables[t]*nstates+membrane_path1[t]];
		double p_m_t = t_ab[membrane_path1[t-1]*nstates+membrane_path1[t]];
		
		double p_g_o = p_obs[observables[t]*nstates+membrane_path2[t]];
		double p_g_t = t_ab[membrane_path2[t-1]*nstates+membrane_path2[t]];

		double fac = (p_g_o * p_g_t) / (p_m_o * p_m_t);

		if( fac > best_increment )
		{
			best_increment = fac;
			p_inc = t;
		}
		if( fac < best_decrement )
		{
			best_decrement = fac;
			p_dec = t;
		}
	}
		
	char recode( int, int );

	char o1 = recode( observables[p_dec-1], 0 ); 
	char o2 = recode( observables[p_dec], 0 ); 
/*
	printf("Best de-enhancement: %le\n", best_decrement);
	printf("OBS1: %c OBS2 %c\n", o1, o2);
	printf("M1: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le.\n", states[membrane_path1[p_dec-1]]->label, states[membrane_path1[p_dec]]->label,
		p_obs[observables[p_dec-1]*nstates+membrane_path1[p_dec-1]],
		p_obs[observables[p_dec]*nstates+membrane_path1[p_dec]],
		t_ab[membrane_path1[p_dec-1]*nstates+membrane_path1[p_dec]]
		);
	printf("M2: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le\n", states[membrane_path2[p_dec-1]]->label, states[membrane_path2[p_dec]]->label,
		p_obs[observables[p_dec-1]*nstates+membrane_path2[p_dec-1]],
		p_obs[observables[p_dec]*nstates+membrane_path2[p_dec]],
		t_ab[membrane_path2[p_dec-1]*nstates+membrane_path2[p_dec]] );
	o1 = recode( observables[p_inc-1], 0 ); 
	o2 = recode( observables[p_inc], 0 ); 
	
	printf("Best enhancement: %le\n", best_increment);
	printf("OBS1: %c OBS2 %c\n", o1, o2);
	printf("M1: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le.\n", states[membrane_path1[p_inc-1]]->label, states[membrane_path1[p_inc]]->label, 
		p_obs[observables[p_inc-1]*nstates+membrane_path1[p_inc-1]],
		p_obs[observables[p_inc]*nstates+membrane_path1[p_inc]],
		t_ab[membrane_path1[p_inc-1]*nstates+membrane_path1[p_inc]]
		);
	printf("M2: state %s to %s, p_obs1: %le p_obs2: %le p_trans: %le\n", states[membrane_path2[p_inc-1]]->label, states[membrane_path2[p_inc]]->label,
		p_obs[observables[p_inc-1]*nstates+membrane_path2[p_inc-1]],
		p_obs[observables[p_inc]*nstates+membrane_path2[p_inc]],
		t_ab[membrane_path2[p_inc-1]*nstates+membrane_path2[p_inc]]
		);
*/




	printf("Membrane path 1: (%le)\n", cur_prob[m1_state] );
	for( int t = 1; t < len; t++ )
	{
		int s1 = membrane_path1[t-1];
		int s2 = membrane_path1[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];

		printf("t=%d state: '%s' p_transition: %le\n", t, states[membrane_path1[t]]->label, p );
			
	}
	
	printf("Membrane path 2: (%le)\n", cur_prob[m2_state] );
	for( int t = 1; t < len; t++ )
	{
		int s1 = membrane_path2[t-1];
		int s2 = membrane_path2[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];

		printf("t=%d state: '%s' p_transition: %le\n", t, states[membrane_path2[t]]->label, p );
			
	}
	printf("The following path was used:\n");
	for( int t = 1; t < len; t++ )
	{
		int s1 = the_path[t-1];
		int s2 = the_path[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];

		printf("t=%d state: '%s' p_transition: %le\n", t, states[the_path[t]]->label, p );
			
	}
	
//	printf("Done.\n");
//	printf("Type return to continue.\n");
//	char tc;
//	scanf("%c", &tc);

	free(back_ptr);

	return 1;
}

int HMMSystem::estimateClassThreeStater( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	double tot_prob[4] = {0,0,0,0};
	const char *path_names[4] = { "MMOO", "OMMO", "OOMM", "MOOM" };

	/*
		There are only four possible states
		
		MMOO
		OMMO
		OOMM
		MOOM

		possible paths:
		
		MMOO:

		1 2 5 0
		3 4 5 0

		OMMO:

		0 1 2 5
		0 3 4 5

		OOMM:

		5 0 1 2
		5 0 3 4

		MOOM:

		2 5 0 1
		2 5 0 3
		4 5 0 1
		4 5 0 3
	*/

	int paths[10][5] = 
		{
			{0, 1, 2, 5, 0 },
			{0, 3, 4, 5, 0 },
			{1, 0, 1, 2, 5 },
			{1, 0, 3, 4, 5 },
			{2, 5, 0, 1, 2 },
			{2, 5, 0, 3, 4 },
			{3, 2, 5, 0, 1 },
			{3, 2, 5, 0, 3 },
			{3, 4, 5, 0, 1 },
			{3, 4, 5, 0, 3 }
		};

	for( int p = 0; p < 10; p++ )
	{
		double path_prob = init_p[paths[p][1]] * p_obs[observables[0]*nstates+paths[p][1]];

		for( int t = 1; t < len; t++ )
			path_prob *= t_ab[paths[p][1+t-1]*nstates+paths[p][1+t]] * p_obs[observables[t]*nstates+paths[p][1+t]];

		tot_prob[paths[p][0]] += path_prob;	
	}

	double sum = 0;
	for( int x = 0; x < 4; x++ )
	{
		sum += tot_prob[x];
	}
	
	for( int x = 0; x < 4; x++ )	
	printf("Path %s prob: %le\n", path_names[x], tot_prob[x]/sum );
	
	for( int p = 0; p < 10; p++ )
	{
		double path_prob = init_p[paths[p][1]] * p_obs[observables[0]*nstates+paths[p][1]];

		for( int t = 1; t < len; t++ )
			path_prob *= t_ab[paths[p][1+t-1]*nstates+paths[p][1+t]] * p_obs[observables[t]*nstates+paths[p][1+t]];

		printf("Path %d%d%d%d has p: %le\n",
			paths[p][1],
			paths[p][2],
			paths[p][3],
			paths[p][4], path_prob / sum );
	}

}



double HMMSystem::rankBasedOnStatePopulations( 
		double *alpha,
		double *beta,
		int len,

		char *state1,
		 char *state2,
		 char *state3,
		 char *state4,
		 char *state5,
		 char *state6,
		 char *state7,
		 char *state8,
		 char *state9,
		 char *state10 )
{
	const char *m_states[10] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
	
	m_states[0] = state1;
	int n_m_states = 1;

	if( state10 )
		m_states[9] = state10;
	if( state9 )
		m_states[8] = state9;
	if( state8 )
		m_states[7] = state8;
	if( state7 )
		m_states[6] = state7;
	if( state6 )
		m_states[5] = state6;
	if( state5 )
		m_states[4] = state5;
	if( state4 )
		m_states[3] = state4;
	if( state3 )
		m_states[2] = state3;
	if( state2 )
		m_states[1] = state2;
	
	for( int x = 1; x < 10; x++ )
	{
		if( !m_states[x] ) break;
		n_m_states++;
	}
		
	double spt = 0;

	for( int t = 0; t < len; t++ )
	{
		for( int s = 0; s < nstates; s++ )
		{
			int s_x = -1;
		
			for( int m = 0; m < n_m_states; m++ )
			{
				if( !strcasecmp( m_states[m], states[s]->label ) )
					s_x = m;
			}

			if( s_x == -1 )
				continue;

			spt += alpha[t*nstates+s] * beta[t*nstates+s];		
		}
	}	
			
	return spt;

}

double HMMSystem::ratio12Viterbi( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *probs)
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.

	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	int do_print = 0;

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );

		if( do_print )
		{
			for( int s = 0; s < nstates; s++ )
				if( cur_prob[s] > 1e-10 )
					printf("Allowed at t=%d state %s\n", t, states[s]->label );
		}
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
		if( t_ab[the_path[t]*nstates+the_path[t+1]] < 1e-20 )
		{
			printf("Fatal error, improbable transition chosen.\n");
//			free(back_ptr);
			return 0;
		} 
	}
	for( int t = 0; t < len; t++ )
	{
		if( ML_state )
			ML_state[t] = the_path[t];

		if( !theClass )
		{
			if( taskid == 0 )
				printf("%c %c\n", 'A' + observables[t], states[the_path[t]]->theClass );
//				printf("%c %c\n", all_res[observables[t]], states[the_path[t]]->theClass );
		}
		else
		{
			if( !class_match( theClass[t], states[the_path[t]]->theClass ) )
			{
				printf("Fatal error in class match!!!!\n");
			//	exit(1);
				free(back_ptr);
				return 0;
			}
			theClass[t] = states[the_path[t]]->theClass;
		}
		states[the_path[t]]->obs_population[observables[t]] += 1;
	}	

	free(back_ptr);

	int m1_state = -1;
	double m1_prob = 0;
	int m2_state = -1;
	double m2_prob = 0;

	memset( probs, 0, sizeof(double) * 8 );
	for( int s = 0; s < nstates; s++ )
	{
		int pt;
		int nr = sscanf( states[s]->label+strlen(states[s]->label)-5, "_%d_SC", &pt );

		if( nr == 1 )
			probs[pt] += cur_prob[s];

	}

	return 0;
}

int HMMSystem::comparePredictedAndTruePaths( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr_wild = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob_wild[nstates]; // the probability.
	
	int *back_ptr_class = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob_class[nstates]; // the probability.

	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob_class[s] = 0;
		}
		else
			cur_prob_class[s] = init_p[s] * p_obs[observables[0]*nstates+s];
				
		cur_prob_wild[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr_class[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr_class[t*nstates+s] = -1;
				continue;
			}


			back_ptr_class[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob_class[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob_class[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob_class[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr_class[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		memcpy( cur_prob_class, new_probs, sizeof(double) * nstates );
	}
	
	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr_wild[t*nstates+s] = -1;
				continue;
			}


			back_ptr_wild[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob_wild[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob_wild[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob_wild[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr_wild[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		int sorter[nstates];

		for( int s = 0; s < nstates; s++ )
			sorter[s] = s;

		int done = 0;
	
		while( !done )
		{
			done = 1;

			for( int x = 0; x < nstates-1; x++ )
			{
				if( new_probs[sorter[x]] < new_probs[sorter[x+1]] )
				{
					int t = sorter[x];
					sorter[x] = sorter[x+1];
					sorter[x+1] = t;
					done = 0;
				}
	
			}
		}
		memcpy( cur_prob_wild, new_probs, sizeof(double) * nstates );
	}


	int the_path_class[len], the_path_wild[len];

	double best_prob_class = 0;
	double best_prob_wild = 0;
	the_path_class[len-1] = 0;
	the_path_wild[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob_class[s] > best_prob_class )
		{
			best_prob_class = cur_prob_class[s];
			the_path_class[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
		the_path_class[t] = back_ptr_class[(t+1)*nstates+the_path_class[t+1]];
	
	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob_wild[s] > best_prob_wild )
		{
			best_prob_wild = cur_prob_wild[s];
			the_path_wild[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
		the_path_wild[t] = back_ptr_wild[(t+1)*nstates+the_path_wild[t+1]];


	
	printf("Restrained path: (%le)\n", best_prob_class );
	double psum = 0;
	double psum_o = 0;
	for( int t = 1; t < len; t++ )
	{
		int s1 = the_path_class[t-1];
		int s2 = the_path_class[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];

		psum += log(p);
		psum_o += log( p_obs[observables[t]*nstates+s2] );

		printf("t=%d state: '%s' p_transition: %le\n", t, states[the_path_class[t]]->label, p );
			
	}
	printf("log Prob: %le probo: %le\n", psum, psum_o );
	
	printf("Chosen path: (%le)\n", best_prob_wild );
	psum = 0;
	psum_o = 0;
	for( int t = 1; t < len; t++ )
	{
		int s1 = the_path_wild[t-1];
		int s2 = the_path_wild[t];

		double p = p_obs[observables[t]*nstates+s2] * t_ab[s1*nstates+s2];
		psum_o += log( p_obs[observables[t]*nstates+s2] );

		printf("t=%d state: '%s' p_transition: %le\n", t, states[the_path_wild[t]]->label, p );
		psum += log(p);	
	}
	printf("log Prob: %le probo: %le\n", psum, psum_o );
	

	free(back_ptr_class);
	free(back_ptr_wild);

	return 1;
}

void calculateAlphaForStringForcePath( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem, int *state_path )
{
	int   nstates = theSystem.getNStates();
	double *p_obs = theSystem.getObservableP();
	double *t_ab  = theSystem.getTransitionP();
	double *init_p = theSystem.getInitP();

	// list of observables.	

	// alpha[0] 

	for( int s = 0; s < nstates; s++ )
			alpha[s] = 0;

	alpha[state_path[0]] = 1.0;

	double ptot= 0;

	for( int x = 0; x < nstates; x++ )
		ptot += alpha[x];
	for( int x = 0; x < nstates; x++ )
		alpha[x] /= ptot;	

	alpha_norm[0] = log(ptot);

	for( int t = 1; t < len; t++ )
	{
		double prev_sum = 0;

		for( int x = 0; x < nstates; x++ )
			alpha[t*nstates+x] = 0;

		alpha[t*nstates+state_path[t]] = t_ab[state_path[t-1]*nstates+state_path[t]] * p_obs[ observables[t] * nstates + state_path[t]];

		ptot =0;

		for( int x = 0; x < nstates; x++ )
			ptot += alpha[t*nstates+x];

		for( int x = 0; x < nstates; x++ )
			alpha[t*nstates+x] /= ptot;	

		alpha_norm[t] = log(ptot) + alpha_norm[t-1];
	}

	// reweight alphas based on 
}


// beta(T), the probability of observing the 
	
void calculateBetaForStringForcePath( double *beta, double *beta_norm, int *observables, char *classes, int len, HMMSystem &theSystem, int *the_path )
{
	int   nstates = theSystem.getNStates();
	double *p_obs = theSystem.getObservableP();
	double *t_ab  = theSystem.getTransitionP();
	double *init_p = theSystem.getInitP();

	// list of observables.	

	for( int s = 0; s < nstates; s++ )
		beta[(len-1)*nstates+s] = 1.0 / nstates;

	double ptot =0;

	for( int x = 0; x < nstates; x++ )
		ptot += beta[(len-1)*nstates+x];

	for( int x = 0; x < nstates; x++ )
		beta[(len-1)*nstates+x] /= ptot;	

	beta_norm[(len-1)] = log(ptot);
	
	for( int t = len-2; t >= 0; t-- )
	{
		for( int s = 0; s < nstates; s++ )
			beta[t*nstates+s] = 0;

		beta[t*nstates+the_path[t]] = t_ab[the_path[t]*nstates+the_path[t+1]] * p_obs[observables[t+1]*nstates+the_path[t+1]];
		
		ptot =0;
	
		for( int x = 0; x < nstates; x++ )
			ptot += beta[t*nstates+x];

		for( int x = 0; x < nstates; x++ )
			beta[t*nstates+x] /= ptot;	
	
		beta_norm[t] = log(ptot) + beta_norm[t+1];
	}

}

double HMMSystem::ProbActual( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.
	
	int do_SC = (len == 40 );

	if( do_SC )
	{
		double n_SC = 0;

		for( int s = 0; s < nstates; s++ )
		{
			char *tname = states[s]->label;
			int tlen = strlen(tname);
			if( tlen >= 3 && tname[tlen-1] == 'C' && tname[tlen-2] == 'S' && tname[tlen-3] == '_' && !strncasecmp( tname, "Cap2StateBackward", strlen("Cap2StateBackward") ) )
				n_SC += 1;
		}

		for( int s = 0; s < nstates; s++ )
		{
			char *tname = states[s]->label;
			int tlen = strlen(tname);
			if( tlen >= 3 && tname[tlen-1] == 'C' && tname[tlen-2] == 'S' && tname[tlen-3] == '_' && !strncasecmp( tname, "Cap2StateBackward", strlen("Cap2StateBackward") ) )
				cur_prob[s] = (1.0/n_SC) * p_obs[ observables[0] * nstates + s];
			else
				cur_prob[s] = 0;
		}
	}
	else
	{
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[0], states[s]->theClass) )
			{
					cur_prob[s] = 0;
			}
			else
				cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
		}	
	}

	double lp_sum = 0;
	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	

		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		lp_sum += log(ptot);

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}

	double psum = 0;

	for( int s = 0; s < nstates; s++ )
		psum += cur_prob[s];	

	free(back_ptr);


	return lp_sum;
}

int HMMSystem::estimateClassViterbiSort( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *h1, double *h2, double *h3, double *h4 )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.
	
	int do_SC = (len == 40 );

	if( do_SC )
	{
		double n_SC = 0;

		for( int s = 0; s < nstates; s++ )
		{
			char *tname = states[s]->label;
			int tlen = strlen(tname);
			if( tlen >= 3 && tname[tlen-1] == 'C' && tname[tlen-2] == 'S' && tname[tlen-3] == '_' && !strncasecmp( tname, "Cap2StateBackward", strlen("Cap2StateBackward") ) )
				n_SC += 1;
		}

		for( int s = 0; s < nstates; s++ )
		{
			char *tname = states[s]->label;
			int tlen = strlen(tname);
			if( tlen >= 3 && tname[tlen-1] == 'C' && tname[tlen-2] == 'S' && tname[tlen-3] == '_' && !strncasecmp( tname, "Cap2StateBackward", strlen("Cap2StateBackward") ) )
				cur_prob[s] = (1.0/n_SC) * p_obs[ observables[0] * nstates + s];
			else
				cur_prob[s] = 0;
		}
	}
	else
	{
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[0], states[s]->theClass) )
			{
					cur_prob[s] = 0;
			}
			else
				cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
		}	
	}


	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	

		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
	/*	if( t_ab[the_path[t]*nstates+the_path[t+1]] < 1e-300 )
		{
//			printf("Fatal error, improbable transition chosen.\n");
//			free(back_ptr);
			printf("Woah, improbable transition chosen!\n");
				
			return 0;
		} */
	}

	*h1 = 0;
	*h2 = 0;
	*h3 = 0;
	*h4 = 0;

	double *hels[4] = { h1, h2, h3, h4 };

	for( int s = 0; s < nstates; s++ )
	{
		char *tpt = states[s]->label;
		
		while( *tpt != '_' && *tpt != '\0' ) tpt += 1;
	
		if( !*tpt) continue;	
	
		char junk[1024];
		int p, h;
		char t1,t2;
		int nr = sscanf( tpt, "_%d_%d_%c%c", &p, &h, &t1, &t2 );
		if( nr != 4 ) continue;
		
		if( t1 != 'S' || t2 != 'C' )
			continue;

		*(hels[h]) += cur_prob[s];	
	}
	
	
	free(back_ptr);

	return 1;
}

double isItProbable( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem )
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	int do_SC = 1;

	if( len != 40 ) do_SC = 0;

	int   nstates = theSystem.getNStates();
	double *p_obs = theSystem.getObservableP();
	double *t_ab  = theSystem.getTransitionP();
	double *init_p = theSystem.getInitP();

	int do_err = 0;
	for( int s = 0; s < nstates; s++ )
	{
		if( !class_match( classes[0], theSystem.states[s]->theClass)  ) // classes[0] != '.' && theSystem.states[s]->theClass != classes[0] )
			alpha[s] = 0;
		else
		{
			alpha[s] = init_p[s] * p_obs[ observables[0] * nstates + s];
		}
	}

	double ptot= 0;

	for( int x = 0; x < nstates; x++ )
		ptot += alpha[x];
	for( int x = 0; x < nstates; x++ )
	{
		alpha[x] /= ptot;	
	}
	alpha_norm[0] = log(ptot);

	FILE *err_file = NULL;
	char err_file_name[256];
	sprintf(err_file_name, "err.%d.file", taskid );

	for( int t = 1; t < len; t++ )
	{
		double prev_sum = 0;

		for( int s = 0; s < nstates; s++ )
		{
			int state_list[MAX_TRANSITIONS];
			int ntrans;

			alpha[t*nstates+s] = 0;

			theSystem.getTransitionsTo( s, state_list, &ntrans, '.' );
			if( !class_match( classes[t], theSystem.states[s]->theClass) ) // classes[t] != '.' && classes[t] != theSystem.states[s]->theClass )
				continue;
			
	
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];
				alpha[t*nstates+s] += alpha[(t-1)*nstates+s2] * t_ab[s2*nstates+s] * p_obs[ observables[t] * nstates + s];
			}
		}
		

			ptot =0;
			for( int x = 0; x < nstates; x++ )
			ptot += alpha[t*nstates+x];

		if( fabs(ptot) < 1e-300 )
			return -1.0;

		for( int x = 0; x < nstates; x++ )
			alpha[t*nstates+x] /= ptot;	

		alpha_norm[t] = log(ptot) + alpha_norm[t-1];
	}

	// reweight alphas based on 

	return alpha_norm[len-1];
}

double HMMSystem::getp7( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *dP_7, int to_proc)
{
	int taskid;
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	// at each state we need a ptr backwards to whichever state led to it

	int *back_ptr = (int *)malloc( sizeof(int) * len * nstates );
	double cur_prob[nstates]; // the probability.

	for( int s = 0; s < nstates; s++ )
	{
		if( theClass && !class_match( theClass[0], states[s]->theClass) )
		{
				cur_prob[s] = 0;
		}
		else
			cur_prob[s] = init_p[s] * p_obs[observables[0]*nstates+s];
	}

	int do_print = 0;

	
	dP_7[0] = 0;

	for( int t = 1; t < len; t++ )
	{
		double new_probs[nstates];
	
		for( int s = 0; s < nstates; s++ )
		{
			if( theClass && !class_match( theClass[t], states[s]->theClass) )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue; 
			}

			new_probs[s] = -1;

			int trans_list[MAX_TRANSITIONS];
			int ntrans;

			getTransitionsTo( s, trans_list, &ntrans, '.' );
			
			if( ntrans == 0 )
			{
				new_probs[s] = 0;
				back_ptr[t*nstates+s] = -1;
				continue;
			}


			back_ptr[t*nstates+s] = trans_list[0];
				int s_prev = trans_list[0];
			new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
		
			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s_prev = trans_list[xs];
				if( cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s] > new_probs[s])
				{
					new_probs[s] = cur_prob[s_prev] * t_ab[s_prev*nstates+s] * p_obs[observables[t]*nstates+s];
					back_ptr[t*nstates+s] = s_prev;
				}
			}		
		}

		for( int s = 0; s < nstates; s++ )
		{
			if( new_probs[s] < 0 )
				new_probs[s] = 0;
		}

		double ptot = 0;
		for( int s = 0; s < nstates; s++ )
			ptot += new_probs[s];
		for( int s = 0; s < nstates; s++ )
			new_probs[s] /= ptot;

		memcpy( cur_prob, new_probs, sizeof(double) * nstates );

		double p7 = 0;

		for( int s = 0; s < nstates; s++ )
		{
			int pt;
			int nr = sscanf( states[s]->label+strlen(states[s]->label)-5, "_%d_SC", &pt );
	
			if( nr == 1 && pt == to_proc )
				p7 += cur_prob[s];
		}

		dP_7[t] = p7;

		if( do_print )
		{
			for( int s = 0; s < nstates; s++ )
				if( cur_prob[s] > 1e-10 )
					printf("Allowed at t=%d state %s\n", t, states[s]->label );
		}
	}


	int the_path[len];

	double best_prob = 0;
	the_path[len-1] = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( cur_prob[s] > best_prob )
		{
			best_prob = cur_prob[s];
			the_path[len-1] = s;
		} 
	}

	for( int t = len-2; t >= 0; t-- )
	{
		if( back_ptr[(t+1)*nstates+the_path[t+1]] == -1 )
		{
			printf("Fatal error, back ptr is uninitialized.\n");
//			exit(1);
			free(back_ptr);
			return 0;
		}
		the_path[t] = back_ptr[(t+1)*nstates+the_path[t+1]];
		if( t_ab[the_path[t]*nstates+the_path[t+1]] < 1e-20 )
		{
			printf("Fatal error, improbable transition chosen.\n");
//			free(back_ptr);
			return 0;
		} 
	}

	for( int t = 0; t < len; t++ )
	{
		if( ML_state )
			ML_state[t] = the_path[t];

		if( !theClass )
		{
		
			if( taskid == 0 )
				printf("%c %c\n", 'A' + observables[t], states[the_path[t]]->theClass );
		}
		else
		{
			if( !class_match( theClass[t], states[the_path[t]]->theClass ) )
			{
				printf("Fatal error in class match!!!!\n");
			//	exit(1);
				free(back_ptr);
				return 0;
			}
			theClass[t] = states[the_path[t]]->theClass;
		}
		states[the_path[t]]->obs_population[observables[t]] += 1;
	}	
	free(back_ptr);


	return 0;
}
