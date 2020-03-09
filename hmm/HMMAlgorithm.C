#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include <sys/time.h>
#define DO_CF
//#define DIVIDE_GS
#define NANCHECK
#define GRAD_NANCHECK
#include "HMM.h"
#include "HMMDirective.h"
//#define DO_PMAT_NORMALIZATION
extern int globalN_OBSERVABLES;
static const double eps = 1e-100;//1e-30;

int my_isnan( double val )
{
	if( !(val<0) && !(val>-1) )
	{
		return 1;
	}
	double ival = val / (1+val);

	if( !(ival<0) && !(ival>-1) )
		return 1;
	return 0;
}

HMMState::HMMState(char classIn, const char *label_in, int possibleInitialStateIn, int emission_group_in, int pseudostate_in, int transition_group_in)
{
	strcpy( label, label_in );
	theClass = classIn;
	pseudostate = pseudostate_in;
	emission_group = emission_group_in;	
	transition_group = transition_group_in;	
	possibleInitialState = possibleInitialStateIn;
	index = -1; // HMMSystem will setup the indexes when it knows what all the states are.
	ntransitions = 0;	
}

void HMMSystem::randomizeTimeSequence( HMMDirective *directive )
{
	for( int d = 0; d < directive->n_data; d ++ )
	{
		int len = directive->translated_len[d];
		int *to_grab = (int *)malloc( sizeof(int) * len );
		int *to_put  = (int *)malloc( sizeof(int) * len );
	
		for( int x = 0; x < len; x++ )
			to_grab[x] = x;
	
		for( int x = 0; x < len; x++ )
		{
			int cur_num = (len - x);
			
			int val = rand() % cur_num;
	
			to_put[x] = to_grab[val];
			to_grab[val] = to_grab[cur_num-1];				
		}
	
	
		int *trans_obs = (int *)malloc( sizeof(int) * len );
	
		for( int x = 0; x < len; x++ )
			trans_obs[x] = directive->translated_observable_string[d][to_put[x]];

		memcpy( directive->translated_observable_string[d], trans_obs, sizeof(int) * len );

		free(to_grab);	
		free(to_put);
		free(trans_obs);
	}
}

void HMMState::addTransition( HMMState *theState )
{
	if( ntransitions >= MAX_TRANSITIONS )
	{	
		printf("Fatal error, attempting to assign %d transitions to a state.\n", ntransitions+1 );
		exit(1);
	}

	p_so[ntransitions].elem = theState;
	p_so[ntransitions].p = 0;

	ntransitions++;
}

HMMSystem::HMMSystem()
{
	nstates = 0;
	nstatesSpace = 10;

	states = (HMMState **)malloc( sizeof(HMMState *) * nstatesSpace );	
}

HMMSystem::~HMMSystem()
{
	for( int s = 0; s < nstates; s++ )
		delete states[s];
	free(states);	

	free(p_obs);
	free(t_ab);
	free(p_obs_BaumWelch);
	free(t_ab_BaumWelch);
	free(p_obs_prev);
	free(t_ab_prev);
	free(p_obs_next_free);
	free(t_ab_next_free);
	free(p_obs_next_clamped);
	free(t_ab_next_clamped);
	free(group_is_fixed);

	free(init_p);
	free(init_p_next_free);
	free(init_p_next_clamped);
	free(init_p_BaumWelch);
	free(pseudocounts);
	free(pseudo_add_clamped);
	free(pseudo_add);
	for( int t = 0; t < nstates; t++ )
	{
		free(transitionListFrom[t]);
		free(transitionListTo[t]);
	}
	free(transitionListFrom);
	free(transitionListTo);
	free(transitionSizeFrom);
	free(transitionSizeTo);
	for( int x = 0; x < nEmissionGroups; x++ )
		free( emission_groups[x] );
	free(emission_groups);
	free(emission_groups_len);
	for( int x = 0; x < nTransitionGroups; x++ )
		free( transition_groups[x] );
	free(transition_groups);
	free(transition_groups_len);

	for( int x = 0; x < MAX_AUX_VARS; x++ )
		free(aux_vars[x]);
}


void HMMSystem::addState( HMMState *theState )
{
	if( nstates == nstatesSpace )
	{
		nstatesSpace *= 2;
		states = (HMMState **)realloc( states, sizeof( HMMState *) * nstatesSpace );
	}

	states[nstates] = theState;
	theState->index = nstates;

	nstates++;
}

void HMMSystem::randomizeObs( void )
{
	for( int s = 0; s < nstates; s++ )
	{
		double sum = 0;
		for( int t = 0; t < globalN_OBSERVABLES; t++ )
		{
			sum += (p_obs[t*nstates+s] = (int)100*(rand()/(double)RAND_MAX));
		}
		for( int t = 0; t < globalN_OBSERVABLES; t++ )
			p_obs[t*nstates+s] /= sum;
	}
}

void HMMSystem::setup( void )
{
	int taskid=0;

#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
	if( taskid == 0 ) printf("Setting up system with %d states.\n", nstates );

	p_obs = (double *)malloc( sizeof(double) * (nstates) * (globalN_OBSERVABLES+1) );
	t_ab =  (double *)malloc( sizeof(double) * nstates * nstates );
	
	p_obs_BaumWelch = (double *)malloc( sizeof(double) * nstates * (globalN_OBSERVABLES+1) );
	t_ab_BaumWelch =  (double *)malloc( sizeof(double) * nstates * nstates );
	
	p_obs_prev = (double *)malloc( sizeof(double) * nstates * (globalN_OBSERVABLES+1) );
	t_ab_prev =  (double *)malloc( sizeof(double) * nstates * nstates );
	
	p_obs_next_free = (double *)malloc( sizeof(double) * nstates * (globalN_OBSERVABLES+1) );
	t_ab_next_free =  (double *)malloc( sizeof(double) * nstates * nstates );
	
	p_obs_next_clamped = (double *)malloc( sizeof(double) * nstates * (globalN_OBSERVABLES+1) );
	t_ab_next_clamped =  (double *)malloc( sizeof(double) * nstates * nstates );
	
	init_p = (double *)malloc( sizeof(double) * nstates );
	init_p_next_clamped = (double *)malloc( sizeof(double) * nstates );	
	init_p_next_free = (double *)malloc( sizeof(double) * nstates );	
	init_p_BaumWelch = (double *)malloc( sizeof(double) * nstates );	

	memset( init_p, 0, sizeof(double) * nstates );
	memset( init_p_next_clamped, 0, sizeof(double) * nstates );
	memset( init_p_next_free, 0, sizeof(double) * nstates );
	memset( init_p_BaumWelch, 0, sizeof(double) * nstates );

	pseudocounts = (double *) malloc( sizeof(double) * 256 * globalN_OBSERVABLES);
	memset( pseudocounts, 0, sizeof(double) * 256 );
	
	pseudo_add_clamped = (double *)malloc( sizeof(double) * globalN_OBSERVABLES * nstates );
	pseudo_add = (double *)malloc( sizeof(double) * globalN_OBSERVABLES );
	memset( pseudo_add_clamped, 0, sizeof(double) * globalN_OBSERVABLES * nstates );
	memset( pseudo_add, 0, sizeof(double) * globalN_OBSERVABLES );

	for( int x = 0; x < MAX_AUX_VARS; x++ )
		aux_vars[x] = (double *)malloc( sizeof(double) * (globalN_OBSERVABLES*nstates+nstates*nstates) );	

	double n_init_states = 0;
	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->possibleInitialState )
			n_init_states += 1;
	}

	for( int s = 0; s < nstates; s++ )
		if( states[s]->possibleInitialState )
			init_p[s] = 1.0 / n_init_states;

	for( int s = 0; s < nstates; s++ )
	{
		double sum = 0;
		for( int t = 0; t < globalN_OBSERVABLES; t++ )
		{
			sum += (p_obs[t*nstates+s] = (int)100*(rand()/(double)RAND_MAX));
		}
		for( int t = 0; t < globalN_OBSERVABLES; t++ )
			p_obs[t*nstates+s] /= sum;
	}
//	for( int s = 0; s < nstates; s++ )
//	for( int t = 0; t < globalN_OBSERVABLES; t++ )
//		p_obs[t*nstates+s] = (1.0 / globalN_OBSERVABLES);
	for( int s = 0; s < nstates; s++ )
		p_obs[globalN_OBSERVABLES*nstates+s] = 1.0;

	memset( t_ab, 0, sizeof(double) * nstates * nstates );

	transitionListFrom = (int **)malloc( sizeof(int *) * nstates );
	transitionSizeFrom =  (int *)malloc( sizeof(int) * nstates );

	
	transitionListTo = (int **)malloc( sizeof(int *) * nstates );
	transitionSizeTo =  (int *)malloc( sizeof(int) * nstates );
	
	memset( transitionSizeFrom, 0, sizeof(int) * nstates);
	memset( transitionSizeTo, 0, sizeof(int) * nstates);

	int temp_n_from[nstates];

	for( int s = 0; s < nstates; s++ )
	{
/*
		if( states[s]->ntransitions == 0 )
		{
			printf("Fatal error. State '%s' has no transitions!\n", states[s]->label );
			exit(1);
		}
*/
		for( int tr = 0; tr < states[s]->ntransitions; tr++ )
		{
			HMMState *to_state = states[s]->p_so[tr].elem;

			

 			int t_from = s;
			int t_to   = to_state->index;
			
			if( t_to == -1 )
			{
				printf("HOLD ON.\n");
				exit(1);
			}


			transitionSizeFrom[t_from] += 1;
			transitionSizeTo[t_to] += 1;
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		int len_alloc = transitionSizeFrom[s];
		if( len_alloc < 1 ) len_alloc = 1;

		transitionListFrom[s] = (int *)malloc( sizeof(int) * len_alloc );
		len_alloc = transitionSizeTo[s];
		if( len_alloc < 1 ) len_alloc = 1;
		transitionListTo[s]   = (int *)malloc( sizeof(int) * len_alloc );

		// from, to
		temp_n_from[s] = transitionSizeFrom[s];

		transitionSizeFrom[s] = 0;
		transitionSizeTo[s] = 0;
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		for( int tr = 0; tr < states[s]->ntransitions; tr++ )
		{
			HMMState *to_state = states[s]->p_so[tr].elem;

 			int t_from = s;
			int t_to   = to_state->index;

			if( t_to == -1 )
			{
				printf("HOLD ON.\n");
				exit(1);
			}

			t_ab[t_from*nstates+t_to] = (1.0 / temp_n_from[t_from] );

			transitionListFrom[t_from][transitionSizeFrom[t_from]] = t_to;
			transitionListTo[t_to][transitionSizeTo[t_to]] = t_from;

			transitionSizeFrom[t_from] += 1;
			transitionSizeTo[t_to] += 1;
		}
	}

	// do Emission groups

	nEmissionGroups = 0;	
	int t_groups[nstates];
	int e_group_lens[nstates];
	memset( e_group_lens, 0, sizeof(int) * nstates );
	int t_n_groups = 0;
	int total_groups = 0;	
	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->emission_group == -1 )
		{
			total_groups++;
			continue;
		}		

		int got_it = 0;

		for( int x = 0; x < t_n_groups; x++ )
		{
			if( t_groups[x] == states[s]->emission_group )
			{
				got_it = 1;
			}
		}

		if( !got_it )
		{
			t_groups[t_n_groups] = states[s]->emission_group;
			t_n_groups++;
			total_groups++;
		}
	}
	nEmissionGroups = total_groups;

	emission_groups = (int **)malloc( sizeof(int *)  * nEmissionGroups );
	emission_groups_len = (int *)malloc( sizeof(int) * nEmissionGroups );

	group_is_fixed = (int *)malloc( sizeof(int) * nEmissionGroups );
	memset( group_is_fixed, 0, sizeof(int) * nEmissionGroups );

	memset( emission_groups_len, 0, sizeof(int) * nEmissionGroups );	

	int temp_emission_group[nstates];

	int g_cntr = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->emission_group == -1 ) // create new group.
		{
			temp_emission_group[s] = g_cntr;
			emission_groups[g_cntr] = (int *)malloc( sizeof(int) * 1 );
			emission_groups[g_cntr][0] = s;
			emission_groups_len[g_cntr] = 1;
			g_cntr++; 
		}
		else
		{
			int got_it = 0;
			for( int g = 0; g < g_cntr ; g++ )
			{
				if( states[emission_groups[g][0]]->emission_group == states[s]->emission_group )
				{
					temp_emission_group[s] = g;
					emission_groups[g][emission_groups_len[g]] = s;
					emission_groups_len[g]++;
					got_it = 1;
					break;
				}
			}

			if( !got_it )
			{
				temp_emission_group[s] = g_cntr;
				emission_groups[g_cntr] = (int *)malloc( sizeof(int) * nstates );
				emission_groups[g_cntr][0] = s;
				emission_groups_len[g_cntr] = 1;
				g_cntr++;
			}
		}
	}

	for( int s = 0; s < nstates; s++ )
		states[s]->emission_group = temp_emission_group[s];

	if( taskid == 0 ) printf("%d total emission groups.\n", nEmissionGroups );
	
	// do Transition groups

	nTransitionGroups = 0;	
	memset( e_group_lens, 0, sizeof(int) * nstates );
	t_n_groups = 0;
	total_groups = 0;	
	int temp_transition_group[nstates];
	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->transition_group == -1 )
		{
			total_groups++;
			continue;
		}		

		int got_it = 0;

		for( int x = 0; x < t_n_groups; x++ )
		{
			if( t_groups[x] == states[s]->transition_group )
			{
				got_it = 1;
			}
		}

		if( !got_it )
		{
			t_groups[t_n_groups] = states[s]->transition_group;
			t_n_groups++;
			total_groups++;
		}
	}
	nTransitionGroups = total_groups;

	transition_groups = (int **)malloc( sizeof(int *)  * nTransitionGroups );
	transition_groups_len = (int *)malloc( sizeof(int) * nTransitionGroups );

	memset( transition_groups_len, 0, sizeof(int) * nTransitionGroups );	

	g_cntr = 0;

	for( int s = 0; s < nstates; s++ )
	{
		if( states[s]->transition_group == -1 ) // create new group.
		{
			temp_transition_group[s] = g_cntr;
			transition_groups[g_cntr] = (int *)malloc( sizeof(int) * 1 );
			transition_groups[g_cntr][0] = s;
			transition_groups_len[g_cntr] = 1;
			g_cntr++; 
		}
		else
		{
			int got_it = 0;
			for( int g = 0; g < g_cntr ; g++ )
			{
				if( states[transition_groups[g][0]]->transition_group == states[s]->transition_group )
				{
					temp_transition_group[s] = g;
					transition_groups[g][transition_groups_len[g]] = s;
					transition_groups_len[g]++;
					got_it = 1;
					break;
				}
			}

			if( !got_it )
			{
				temp_transition_group[s] = g_cntr;
				transition_groups[g_cntr] = (int *)malloc( sizeof(int) * nstates );
				transition_groups[g_cntr][0] = s;
				transition_groups_len[g_cntr] = 1;
				g_cntr++;
			}
		}
	}

	for( int s = 0; s < nstates; s++ )
		states[s]->transition_group = temp_transition_group[s];

	TransitionGroupCheck();
}


void HMMSystem::getTransitionsFrom( int s, int *state_list, int *ntransitions, char class_restriction )
{
	if( class_restriction == '.' )
	{
		memcpy( state_list, transitionListFrom[s], sizeof(int) * transitionSizeFrom[s] );
		*ntransitions = transitionSizeFrom[s];
		return;
	}

	*ntransitions = 0;

	for( int xs = 0; xs < transitionSizeFrom[s]; xs++ )
	{
		int s_abs = transitionListFrom[s][xs];

		if( states[s_abs]->theClass == class_restriction || class_restriction == '.' )
		{	
			state_list[*ntransitions] = s_abs;
			(*ntransitions)++;
		}
	} 
}

void HMMSystem::getTransitionsTo( int s, int *state_list, int *ntransitions, char class_restriction )
{
	if( class_restriction == '.' )
	{
		memcpy( state_list, transitionListTo[s], sizeof(int) * transitionSizeTo[s] );
		*ntransitions = transitionSizeTo[s];
		return;
	}

	*ntransitions = 0;

	for( int xs = 0; xs < transitionSizeTo[s]; xs++ )
	{
		int s_abs = transitionListTo[s][xs];

		if( states[s_abs]->theClass == class_restriction || class_restriction == '.' )
		{	
			state_list[*ntransitions] = s_abs;
			(*ntransitions)++;
		}
	} 
}



int HMMSystem::getNStates(void) { return nstates; }

double *HMMSystem::getObservableP( void ) { return p_obs; } 
double *HMMSystem::getTransitionP( void ) { return t_ab; }
double *HMMSystem::getInitP( void ) { return init_p; }


void HMMSystem::zeroNextPMat( void )
{
	memset( t_ab_next_free, 0, sizeof(double) * nstates * nstates );
	memset( p_obs_next_free, 0, sizeof(double) * nstates * globalN_OBSERVABLES );
		
	memset( t_ab_next_clamped, 0, sizeof(double) * nstates * nstates );
	memset( p_obs_next_clamped, 0, sizeof(double) * nstates * globalN_OBSERVABLES );
	
	memset( init_p_next_clamped, 0, sizeof(double) * nstates );
	memset( init_p_next_free, 0, sizeof(double) * nstates );
}

void HMMSystem::incrementNextPMat_Free( double *alphas, double *betas, int *observables,  int len, double weight, double scale  )
{
	incrementNextPMat_GEN( alphas, betas, observables, len, t_ab_next_free, p_obs_next_free, init_p_next_free, weight, scale );
}

void HMMSystem::incrementNextPMat_Clamped( double *alphas, double *betas, int *observables, int len, double weight, double scale  )
{
	incrementNextPMat_GEN( alphas, betas, observables, len, t_ab_next_clamped, p_obs_next_clamped, init_p_next_clamped, weight, scale );
}

void HMMSystem::incrementNextPMat_GEN( double *alphas,  double *betas, int *observables, int len, double *t_ab_next, double *p_obs_next, double *init_p_next, double weight, double scale )
{
	double normal_t[len];
	memset( normal_t, 0, sizeof(double) * len );
	
	double normal_p[len];
	memset( normal_p, 0, sizeof(double) * len );

	for( int s = 0; s < nstates; s++ )
		if( alphas[s] * betas[s] > 1e-50 )
		{
//			printf("%s was used as a start state (%le).\n",
//				states[s]->label, alphas[s]*betas[s] );
		}

	for( int t = 0; t < len-1; t++ )
	{
		normal_t[t] = 0;

		for( int s1 = 0; s1 < nstates; s1++ )
		{
			int state_list[MAX_TRANSITIONS];
			double pvals[MAX_TRANSITIONS];

			int ntrans;

			getTransitionsFrom( s1, state_list, &ntrans, '.' );

			for( int xs = 0; xs < ntrans; xs++ )
			{
				int s2 = state_list[xs];

					normal_t[t] += alphas[t*nstates+s1] * betas[(t+1)*nstates+s2] * 
					p_obs[observables[t+1]*nstates+s2]  *t_ab[s1*nstates+s2];
/*
				if( normal_t[t] < 0 )
				{
					printf("alpha: %lf betas: %lf p_obs: %lf t_ab: %lf\n",	
						alphas[t*nstates+s1], betas[(t+1)*nstates+s2], p_obs[observables[t+1]*nstates+s2], t_ab[s1*nstates+s2] );
					exit(1);
				} */
			}
		}
	}

	for( int t = 0; t < len; t++ )
	{
		normal_p[t] = 0;

		for( int s1 = 0; s1 < nstates; s1++ )
			normal_p[t] +=alphas[t*nstates+s1] * betas[t*nstates+s1];
	}

	for( int s = 0; s < nstates; s++ )
	{
		init_p_next[s] += scale * weight * alphas[0*nstates+s] * betas[0*nstates+s] / normal_p[0];
	}
	double state_p[nstates];
	memset( state_p, 0, sizeof(double) * nstates );
	
	double den = 0;

	for( int s = 0; s < nstates; s++ )
	{
		state_p[s] = 1e-100;
		for( int t = 0; t < len; t++ )
		{
			state_p[s] += alphas[t*nstates+s] * betas[t*nstates+s];
		}	
	}

	for( int s = 0; s < nstates; s++ )
	{
		int state_list[MAX_TRANSITIONS];
		double pvals[MAX_TRANSITIONS];

		int ntrans;

		getTransitionsFrom( s, state_list, &ntrans, '.' );

		for( int xs = 0; xs < ntrans; xs++ )
		{
			pvals[xs] = 0;

			int to_state = state_list[xs];

			for( int t = 0; t < len-1; t++ )
			{
#ifdef NANCHECK
				if( my_isnan( alphas[t*nstates+s] ) || my_isnan( betas[(t+1)*nstates+to_state]  ) || my_isnan(normal_t[t])  )
				{
					printf("Nan problem here.\n");
					exit(1);
				}
#endif

				pvals[xs] += alphas[t*nstates+s] * t_ab[s*nstates+to_state] * betas[(t+1)*nstates+to_state] * p_obs[observables[t+1]*nstates+to_state] / normal_t[t];
		
			}

#ifdef NANCHECK
			if( !(pvals[xs]<0) && !(pvals[xs]>-1))
			{
				printf("nan pval!!\n");
				exit(1);
			}		
#endif
			t_ab_next[s*nstates+state_list[xs]] += scale * pvals[xs] * weight;
		}
	}

	for( int t = 0; t < len; t++ )
	{
		for( int s = 0; s < nstates; s++ )
		{
			if( normal_p[t] < 1e-50 ) continue;
/*			{
				printf("bailing, check this.\n");
				exit(1);
			} */

			if( s == 0 && observables[t] == 0 )
			{
//				printf("at zero-zero increment step, alpha: %le beta: %le normal_p: %le\n",
//						alphas[t*nstates+s], betas[t*nstates+s], normal_p[t] );
			}

#ifdef NANCHECK
			double tval = scale * weight * alphas[t*nstates+s] * betas[t*nstates+s] / normal_p[t];
			if( !(tval<0) && !(tval>-1))
			{
				printf("eh?\n");
				exit(1);
			}
#endif
		
			p_obs_next[observables[t]*nstates+s] += scale * weight * alphas[t*nstates+s] * betas[t*nstates+s] / normal_p[t]; // / state_p[s];
		}
	}
}

void HMMSystem::finishNextPMat( void )
{
//#ifndef DO_PMAT_NORMALIZATION
//	return;
//#endif

	for( int s = 0; s < nstates; s++ )
	{
		int state_list[MAX_TRANSITIONS];
		double pvals[MAX_TRANSITIONS];

		int ntrans;

		getTransitionsFrom( s, state_list, &ntrans, '.' );

		double ptot_clamped = 1e-30;
		double ptot_free = 1e-30;

		for( int xs = 0; xs < ntrans; xs++ )
		{
			int to_state = state_list[xs];

//			ptot += t_ab_next_clamped[s*nstates+state_list[xs]] - t_ab_next_free[s*nstates+state_list[xs]];
			ptot_clamped += t_ab_next_clamped[s*nstates+state_list[xs]];
			ptot_free += t_ab_next_free[s*nstates+state_list[xs]];
		}

		if( ptot_clamped < 1e-30 || (!(ptot_clamped<0) && !(ptot_clamped >-1)))
		{
			printf("Invalid ptot %le for computing transitions from state %d.\n", ptot_clamped, s );
//			exit(1);
		}
		
		for( int xs = 0; xs < ntrans; xs++ )
		{
			int to_state = state_list[xs];

			t_ab_next_clamped[s*nstates+state_list[xs]] /= ptot_clamped;
			t_ab_next_free[s*nstates+state_list[xs]] /= ptot_free;

//			if( !(t_ab_next[s*nstates+state_list[xs]] < 0) && !( t_ab_next[s*nstates+state_list[xs]] > -1 ) )
//			{
//				printf("NAN error.\n");
//				exit(1);
//			}
		}
	}


/*
	RIGHT HERE WE NEED TO HAVE NORMALIZATION FOR EACH GROUP!!!

	every observation accmulated should be combined?

	*/
	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		int GS = emission_groups_len[EG];
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			double val_c = 0;
			double val_f = 0;

			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
	
				val_c += p_obs_next_clamped[o*nstates+s];
				val_f += p_obs_next_free[o*nstates+s];
			}
			
			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				p_obs_next_clamped[o*nstates+s] = val_c;
				p_obs_next_free[o*nstates+s] = val_f;
			}
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		double ptot_clamped = 1e-30;
		double ptot_free = 1e-30;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			ptot_clamped += p_obs_next_clamped[o*nstates+s];
			ptot_free += p_obs_next_free[o*nstates+s];
		}
		if( ptot_clamped < 1e-30 || ( !(ptot_clamped<0) && !(ptot_clamped>-1)) )
		{
			printf("Invalid ptot %le.\n", ptot_clamped );
		}

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			p_obs_next_clamped[o*nstates+s] /= ptot_clamped;
			p_obs_next_free[o*nstates+s] /= ptot_free;
//			if( !(p_obs_next[o*nstates+s]<0) && !(p_obs_next[o*nstates+s]>-1))
//			{
//				printf("NAN error.\n");
//				exit(1);
//			}
		}
	}

}

double noise( void )
{
	return (double)rand() / (double)RAND_MAX;
}

double HMMSystem::setNextPMat( double factor, int doBaumWelch, double free_multiplier, double A )
{
	int taskid=0;
#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
/*
//	memcpy( t_ab_prev, t_ab, sizeof(double) * nstates * nstates );
	
	for( int t = 0; t < nstates; t++ )
	{
		int ntrans;
		int trans_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
		
		for( int tx = 0; tx < ntrans; tx++ )
			t_ab_prev[t*nstates+trans_list[tx]] = t_ab[t*nstates+trans_list[tx]];

	}
*/
	memcpy( p_obs_prev, p_obs, sizeof(double) * nstates * globalN_OBSERVABLES );

	double min_alpha = 100.0;
	double squared_change = 0;
	double max_change_p = 0;

	double init_sum_c = 0;
	double init_sum_f = 0;

	for( int s = 0; s < nstates; s++ )
	{
		init_sum_f += init_p_next_free[s];
		init_sum_c += init_p_next_clamped[s];
		// really, these should be the same.
	}

	// Baum-Welch, just set it.
	for( int s = 0; s < nstates; s++ )
	{
		init_p[s] = init_p_next_clamped[s] / init_sum_c;
	}
#ifndef DO_TRANSITION_GROUPS
	for( int s = 0; s < nstates; s++ )
	{
		double expec_c = eps;
		double expec_f = eps;

		int ntrans;
		int trans_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, trans_list, &ntrans, '.' );

		for( int tx = 0; tx < ntrans; tx++ )
			expec_c += t_ab_next_clamped[s*nstates+trans_list[tx]];
		for( int tx = 0; tx < ntrans; tx++ )
			expec_f += t_ab_next_free[s*nstates+trans_list[tx]];
		
		for( int tx = 0; tx < ntrans; tx++ )
		{
			int t = trans_list[tx];

			double old_p = t_ab[s*nstates+t];
			double dn_c = (t_ab_next_clamped[s*nstates+t] - expec_c * t_ab[s*nstates+t])/expec_c;
			double dn_f = (t_ab_next_free[s*nstates+t] - expec_f * t_ab[s*nstates+t])/expec_f;
	
			double del_w = dn_c - dn_f;
			double cur_w = log( t_ab[s*nstates+t] );
			
			double new_w = cur_w + factor * del_w;

			if( doBaumWelch )
				t_ab[s*nstates+t] = exp( del_w );
			else
				t_ab[s*nstates+t] = exp( cur_w ); 
		}
	}
#endif


#ifdef DO_TRANSITION_GROUPS	
	for( int TG = 0; TG < nTransitionGroups; TG++ )
	{
		int a_state = 0;

		int GS = transition_groups_len[TG];
		int nt;
		int trans_list[MAX_TRANSITIONS];
		getTransitionsFrom( transition_groups[TG][0], trans_list, &nt, '.' );				

		double dn_c[nt];
		double dn_f[nt];

		memset( dn_c, 0, sizeof(double) * nt );
		memset( dn_f, 0, sizeof(double) * nt );
			
		for( int GE = 0; GE < GS; GE++ )
		{
			double expec_c = eps;
			double expec_f = eps;
			int nt;
			int trans_list[MAX_TRANSITIONS];
			getTransitionsFrom( transition_groups[TG][GE], trans_list, &nt, '.'  );				
			int s_from = transition_groups[TG][GE];

			for( int o = 0; o < nt; o++ )
			{	
				int s_to = trans_list[o];

				expec_c += t_ab_next_clamped[s_from*nstates+s_to];
				expec_f += t_ab_next_free[s_from*nstates+s_to];			
			}
	
			for( int o = 0; o < nt; o++ )
			{	
				int s_to = trans_list[o];
				
				dn_c[o] += (t_ab_next_clamped[s_from*nstates+s_to]);
			}
		}
	
		
		for( int GE = 0; GE < GS; GE++ )
		{
			int nt;
			int trans_list[MAX_TRANSITIONS];
			getTransitionsFrom( transition_groups[TG][GE], trans_list, &nt, '.' );				
			int s_from = transition_groups[TG][GE];

			for( int o = 0; o < nt; o++ )
			{
				int s_to = trans_list[o];

				t_ab[s_from*nstates+s_to] = dn_c[o];
			}
		}
	}
#endif

	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		int a_state = 0;

		int GS = emission_groups_len[EG];
				
		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];
		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
			
		for( int GE = 0; GE < GS; GE++ )
		{
			double expec_c = eps;
			double expec_f = eps;

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{	
				int s = emission_groups[EG][GE];
				a_state = s;

				expec_c += p_obs_next_clamped[o*nstates+s];
				expec_f += p_obs_next_free[o*nstates+s];			
			}
	
			int a_state = 0;

			int s = emission_groups[EG][GE];
		
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{	
				if( doBaumWelch )
				{
					dn_c[o] += (p_obs_next_clamped[o*nstates+s]);
				}
				else
				{
					dn_c[o] += (p_obs_next_clamped[o*nstates+s] - expec_c * p_obs[o*nstates+s])/expec_c;
					dn_f[o] += (p_obs_next_free[o*nstates+s]    - expec_f * p_obs[o*nstates+s])/expec_f;
				}
			}
		}
	
		
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
				p_obs[o*nstates+s] = dn_c[o];
		}
	}


	Normalize();
	addNoise(A);	
	Normalize();

	squared_change = 0;

	double biggest_change = 0;
	int b_s = 0;
	int b_o = 0;
	double p_p = 0, n_p = 0;
	for( int s = 0; s < nstates; s++ )
	{
//	for( int t = 0; t < nstates; t++ )
//	{
//		squared_change += (t_ab[s*nstates+t] - t_ab_prev[s*nstates+t])*(t_ab[s*nstates+t] - t_ab_prev[s*nstates+t]);
//	}
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
	{
		if( (p_obs[o*nstates+s] - p_obs_prev[o*nstates+s])*(p_obs[o*nstates+s] - p_obs_prev[o*nstates+s]) > biggest_change )
		{
			n_p = p_obs[o*nstates+s];	
			p_p = p_obs_prev[o*nstates+s];
			biggest_change = (p_obs[o*nstates+s] - p_obs_prev[o*nstates+s])*(p_obs[o*nstates+s] - p_obs_prev[o*nstates+s]);
			b_s = s;
			b_o = o;
		}
		squared_change += (p_obs[o*nstates+s] - p_obs_prev[o*nstates+s])*(p_obs[o*nstates+s] - p_obs_prev[o*nstates+s]);
	}
	}

//	printf("biggest change: %le %s %d (%le %le)\n", biggest_change, states[b_s]->label, b_o, n_p, p_p );
//	if(taskid == 0) printf("Squared change: %.14le\n", squared_change );

	return squared_change;
}

double HMMSystem::setNextPMatCF( double factor, double A, double alpha )
{	/* set the next probabilities, trying to maximize clamped over free */
	int taskid=0;
#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
/*
//	memcpy( t_ab_prev, t_ab, sizeof(double) * nstates * nstates );
	
	for( int t = 0; t < nstates; t++ )
	{
		int ntrans;
		int trans_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
		
		for( int tx = 0; tx < ntrans; tx++ )
			t_ab_prev[t*nstates+trans_list[tx]] = t_ab[t*nstates+trans_list[tx]];

	}
*/
	memcpy( p_obs_prev, p_obs, sizeof(double) * nstates * globalN_OBSERVABLES );

	double min_alpha = 100.0;
	double squared_change = 0;
	double max_change_p = 0;

	double init_sum_c = 0;
	double init_sum_f = 0;

	for( int s = 0; s < nstates; s++ )
	{
		init_sum_f += init_p_next_free[s];
		init_sum_c += init_p_next_clamped[s];
		// really, these should be the same.
	}

	// Baum-Welch, just set it.
	for( int s = 0; s < nstates; s++ )
	{
		init_p[s] = init_p_next_clamped[s] / init_sum_c;
	}
#ifndef DO_TRANSITION_GROUPS
	for( int s = 0; s < nstates; s++ )
	{
		double expec_c = 1e-20;
		double expec_f = 1e-20;

		int ntrans;
		int trans_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, trans_list, &ntrans, '.' );

		for( int tx = 0; tx < ntrans; tx++ )
			expec_c += t_ab_next_clamped[s*nstates+trans_list[tx]];
		for( int tx = 0; tx < ntrans; tx++ )
			expec_f += t_ab_next_free[s*nstates+trans_list[tx]];
		
		for( int tx = 0; tx < ntrans; tx++ )
		{
			int t = trans_list[tx];

			double old_p = t_ab[s*nstates+t];

			double dn_c = (t_ab_next_clamped[s*nstates+t] - expec_c * t_ab[s*nstates+t]);
			double dn_f = (t_ab_next_free[s*nstates+t] - expec_f * t_ab[s*nstates+t]);
	
			double del_w = dn_c - dn_f;
			double tval = expec_c * t_ab[s*nstates+t] + alpha * del_w;	
				if( !( tval >= 0) )
					tval = 1e-20;
			t_ab[s*nstates+t] = tval;
		}
	}
#endif


#ifdef DO_TRANSITION_GROUPS	
	for( int TG = 0; TG < nTransitionGroups; TG++ )
	{
		int a_state = 0;

		int GS = transition_groups_len[TG];
		int nt;
		int trans_list[MAX_TRANSITIONS];
		getTransitionsFrom( transition_groups[TG][0], trans_list, &nt, '.' );				

		double dn_c[nt];
		double dn_f[nt];
		double expec_c_save = 1e-20;

		memset( dn_c, 0, sizeof(double) * nt );
		memset( dn_f, 0, sizeof(double) * nt );
		
			
		for( int GE = 0; GE < GS; GE++ )
		{
			double expec_c = eps;
			double expec_f = eps;
			int nt;
			int trans_list[MAX_TRANSITIONS];
			getTransitionsFrom( transition_groups[TG][GE], trans_list, &nt, '.'  );				
			int s_from = transition_groups[TG][GE];

			for( int o = 0; o < nt; o++ )
			{	
				int s_to = trans_list[o];

				expec_c_save += t_ab_next_clamped[s_from*nstates+s_to];

				expec_c += t_ab_next_clamped[s_from*nstates+s_to];
				expec_f += t_ab_next_free[s_from*nstates+s_to];			
			}
	
			for( int o = 0; o < nt; o++ )
			{	
				int s_to = trans_list[o];
				
				dn_c[o] += (t_ab_next_clamped[s_from*nstates+s_to] - expec_c * t_ab[s_from*nstates+s_to]);
				dn_f[o] += (t_ab_next_free[s_from*nstates+s_to] - expec_f * t_ab[s_from*nstates+s_to]);
			}
		}
	
		
		for( int GE = 0; GE < GS; GE++ )
		{
			int nt;
			int trans_list[MAX_TRANSITIONS];
			getTransitionsFrom( transition_groups[TG][GE], trans_list, &nt, '.' );				
			int s_from = transition_groups[TG][GE];

			for( int o = 0; o < nt; o++ )
			{
				int s_to = trans_list[o];
				double tval = expec_c_save * t_ab[s_from*nstates+s_to] + alpha * (dn_c[o] - dn_f[o]);
				if( !( tval >= 0) )
					tval = 1e-20;
				t_ab[s_from*nstates+s_to] = tval;
			}
		}
	}
#endif

	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		int a_state = 0;

		int GS = emission_groups_len[EG];
				
		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];

		double expec_c_save = 1e-20;

		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
		
		for( int GE = 0; GE < GS; GE++ )
		{
			double expec_c = eps;
			double expec_f = eps;

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{	
				int s = emission_groups[EG][GE];
				a_state = s;

				expec_c_save += p_obs_next_clamped[o*nstates+s];
				
				expec_c += p_obs_next_clamped[o*nstates+s];
				expec_f += p_obs_next_free[o*nstates+s];			
			}
	
			int a_state = 0;

			int s = emission_groups[EG][GE];
		
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{	
				dn_c[o] += (p_obs_next_clamped[o*nstates+s] - expec_c * p_obs[o*nstates+s]);
				dn_f[o] += (p_obs_next_free[o*nstates+s]    - expec_f * p_obs[o*nstates+s]);
			}
		}
	
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				double tval = expec_c_save * p_obs[o*nstates+s] + alpha * (dn_c[o] - dn_f[o]);
				if( !( tval >= 0 ) )
					tval = 1e-20;
				p_obs[o*nstates+s] = tval;
			}
		}
	}


	Normalize();
	addNoise(A);	
	Normalize();

	squared_change = 0;

	for( int s = 0; s < nstates; s++ )
	{
//	for( int t = 0; t < nstates; t++ )
//	{
//		squared_change += (t_ab[s*nstates+t] - t_ab_prev[s*nstates+t])*(t_ab[s*nstates+t] - t_ab_prev[s*nstates+t]);
//	}
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
		squared_change += (p_obs[o*nstates+s] - p_obs_prev[o*nstates+s])*(p_obs[o*nstates+s] - p_obs_prev[o*nstates+s]);
	}

//	if(taskid == 0) printf("Squared change: %.14le\n", squared_change );

	return squared_change;
}

void HMMSystem::addNoise( double A )
{
	double t_sum = 0;
	double o_sum = 0;

	for( int s = 0; s < nstates; s++ )
	{
		int ntrans;
		int state_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, state_list, &ntrans, '.' );
		double ptot = 0;

		for( int xs = 0; xs < ntrans; xs++ )
			t_sum += t_ab_next_free[s*nstates+state_list[xs]];

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			o_sum += p_obs_next_free[s+nstates*o];
	}

	for( int s = 0; s < nstates; s++ )
	{	
		int ntrans;
		int state_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, state_list, &ntrans, '.' );
		double ptot = 0;

		for( int xs = 0; xs < ntrans; xs++ )
			t_ab_next_free[s*nstates+state_list[xs]] += A * noise() * t_sum / nstates;
	}
	
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
	{
		for( int EG = 0; EG < nEmissionGroups; EG++ )
		{
			double del = 0;
			int a_state = 0;

			int GS = emission_groups_len[EG];
	
			double n = A * noise() * o_sum / (globalN_OBSERVABLES*nstates);

			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				a_state = s;
			
				p_obs_next_free[o*nstates+s] += n;
			}
		}
	}	
}

void HMMSystem::Normalize( void )
{
	double p_init_tot = 0;

	for( int s = 0; s < nstates; s++ )
		p_init_tot += init_p[s];
	for( int s = 0; s < nstates; s++ )
		init_p[s] /= p_init_tot;

	for( int s = 0; s < nstates; s++ )
	{
		int ntrans;
		int state_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, state_list, &ntrans, '.' );
		double ptot = 0;
		
		for( int xs = 0; xs < ntrans; xs++ )
			t_ab[s*nstates+state_list[xs]] += eps;

		for( int xs = 0; xs < ntrans; xs++ )
			ptot += t_ab[s*nstates+state_list[xs]];

		if( my_isnan( ptot ) )
		{
//			printf("Trouble, ptot is nan for state %s.\n", states[s]->label );
		}

		if( fabs(ptot) < 1e-30 )
		{
//			printf("Warning, state '%s' never accessed.\n",
//				states[s]->label );
			ptot = 1.0;
		}

		for( int xs = 0; xs < ntrans; xs++ )
			t_ab[s*nstates+state_list[xs]] /= ptot;

		ptot =0;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			p_obs[o*nstates+s] += eps;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			ptot += p_obs[o*nstates+s];
		if( my_isnan( ptot ) )
		{
//			printf("Trouble, ptot is nan for state %s.\n", states[s]->label );
		}

		if( fabs(ptot) < 1e-30 )
			ptot = 1.0;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			p_obs[o*nstates+s] /= ptot;

	}

}

void HMMSystem::checkNormalization( void )
{
	int failure = 0;
	for( int s = 0; s < nstates; s++ )
	{
		int ntrans;
		int state_list[MAX_TRANSITIONS];

		getTransitionsFrom( s, state_list, &ntrans, '.' );
		double ptot = 0;
		for( int xs = 0; xs < ntrans; xs++ )
			ptot += t_ab[s*nstates+state_list[xs]];

		if( fabs(ptot-1.0) > 1e-9 || my_isnan( ptot)  )
		{
			printf("Warning, state '%s' never accessed.\n", states[s]->label );
			printf("NORMALIZATION TRANSITION FAILURE: %.14le from state %s.\n", ptot, states[s]->label );
			failure = 1;
			
		}
		ptot =0;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			ptot += p_obs[o*nstates+s];

		if( fabs(ptot-1.0) > 1e-9 || my_isnan(ptot) )
		{
			printf("NORMALIZATION EMISSION FAILURE: %.14le from state %s.\n", ptot, states[s]->label );
			failure = 1;
		}

	}
}


#define MAX_FILE_TRANSITIONS	100
#define MAX_FILE_OBSERVABLES	100

struct HMMSystemFileHeader
{
	int nstates;
	int interpreter;
	int model_id;
};

struct HMMStateFileEntry
{
	int file_id;
	char label[256];
	double transmission[MAX_FILE_TRANSITIONS];
	double emission[MAX_FILE_OBSERVABLES];
};




void HMMSystem::SaveModel( const char *fileName )
{	
	int taskid=0;
	int ierr;	

#ifdef PARALLEL
        if (ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid)) {
          exit(1);
        }
#endif
	if( taskid != 0 )
		return;

	FILE *theFile = fopen(fileName, "w");

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", fileName	);
		exit(1);
	}

	HMMSystemFileHeader theHeader;

	theHeader.nstates = nstates;
	theHeader.interpreter = 0;
	theHeader.model_id = 0;

	fwrite( &theHeader, sizeof(HMMSystemFileHeader), 1, theFile);

	for( int s = 0; s < nstates; s++ )
	{
		HMMStateFileEntry theEntry;

		strcpy( theEntry.label, states[s]->label );
		
		int trans_list[MAX_TRANSITIONS];
		int ntrans;
	
		getTransitionsFrom( s, trans_list, &ntrans, '.' );

		memset( theEntry.transmission, 0, sizeof(double) * MAX_FILE_TRANSITIONS ); 
		memset( theEntry.emission, 0, sizeof(double) * MAX_FILE_OBSERVABLES ); 

		for( int xs = 0; xs < ntrans; xs++ )
			theEntry.transmission[xs] = t_ab[s*nstates+trans_list[xs]];
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			theEntry.emission[o] = p_obs[o*nstates+s];

		theEntry.file_id = s;

		fwrite( &theEntry, sizeof(HMMStateFileEntry), 1, theFile );
	}

	fclose(theFile);
}

void HMMSystem::ReadModel( const char *fileName )
{	
	FILE *theFile = fopen(fileName, "r");

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", fileName	);
		exit(1);
	}

	HMMSystemFileHeader theHeader;

	fread( &theHeader,  sizeof(HMMSystemFileHeader), 1, theFile );

	if( theHeader.nstates != nstates )
	{
		printf("Number of states in saved parameters does not match states in model.\n");
		exit(1);
	}

	for( int ts = 0; ts < nstates; ts++ )
	{
		HMMStateFileEntry theEntry;
		fread( &theEntry, sizeof(HMMStateFileEntry), 1, theFile);

		int trans_list[MAX_TRANSITIONS];
		int ntrans;
		
		int s = theEntry.file_id;
	
		if( strcmp( theEntry.label, states[s]->label ) )
		{
			printf("State '%s' does not match state '%s'.\n", theEntry.label, states[s]->label );
			exit(1);
		}

		getTransitionsFrom( s, trans_list, &ntrans, '.' );

		for( int xs = 0; xs < ntrans; xs++ )
			t_ab[s*nstates+trans_list[xs]] = theEntry.transmission[xs] ;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			p_obs[o*nstates+s] = theEntry.emission[o];
	}

	fclose(theFile);
}

void HMMSystem::checkEGConsistency( void )
{
	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		int GS = emission_groups_len[EG];
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			double val_c = 0;
			double val_f = 0;

			double p_zero_c = p_obs[o*nstates+emission_groups[EG][0]];
			double p_zero_f = p_obs[o*nstates+emission_groups[EG][0]];

			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
	
				if( fabs(p_obs[o*nstates+s]-p_zero_c) > 1e-14 )
				{
					printf("EMISSION GROUP IS INCONSISTENT, BAILING %le.\n",p_obs_next_clamped[o*nstates+s]-p_zero_c);
					exit(1);
				}
				if( fabs(p_obs[o*nstates+s] - p_zero_f ) > 1e-14 )
				{
					printf("EMISSION GROUP IS INCONSISTENT, BAILING %le.\n",p_obs_next_free[o*nstates+s]-p_zero_f);
					exit(1);
				}
			}
		}
	}
}

/*
void HMMSystem::setGradient( double *dw )
{
#ifdef DO_TRANSITION_GROUPS
	printf("setGradient hasn't been coded for transition groups.\n");
	exit(1);
#endif

	int taskid;
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		int GS = emission_groups_len[EG];
			
		if( group_is_fixed[EG] ) 
		{
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				dw[1+s*globalN_OBSERVABLES+o] = 0;
			}

			continue;
		}


		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];

		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
	
		for( int GE = 0; GE < GS; GE++ )
		{
			double expec_c = eps;
			double expec_f = eps;
				
			int s = emission_groups[EG][GE];

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				expec_c += p_obs_next_clamped[o*nstates+s];
				expec_f += p_obs_next_free[o*nstates+s];			
			}

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				int a_state = 0;
	
				if( s == 0 && o == 0 )
				{
					printf("clamped got: %lf clamped expected: %lf\n", p_obs_next_clamped[o*nstates+s],
											  expec_c * p_obs[o*nstates+s] );
					printf("free got: %lf free expected: %lf\n", p_obs_next_free[o*nstates+s],
											  expec_f * p_obs[o*nstates+s] );
				}
			
				dn_c[o] += (p_obs_next_clamped[o*nstates+s] - expec_c * p_obs[o*nstates+s]);
				dn_f[o] += (p_obs_next_free[o*nstates+s] - expec_f * p_obs[o*nstates+s]);
			}
		}

			
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			double old_w = log( p_obs[o*nstates+s] );
			dw[1+s*globalN_OBSERVABLES+o] = dn_c[o] - dn_f[o];
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		double expec_c = eps;
		double expec_f = eps;

		for( int t = 0; t < nstates; t++ )
			expec_c += t_ab_next_clamped[s*nstates+t];
		for( int t = 0; t < nstates; t++ )
			expec_f += t_ab_next_free[s*nstates+t];

		for( int t = 0; t < nstates; t++ )
		{
			double old_p = t_ab[s*nstates+t];
			double dn_c = (t_ab_next_clamped[s*nstates+t] - expec_c * t_ab[s*nstates+t]);
			double dn_f = (t_ab_next_free[s*nstates+t] - expec_f * t_ab[s*nstates+t]);
	
			double del_w = dn_c - dn_f;

			dw[1+nstates*globalN_OBSERVABLES+s*nstates+t] = del_w;
		}
	}
}
*/

void HMMSystem::CommunicateGradients( double *dw )
{
#ifdef DO_TRANSITION_GROUPS
	printf("Communicate gradients doesn't work for groups.\n");
#endif

	int nprocs=1;

#ifdef PARALLEL
       	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
	double *from_vec = (double *)malloc( sizeof(double) * nstates * (nstates + globalN_OBSERVABLES) );
	int tp = 0;

	for( int x = 0; x < nstates; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
			from_vec[tp++] = dw[1+emission_groups[x][0]*globalN_OBSERVABLES+y];
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
			from_vec[tp++] = dw[1+globalN_OBSERVABLES*nstates+s*nstates+trans_list[x]];
	}

	double *to_vec = (double *)malloc( sizeof(double) * tp );
	memset( to_vec, 0, sizeof(double) * tp );

#ifdef PARALLEL
	MPI_Allreduce( from_vec, to_vec, tp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	memcpy( to_vec, from_vec, tp * sizeof(double) );
#endif
	tp = 0;

	memset( dw, 0, sizeof(double) * nstates * (nstates + globalN_OBSERVABLES) );

	for( int x = 0; x < nEmissionGroups; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
			for( int z = 0; z < emission_groups_len[x]; z++ )
				dw[1+y+globalN_OBSERVABLES*emission_groups[x][z]] = to_vec[tp];
			tp++;
		}
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
			dw[1+globalN_OBSERVABLES*nstates+s*nstates+trans_list[x]] = to_vec[tp++];
	}
	
	free(from_vec);
	free(to_vec);
}

/*
void HMMSystem::CommunicatePMats( void )
{
	return;
	int nprocs;
       	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	double *from_vec = (double *)malloc( sizeof(double) * nstates * (nstates + globalN_OBSERVABLES) );
	int tp = 0;

	for( int x = 0; x < nEmissionGroups; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
			from_vec[tp] = 0;
			for( int z = 0; z < emission_groups_len[x]; z++ )
				from_vec[tp] += p_obs_next_clamped[y*nstates+emission_groups[x][z]];
			from_vec[tp] /= emission_groups_len[x];
			tp++;
			from_vec[tp] = 0;
			for( int z = 0; z < emission_groups_len[x]; z++ )
				from_vec[tp] += p_obs_next_free[y*nstates+emission_groups[x][z]];
			from_vec[tp] /= emission_groups_len[x];
			tp++;
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
		{
			from_vec[tp++] = t_ab_next_clamped[s*nstates+trans_list[x]];
			from_vec[tp++] = t_ab_next_free[s*nstates+trans_list[x]];
		}
	}

	double *to_vec = (double *)malloc( sizeof(double) * tp );
	memset( to_vec, 0, sizeof(double) * tp );

	MPI_Allreduce( from_vec, to_vec, tp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	tp = 0;

	for( int x = 0; x < nEmissionGroups; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
			for( int z = 0; z < emission_groups_len[x]; z++ )
			{
				p_obs_next_clamped[y*nstates+emission_groups[x][z]] = to_vec[tp];
			}
			tp++;
			for( int z = 0; z < emission_groups_len[x]; z++ )
			{
				p_obs_next_free[y*nstates+emission_groups[x][z]] = to_vec[tp];
			}
			tp++;
		}
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
		{
			t_ab_next_clamped[s*nstates+trans_list[x]] = to_vec[tp++];
			t_ab_next_free[s*nstates+trans_list[x]] = to_vec[tp++];
		}
	}
	
	free(from_vec);
	free(to_vec);

}
*/

void HMMSystem::CommunicatePMats( void )
{
	timeval thet1;
	timeval thet2;

	gettimeofday( &thet1, NULL );	

	double *recv_init_p_next = (double *)malloc( sizeof(double) * nstates );

	memset( recv_init_p_next, 0, sizeof(double) * nstates);
#ifdef PARALLEL
	MPI_Allreduce( init_p_next_free, recv_init_p_next, nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	memcpy( recv_init_p_next, init_p_next_free, sizeof(double) * nstates );
#endif
	memcpy( init_p_next_free, recv_init_p_next, sizeof(double) * nstates );
	
	memset( recv_init_p_next, 0, sizeof(double) * nstates);
#ifdef PARALLEL
	MPI_Allreduce( init_p_next_clamped, recv_init_p_next, nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	memcpy( recv_init_p_next, init_p_next_clamped, sizeof(double) * nstates );
#endif
	memcpy( init_p_next_clamped, recv_init_p_next, sizeof(double) * nstates );

	free( recv_init_p_next );
				double *recv_t_ab_next = (double *)malloc( sizeof(double) * nstates * nstates ); 
				double *recv_p_obs_next = (double *)malloc( sizeof(double) * nstates * globalN_OBSERVABLES ); 
				memset( recv_t_ab_next, 0, sizeof(double) * nstates * nstates );
				memset( recv_p_obs_next, 0, sizeof(double) * nstates * globalN_OBSERVABLES );

				int ar_err;

#ifdef PARALLEL
				MPI_Allreduce( t_ab_next_free, recv_t_ab_next, nstates*nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( p_obs_next_free, recv_p_obs_next, nstates*globalN_OBSERVABLES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
				memcpy( recv_t_ab_next, t_ab_next_free, nstates * nstates * sizeof(double) );
				memcpy( recv_p_obs_next, p_obs_next_free, globalN_OBSERVABLES * nstates * sizeof(double) );
#endif			
				memcpy( t_ab_next_free, recv_t_ab_next, nstates*nstates*sizeof(double) );
				memcpy( p_obs_next_free, recv_p_obs_next, globalN_OBSERVABLES*nstates*sizeof(double) );
				
				memset( recv_t_ab_next, 0, sizeof(double) * nstates * nstates );
				memset( recv_p_obs_next, 0, sizeof(double) * nstates * globalN_OBSERVABLES );
				
#ifdef PARALLEL
				MPI_Allreduce( t_ab_next_clamped, recv_t_ab_next, nstates*nstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
				MPI_Allreduce( p_obs_next_clamped, recv_p_obs_next, nstates*globalN_OBSERVABLES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
				memcpy( recv_t_ab_next, t_ab_next_clamped, nstates * nstates * sizeof(double) );
				memcpy( recv_p_obs_next, p_obs_next_clamped, globalN_OBSERVABLES * nstates * sizeof(double) );
#endif				
				memcpy( t_ab_next_clamped, recv_t_ab_next, nstates*nstates*sizeof(double) );
				memcpy( p_obs_next_clamped, recv_p_obs_next, globalN_OBSERVABLES*nstates*sizeof(double) );

				free( recv_t_ab_next);
				free( recv_p_obs_next);
	
	gettimeofday( &thet2, NULL );	

	double secs = thet2.tv_sec - thet1.tv_sec + (1e-6) * (thet2.tv_usec - thet1.tv_usec);

	static int mpi_it = 0;

//	if( mpi_it % 100 == 0 )
//		printf("MPI communication wall time: %le seconds.\n",secs );
	mpi_it++;

}

void HMMSystem::CommunicateShortPMats( void )
{
	int nprocs=1;
#ifdef PARALLEL
       	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
	double *from_vec = (double *)malloc( sizeof(double) * (2 * nstates * (nstates + globalN_OBSERVABLES)  + nstates * 2) );
	int tp = 0;

	for( int x = 0; x < nstates; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
			from_vec[tp] = 0;
				from_vec[tp] += p_obs_next_clamped[y*nstates+x];
			tp++;
			from_vec[tp] = 0;
				from_vec[tp] += p_obs_next_free[y*nstates+x];
			tp++;
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
		{
			from_vec[tp++] = t_ab_next_clamped[s*nstates+trans_list[x]];
			from_vec[tp++] = t_ab_next_free[s*nstates+trans_list[x]];
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		from_vec[tp++] = init_p_next_clamped[s];
		from_vec[tp++] = init_p_next_free[s];
	}
	double *to_vec = (double *)malloc( sizeof(double) * tp );
	memset( to_vec, 0, sizeof(double) * tp );

	timeval thet1;
	timeval thet2;

	gettimeofday( &thet1, NULL );	
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Allreduce( from_vec, to_vec, tp, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	memcpy( to_vec, from_vec, sizeof(double) * tp );
#endif	
	gettimeofday( &thet2, NULL );	

	double secs = thet2.tv_sec - thet1.tv_sec + (1e-6) * (thet2.tv_usec - thet1.tv_usec);

//	printf("MPI communication wall time: %le seconds.\n",secs );

	tp = 0;

	for( int x = 0; x < nstates; x++ )
	{
		for( int y = 0; y < globalN_OBSERVABLES; y++ )
		{
				p_obs_next_clamped[y*nstates+x] = to_vec[tp];
			tp++;
				p_obs_next_free[y*nstates+x] = to_vec[tp];
			tp++;
		}
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		getTransitionsFrom( s, trans_list, &ntrans, '.' );
			
		for( int x = 0; x < ntrans; x++ )
		{
			t_ab_next_clamped[s*nstates+trans_list[x]] = to_vec[tp++];
			t_ab_next_free[s*nstates+trans_list[x]] = to_vec[tp++];
		}
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		init_p_next_clamped[s] = to_vec[tp++];
		init_p_next_free[s] = to_vec[tp++];
	}
	
	free(from_vec);
	free(to_vec);

}

void HMMSystem::ReportEmissionGroups( int trans)
{
/*
        for( int c = 0; c < 30; c++ )
                printf(" ");
        for( int o = 0; o < globalN_OBSERVABLES; o++ )
        {
                if( trans == INTERPRETER_AA )
                {
                        printf("%c", all_res[o] );
                        for( int c = 1; c < 16; c++ )
                                printf(" ");
                }
                else if ( trans == INTERPRETER_SUPER )
                        printSTGroupName(o);
		else if( trans >= INTERPRETER_SCHEMES && trans <= INTERPRETER_SCHEMES+9 )
		{
			printSchemeGroupName(o, trans-INTERPRETER_SCHEMES);
		}
        }
        printf("\n");
        for( int EG = 0; EG < nEmissionGroups; EG++ )
        {
                int astate = emission_groups[EG][0];
                char tbuf[30];
                sprintf(tbuf, "%s ", states[astate]->label );
                printf("%s", tbuf );
                for( int x = strlen(tbuf); x < 30; x++ )
                        printf(" ");

                for( int o = 0; o < globalN_OBSERVABLES; o++ )
                {
                        sprintf(tbuf, "%le ", p_obs[o*nstates+astate] );
                        printf("%s", tbuf );
                        for( int x = strlen(tbuf); x < 16; x++ )
                                printf(" ");
                }
                printf("\n");
        }

	for( int s = 0; s < nstates; s++ )
	{
		printf("State: %s\n", states[s]->label );
	
			int state_list[MAX_TRANSITIONS];

			int ntrans;

			getTransitionsFrom( s, state_list, &ntrans, '.' );
		
		for( int xs = 0; xs < ntrans; xs++ )
		{
			int	s2 = state_list[xs];
			printf("\tto state %s, p: %le\n", states[state_list[xs]]->label, t_ab[s*nstates+s2] );
			
		}
	}
*/
}

void HMMSystem::setGradient( double *dw )
{
	double f_obs[nstates];
	double f_t_ab[nstates];

//	setFArray( f_in, f_obs, f_t_ab );

	int taskid=0;
#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
	double init_sum_f = 0;

	for( int s = 0; s < nstates; s++ )
	{
		init_sum_f += init_p_next_free[s];
		// really, these should be the same.
	}

	for( int s = 0; s < nstates; s++ )
	{
		double del_f = init_p_next_free[s] - init_sum_f * init_p[s];

		dw[1+globalN_OBSERVABLES*nstates+nstates*nstates+s] = del_f;
	}

	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{

		//double fp = 1.0 - f_obs[emission_groups[EG][0]];
		//double f = f_obs[emission_groups[EG][0]];

		double f =1;
		double fp = 0;

		int GS = emission_groups_len[EG];
		
		if( group_is_fixed[EG] ) 
		{
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				dw[1+s*globalN_OBSERVABLES+o] = 0;
			}

			continue;
		}

		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];

		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
	
		for( int GE = 0; GE < GS; GE++ )
		{
			/* intermediates we need */

			double n_c = eps;
			double n_f = eps;
			
			double n_p_t_c = eps;
			double n_p_t_f = eps;				

			double n_t_c = eps;
			double n_t_f = eps;

			double n_p_c = eps;
			double n_p_f = eps;

			int s = emission_groups[EG][GE];

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				n_c += p_obs_next_clamped[o*nstates+s];
				n_f += p_obs_next_free[o*nstates+s];

				n_p_t_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_p_t_f +=    p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);

				n_t_c += p_obs_next_clamped[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_t_f += p_obs_next_free[o*nstates+s] / (p_obs[o*nstates+s]+eps);
							
				n_p_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
				n_p_f += p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
			}

			/* for fp-> zero, we must get the other gradient, clearly true... f goes to zero is dumb */

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				int a_state = 0;
			
				double tval =

					      p_obs_next_clamped[o*nstates+s] 
					    - n_c * p_obs[o*nstates+s]; 
				dn_c[o] += tval;				
	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_c[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif

				tval = p_obs_next_free[o*nstates+s] 
					    - n_f * p_obs[o*nstates+s]; 
				dn_f[o] += tval;
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_f[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			}
		}

/*	
#ifdef DIVIDE_GS
			dn_c /= GS;
			dn_f /= GS;
#endif		
*/	
			
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			double old_w = log( p_obs[o*nstates+s] );
			dw[1+s*globalN_OBSERVABLES+o] = dn_f[o];

		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		double f_t = 1.0;
		double fp_t = 0;

		double n_c = eps;
		double n_f = eps;
		
		double n_p_t_c = eps;
		double n_p_t_f = eps;				

		double n_t_c = eps;
		double n_t_f = eps;

		double n_p_c = eps;
		double n_p_f = eps;


		for( int t = 0; t < nstates; t++ )
		{
			n_c += t_ab_next_clamped[s*nstates+t];
			n_f += t_ab_next_free[s*nstates+t];

			n_p_t_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_p_t_f +=    t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);

			n_t_c += t_ab_next_clamped[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_t_f += t_ab_next_free[s*nstates+t] / (t_ab[s*nstates+t]+eps);
						
			n_p_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
			n_p_f += t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
		}

		double dn_c[nstates];
		double dn_f[nstates];

		memset( dn_c, 0, sizeof(double) * nstates );
		memset( dn_f, 0, sizeof(double) * nstates );

		for( int t = 0; t < nstates; t++ )
		{
			int a_state = 0;
		
			double tval =
			              t_ab_next_clamped[s*nstates+t] 
				    - n_c * t_ab[s*nstates+t]; 
			dn_c[t] += tval;	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_c[t] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			tval =
				      t_ab_next_free[s*nstates+t] 
				    - n_f * t_ab[s*nstates+t]; 
			dn_f[t] += tval;	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_f[t] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
		}

		for( int t = 0; t < nstates; t++ )
		{
			double del_w = dn_f[t];

			dw[1+nstates*globalN_OBSERVABLES+s*nstates+t] = del_w;
		}
	}
}


void HMMSystem::setFractionalGradientSplit( double *dw, double f_in)
{
	double f_obs[nstates];
	double f_t_ab[nstates];

	setFArray( f_in, f_obs, f_t_ab );

	int taskid=0;
#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
	double init_sum_c = 0;
	double init_sum_f = 0;

	for( int s = 0; s < nstates; s++ )
	{
		init_sum_f += init_p_next_free[s];
		init_sum_c += init_p_next_clamped[s];
		// really, these should be the same.
	}

	for( int s = 0; s < nstates; s++ )
	{
		double del_c = init_p_next_clamped[s] - init_sum_c * init_p[s];
		double del_f = init_p_next_free[s] - init_sum_f * init_p[s];

		dw[1+globalN_OBSERVABLES*nstates+nstates*nstates+s] = (del_c - del_f);
	}

	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{

		double fp = 1.0 - f_obs[emission_groups[EG][0]];
		double f = f_obs[emission_groups[EG][0]];

		int GS = emission_groups_len[EG];
		
		if( group_is_fixed[EG] ) 
		{
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				dw[1+s*globalN_OBSERVABLES+o] = 0;
			}

			continue;
		}

		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];

		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
	
		for( int GE = 0; GE < GS; GE++ )
		{
			/* intermediates we need */

			double n_c = eps;
			double n_f = eps;
			
			double n_p_t_c = eps;
			double n_p_t_f = eps;				

			double n_t_c = eps;
			double n_t_f = eps;

			double n_p_c = eps;
			double n_p_f = eps;

			int s = emission_groups[EG][GE];

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				n_c += p_obs_next_clamped[o*nstates+s];
				n_f += p_obs_next_free[o*nstates+s];

				n_p_t_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_p_t_f +=    p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);

				n_t_c += p_obs_next_clamped[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_t_f += p_obs_next_free[o*nstates+s] / (p_obs[o*nstates+s]+eps);
							
				n_p_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
				n_p_f += p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
			}

			/* for fp-> zero, we must get the other gradient, clearly true... f goes to zero is dumb */

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				int a_state = 0;
			
				double tval =

					      p_obs_next_clamped[o*nstates+s] 
					    - p_obs_next_clamped[o*nstates+s] * (p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps)) * fp 
					    - n_c * p_obs[o*nstates+s] / f 
					    - (fp*n_p_t_c) * (fp*p_obs_BaumWelch[o*nstates+s]) * (1.0 / f)
//					    + n_p_c * ( fp / f)
//					    + n_c * p_obs_BaumWelch[o*nstates+s] * ( fp / f)
					    + n_c * p_obs_BaumWelch[o*nstates+s] * (fp / f)
					    + n_p_t_c * p_obs[o*nstates+s] * ( fp / f);
				dn_c[o] += tval;				
	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_c[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif

				tval = p_obs_next_free[o*nstates+s] 
					    - p_obs_next_free[o*nstates+s] * (p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps)) * fp 
					    - n_f * p_obs[o*nstates+s] / f 
					    - (fp*n_p_t_f) * (fp*p_obs_BaumWelch[o*nstates+s]) * (1.0 / f)
					    //+ n_p_f * ( fp / f)
					    //+ n_f * p_obs_BaumWelch[o*nstates+s] * ( fp / f)
					    + n_f * p_obs_BaumWelch[o*nstates+s] * (fp / f)
					    + n_p_t_f * p_obs[o*nstates+s] * ( fp / f);
				dn_f[o] += tval;
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_f[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			}
		}

/*	
#ifdef DIVIDE_GS
			dn_c /= GS;
			dn_f /= GS;
#endif		
*/	
			
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			double old_w = log( p_obs[o*nstates+s] );
			dw[1+s*globalN_OBSERVABLES+o] = dn_c[o] - dn_f[o];

		}
	}

#ifdef DO_TRANSITION_GROUPS	
	for( int TG = 0; TG < nTransitionGroups; TG++ )
	{
		double fp = 1.0 - f_obs[transition_groups[TG][0]];
		double f = f_obs[transition_groups[TG][0]];
		int GS = transition_groups_len[TG];
		
		int nt;
		int trans_list[MAX_TRANSITIONS];
	
		getTransitionsFrom( transition_groups[TG][0], trans_list, &nt, '.' );
		
		double dn_c[nt];
		double dn_f[nt];

		memset( dn_c, 0, sizeof(double) * nt );
		memset( dn_f, 0, sizeof(double) * nt );
	
		for( int GE = 0; GE < GS; GE++ )
		{
			/* intermediates we need */

			double n_c = eps;
			double n_f = eps;
			
			double n_p_t_c = eps;
			double n_p_t_f = eps;				

			double n_t_c = eps;
			double n_t_f = eps;

			double n_p_c = eps;
			double n_p_f = eps;

			int s_from = transition_groups[TG][GE];

			int nt;
			int trans_list[MAX_TRANSITIONS];

			getTransitionsFrom( s_from, trans_list, &nt, '.' );

			for( int o = 0; o < nt; o++ )
			{
				int s_to = trans_list[o];

				n_c += t_ab_next_clamped[s_from*nstates+s_to];
				n_f += t_ab_next_free[s_from*nstates+s_to];

				n_p_t_c += t_ab_next_clamped[s_from*nstates+s_to] * t_ab_BaumWelch[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps);
				n_p_t_f +=    t_ab_next_free[s_from*nstates+s_to] * t_ab_BaumWelch[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps);

				n_t_c += t_ab_next_clamped[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps);
				n_t_f += t_ab_next_free[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps);
							
				n_p_c += t_ab_next_clamped[s_from*nstates+s_to] * t_ab_BaumWelch[s_from*nstates+s_to];			
				n_p_f += t_ab_next_free[s_from*nstates+s_to] * t_ab_BaumWelch[s_from*nstates+s_to];			
			}

			/* for fp-> zero, we must get the other gradient, clearly true... f goes to zero is dumb */
			
			for( int o = 0; o < nt; o++ )
			{
				int a_state = 0;
			
				double tval ;
				int s_to = trans_list[o];
				tval = t_ab_next_clamped[s_from*nstates+s_to] 
					    - t_ab_next_clamped[s_from*nstates+s_to] * (t_ab_BaumWelch[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps)) * fp 
					    - n_c * t_ab[s_from*nstates+s_to] / f 
					    - (fp*n_p_t_c) * (fp*t_ab_BaumWelch[s_from*nstates+s_to]) * (1.0 / f)
					    //+ n_p_f * ( fp / f)
					    //+ n_f * t_ab_BaumWelch[s_from*nstates+s_to] * ( fp / f)
					    + n_c * t_ab_BaumWelch[s_from*nstates+s_to] * (fp / f)
					    + n_p_t_c * t_ab[s_from*nstates+s_to] * ( fp / f);
				dn_c[o] += tval;
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_c[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			}

			for( int o = 0; o < nt; o++ )
			{
				int a_state = 0;
			
				double tval ;
				int s_to = trans_list[o];
				tval = t_ab_next_free[s_from*nstates+s_to] 
					    - t_ab_next_free[s_from*nstates+s_to] * (t_ab_BaumWelch[s_from*nstates+s_to] / (t_ab[s_from*nstates+s_to]+eps)) * fp 
					    - n_f * t_ab[s_from*nstates+s_to] / f 
					    - (fp*n_p_t_f) * (fp*t_ab_BaumWelch[s_from*nstates+s_to]) * (1.0 / f)
					    //+ n_p_f * ( fp / f)
					    //+ n_f * t_ab_BaumWelch[s_from*nstates+s_to] * ( fp / f)
					    + n_f * t_ab_BaumWelch[s_from*nstates+s_to] * (fp / f)
					    + n_p_t_f * t_ab[s_from*nstates+s_to] * ( fp / f);
				dn_f[o] += tval;
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_f[o] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			}
		}

/*	
#ifdef DIVIDE_GS
			dn_c /= GS;
			dn_f /= GS;
#endif		
*/	
			
		for( int o = 0; o < nt; o++ )
		for( int GE = 0; GE < GS; GE++ )
		{
			int ntt;
			int trans_list[MAX_TRANSITIONS];
			int s_from = transition_groups[TG][GE];
			getTransitionsFrom( s_from, trans_list, &ntt, '.' );
			int s_to = trans_list[o];
			dw[1+globalN_OBSERVABLES*nstates+s_from*nstates+s_to] = dn_c[o] - dn_f[o];
		}
	}

#else // Don't do transition groups
	for( int s = 0; s < nstates; s++ )
	{
		double fp_t = 1.0 - f_t_ab[s];
		double f_t = f_t_ab[s];

		double n_c = eps;
		double n_f = eps;
		
		double n_p_t_c = eps;
		double n_p_t_f = eps;				

		double n_t_c = eps;
		double n_t_f = eps;

		double n_p_c = eps;
		double n_p_f = eps;


		for( int t = 0; t < nstates; t++ )
		{
			n_c += t_ab_next_clamped[s*nstates+t];
			n_f += t_ab_next_free[s*nstates+t];

			n_p_t_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_p_t_f +=    t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);

			n_t_c += t_ab_next_clamped[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_t_f += t_ab_next_free[s*nstates+t] / (t_ab[s*nstates+t]+eps);
						
			n_p_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
			n_p_f += t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
		}

		double dn_c[nstates];
		double dn_f[nstates];

		memset( dn_c, 0, sizeof(double) * nstates );
		memset( dn_f, 0, sizeof(double) * nstates );

		for( int t = 0; t < nstates; t++ )
		{
			int a_state = 0;
		
			double tval =
			              t_ab_next_clamped[s*nstates+t] 
			            - t_ab_next_clamped[s*nstates+t] * (t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps)) * fp_t 
				    - n_c * t_ab[s*nstates+t] / f_t 
				    - (fp_t*n_p_t_c) * (fp_t*t_ab_BaumWelch[s*nstates+t]) * (1.0 / f_t)
//				    + n_p_c * ( fp / f)
				   // + n_c * t_ab_BaumWelch[s*nstates+t] * ( fp / f)
			            + n_c * t_ab_BaumWelch[s*nstates+t] * (fp_t / f_t)
				    + n_p_t_c * t_ab[s*nstates+t] * ( fp_t / f_t);
			dn_c[t] += tval;	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_c[t] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
			tval =
				      t_ab_next_free[s*nstates+t] 
			            - t_ab_next_free[s*nstates+t] * (t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps)) * fp_t 
				    - n_f * t_ab[s*nstates+t] / f_t 
				    - (fp_t*n_p_t_f) * (fp_t*t_ab_BaumWelch[s*nstates+t]) * (1.0 / f_t)
//				    + n_p_f * ( fp / f)
//				    + n_f * t_ab_BaumWelch[s*nstates+t] * ( fp / f)
			            + n_f * t_ab_BaumWelch[s*nstates+t] * (fp_t / f_t)
				    + n_p_t_f * t_ab[s*nstates+t] * ( fp_t / f_t);
			dn_f[t] += tval;	
#ifdef GRAD_NANCHECK
				if( my_isnan( dn_f[t] ) )
				{
					printf("GRAD NAN.\n");
				}
#endif
		}

		for( int t = 0; t < nstates; t++ )
		{
			double del_w = dn_c[t] - dn_f[t];

			dw[1+nstates*globalN_OBSERVABLES+s*nstates+t] = del_w;
		}
	}
#endif
}

void HMMSystem::setFractionalGradient( double *dw, double f_in )
{
#ifdef DO_TRANSITION_GROUPS
	printf("setFractionalGradient not yet implemented with transition groups.\n");
	exit(1);
#endif

	int taskid=0;
#ifdef PARALLEL
       	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif
	double f_obs[nstates];
	double f_t_ab[nstates];

	setFArray( f_in, f_obs, f_t_ab );
	
	for( int EG = 0; EG < nEmissionGroups; EG++ )
	{
		double f = f_obs[emission_groups[EG][0]];
		double fp = (1.0-f);

		int GS = emission_groups_len[EG];
		
		if( group_is_fixed[EG] ) 
		{
			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			for( int GE = 0; GE < GS; GE++ )
			{
				int s = emission_groups[EG][GE];
				dw[1+s*globalN_OBSERVABLES+o] = 0;
			}

			continue;
		}

		double dn_c[globalN_OBSERVABLES];
		double dn_f[globalN_OBSERVABLES];

		memset( dn_c, 0, sizeof(double) * globalN_OBSERVABLES );
		memset( dn_f, 0, sizeof(double) * globalN_OBSERVABLES );
	
		for( int GE = 0; GE < GS; GE++ )
		{
			/* intermediates we need */

			double n_c = eps;
			double n_f = eps;
			
			double n_p_t_c = eps;
			double n_p_t_f = eps;				

			double n_t_c = eps;
			double n_t_f = eps;

			double n_p_c = eps;
			double n_p_f = eps;

			int s = emission_groups[EG][GE];

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				n_c += p_obs_next_clamped[o*nstates+s];
				n_f += p_obs_next_free[o*nstates+s];

				n_p_t_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_p_t_f +=    p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps);

				n_t_c += p_obs_next_clamped[o*nstates+s] / (p_obs[o*nstates+s]+eps);
				n_t_f += p_obs_next_free[o*nstates+s] / (p_obs[o*nstates+s]+eps);
							
				n_p_c += p_obs_next_clamped[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
				n_p_f += p_obs_next_free[o*nstates+s] * p_obs_BaumWelch[o*nstates+s];			
			}

			/* for fp-> zero, we must get the other gradient, clearly true... f goes to zero is dumb */

			for( int o = 0; o < globalN_OBSERVABLES; o++ )
			{
				int a_state = 0;
			
				double tval =

					      p_obs_next_clamped[o*nstates+s] 
					    - p_obs_next_clamped[o*nstates+s] * (p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps)) * fp 
					    - n_c * p_obs[o*nstates+s] / f 
					    - (fp*n_p_t_c) * (fp*p_obs_BaumWelch[o*nstates+s]) * (1.0 / f)
//					    + n_p_c * ( fp / f)
//					    + n_c * p_obs_BaumWelch[o*nstates+s] * ( fp / f)
					    + n_c * p_obs_BaumWelch[o*nstates+s] * (fp / f)
					    + n_p_t_c * p_obs[o*nstates+s] * ( fp / f);
				dn_c[o] += tval;				
				tval = p_obs_next_free[o*nstates+s] 
					    - p_obs_next_free[o*nstates+s] * (p_obs_BaumWelch[o*nstates+s] / (p_obs[o*nstates+s]+eps)) * fp 
					    - n_f * p_obs[o*nstates+s] / f 
					    - (fp*n_p_t_f) * (fp*p_obs_BaumWelch[o*nstates+s]) * (1.0 / f)
					    //+ n_p_f * ( fp / f)
					    //+ n_f * p_obs_BaumWelch[o*nstates+s] * ( fp / f)
					    + n_f * p_obs_BaumWelch[o*nstates+s] * (fp / f)
					    + n_p_t_f * p_obs[o*nstates+s] * ( fp / f);
				dn_f[o] += tval;
			}
		}

/*	
#ifdef DIVIDE_GS
			dn_c /= GS;
			dn_f /= GS;
#endif		
*/	
			
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		for( int GE = 0; GE < GS; GE++ )
		{
			int s = emission_groups[EG][GE];
			double old_w = log( p_obs[o*nstates+s] );
			dw[1+s*globalN_OBSERVABLES+o] = dn_c[o] - dn_f[o];
		}
	}

	for( int s = 0; s < nstates; s++ )
	{
		double f = f_t_ab[s];
		double fp = (1.0-f);

		double n_c = eps;
		double n_f = eps;
		
		double n_p_t_c = eps;
		double n_p_t_f = eps;				

		double n_t_c = eps;
		double n_t_f = eps;

		double n_p_c = eps;
		double n_p_f = eps;


		for( int t = 0; t < nstates; t++ )
		{
			n_c += t_ab_next_clamped[s*nstates+t];
			n_f += t_ab_next_free[s*nstates+t];

			n_p_t_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_p_t_f +=    t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps);

			n_t_c += t_ab_next_clamped[s*nstates+t] / (t_ab[s*nstates+t]+eps);
			n_t_f += t_ab_next_free[s*nstates+t] / (t_ab[s*nstates+t]+eps);
						
			n_p_c += t_ab_next_clamped[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
			n_p_f += t_ab_next_free[s*nstates+t] * t_ab_BaumWelch[s*nstates+t];			
		}

		double dn_c[nstates];
		double dn_f[nstates];

		memset( dn_c, 0, sizeof(double) * nstates );
		memset( dn_f, 0, sizeof(double) * nstates );

		for( int t = 0; t < nstates; t++ )
		{
			int a_state = 0;
		
			double tval =
			              t_ab_next_clamped[s*nstates+t] 
			            - t_ab_next_clamped[s*nstates+t] * (t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps)) * fp 
				    - n_c * t_ab[s*nstates+t] / f 
				    - (fp*n_p_t_c) * (fp*t_ab_BaumWelch[s*nstates+t]) * (1.0 / f)
//				    + n_p_c * ( fp / f)
				   // + n_c * t_ab_BaumWelch[s*nstates+t] * ( fp / f)
			            + n_c * t_ab_BaumWelch[s*nstates+t] * (fp / f)
				    + n_p_t_c * t_ab[s*nstates+t] * ( fp / f);
			dn_c[t] += tval;	
			tval =
				      t_ab_next_free[s*nstates+t] 
			            - t_ab_next_free[s*nstates+t] * (t_ab_BaumWelch[s*nstates+t] / (t_ab[s*nstates+t]+eps)) * fp 
				    - n_f * t_ab[s*nstates+t] / f 
				    - (fp*n_p_t_f) * (fp*t_ab_BaumWelch[s*nstates+t]) * (1.0 / f)
//				    + n_p_f * ( fp / f)
//				    + n_f * t_ab_BaumWelch[s*nstates+t] * ( fp / f)
			            + n_f * t_ab_BaumWelch[s*nstates+t] * (fp / f)
				    + n_p_t_f * t_ab[s*nstates+t] * ( fp / f);
			dn_f[t] += tval;	
		}

		for( int t = 0; t < nstates; t++ )
		{
			double del_w = dn_c[t] - dn_f[t];

			dw[1+nstates*globalN_OBSERVABLES+s*nstates+t] = del_w;
		}
	}
}

void HMMSystem::resetStatePopulations( int do_all )
{
	for( int s = 0; s < nstates; s++ )
		memset( states[s]->obs_population, 0, sizeof(double) * MAX_OBSERVABLES );
	if( do_all )
	{
		for( int s = 0; s < nstates; s++ )
			memset( states[s]->obs_population_copy, 0, sizeof(double) * MAX_OBSERVABLES );
	}
}

void HMMSystem::swapStatePopulations( void )
{
	for( int s = 0; s < nstates; s++ )
	{
		double copy[MAX_OBSERVABLES];
		memcpy( copy, states[s]->obs_population, sizeof(double) * MAX_OBSERVABLES );
		memcpy( states[s]->obs_population, states[s]->obs_population_copy, sizeof(double) * MAX_OBSERVABLES );
		memcpy( states[s]->obs_population_copy, copy, sizeof(double) * MAX_OBSERVABLES );
	}
}

void HMMSystem::addStatePopulations( void )
{
	for( int s = 0; s < nstates; s++ )
	{
		for( int p = 0; p < MAX_OBSERVABLES; p++ )
			states[s]->obs_population_copy[p] += states[s]->obs_population[p];
	}
}


void HMMSystem::printStatePopulations( int trans)
{
/*
        for( int c = 0; c < 30; c++ )
                printf(" ");
        for( int o = 0; o < globalN_OBSERVABLES; o++ )
        {
                if( trans == INTERPRETER_AA )
                {
                        printf("%c", all_res[o] );
                        for( int c = 1; c < 16; c++ )
                                printf(" ");
                }
                else if ( trans == INTERPRETER_SUPER )
                        printSTGroupName(o);
		else if( trans >= INTERPRETER_SCHEMES && trans <= INTERPRETER_SCHEMES+9 )
		{
			printSchemeGroupName(o, trans-INTERPRETER_SCHEMES);
		}
        }
	printf(" Ignored");
        printf("\n");
	
	double av_pop = 0;
	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o <= globalN_OBSERVABLES; o++ )
			av_pop += states[s]->obs_population[o];
	}	

	av_pop /= nstates;

	for( int s = 0; s < nstates; s++ )
	{
		double tot_pop = 0;

		printf("%s ", states[s]->label );
		
		for( int o = 0; o <= globalN_OBSERVABLES; o++ )
		{
			tot_pop += states[s]->obs_population[o];
			printf("%lf ", states[s]->obs_population[o] );
		}
		
		printf("\n");
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		double tot_pop = 0;
		for( int o = 0; o <= globalN_OBSERVABLES; o++ )
			tot_pop += states[s]->obs_population[o];
		if( tot_pop / av_pop < 0.1 )
			printf("State %s is underpopulated %lf/%lf.\n", states[s]->label, tot_pop, av_pop );	
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		double ilv_pop = states[s]->obs_population[1] + states[s]->obs_population[2] + states[s]->obs_population[3];
		double gast_pop = states[s]->obs_population[18] + states[s]->obs_population[17] + states[s]->obs_population[16] + states[s]->obs_population[15];

		if( gast_pop - ilv_pop > 300 && states[s]->theClass == 'M' || states[s]->theClass == 'I' )
			printf("State %s of type %c has GAST %lf and ILV only %lf.\n",				
				states[s]->label, states[s]->theClass, gast_pop, ilv_pop );
		if( ilv_pop - gast_pop > 300 && states[s]->theClass == 'M' || states[s]->theClass == 'I' )
			printf("State %s of type %c has GAST %lf and ILV only %lf.\n",				
				states[s]->label, states[s]->theClass, gast_pop, ilv_pop );
	}
*/

}

void HMMSystem::FreezeGlobular( void )
{
	for( int s = 0; s < nstates; s++ )	
	{
		if( !strncasecmp( states[s]->label, "glob", 4 ) )
			group_is_fixed[states[s]->emission_group] = 1;
	}
}

double HMMSystem::MapP0ObsFactor( double f, int s )
{
	double ptot = 0;
		
	for( int x = 0; x < globalN_OBSERVABLES; x++ )
	{	
		ptot += p_obs_BaumWelch[x*nstates+s];
	}

	if( fabs(ptot-1.0) < 1e-5 )
	{
//		printf("State %s is using p0 for its emissions (%le).\n", states[s]->label, f);
		return f;
	}
	return 1.0;
}	


double HMMSystem::MapP0TFactor( double f, int s )
{
	double ptot = 0;
		
	for( int x = 0; x < nstates; x++ )
		ptot += t_ab_BaumWelch[s*nstates+x];

	if( fabs(ptot-1.0) < 1e-5 )
		return f;
	return 1.0;
}	

void HMMSystem::setFArray( double f, double *obs_f, double *t_f )
{
	for( int s = 0; s < nstates; s++ )
	{
		obs_f[s] = MapP0ObsFactor( f, s );
//		printf("state %d f mapped to %lf.\n", s, obs_f[s] );
	}
	for( int s = 0; s < nstates; s++ )
	{
		t_f[s] = MapP0TFactor( f, s );
		if( t_f[s] != 1.0 )	
		{
//			printf("State %s is using p0 for its transitions (%le).\n", states[s]->label, t_f[s] );
		}
	}
} 


void HMMSystem::printBorderStatePopulations( int trans)
{
/*
        for( int c = 0; c < 30; c++ )
                printf(" ");
        for( int o = 0; o < globalN_OBSERVABLES; o++ )
        {
                if( trans == INTERPRETER_AA )
                {
                        printf("%c", all_res[o] );
                        for( int c = 1; c < 16; c++ )
                                printf(" ");
                }
                else if ( trans == INTERPRETER_SUPER )
                        printSTGroupName(o);
		else if( trans >= INTERPRETER_SCHEMES && trans <= INTERPRETER_SCHEMES+9 )
		{
			printSchemeGroupName(o, trans-INTERPRETER_SCHEMES);
		}
        }
        printf("\n");
	
	double av_pop = 0;
	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			av_pop += states[s]->obs_population[o];
	}	

	av_pop /= nstates;

	for( int s = 0; s < nstates; s++ )
	{
		if( strncasecmp( states[s]->label, "Border", 6) ) 
			continue;

		double tot_pop = 0;

		printf("%s ", states[s]->label );
		
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			tot_pop += states[s]->obs_population[o];
			printf("%lf ", states[s]->obs_population[o] );
		}
		
		printf("\n");
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		if( strncasecmp( states[s]->label, "Border", 6) ) 
			continue;
		double tot_pop = 0;
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			tot_pop += states[s]->obs_population[o];
		if( tot_pop / av_pop < 0.1 )
			printf("State %s is underpopulated %lf/%lf.\n", states[s]->label, tot_pop, av_pop );	
	}
	
	for( int s = 0; s < nstates; s++ )
	{
		if( strncasecmp( states[s]->label, "Border", 6) ) 
			continue;
		double ilv_pop = states[s]->obs_population[1] + states[s]->obs_population[2] + states[s]->obs_population[3];
		double gast_pop = states[s]->obs_population[18] + states[s]->obs_population[17] + states[s]->obs_population[16] + states[s]->obs_population[15];

		if( gast_pop - ilv_pop > 300 && states[s]->theClass == 'M' || states[s]->theClass == 'I' )
			printf("State %s of type %c has GAST %lf and ILV only %lf.\n",				
				states[s]->label, states[s]->theClass, gast_pop, ilv_pop );
		if( ilv_pop - gast_pop > 300 && states[s]->theClass == 'M' || states[s]->theClass == 'I' )
			printf("State %s of type %c has GAST %lf and ILV only %lf.\n",				
				states[s]->label, states[s]->theClass, gast_pop, ilv_pop );
	}
*/
}

void HMMSystem::TransitionGroupCheck( void )
{
	for( int TG = 0; TG < nTransitionGroups; TG++ )
	{
		int a_state = 0;

		int GS = transition_groups_len[TG];
		int nt;
		int trans_list[MAX_TRANSITIONS];
		int cur_len = -1;

		for( int GE = 0; GE < GS; GE++ )
		{
			getTransitionsFrom( transition_groups[TG][GE], trans_list, &nt, '.' );

			if( cur_len != -1 && nt != cur_len )	
			{
				printf("Number of transitions for state %s does not match group's.\n",
					states[transition_groups[TG][GE]]->label );
				exit(1);
			}
			
			for( int el1 = 0; el1 < nt; el1++)
			{
				for( int el2 = 0; el2 < nt; el2++ )
				{
					if( el1 == el2 ) continue;

					if( trans_list[el1] == trans_list[el2] )
					{
						printf("Duplicate members %s and %s in transition group for state %s.\n", states[trans_list[el1]]->label, states[trans_list[el2]]->label, states[transition_groups[TG][GE]]->label );
						exit(1);
					}
				}
			}
		}	
	}
}

void HMMSystem::setPseudocounts( HMMDirective *directive  )
{
	memset( pseudo_add, 0, sizeof(double) * globalN_OBSERVABLES );
	for( int d = 0; d < directive->n_data; d ++ )
	{
		int *obs = directive->translated_observable_string[d];	
		int len = directive->translated_len[d];

		for( int t = 0; t < len; t++ )
			pseudo_add[obs[t]] += 1;
	}
/*
	double tsum = 0;
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
		tsum += pseudo_add[o];
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
		pseudo_add[o] /= tsum;
*/
}

void HMMSystem::addPseudocounts( double frac )
{
	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
//			printf("s: %d o: %d was: %lf adding: %lf\n", s, o, p_obs_next_free[s+nstates*o], pseudo_add[o] * frac / nstates ); 
			p_obs_next_free[s+nstates*o] += pseudo_add[o] * frac / nstates;
		}
	}
}





