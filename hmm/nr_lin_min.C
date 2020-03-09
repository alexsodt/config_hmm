#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define USE_K	20.0

double USE_K = 0.0;

//#define DO_DEBUG_SURFACE
#include "HMM.h"
#include "HMMDirective.h"
#include <math.h>

#ifdef PARALLEL
#include "mpi.h"
#endif
#include "GSL_interface.h"
#include "Fuzzy.h"
#define DISABLE_GRADIENT_INIT
#define DO_NR

int do_print_sigmas = 0;
int do_detailed_balance = 1;

#define ADD_EPS (0.0)

double mmin( double a, double b )
{
	if( a > b ) return b;

	return a;
}

int my_isnan( double );

void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
	void (*dfunc)(double [], double []));
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));

HMMSystem *the_local_system = NULL;
HMMDirective *the_directive = NULL;
extern int globalN_OBSERVABLES;
	
extern "C" void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
        extern "C" void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []));

double global_f = 0.0;

void test_gradient( HMMSystem &theSystem, HMMDirective &theDirective, double *w_ref=NULL );
void AnalyzeGradient( double *dw, double *Grad  );

int do_crash_on_bad = 0;
int do_bwgrad = 0;
	
int ConvertSystemParametersFrom( double *w, double fp  ) // use "fp" fraction BaumWelch.
{
	int nstates = the_local_system->nstates;
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( (1-fp), f_obs, f_t_ab );

	
	int tp = 1;


	for( int EG = 0; EG < the_local_system->nEmissionGroups; EG++ )
	{
		double fp = 1 - f_obs[the_local_system->emission_groups[EG][0]];

		double p[globalN_OBSERVABLES];

		double n_sum = 1e-20;

		int ttp = tp;

		for( int o = 0; o < globalN_OBSERVABLES; o++, ttp++ )
		{
/*			if( w[ttp] > 50.0 ) */
			if( my_isnan( exp(mmin(w[ttp],50.0)) ) )
			{
				printf("%s %d\n", the_local_system->states[the_local_system->emission_groups[EG][0]]->label, o );
				printf("Not okay: w[%d]: %le.\n", ttp, w[ttp] );
				return 0;
			}
			n_sum += exp(mmin(w[ttp],50.0) );
		}
		
		ttp = tp;
			
		for( int o = 0; o < globalN_OBSERVABLES; o++, ttp++ )
		{
			double tval = w[ttp];
	
			p[o] = exp( mmin(tval,50.0) ) / n_sum;
		}
		
		// this is the "w" probability.

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		for( int GE = 0; GE < the_local_system->emission_groups_len[EG]; GE++ )
			the_local_system->p_obs[o*nstates+the_local_system->emission_groups[EG][GE]] = p[o];// + fp * the_local_system->p_obs_BaumWelch[o*nstates+the_local_system->emission_groups[EG][GE]];

		tp = ttp;
	}
	
#ifdef DO_TRANSITION_GROUPS
	for( int TG = 0; TG < the_local_system->nTransitionGroups; TG++ )
	{
		double fp = 1 - f_obs[the_local_system->transition_groups[TG][0]];

		int trans_list[MAX_TRANSITIONS];
		int nt;

		the_local_system->getTransitionsFrom( the_local_system->transition_groups[TG][0], trans_list, &nt, '.' );
	
		double p[nt];

		double n_sum = 1e-20;

		int ttp = tp;

		for( int o = 0; o < nt; o++, ttp++ )
			n_sum += exp(mmin(w[ttp],50.0) );
		
		ttp = tp;
			
		for( int o = 0; o < nt; o++, ttp++ )
		{
			double tval = w[ttp];
	
			p[o] = exp( mmin(tval,50.0) ) / n_sum;
		}
		
		// this is the "w" probability.

		for( int GE = 0; GE < the_local_system->transition_groups_len[TG]; GE++ )
		{
			int trans_list[MAX_TRANSITIONS];
			int nt;

			the_local_system->getTransitionsFrom( the_local_system->transition_groups[TG][GE], trans_list, &nt, '.' );
			int s_from = the_local_system->transition_groups[TG][GE];

			for( int o = 0; o < nt; o++ )
			{
				the_local_system->t_ab[s_from*nstates+trans_list[o]] = p[o]; //(1-fp)*p[o] + fp * 
//				the_local_system->t_ab_BaumWelch[s_from*nstates+trans_list[o]];
			}
		}

		tp = ttp;
	}

#else
	if( do_detailed_balance )
	{
		// we are given the equilibration populations, w, and symmetric rate generator, Q_s,	
		// t_ab = W^{-1} . Q_s
		// where W is diag(w).
		double *wi = (double *)malloc( sizeof(double) * nstates );
		double *Qs = (double *)malloc( sizeof(double) * nstates*nstates);
		memset( Qs, 0, sizeof(double) * nstates * nstates );
		memset( wi, 0, sizeof(double) * nstates );
	}
	else
	{
		for( int o = 0; o < nstates; o++ )
		{
			double fp_t = 1 - f_t_ab[o];
			double n_sum = 1e-20;
			int ttp = tp;
	
			int trans_list[MAX_TRANSITIONS];
			int ntrans;
	
			the_local_system->getTransitionsFrom( o, trans_list, &ntrans, '.' );
		
			for( int sx = 0; sx < ntrans; sx++, ttp++ )
				n_sum += exp( mmin(w[ttp],50.0) );
			
			if( n_sum > 1e100 )
				n_sum = 1e100;
	
			ttp = tp;
	
			for( int sx = 0; sx < ntrans; sx++, ttp++ )
			{
				double tval = w[ttp];
	
				//if( tval > 50.0 ) 
				if( my_isnan( exp(mmin(w[ttp],50.0)) ) )
				{
					printf("%s %s\n", the_local_system->states[o]->label, the_local_system->states[trans_list[sx]]->label );
					printf("Not okay w[%d]: %le.\n", ttp, w[ttp] );
					return 0;
				}
	
				int s2 = trans_list[sx];
	
				the_local_system->t_ab[o*nstates+s2] =  /*fp_t * the_local_system->t_ab_BaumWelch[o*nstates+s2] + (1.0-fp_t) **/ exp( mmin(tval,50.0) ) / n_sum;
	
				if( my_isnan( the_local_system->t_ab[o*nstates+s2]) )
				{
					printf("w: %le n_sum: %le\n", w[ttp], n_sum );
					printf("NAN error in nr_lin_min.\n");	
					return 0;	
		//			exit(1);
				}
			}
			tp = ttp;
		}
	}
#endif
	
	double init_sum = 1e-50;
	int ttp = tp;
	for( int o = 0; o < nstates; o++, ttp++ )
		init_sum += exp( w[ttp] );
	ttp = tp;
	for( int o = 0; o < nstates; o++, ttp++ )
		the_local_system->init_p[o] = exp( w[ttp] ) / init_sum;
	tp = ttp;

	return 1; // okay
}

int GetNParams( void )
{
	int nstates = the_local_system->nstates;
	int np = globalN_OBSERVABLES * the_local_system->nEmissionGroups + nstates;

#ifdef DO_TRANSITION_GROUPS
	for( int o = 0; o < the_local_system->nTransitionGroups; o++ )
	{
		int s_from = the_local_system->transition_groups[o][0];

		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		the_local_system->getTransitionsFrom( s_from, trans_list, &ntrans, '.' );

		np += ntrans;
	}
#else
	for( int o = 0; o < nstates; o++ )
	{
		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		the_local_system->getTransitionsFrom( o, trans_list, &ntrans, '.' );

		np += ntrans;
	}
#endif
	return np;
}

int ParamIsEmission( int p )
{
	int np = globalN_OBSERVABLES * the_local_system->nEmissionGroups;

	if( p < np )
		return 1;
	else
		return 0;
}

void printParameterInfo( int p, char *buffer )
{
	int ntrans = 0;
	for( int EG = 0; EG < the_local_system->nTransitionGroups; EG++ )
	{

		int trans_list[MAX_TRANSITIONS];
		int nt;

		the_local_system->getTransitionsFrom( the_local_system->transition_groups[EG][0], trans_list, &nt, '.' );
		ntrans += nt;
	}
	int nstates = the_local_system->nstates;
	if( p < globalN_OBSERVABLES * the_local_system->nEmissionGroups)
	{
		int EG = p / globalN_OBSERVABLES;
		int o = p % globalN_OBSERVABLES;

		printf("Emission of o: %d from states:\n", o );
		double n_c = 0, n_f = 0;
		for( int x = 0; x < the_local_system->emission_groups_len[EG]; x++ )
		{
			n_c = the_local_system->p_obs_next_clamped[o*nstates+the_local_system->emission_groups[EG][x]];
			n_f = the_local_system->p_obs_next_free[o*nstates+the_local_system->emission_groups[EG][x]];
			printf("\t%s (%le, %le) p: %le\n", the_local_system->states[the_local_system->emission_groups[EG][x]]->label, n_c, n_f, the_local_system->p_obs[o*nstates+the_local_system->emission_groups[EG][x]] );
		}
	}
	else if( p < globalN_OBSERVABLES * the_local_system->nEmissionGroups + ntrans )
	{
	
		int tp = globalN_OBSERVABLES * the_local_system->nEmissionGroups;
	
		for( int EG = 0; EG < the_local_system->nTransitionGroups; EG++ )
		{
			double n_sum = 1e-20;
	
			int trans_list[MAX_TRANSITIONS];
			int nt;
	
			the_local_system->getTransitionsFrom( the_local_system->transition_groups[EG][0], trans_list, &nt, '.' );
	
			for( int o = 0; o < nt; o++,tp++ )
			{
				if( tp == p )
				{
					for( int GE = 0; GE < the_local_system->transition_groups_len[EG]; GE++ )
					{
		
						int trans_list2[MAX_TRANSITIONS];
						int nt2;
	
						the_local_system->getTransitionsFrom( the_local_system->transition_groups[EG][GE], trans_list2, &nt2, '.' );
	
						printf("transition group %d Transition from %s to %s.\n", o,
								the_local_system->states[the_local_system->transition_groups[EG][GE]]->label,
								the_local_system->states[trans_list2[o]]->label );
					}	
	
					return;
				}
			}
		}
	}
	else
	{
		int s = p - (globalN_OBSERVABLES * the_local_system->nEmissionGroups + ntrans);

		printf("Init prob for state '%s'.\n", the_local_system->states[s]->label ); 
	}

}

int ConvertSystemParametersTo( double *w, double fp )
{
	int nstates = the_local_system->nstates;
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( (1-fp), f_obs, f_t_ab );
	int tp = 1;

	for( int EG = 0; EG < the_local_system->nEmissionGroups; EG++ )
	{
		double fp = 1.0 - f_obs[the_local_system->emission_groups[EG][0]];

		double p[globalN_OBSERVABLES];

		double n_sum = 1e-20;

		for( int o = 0; o < globalN_OBSERVABLES; o++,tp++ )
		{
			double p_temp =  (the_local_system->p_obs[o*nstates+the_local_system->emission_groups[EG][0]] - fp*the_local_system->p_obs_BaumWelch[o*nstates+the_local_system->emission_groups[EG][0]]) / (1-fp);
			w[tp] = log( p_temp );

			if( !(w[tp] < 1) && !(w[tp]>0) )
			{
				printf("p_obs: %le p_obs_BW: %le\n", the_local_system->p_obs[o*nstates+the_local_system->emission_groups[EG][0]], the_local_system->p_obs_BaumWelch[o*nstates+the_local_system->emission_groups[EG][0]] );
//				printf("Parameter EG: %d obs: %d in error, p_temp: %le fp: %le\n", p_temp, fp );
				if( do_crash_on_bad )
					exit(1);
				return 0;
			}
		}
	}
#ifdef DO_TRANSITION_GROUPS	
	for( int EG = 0; EG < the_local_system->nTransitionGroups; EG++ )
	{
		double fp = 1.0 - f_obs[the_local_system->transition_groups[EG][0]];
		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int nt;

		the_local_system->getTransitionsFrom( the_local_system->transition_groups[EG][0], trans_list, &nt, '.' );

		double p[nt];

		for( int o = 0; o < nt; o++,tp++ )
		{
			int s_from = the_local_system->transition_groups[EG][0];

			double p_temp =  (the_local_system->t_ab[s_from*nstates+trans_list[o]] - 
			               fp*the_local_system->t_ab_BaumWelch[s_from*nstates+trans_list[o]]) / (1-fp);
			w[tp] = log( p_temp );

			if( !(w[tp] < 1) && !(w[tp]>0) )
			{
				printf("OOF.\n");
				if( do_crash_on_bad )
					exit(1);
				return 0;
			}
		}
	}
#else
	for( int o = 0; o < nstates; o++ )
	{
		double fp_t = 1.0 - f_t_ab[o];

		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		the_local_system->getTransitionsFrom( o, trans_list, &ntrans, '.' );
	
		for( int sx = 0; sx < ntrans; sx++, tp++ )
		{
			int s2 = trans_list[sx];
			
			double p_temp =  (the_local_system->t_ab[o*nstates+s2] - fp_t*the_local_system->t_ab_BaumWelch[o*nstates+s2]) / (1-fp_t);
			w[tp] = log( p_temp );
			if( !(w[tp] < 1) && !(w[tp]>0) )
			{
				if( do_crash_on_bad )
					exit(1);
				return 0;
			}
		}
	}
#endif
	for( int o = 0; o < nstates; o++, tp++ )
		w[tp] = log( the_local_system->init_p[o] + 1e-50 );
	return 1;
}

void ConvertGradientTo( double *dw, double *Grad  )
{
	int tp = 1;
	int nstates = the_local_system->nstates;

	for( int EG = 0; EG < the_local_system->nEmissionGroups; EG++ )
	{
		double p[globalN_OBSERVABLES];

		double n_sum = 1e-20;

		int astate = the_local_system->emission_groups[EG][0];

		for( int o = 0; o < globalN_OBSERVABLES; o++, tp++ )
			dw[tp] = Grad[1+astate*globalN_OBSERVABLES+o];
	}
#ifdef DO_TRANSITION_GROUPS	
	for( int EG = 0; EG < the_local_system->nTransitionGroups; EG++ )
	{
		int trans_list[MAX_TRANSITIONS];
		int nt;

		the_local_system->getTransitionsFrom( the_local_system->transition_groups[EG][0], trans_list, &nt, '.' );
		double p[nt];

		double n_sum = 1e-20;

		int astate = the_local_system->transition_groups[EG][0];

		for( int o = 0; o < nt; o++, tp++ )
			dw[tp] = Grad[1+nstates*globalN_OBSERVABLES+astate*nstates+trans_list[o]];
	}
#else
	for( int o = 0; o < nstates; o++ )
	{
		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		the_local_system->getTransitionsFrom( o, trans_list, &ntrans, '.' );
	
		for( int sx = 0; sx < ntrans; sx++, tp++ )
		{
			int s2 = trans_list[sx];

			dw[tp] = Grad[1+nstates*globalN_OBSERVABLES+o*nstates+s2];
		}
	}
#endif
	for( int o = 0; o < nstates; o++, tp++ )
	{
#ifndef DISABLE_GRADIENT_INIT
		dw[tp] = Grad[1+nstates*globalN_OBSERVABLES+nstates*nstates+o];
#else
		dw[tp] = 0;
#endif
	}
}

void AnalyzeGradient( double *dw, double *Grad  )
{
	int tp = 1;
	int nstates = the_local_system->nstates;

	int neg = the_local_system->nEmissionGroups;
	int gsorter[neg];
	double gde[neg];

	for( int s = 0; s < the_local_system->nEmissionGroups; s++ )
		gsorter[s] = s;

	for( int EG = 0; EG < the_local_system->nEmissionGroups; EG++ )
	{
		gde[EG] = 0;

		int astate = the_local_system->emission_groups[EG][0];
		for( int o = 0; o < globalN_OBSERVABLES; o++, tp++ )
			gde[EG] += Grad[1+astate*globalN_OBSERVABLES+o] * Grad[1+astate*globalN_OBSERVABLES+o];
	}

	int done = 0;

	while( !done )
	{
		done = 1;

		for( int e = 0; e < neg-1; e++ )
		{
			if( gde[gsorter[e]] < gde[gsorter[e+1]] )
			{
				int t= gsorter[e];
				gsorter[e] = gsorter[e+1];
				gsorter[e+1] = t;
				done = 0;
			}
		}
	}

	for( int s = 0; s < neg; s++ )
	{
		printf("Group with member '%s' has dE: %le:\n",
			the_local_system->states[the_local_system->emission_groups[gsorter[s]][0]]->label,
			gde[gsorter[s]] );
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			printf("obs %d value %le has gradient %le.\n",
				o, the_local_system->p_obs[o*nstates+the_local_system->emission_groups[gsorter[s]][0]], Grad[1+the_local_system->emission_groups[gsorter[s]][0]*globalN_OBSERVABLES+o] );
		}
	}
	

}

double f( double *w )
{

	int nstates = the_local_system->getNStates();
	int n_params = GetNParams(); //nstates*nstates+nstates*globalN_OBSERVABLES;
	

	int okay = ConvertSystemParametersFrom( w, 1 - global_f );	
	
	if( !okay) 
	{
		printf("returning straight out with not okay!!\n");

		return 1e30;
	}
//	else
//		printf("Parameters converted fine at start of f.\n");
	double fval = 0;

	int taskid=0, nprocs=1;

#ifdef PARALLEL
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
	double log_sum_free = 0;
	double log_sum_clamped = 0;
	
	double surface_f = 0;
	double surface_k = USE_K;

	int prot_sys_off[the_directive->n_data];
	int prot_sys_len[the_directive->n_data];

	int prot_n_sites[the_directive->n_data*2];
	double prot_avz    [the_directive->n_data*2];

	double par_f = 0;
	
	for( int d = taskid; d < the_directive->n_data; d += nprocs )
	{
		int len = the_directive->translated_len[d];
		char *wild_string = (char *)malloc( sizeof(char) * (len+1) );
		for( int s = 0; s < len; s++ )
			wild_string[s] = '.';
		wild_string[len] = '\0';
					
		double *alpha_constrained = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha_norm = (double *)malloc( sizeof(double) * len * nstates );

		calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, the_directive->translated_observable_string[d], wild_string, len, *the_local_system );
		log_sum_clamped += alpha_constrained_norm[len-1];

		fval += ( alpha_constrained_norm[len-1]);

		free(wild_string);

		free( alpha_constrained );
		free( alpha);

		free( alpha_constrained_norm );
		free( alpha_norm );
	}
	
	double f_final = 0;
#ifdef PARALLEL
	MPI_Allreduce( &fval, &f_final, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	f_final = fval;
#endif
	printf("prob_f: %le\n", f_final );

	fval = f_final;

	if( my_isnan(fval) ) 
	{
		printf("Parameter error at the end of f.\n");
		return 1e30;	
	}
	return -fval;
}


int in_gradient_debug = 0;


double fdf( double *w, double *dw )
{
	int taskid=0, nprocs=1;
#ifdef PARALLEL
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif	
	int nstates = the_local_system->getNStates();
	int n_params = GetNParams(); //nstates*nstates+nstates*globalN_OBSERVABLES;

	int err = ConvertSystemParametersFrom( w, 1 - global_f );	

	if( !err )
	{
		printf("Gradient parameter error.\n");
	}

	double fval = 0;
	
	the_local_system->zeroNextPMat();
	
	/*
		for d/dp expec(z) we just need to scale the fixed derivatives by z / nsites
		for (av z)^2, we need 2 * av z * z / nsites.. for each protein system and side we need nsites and av z.
	
		for each site, we will constrain the topology to pass a surface through a specific point.
	*/
					
	
	for( int d = taskid; d < the_directive->n_data; d += nprocs )
	{
		int len = the_directive->translated_len[d];
		char *wild_string = (char *)malloc( sizeof(char) * (len+1) );
		for( int s = 0; s < len; s++ )
			wild_string[s] = '.';
		wild_string[len] = '\0';
					
		double *alpha_constrained = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
		double *alpha_norm = (double *)malloc( sizeof(double) * len * nstates );
		
		double *beta_constrained = (double *)malloc( sizeof(double) * len * nstates );
		double *beta = (double *)malloc( sizeof(double) * len * nstates );
		double *beta_constrained_norm = (double *)malloc( sizeof(double) * len * nstates );
		double *beta_norm = (double *)malloc( sizeof(double) * len * nstates );

		 
		calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, the_directive->translated_observable_string[d], the_directive->raw_class_string[d],  len, *the_local_system );
		calculateBetaForStringConstrained( beta_constrained, beta_constrained_norm, the_directive->translated_observable_string[d], the_directive->raw_class_string[d],  len, *the_local_system );
						
		the_local_system->incrementNextPMat_Clamped( alpha_constrained, beta_constrained, the_directive->translated_observable_string[d],  len,  the_directive->weights[d] );				

		fval += alpha_constrained_norm[len-1] * the_directive->weights[d];

		if( my_isnan(alpha_constrained_norm[len-1]) )
		{
			printf("NAN for %s %s\n", the_directive->raw_observable_string[d], the_directive->raw_class_string[d] );
#ifdef PARALLEL
			MPI_Finalize();
#endif
			exit(1);
		}


		free( alpha_constrained );
		free( alpha);

		free( alpha_constrained_norm );
		free( alpha_norm );
		
		free( beta_constrained );
		free( beta);

		free( beta_constrained_norm );
		free( beta_norm );
			
		free(wild_string);
	}

	the_local_system->CommunicateShortPMats();
	
	double pseudo_free, pseudo_clamped;
	
	double f_final = 0;

#ifdef PARALLEL
	MPI_Allreduce( &fval, &f_final, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#endif
	double *Grad = (double *)malloc( sizeof(double) * ( 1 + nstates*globalN_OBSERVABLES + nstates*nstates + nstates) );
	memset( Grad, 0, sizeof(double)* ( 1 + nstates*globalN_OBSERVABLES + nstates*nstates) );

	the_local_system->setFractionalGradientSplit( Grad, global_f );

	ConvertGradientTo( dw, Grad );
	
	double G_sum = 0;

	for( int x = 0; x < n_params; x++ )
		G_sum += fabs(dw[1+x]);// * dw[1+x];
	
	printf("G_sum: %le %le\n", G_sum, f_final);
	fflush(stdout);
	
	for( int x = 0; x < n_params; x++ )
		dw[1+x] *= -1;

	free(Grad);

	return -f_final ;
}

void df( double *w, double *dw )
{
	fdf( w, dw );
}

double SA_opt( HMMSystem &theSystem, HMMDirective &theDirective )
{
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );

	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			double tc = log( the_local_system->p_obs[o*nstates+s] );
			if( my_isnan(tc) || tc > 50 )
				printf("Trouble, parameter is %le.\n", tc );
			the_local_system->p_obs[o*nstates+s] = the_local_system->p_obs[o*nstates+s] * f_obs[s] + (1-f_obs[s]) *the_local_system->p_obs_BaumWelch[o*nstates+s];
		}
	}
	for( int s = 0; s < nstates; s++ )
	{
		int trans_list[MAX_TRANSITIONS];	
		int ntrans;

		the_local_system->getTransitionsFrom( s, trans_list, &ntrans, '.' );

		for( int sx = 0; sx < ntrans; sx++ )
		{
			double tc = log( the_local_system->t_ab[s*nstates+trans_list[sx]] );
			if( my_isnan(tc) || tc > 50 )
				printf("Trouble (2), parameter is %le.\n", tc );

		}	
	}

	printf("Done with trouble check.\n");
	
	if( !ConvertSystemParametersTo( p, 1 - global_f ) )
	{
		printf("Bad input to grad_opt.\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
		return 1e30;
	}
	else
		printf("The parameters were fine.\n");
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 1");

	int popt[n_params];
	int n_2_opt = 0;

	for( int x = 0; x < n_params; x++ )
	{
		if( fabs(xi[1+x]) > 1e-7 )
			popt[n_2_opt++] = x; 
	}
		
	double *p0 = (double *)malloc( sizeof(double) * (1+n_params) );
	memcpy( p0, p, sizeof(double) * (1+n_params) );

	double F0 = f(p);
	double del_F = 0;
	int ndf = 0;
	// calibrate temperature..

	for( int s = 0; s < 10; s++ )
	{
		int p_c = 1 + rand() % n_params;
		double p_v = 0.95 + 0.10 * (rand () % 100) / 100.0;

		p[p_c] *= p_v;

		double fp = f(p);

		del_F += ((fp-F0)*(fp-F0));
		ndf += 1;		
	}	

	del_F /= ndf;
	del_F = sqrt(del_F);

	double T = del_F;
	double cur_F = 0;

	memcpy( p, p0, sizeof(double) * (1+n_params) );

	for( int ot = 0; ot < 10; ot++, T /= 2 )
	{
		for( int s = 0; s < 100; s++ )
		{
			memcpy( p, p0, sizeof(double) * (1+n_params) );

			int p_c_t = rand() % n_2_opt;
			int p_c = 1+ popt[p_c_t];
			double p_v = 0.95 + 0.10 * (rand () % 100) / 100.0;
						
			p[p_c] *= p_v;

			double f1 = f(p);

			double pr = exp( -(f1-cur_F) / T );

			double rn = (double)rand() / (double)RAND_MAX;

			if( rn < pr )
			{
				printf("Accepted F=%le, delF: %le T: %le p: %le\n", f1, f1-cur_F, T, (pr>1.0?1.0:pr) );
				memcpy( p0, p, sizeof(double) * (1+n_params) );
				cur_F = f1;	
			}
		}
	}

	return f(p0);
}

double grad_opt( HMMSystem &theSystem, HMMDirective &theDirective )
{	
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );
	
	if( !ConvertSystemParametersTo( p, 1 - global_f ) )
	{
		printf("Bad input to grad_opt.\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
		return 1e30;
	}
	else
		printf("The parameters were fine.\n");
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 1");

	

#ifdef DO_NR

	double fstart = f(p);
	printf("fstart: %le\n", fstart );	
	df( p, xi );	

	double fret;

	//dlinmin( p, xi, n_params, &fret, f, df );
//	linmin( p, xi, n_params, &fret, f );


	int iter;
	dfpmin( p, n_params, 1e-14, &iter, &fret, f, df ); 
	//frprmn( p, n_params, 1e-14, &iter, &fret, f, df ); 

#else

	double fret = 0;	
	do_GSL_minimization( USE_GSL_BFGS);	
#endif


	int okay = ConvertSystemParametersFrom(p , 1 - global_f  );	

	if( !okay )
	{
		printf("Grad opt, garbage out!\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
	}
//	printf("Crash barrier 2.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 2");

	free(p);
	free(xi);

	return fret;
		
}


double bwgrad_opt( HMMSystem &theSystem, HMMDirective &theDirective )
{	
	do_bwgrad = 1;

	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );

	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			the_local_system->p_obs[o*nstates+s] = the_local_system->p_obs[o*nstates+s] * f_obs[s] + (1-f_obs[s]) *the_local_system->p_obs_BaumWelch[o*nstates+s];
	}

	
	if( !ConvertSystemParametersTo( p, 1 - global_f ) )
	{
		printf("Bad input to grad_opt.\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
		return 1e30;
	}
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 1");

	double fstart = f(p);
//	printf("fstart: %le\n", fstart );	
	df( p, xi );	

	double fret;

	//dlinmin( p, xi, n_params, &fret, f, df );
//	linmin( p, xi, n_params, &fret, f );


	int iter;
	dfpmin( p, n_params, 1e-14, &iter, &fret, f, df ); 
	//frprmn( p, n_params, 1e-14, &iter, &fret, f, df ); 

	int okay = ConvertSystemParametersFrom(p , 1 - global_f  );	

	if( !okay )
	{
		printf("Grad opt, garbage out!\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
	}
//	printf("Crash barrier 2.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 2");

	free(p);
	free(xi);


	return fret;
		
}

void test_gradient( HMMSystem &theSystem, HMMDirective &theDirective, double *w_ref )
{	
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	the_local_system->checkNormalization();
	the_local_system->checkEGConsistency();
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );
/*
	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			the_local_system->p_obs[o*nstates+s] = the_local_system->p_obs[o*nstates+s] * f_obs[s] + (1-f_obs[s]) *the_local_system->p_obs_BaumWelch[o*nstates+s];
	}

*/
	
	ConvertSystemParametersTo( p, 1 - global_f );

	double fstart = f(p);
	printf("fstart: %le\n", fstart );	
	df( p, xi );	
	printf("done with gradient.\n");

	
	int iters = 0;
	for( double dG = 1e-15; dG < 10.0 && iters < 1000; dG *= 1.25, iters++ )
	{
		double *new_p = (double *)malloc( sizeof(double) * ( 1 + n_params) );

		double new_E = fstart;

		for( int i = 1; i <= n_params; i++ )
			new_p[i] = p[i];

		for( int i = 1; i <= n_params; i++ )
		{
			new_E -= dG * xi[i] * xi[i];
			new_p[i] = p[i] - dG * xi[i];
		}
	
		double f_new = f(new_p);

		printf("%le %le %le ratio: %le\n", dG, f_new - fstart, new_E-fstart, (f_new - fstart)/(new_E-fstart) ); 

		free(new_p);

	}
}


void test_gradient_element_by_element( HMMSystem &theSystem, HMMDirective &theDirective, double *w_ref )
{	
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	the_local_system->checkNormalization();
	the_local_system->checkEGConsistency();
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );
/*
	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			the_local_system->p_obs[o*nstates+s] = the_local_system->p_obs[o*nstates+s] * f_obs[s] + (1-f_obs[s]) *the_local_system->p_obs_BaumWelch[o*nstates+s];
	}

*/
	
	ConvertSystemParametersTo( p, 1 - global_f );

	double fstart = f(p);
	printf("fstart: %le\n", fstart );	
	df( p, xi );	
	printf("done with gradient.\n");

	
	for( int elem = 1/* + n_params - nstates*/; elem <= n_params; elem++ )
	{
//		if( elem != 235 && elem != 252 && elem != 254 && elem != 261 ) continue;

		int iters = 0;

		if( fabs(xi[elem]) < 1e-7 )
		{
			printParameterInfo( elem-1, NULL);
			printf("Zero gradient.\n");
			continue;
		}
		double best_rat = 1e10;

		for( double dG = 1e-5; dG < 1000.0 && iters < 10; dG *= 3, iters++ )
		{
			double *new_p = (double *)malloc( sizeof(double) * ( 1 + n_params) );
	
			double new_E = fstart;
	
			for( int i = 1; i <= n_params; i++ )
				new_p[i] = p[i];
	
			for( int i = 1; i <= n_params; i++ )
			{
				if( i == elem )
				{
					new_E -= dG * xi[i] * xi[i];
					new_p[i] = p[i] - dG * xi[i];
				}
				else
					new_p[i] = p[i];
			}
		
			double f_new = f(new_p);
	
			double cur_rat = (f_new - fstart)/(new_E-fstart);

			printf("dG: %le f: %le cur_rat: %le f_new-f_start: %le, new_E-fstart:%le\n", dG, f_new-fstart, cur_rat, f_new-fstart, new_E-fstart );
			if( fabs(1.0-cur_rat) < fabs(1.0-best_rat) )
				best_rat = cur_rat;

	
			free(new_p);

		}
		
		
		printParameterInfo( elem-1, NULL);
	
		printf("elem: %d best ratio: %le dF: %le FD dF: %le\n", elem, best_rat, xi[elem], best_rat * xi[elem] ); 
	}
}



double print_sigmas( HMMSystem &theSystem, HMMDirective &theDirective )
{	
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( global_f, f_obs, f_t_ab );
	printf("Done with trouble check.\n");
	
	if( !ConvertSystemParametersTo( p, 1 - global_f ) )
	{
		printf("Bad input to grad_opt.\n");
#ifdef PARALLEL
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		exit(1);
		return 1e30;
	}
	else
		printf("The parameters were fine.\n");
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 1");

	
	do_print_sigmas = 1;

	return f(p);
		
}
