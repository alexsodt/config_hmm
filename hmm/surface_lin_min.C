#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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

static int do_print_sigmas = 0;

#define ADD_EPS (0.0)

double mmin( double a, double b );
int my_isnan( double );

void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double []),
	void (*dfunc)(double [], double []));
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []));

extern HMMSystem *the_local_system;
extern HMMDirective *the_directive;
extern int globalN_OBSERVABLES;
	
extern "C" void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));
        extern "C" void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []));

extern double global_f;

void surface_test_gradient( HMMSystem &theSystem, HMMDirective &theDirective, double *w_ref=NULL );
void AnalyzeGradient( double *dw, double *Grad  );
void SurfaceReport( HMMSystem &theSystem, HMMDirective &theDirective );

extern int do_crash_on_bad;
extern int do_bwgrad;
	
int surface_ConvertSystemParametersFrom( double *w ) // use "fp" fraction BaumWelch.
{
	int nstates = the_local_system->nstates;
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( 1.0, f_obs, f_t_ab );

	// all observables fixed.

	
	int tp = 1;


	for( int s = 0; s < nstates; s++ )
	{
		// fix these emission probabilities: border states 50/50, main states fixed at 100%
		
		for( int o = 0; o < globalN_OBSERVABLES; o++ )	
			the_local_system->p_obs[o*nstates+s] = 0;

		if( !strncasecmp( the_local_system->states[s]->label, "Border", 6) )
		{
			// a border state. format: Border_A_B
			int o1 = the_local_system->states[s]->label[7] - 'A';
			int o2 = the_local_system->states[s]->label[9] - 'A';

			if( the_local_system->states[s]->label[7] >= 'a' )
				o1 = globalN_OBSERVABLES/2 + the_local_system->states[s]->label[7] - 'a';
			if( the_local_system->states[s]->label[9] >= 'a' )
				o2 = globalN_OBSERVABLES/2 + the_local_system->states[s]->label[9] - 'a';

			the_local_system->p_obs[o1*nstates+s] = 0.5;
			the_local_system->p_obs[o2*nstates+s] = 0.5;
		} 
		else
		{
			int o1 = the_local_system->states[s]->label[0] - 'A';
			
			if( the_local_system->states[s]->label[0] >= 'a' )
				o1 = globalN_OBSERVABLES/2 + the_local_system->states[s]->label[0] - 'a';
			
			the_local_system->p_obs[o1*nstates+s] = 1.0;
		}
	}
	
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

void check_surface_consistency( HMMSystem &theSystem, HMMDirective &theDirective )
{
	// make sure our surfaces don't violate the network.

	// do we always make an allowed transition?

	int nstates = theSystem.nstates;

	int reported[nstates*nstates];
	int allowed_transitions[nstates*nstates];
	memset( allowed_transitions, 0, sizeof(int) * nstates * nstates );
	memset( reported, 0, sizeof(int) * nstates*nstates);

	for( int s = 0; s < nstates; s++ )
	{
		if( !strncasecmp( theSystem.states[s]->label, "Border", 6) )
		{
			// a border state. format: Border_A_B
			int o1 = theSystem.states[s]->label[7] - 'A';
			int o2 = theSystem.states[s]->label[9] - 'A';

			if( theSystem.states[s]->label[7] >= 'a' )
				o1 = globalN_OBSERVABLES/2 + theSystem.states[s]->label[7] - 'a';
			if( theSystem.states[s]->label[9] >= 'a' )
				o2 = globalN_OBSERVABLES/2 + theSystem.states[s]->label[9] - 'a';

			allowed_transitions[o1*nstates+o2] = 1;
			allowed_transitions[o2*nstates+o1] = 1;
		}
	}

	HMMDirective *the_directive = &theDirective;

	int error = 0;
	for( int d = 0; d < the_directive->n_data; d ++ )
	{
		int len = the_directive->translated_len[d];
		
		int pstate = the_directive->translated_observable_string[d][0];
		
		for( int t = 1; t < len; t++ )
		{
			int nstate = the_directive->translated_observable_string[d][t];

			if( pstate == nstate ) continue;

			if( !allowed_transitions[pstate*nstates+nstate] && !reported[pstate*nstates+nstate] )
			{
				char code1 = 'A'+pstate;
				char code2 = 'A'+nstate;
				if( pstate >= globalN_OBSERVABLES/2 )
					code1 = 'a' + (pstate - globalN_OBSERVABLES/2);
				if( nstate >= globalN_OBSERVABLES/2 )
					code2 = 'a' + (nstate - globalN_OBSERVABLES/2);
 
				printf("ERROR transition from %c to %c not allowed at position %d in string '%s'.\n",
					code1, code2, t, the_directive->raw_observable_string[d] );
				reported[pstate*nstates+nstate] = 1;
				error = 1;		
			}

			pstate = nstate;
		}					
	}

	if( error )
	{
		printf("System is not consistent.\n");
		exit(1);
	}
	
	printf("System OK.\n");
}

int surface_GetNParams( void )
{
	int nstates = the_local_system->nstates;
	int np = nstates; // initial probabilities
	int tp = 1;
	
	for( int o = 0; o < nstates; o++ )
	{
		double n_sum = 1e-20;

		int trans_list[MAX_TRANSITIONS];
		int ntrans;

		the_local_system->getTransitionsFrom( o, trans_list, &ntrans, '.' );

		np += ntrans;
	
	}
	
	return np; // okay
}

int surface_ParamIsEmission( int p )
{
	return 0;
}

void surface_printParameterInfo( int p, char *buffer )
{
	

}

int surface_ConvertSystemParametersTo( double *w  )
{
	int nstates = the_local_system->nstates;
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( 1.0, f_obs, f_t_ab );
	int tp = 1;

	/* Currenty there are no emission probabilities used */

	/* a single transition probability, to leave the state */

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
			
			double p_temp =  the_local_system->t_ab[o*nstates+s2];
			w[tp] = log( p_temp );
			if( !(w[tp] < 1) && !(w[tp]>0) )
			{
				if( do_crash_on_bad )
					exit(1);
				return 0;
			}
		}
	}
	for( int o = 0; o < nstates; o++, tp++ )
		w[tp] = log( the_local_system->init_p[o] + 1e-50 );
	return 1;
}

void surface_ConvertGradientTo( double *dw, double *Grad  )
{
	int tp = 1;
	int nstates = the_local_system->nstates;

	// currently nothing done with emission.

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
	for( int o = 0; o < nstates; o++, tp++ )
	{
#ifndef DISABLE_GRADIENT_INIT
		dw[tp] = Grad[1+nstates*globalN_OBSERVABLES+nstates*nstates+o];
#else
		dw[tp] = 0;
#endif
	}
}

#if 0
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
#endif

double surface_f( double *w )
{

	int nstates = the_local_system->getNStates();
	int n_params = surface_GetNParams(); //nstates*nstates+nstates*globalN_OBSERVABLES;
	

	int okay = surface_ConvertSystemParametersFrom( w );	
	
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
	double surface_k = 0;

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


extern int in_gradient_debug;


double surface_fdf( double *w, double *dw )
{
	int taskid=0, nprocs=1;
#ifdef PARALLEL
        MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif	
	int nstates = the_local_system->getNStates();
	int n_params = surface_GetNParams(); //nstates*nstates+nstates*globalN_OBSERVABLES;

	int err = surface_ConvertSystemParametersFrom( w );	

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
#else
	f_final = fval;
#endif
	double *Grad = (double *)malloc( sizeof(double) * ( 1 + nstates*globalN_OBSERVABLES + nstates*nstates + nstates) );
	memset( Grad, 0, sizeof(double)* ( 1 + nstates*globalN_OBSERVABLES + nstates*nstates) );

	the_local_system->setFractionalGradientSplit( Grad, 1.0);

	surface_ConvertGradientTo( dw, Grad );
	
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

void surface_df( double *w, double *dw )
{
	surface_fdf( w, dw );
}

#if 0
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
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
		return 1e30;
	}
	else
		printf("The parameters were fine.\n");
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
	MPI_Barrier(MPI_COMM_WORLD);
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

	f(p0);
}
#endif
double surface_grad_opt( HMMSystem &theSystem, HMMDirective &theDirective )
{	
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = surface_GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( 1.0, f_obs, f_t_ab );
/*
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
*/
	printf("Done with trouble check.\n");
	
	if( !surface_ConvertSystemParametersTo( p ) )
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

	//SurfaceReport(theSystem,theDirective);
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif
//	printf("After crash barrier 1");

//	surface_test_gradient( theSystem, theDirective, NULL ); 
	

#ifdef DO_NR

	double fstart =surface_f(p);
	printf("fstart: %le\n", fstart );	
	surface_df( p, xi );	

	double fret;

	//dlinmin( p, xi, n_params, &fret, f, df );
//	linmin( p, xi, n_params, &fret, f );


	int iter;
	dfpmin( p, n_params, 1e-14, &iter, &fret, surface_f, surface_df ); 
	//frprmn( p, n_params, 1e-14, &iter, &fret, f, df ); 

#else

	double fret = 0;	
	do_GSL_minimization( USE_GSL_BFGS);	
#endif


	int okay = surface_ConvertSystemParametersFrom(p );	

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
	
	SurfaceReport(theSystem,theDirective);

	return fret;
		
}


double surface_bwgrad_opt( HMMSystem &theSystem, HMMDirective &theDirective )
{	
	do_bwgrad = 1;

	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = surface_GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

	double *p = (double *)malloc( sizeof(double) * (n_params+1) );
	double *xi = (double *)malloc( sizeof(double) * (n_params+1) );

	do_crash_on_bad = 1;
	
	double f_obs[nstates];
	double f_t_ab[nstates];

	the_local_system->setFArray( 1.0, f_obs, f_t_ab );

	for( int s = 0; s < nstates; s++ )
	{
		for( int o = 0; o < globalN_OBSERVABLES; o++ )
			the_local_system->p_obs[o*nstates+s] = the_local_system->p_obs[o*nstates+s] * f_obs[s] + (1-f_obs[s]) *the_local_system->p_obs_BaumWelch[o*nstates+s];
	}

	
	if( !surface_ConvertSystemParametersTo( p ) )
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

	double fstart = surface_f(p);
//	printf("fstart: %le\n", fstart );	
	surface_df( p, xi );	

	double fret;

	//dlinmin( p, xi, n_params, &fret, f, df );
//	linmin( p, xi, n_params, &fret, f );


	int iter;
	dfpmin( p, n_params, 1e-14, &iter, &fret, surface_f, surface_df ); 
	//frprmn( p, n_params, 1e-14, &iter, &fret, f, df ); 

	int okay = surface_ConvertSystemParametersFrom(p);	

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

#if 1 
void surface_test_gradient( HMMSystem &theSystem, HMMDirective &theDirective, double *w_ref )
{	
	printf("Testing gradient.\n");
	int nstates = theSystem.getNStates();
	the_local_system = &theSystem;
	the_directive =  &theDirective;
	int n_params = surface_GetNParams(); //nstates*globalN_OBSERVABLES + nstates*nstates;	

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
	
	surface_ConvertSystemParametersTo( p );

	double fstart = surface_f(p);
	printf("fstart: %le\n", fstart );	
	surface_df( p, xi );	
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
	
		double f_new = surface_f(new_p);

		printf("%le %le %le ratio: %le\n", dG, f_new - fstart, new_E-fstart, (f_new - fstart)/(new_E-fstart) ); 

		free(new_p);

	}
}
#endif

#if 0 
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

#endif



#if 0
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
		MPI_Barrier(MPI_COMM_WORLD);
		exit(1);
		return 1e30;
	}
	else
		printf("The parameters were fine.\n");
	do_crash_on_bad = 0;
	
//	printf("Crash barrier 1.\n");
	MPI_Barrier(MPI_COMM_WORLD);
//	printf("After crash barrier 1");

	
	do_print_sigmas = 1;

	f(p);
		
}
#endif

void TrialRun( HMMSystem &theSystem, HMMDirective &theDirective, int time);
void getExpected( HMMSystem &theSystem, HMMDirective &theDirective );

void SurfaceReport( HMMSystem &theSystem, HMMDirective &theDirective )
{

	int nstates = theSystem.getNStates();
	for( int s = 0; s < nstates; s++ )
	{
		printf("State %s\n", theSystem.states[s]->label );

		double p_stay = theSystem.t_ab[s*nstates+s];
				

		printf("Stays in state p: %lf (1-p: %lf)\n", p_stay, 1-p_stay );
	}	

	
	TrialRun( theSystem, theDirective, 2477 ); 
	getExpected( theSystem, theDirective );
}

 


void TrialRun( HMMSystem &theSystem, HMMDirective &theDirective, int time)
{
	int nstates = theSystem.getNStates();
	int cur_state = rand() % nstates;

	double *t_ab = theSystem.t_ab;
	double *p_obs = theSystem.p_obs;
	printf("RUN: ");
	for( int t = 0; t < time; t++ )
	{
		double cum_p[globalN_OBSERVABLES];
		double prev = 0;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			cum_p[o] = p_obs[o*nstates+cur_state] + prev;
			prev = cum_p[o];
		}
		
		double rn = rand() / (double)RAND_MAX;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			if( rn <= cum_p[o] )
			{
				printf("%c", 'A' + o );
				break;
			}
		}

		double cum_t[nstates];
		prev = 0;
		for( int o = 0; o < nstates; o++ )
		{
			cum_t[o] = t_ab[cur_state*nstates+o] + prev;
			prev = cum_t[o];
		}
		
		rn = rand() / (double)RAND_MAX;

		for( int o = 0; o < nstates; o++ )
		{
			if( rn <= cum_t[o] )
			{
				cur_state = o;
				break;
			}
		}
		
	}	

	printf("\n");

}


void getExpected( HMMSystem &theSystem, HMMDirective &theDirective )
{
	// Run viterbi on the directives to count the number of real transitions

	int nstates = theSystem.getNStates();
	int nreal = 0;
	int is_real[nstates];
	memset( is_real, 0, sizeof(int) * nstates );
	int real_off[nstates];
	int is_hbonded[nstates];
	memset( is_hbonded, 0, sizeof(int) * nstates );

	for( int s = 0; s < theSystem.nstates; s++ )
	{
		real_off[s] = -1;
		if( !strncasecmp( theSystem.states[s]->label, "Border", 6 ) )
		{
			char code1 = theSystem.states[s]->label[7];
			char code2 = theSystem.states[s]->label[9];

			if( code1 >= 'a' && code1 <= 'z' ) is_hbonded[s] = 1;
			if( code2 >= 'a' && code2 <= 'z' ) is_hbonded[s] = 1;

			continue;
		}

		
		if( s >= globalN_OBSERVABLES/2 && s < globalN_OBSERVABLES )
			is_hbonded[s] = 1;

		real_off[s] = nreal;
		is_real[s] = 1;
		nreal++;
	}
	
	memset( is_hbonded, 0, sizeof(int) * nstates ); 

	
	char *real_codes = (char *)malloc( sizeof(char) * nreal );
	
	for( int x = 0; x < globalN_OBSERVABLES/2; x++ )
		real_codes[x] = 'A' + x;
	for( int x = 0; x < globalN_OBSERVABLES/2; x++ )
		real_codes[x+globalN_OBSERVABLES/2] = 'a' + x;

	double transitions[nreal*nreal];
	memset( transitions, 0, sizeof(double) * nreal*nreal );


	for( int d = 0; d < theDirective.n_data; d ++ )
	{
		int len = theDirective.translated_len[d];
		
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
	
	
		calculateAlphaForStringConstrained( alpha_constrained, alpha_constrained_norm, theDirective.translated_observable_string[d], theDirective.raw_class_string[d],  len, theSystem );
		calculateBetaForStringConstrained( beta_constrained, beta_constrained_norm,    theDirective.translated_observable_string[d],  theDirective.raw_class_string[d],  len, theSystem );
		
		int *predicted_state = (int *)malloc( sizeof(int) * len );
		
		theSystem.estimateClassViterbi( alpha_constrained, beta_constrained, theDirective.translated_observable_string[d], theDirective.raw_class_string[d], len, predicted_state  );
		
		int cur_state_real = -1; 

		printf("d: %d ", d   );

		// map hydrogen bonded states back to the original.
		int state_map[nreal];
		for( int i = 0; i < globalN_OBSERVABLES/2; i++ )
		{		
			state_map[i] = i;
//			state_map[i+globalN_OBSERVABLES/2] = i;
			state_map[i+globalN_OBSERVABLES/2] = i+globalN_OBSERVABLES/2;
		}

		/*
			transitions --
				I am trying to figure out the curvature-dependent energy on the state.
				I do not want hydrogen bonding to be an influence.
				I will leave out all h-bonded border and h-bonded states from the transitions.
		

		*/

		for( int t = 0; t < len; t++ )
		{
			if( is_real[predicted_state[t]] )
			{
				if( is_hbonded[predicted_state[t]] || is_hbonded[cur_state_real] )
				{
				}
				else
				{
					if( cur_state_real >= 0 && cur_state_real != predicted_state[t] )
						transitions[real_off[cur_state_real]*nreal+real_off[predicted_state[t]]] += 1;
					else if( cur_state_real >= 0 )
						transitions[real_off[cur_state_real]*nreal+real_off[cur_state_real]] += 1;
				}
	
				cur_state_real = predicted_state[t];
				printf("%c", real_codes[real_off[cur_state_real]] );
			}
			else if( cur_state_real >= 0 )
			{
				if( is_hbonded[predicted_state[t]] || is_hbonded[cur_state_real] )
				{
				}
				else
				{
					transitions[real_off[cur_state_real]*nreal+real_off[cur_state_real]] += 1;
				}
				printf("." );
			}
		}
		printf("\n");

		free(predicted_state);
		free(alpha_norm);
		free(beta_norm);
		free(alpha_constrained_norm);
		free(beta_constrained_norm);
		free(alpha);
		free(beta);
		free(alpha_constrained);
		free(beta_constrained);
		free(wild_string);
	}

	double time_per_click = 0.5; // nanoseconds



	printf("SRC");
	for( int r = 0; r < nreal; r++ )
		printf(" %c", real_codes[r] );
	printf("\n");
	for( int r = 0; r < nreal; r++ )
	{
		printf("%c", real_codes[r] );
	
		double ntrans = 0;
		double nstay = transitions[r*nreal+r];

		for( int r2 = 0; r2 < nreal; r2++ )
		{
			printf(" %lf", transitions[r*nreal+r2] );		
			if( r2 >0 && r2 != r )
				ntrans += transitions[r*nreal+r2]; 
		}
		printf(" lifetime %lf\n", time_per_click / ( (ntrans/(ntrans+nstay)) ) );
	}

	double expected[globalN_OBSERVABLES];
	memset(expected, 0, sizeof(double) * globalN_OBSERVABLES );

	int cur_state = rand() % nstates;

	double *t_ab = theSystem.t_ab;
	double *p_obs = theSystem.p_obs;

	int max_time = 1000000;
	
	double certain_occ[nreal];
	memset( certain_occ, 0, sizeof(double) * nreal );

	for( int t = 0; t < max_time; t++ )
	{
		double cum_p[globalN_OBSERVABLES];
		double prev = 0;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			cum_p[o] = p_obs[o*nstates+cur_state] + prev;
			prev = cum_p[o];
		}
		
		double rn = rand() / (double)RAND_MAX;

		for( int o = 0; o < globalN_OBSERVABLES; o++ )
		{
			if( rn <= cum_p[o] )
			{
				expected[o] += 1.0;	
				break;
			}
		}

		double cum_t[nstates];
		prev = 0;
		for( int o = 0; o < nstates; o++ )
		{
			cum_t[o] = t_ab[cur_state*nstates+o] + prev;
			prev = cum_t[o];
		}
		
		rn = rand() / (double)RAND_MAX;

		for( int o = 0; o < nstates; o++ )
		{
			if( rn <= cum_t[o] )
			{
				cur_state = o;
				break;
			}
		}

		//if( cur_state < globalN_OBSERVABLES/2 )
		if( cur_state < globalN_OBSERVABLES )
			certain_occ[cur_state] += 1;
	}	

	double fsum = 0;

	double found[globalN_OBSERVABLES];
	memset( found, 0, sizeof(double) * globalN_OBSERVABLES );
	
	for( int d = 0; d < theDirective.n_data; d ++ )
	{
		int len = theDirective.translated_len[d];
		
		
		for( int t = 0; t < len; t++ )
		{
			int nstate = theDirective.translated_observable_string[d][t];
			found[nstate] += 1;
		}					
	}


	for( int o = 0; o < globalN_OBSERVABLES; o++ )
		fsum += found[o]; 
	
	double sum = 0;
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
		sum += expected[o];
	
	double csum = 0;
	for( int o = 0; o < globalN_OBSERVABLES/2; o++ )
		csum += certain_occ[o];
	
	for( int o = 0; o < globalN_OBSERVABLES; o++ )
	{	
		if( found[o] > 0 || expected[o] > 0 )
		printf("%c obs %lf modeled %lf", real_codes[o], found[o] / fsum, expected[o] / sum );
		if( o < globalN_OBSERVABLES/2 )
			printf(" certain occ %lf", certain_occ[o] / csum );
		printf("\n");
	}

}


