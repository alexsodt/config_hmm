
/*
	uses the GSL minimizers.
*/

#include "mpi.h"
#include "gsl/gsl_multimin.h"
#include "GSL_interface.h"

int GetNParams( void );
double f( double * );
double fdf( double *, double *);
extern double global_f;
int GSL_ITMAX = 500;
int ConvertSystemParametersTo( double *w, double fp );
void do_NR_BFGS( void )
{

}

double GSL_f_wrapper( const gsl_vector *v, void *params )
{
	int NP = GetNParams();

	double input_vals[1+NP];

	for( int x = 0; x < NP; x++ )	
		input_vals[1+x] = gsl_vector_get( v, x );

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Here at f.\n");
	MPI_Barrier(MPI_COMM_WORLD);
	double ret = f(input_vals);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Done with f.\n");
	MPI_Barrier(MPI_COMM_WORLD);
	
	return ret;
}

void GSL_df_wrapper( const gsl_vector *v, void *params, gsl_vector *g )
{
	int NP = GetNParams();

	double input_vals[1+NP];
	double output_vals[1+NP];

	for( int x = 0; x < NP; x++ )	
		input_vals[1+x] = gsl_vector_get( v, x );
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Here at df.\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	fdf( input_vals, output_vals );
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Done with df.\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	for( int x = 0; x < NP; x++ )
		gsl_vector_set( g, x, output_vals[1+x] );
}

void GSL_fdf_wrapper( const gsl_vector *v, void *params, double *fret, gsl_vector *g )
{
	int NP = GetNParams();

	double input_vals[1+NP];
	double output_vals[1+NP];

	for( int x = 0; x < NP; x++ )	
		input_vals[1+x] = gsl_vector_get( v, x );
	
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Here at fdf.\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	*fret = fdf( input_vals, output_vals );
	MPI_Barrier(MPI_COMM_WORLD);
	printf("Done with fdf.\n");
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	for( int x = 0; x < NP; x++ )
		gsl_vector_set( g, x, output_vals[1+x] );
}

void do_GSL_minimization( int use_GSL_Minimizer )
{
	int np = GetNParams();

	double start_pt[1+np];
	gsl_vector *vector_start_pt = gsl_vector_alloc( np );

	ConvertSystemParametersTo( start_pt, 1-global_f );
	
	for( int x = 0; x < np; x++ )
		gsl_vector_set( vector_start_pt, x, start_pt[1+x] );

	if( use_GSL_Minimizer == USE_GSL_SIMPLEX )
	{
		gsl_multimin_function my_func;

		my_func.n = np;
		my_func.f = &GSL_f_wrapper;
		my_func.params = NULL;
		
	       	gsl_vector *step_size_vector = gsl_vector_alloc (np);
       		gsl_vector_set_all(step_size_vector, 1);

		gsl_multimin_fminimizer * minimizer = gsl_multimin_fminimizer_alloc( gsl_multimin_fminimizer_nmsimplex, np); 
		gsl_multimin_fminimizer_set( minimizer, &my_func, vector_start_pt, step_size_vector );
		
		int iter = 0;
		int status;
		double size;
		do {
			iter++;
			status = gsl_multimin_fminimizer_iterate(minimizer);

			if( status) break;

			size = gsl_multimin_fminimizer_size( minimizer );
			status = gsl_multimin_test_size(size,  1e-2 );
		
			if( status = GSL_SUCCESS )
				printf("GSL minimum reached.\n");

		} while( status == GSL_CONTINUE && iter < GSL_ITMAX );

		gsl_vector_free( step_size_vector );
		gsl_multimin_fminimizer_free( minimizer );
	}
	else
	{	
		gsl_multimin_function_fdf my_func;

		my_func.n = np;
		my_func.f = &GSL_f_wrapper;
		my_func.df = &GSL_df_wrapper;
		my_func.fdf = &GSL_fdf_wrapper;
		my_func.params = NULL;
		gsl_multimin_fdfminimizer* minimizer = NULL;
		switch( use_GSL_Minimizer )
		{
			case USE_GSL_FR:
				minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr , np );
				break;
			case USE_GSL_PR:
				minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr , np );
				break;
			case USE_GSL_BFGS:
				minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2 , np );
				break;
			case USE_GSL_SD:
				minimizer = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent , np );
				break;
			default:
				printf("Unknown GSL minimization routine requested.\n");
				exit(1);
				break;
		}

		int status;
//		int outer_iter = 0;
//		do {
//			outer_iter++;
			gsl_multimin_fdfminimizer_set( minimizer, &my_func, vector_start_pt, 100, 0.1 );

			int iter = 0;
			do {
				iter++;
				printf("Doing optimization round.\n");
				status = gsl_multimin_fdfminimizer_iterate(minimizer);
				MPI_Barrier(MPI_COMM_WORLD);
				printf("Done.\n");
				fflush(stdout);
				MPI_Barrier(MPI_COMM_WORLD);
	//			printf("Stopping at %le.\n", gsl_multimin_fdfminimizer_minimum(minimizer) );
				printf("Iteration Status: %d\n", status );
				if( status) break;
	
				status = gsl_multimin_test_gradient( minimizer->gradient, 1e-3 );
				printf("Grad Status: %d\n", status );
			
				if( status == GSL_SUCCESS )
					printf("GSL minimum reached.\n");
	
			} while( status == GSL_CONTINUE && iter < GSL_ITMAX );

		//	gsl_multimin_fdfminimizer_restart(minimizer);	
		
//			ConvertSystemParametersTo( start_pt, 1-global_f );
		
//			for( int x = 0; x < np; x++ )
//				gsl_vector_set( vector_start_pt, x, start_pt[1+x] );
				
//		} while( status == GSL_ENOPROG && outer_iter < 10);

		gsl_multimin_fdfminimizer_free( minimizer );
	}

	gsl_vector_free(vector_start_pt);

}

