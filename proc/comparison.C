#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "alignSet.h"

#define __comparisonc__

#include "comparison.h"

int debug_trigger = 0;
int current_letter = 0;

double getChi2Symm( double *xyz1, double *xyz2, int nat )
{
	double tmp1[3*nat];
	double tmp2[3*nat];

	memcpy( tmp1, xyz1, 3 * nat * sizeof(double) );
	memcpy( tmp2, xyz2, 3 * nat * sizeof(double) );
	

	double chi2_a = getChi2( tmp1, tmp2, nat );
	
	if( debug_trigger  && (current_letter == 'Y'-'A' || current_letter == 26 + 'i'-'a'  ) )
	{
		int all_atoms[nat];
		for( int i = 0; i < nat; i++ )
			all_atoms[i] = i;

		alignStructuresOnAtomSet( tmp1, all_atoms, tmp2, all_atoms, nat, nat );

		printf("%d\n", nat );
		printf("comment\n");
		for( int x = 0; x < nat; x++ )
			printf("%c %lf %lf %lf\n", (x < nat/2 ? 'C' : 'O'), tmp1[3*x+0], tmp1[3*x+1], tmp1[3*x+2] );
		printf("%d\n", nat );
		printf("comment\n");
		for( int x = 0; x < nat; x++ )
			printf("%c %lf %lf %lf\n", (x < nat/2 ? 'N' : 'F'), tmp2[3*x+0], tmp2[3*x+1], tmp2[3*x+2] );
	}
	
	memcpy( tmp1, xyz1, 3 * nat * sizeof(double) );
	memcpy( tmp2, xyz2, 3 * nat * sizeof(double) );


	for( int x = 0; x < nat/2; x++ )
	{
		tmp1[(nat - nat/2+x)*3+0] = xyz1[3*x+0];
		tmp1[(nat - nat/2+x)*3+1] = xyz1[3*x+1];
		tmp1[(nat - nat/2+x)*3+2] = xyz1[3*x+2];
		
		tmp1[x*3+0] = xyz1[3*(nat - nat/2+x)+0];
		tmp1[x*3+1] = xyz1[3*(nat - nat/2+x)+1];
		tmp1[x*3+2] = xyz1[3*(nat - nat/2+x)+2];
	}
	
	double chi2_b = getChi2( tmp1, tmp2, nat );
	
	if( debug_trigger  && (current_letter == 'Y'-'A' || current_letter == 26 + 'i'-'a' ) )
	{
		int all_atoms[nat];
		for( int i = 0; i < nat; i++ )
			all_atoms[i] = i;

		alignStructuresOnAtomSet( tmp1, all_atoms, tmp2, all_atoms, nat, nat );
		printf("and..\n");
		printf("%d\n", nat );
		printf("comment\n");
		for( int x = 0; x < nat; x++ )
			printf("%c %lf %lf %lf\n", (x < nat/2 ? 'C' : 'O'), tmp1[3*x+0], tmp1[3*x+1], tmp1[3*x+2] );
		printf("%d\n", nat );
		printf("comment\n");
		for( int x = 0; x < nat; x++ )
			printf("%c %lf %lf %lf\n", (x < nat/2 ? 'N' : 'F'), tmp2[3*x+0], tmp2[3*x+1], tmp2[3*x+2] );
	}

	if( debug_trigger  && (current_letter == 'Y'-'A' || current_letter == 26 + 'i'-'a' ) )
	{
		printf("chi2: %le %le\n", chi2_a, chi2_b );
		if( current_letter == 26 + 'i'-'a' )
			exit(1);
	}
	if( chi2_a < chi2_b ) return chi2_a;

	return chi2_b;
}

double value_comparison( aStructure *data1, aStructure *data2, int do_swap_config )
{
	double tmp1[3*use_atoms];
	double tmp2[3*use_atoms];

	memcpy( tmp1, data1->atoms, sizeof(double) * use_atoms * 3 );
	memcpy( tmp2, data2->atoms, sizeof(double) * use_atoms * 3 );

	double value = 0;

	for( int i = 0; i < comparison_n_binary; i++ )
	{
		if( (data1->binary_data & (1<<i)) && (data2->binary_data & (1<<i)) )
			value += ( i == 0 ? hbond0_benefit : binary_benefit);
		else if( (data1->binary_data & (1<<i)) || (data2->binary_data & (1<<i)) )  
			value += ( i == 0 ? hbond0_penalty : binary_penalty );		
	}
	
	double chi2 = 0;
	if( do_swap_config )
		chi2 = getChi2Symm( tmp1, tmp2, use_atoms );
	else
		chi2 = getChi2( tmp1, tmp2, use_atoms );

	value += chi2;

	return value;
}
 
