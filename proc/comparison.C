#define __comparisonc__

#include "comparison.h"

double getChi2Symm( double *xyz1, double *xyz2, int nat )
{
	double tmp1[3*nat];

	double chi2_a = getChi2( xyz1, xyz2, nat );

	for( int x = 0; x < nat/2; x++ )
	{
		tmp1[(nat/2+x)*3+0] = xyz1[3*x+0];
		tmp1[(nat/2+x)*3+1] = xyz1[3*x+1];
		tmp1[(nat/2+x)*3+2] = xyz1[3*x+2];
		
		tmp1[x*3+0] = xyz1[3*(nat/2+x)+0];
		tmp1[x*3+1] = xyz1[3*(nat/2+x)+1];
		tmp1[x*3+2] = xyz1[3*(nat/2+x)+2];
	}
	double chi2_b = getChi2( tmp1, xyz2, nat );

	if( chi2_a < chi2_b ) return chi2_a;

	return chi2_b;
}

double value_comparison( aStructure *data1, aStructure *data2, int do_swap_config )
{
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
		chi2 = getChi2Symm( data1->atoms, data2->atoms, use_atoms );
	else
		chi2 = getChi2( data1->atoms, data2->atoms, use_atoms );

	value += chi2;

	return value;
}
 
