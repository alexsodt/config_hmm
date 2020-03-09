#ifndef __comparisonh__

#ifdef __comparisonc__
double binary_penalty = 0.1;
double binary_benefit = -0.5;
double hbond0_penalty = 100;
double hbond0_benefit = -100;
int comparison_n_binary = 31;
int use_atoms=0;
#else
extern int comparison_n_binary;
extern double hbond0_penalty;
extern double hbond0_benefit;
extern double binary_penalty;
extern double binary_benefit;
extern int use_atoms;
#endif

struct aStructure {
	long binary_data;
	double *atoms;
	int cluster;
};

void doKMedoidClusteringLR( aStructure *xyz, int nxyz, int nat, int KMED, int *outputMedoids, double r, double theta, int kmax);

double value_comparison( aStructure *data1, aStructure *data2, int do_swap_config );
double getChi2( double *xyz1, 
			double *xyz2,
			int nat );
double getChi2Symm( double *, double *, int nat );

int useAtomInChi2( struct atom_rec *at );
int binaryHBonder( struct atom_rec *at1, struct atom_rec *at2); // is this a binary hbond motif?

struct bond_data
{
	const char *atName1; 
	const char *atName2;
	int res1;
	int res2;
};

#ifndef __spec_comparisonc__
extern bond_data allBonders[];
#endif
#endif

