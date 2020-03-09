
#ifndef __HMM_H__
#define __HMM_H__

//#define PROTEIN_STATE
struct HMMDirective;

#define MAX_OBSERVABLES		50
#define MAX_TRANSITIONS		50
//#define DO_TRANSITION_GROUPS
#define MAX_AUX_VARS		10

struct HMMState;

struct prob
{
	HMMState *elem;
	double p;
};

#define MAX_SEGMENTS 128

struct border_assignment
{
	int nsegs;

	int start_pos[MAX_SEGMENTS];
	int stop_pos[MAX_SEGMENTS];	
	int topology;
};

struct HMMState
{
	char label[256];
	char theClass;
	int pseudostate;
	int emission_group;
	int transition_group;
	int index;
	int possibleInitialState;
	int ntransitions; // number of transitions present
	prob p_so[MAX_TRANSITIONS]; // transition probabilities, probability of transitioning to this state.
	double f_ys[MAX_OBSERVABLES]; // emission probabilities, probability of observing y in this state.	
	double obs_population[MAX_OBSERVABLES];
	double obs_population_copy[MAX_OBSERVABLES];

	HMMState(char theClass, const char *label, int possibleInitialState, int emission_group, int pseudostate=1, int transition_group=-1 );

	void addTransition( HMMState *theState );
};

struct HMMSystem
{
	int nstates;
	int nstatesSpace;

	int **transitionListFrom;
	int *transitionSizeFrom;
	
	int **transitionListTo;
	int *transitionSizeTo;

	int nEmissionGroups;
	int **emission_groups;
	int *emission_groups_len;
	
	int nTransitionGroups;
	int **transition_groups;
	int *transition_groups_len;

	int *group_is_fixed;

	double *p_obs_prev;
	double *t_ab_prev;

	double *p_obs;
	double *t_ab;
	
	double *p_obs_BaumWelch;
	double *t_ab_BaumWelch;
	
	double *p_obs_next_free;
	double *t_ab_next_free;
	double *aux_vars[MAX_AUX_VARS];
	
	double *p_obs_next_clamped;
	double *t_ab_next_clamped;
	
	double *p_obs_next_free_log_norm;
	double *t_ab_next_free_log_norm;
	
	double *p_obs_next_clamped_log_norm;
	double *t_ab_next_clamped_log_norm;
	
	double *init_p;
	double *init_p_next_clamped;
	double *init_p_next_free;
	double *init_p_BaumWelch;

	double *pseudocounts;
	double *pseudo_add_clamped;
	double *pseudo_add;

	HMMState **states;

	HMMSystem();
	~HMMSystem();
	void addState( HMMState *theState );

	void getTransitionsFrom( int s, int *state_list, int *ntransitions, char theClass );
	void getTransitionsTo( int s, int *state_list, int *ntransitions, char theClass );

	int getNStates(void);
	double *getObservableP( void );
	double *getTransitionP( void );
	double *getInitP( void );

	void zeroNextPMat( void );
	void incrementNextPMat_GEN( double *alphas, double *betas, int *observables, int len, double *t_ab, double *p_obs, double *init_p, double weight, double scale=1.0 );
	void incrementNextPMat_Free( double *alphas, double *betas, int *observables, int len, double weight, double scale=1.0 );
	void incrementNextPMat_Clamped( double *alphas, double *betas, int *observables, int len, double weight, double scale=1.0 );
	void finishNextPMat( void );
	double setNextPMat( double factor, int doBaumWelch, double free_multiplier, double A );
	void setGradient( double *dw);
	void setFractionalGradient( double *dw, double f);
	void setFractionalGradientSplit( double *dw, double f_in);
	void checkNormalization( void );
	void setup();
	void Normalize(void);
	void addNoise( double A);
	int setNextPMatByGradient( double factor  );
	void SaveModel(const char *fileName);
	void ReadModel(const char *fileName);
	double estimateClass1Best( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	int estimateClassViterbi( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	int estimateClassViterbiSort( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *h1, double *h2, double *h3, double *h4 );
	int comparePredictedAndTruePaths( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	double ratio12Viterbi( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *p );
	double getp7( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state, double *p7, int to_proc );
	int estimateClassPlayground( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	int estimateClassPlayground2( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	void checkEGConsistency( void );
	void CommunicatePMats(void);
	void CommunicateShortPMats( void );
	void CommunicateGradients( double *);
	void ReportEmissionGroups( int translator );
	void measurePseudocounts( int *observables, char *classes, int len );
	void addPseudocounts( double *pseudo_clamped, double *pseudo_free );
	void addPseudocounts( double frac );
	void randomizeObs( void );
	void zeroPseudocounts( void );
	void fixPseudocounts( void );
	void generateSequence( int *observables, char *classes, int len );
	void swapStatePopulations( void );
	void addStatePopulations( void );
	void resetStatePopulations( int do_all=0);
	void printStatePopulations( int trans);
	void printBorderStatePopulations( int trans );
	void FreezeGlobular( void );
	void SetHeuristicP0( void );
	void SetGenericParameters( void );
	void setPseudocounts( HMMDirective * theDirective );
	double MapP0ObsFactor( double f, int s );
	double MapP0TFactor( double f, int s );
	void setFArray( double f, double *obs_f, double *t_f);
	int estimateClassThreeStater( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	void ProbabilityCompress( void);
double rankBasedOnStatePopulations( 
		double *alpha,
		double *beta,
		int len,
		char *state1,
		 char *state2 = NULL,
		 char *state3 = NULL,
		 char *state4 = NULL,
		 char *state5 = NULL,
		 char *state6 = NULL,
		 char *state7 = NULL,
		 char *state8 = NULL,
		 char *state9 = NULL,
		 char *state10 = NULL );
	void LoadFixASCII( const char *fileName );
	double ProbActual( double *alphas, double *betas, int *observables, char *theClass, int len, int *ML_state );
	double setNextPMatCF( double factor, double A, double alpha );
	void SAPerturb(double T);
	double BaumWelchIteration( HMMDirective *theDirective, int doUpdate ); 
	void TransitionGroupCheck(void);
	double outerModuleProbability( void );
	void printSurfaceInfo( int *sequence, char *boundaries, int len );
	void SurfaceSimulatedAnnealing( int *sequence, int len, border_assignment *border, int do_train );
	void SurfaceShortLoop( int *sequence, char *topo, int len, border_assignment *border );
	double SequenceProbability( int *sequence, border_assignment *borders, int len );
	void writeBorderProbabilities( double *borderProb, int *sequence, int len );
	void CopyViterbiBorders( border_assignment *borders, int *sequence, int len ); 
	void randomizeTimeSequence( HMMDirective *directive );
};


void calculateAlphaForString( double *alpha, double *alpha_norm, int *observables, int len, HMMSystem &theSystem );
void calculateBetaForString( double *beta, double *alpha_norm, int *observables, int len, HMMSystem &theSystem );
void calculateAlphaForStringConstrained( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem );
double isItProbable( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem );
void calculateBetaForStringConstrained( double *beta, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem );
void calculateAlphaForStringForcePath( double *alpha, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem, int *the_path );
void calculateBetaForStringForcePath( double *beta, double *alpha_norm, int *observables, char *classes, int len, HMMSystem &theSystem, int *the_path );

#endif

