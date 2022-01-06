#ifndef __proc_definitionh__
#define __proc_definitionh__

#define USE_ANYWHERE	0
#define USE_CHI2	1

int info_mode( void );
int is_sym( void);
int is_seg( int pair );
int nextra( void );
int force_hbonds( void );
int use_atom( int use_type, int seg, int res, const char *atname, const char *resname, const char *segname );
int binaryHBonder( struct atom_rec *at1, struct atom_rec *at2); // is this a binary hbond motif?
void load_pair_definition( const char *pairDefinitionFileName );
int used_as_ion( int res, const char *resname, const char *segname, double *check_rval );
double pair_binary_penalty(void);
double pair_binary_benefit(void);
double pair_hbond0_penalty(void);
int allowNoIons( void );

double getCutoff( void );
int getPDefIndex( const char *atName, const char *resName, const char *segName, int res );
void clearPDefPSF( void );
int loadPDefPSF( void );


#endif
