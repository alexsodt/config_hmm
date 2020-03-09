#ifndef __proc_definitionh__
#define __proc_definitionh__

#define USE_ANYWHERE	0
#define USE_CHI2	1

int info_mode( void );
int is_sym( void);
int is_seg( int pair );
int use_atom( int use_type, int seg, int res, const char *atname, const char *resname, const char *segname );
int binaryHBonder( struct atom_rec *at1, struct atom_rec *at2); // is this a binary hbond motif?
void load_pair_definition( const char *pairDefinitionFileName );
double pair_binary_penalty(void);
double pair_binary_benefit(void);
double pair_hbond0_penalty(void);

#endif
