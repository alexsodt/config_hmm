#ifndef __INPUT_FILE_H__
#define __INPUT_FILE_H__

#include "HMM.h"

#define MAX_DIRECTIVES 1000

#define COMMAND_BIT	(1<<12)
#define NETWORK_BIT	(1<<13)
#define INTERPRETER_BIT	(1<<14)

#define COMMAND_TRAIN	    (COMMAND_BIT+0)
#define COMMAND_DECODE	    (COMMAND_BIT+1)
#define COMMAND_TRAIN_FUZZY (COMMAND_BIT+2)
#define COMMAND_TRAIN_GRADIENT (COMMAND_BIT+3)
#define COMMAND_LOAD	    (COMMAND_BIT+4)
#define COMMAND_SAVE	    (COMMAND_BIT+5)
#define COMMAND_COPY	    (COMMAND_BIT+6)
#define COMMAND_ASSESS	    (COMMAND_BIT+7)
#define COMMAND_CHECK	    (COMMAND_BIT+8)
#define COMMAND_TRAIN_EXT   (COMMAND_BIT+9)
#define COMMAND_PSEUDO	    (COMMAND_BIT+10)
#define COMMAND_TESTG	    (COMMAND_BIT+11)
#define COMMAND_GEN	    (COMMAND_BIT+12)
#define COMMAND_ALPHA	    (COMMAND_BIT+13)
#define COMMAND_NHELICES    (COMMAND_BIT+14)
#define COMMAND_TRAIN_BWGRAD	    (COMMAND_BIT+15)
#define COMMAND_ASSIGN_IO   (COMMAND_BIT+16)
#define COMMAND_ACTIVATE_N  (COMMAND_BIT+17)
#define COMMAND_DEACTIVATE_N (COMMAND_BIT+18)
#define COMMAND_SEARCH	    (COMMAND_BIT+19)
#define COMMAND_FIX	    (COMMAND_BIT+20)
#define COMMAND_SC_WEIGHT   (COMMAND_BIT+21)
#define COMMAND_MISC1	    (COMMAND_BIT+22)
#define COMMAND_MISC2	    (COMMAND_BIT+23)
#define COMMAND_MISC3	    (COMMAND_BIT+24)
#define COMMAND_MISC4	    (COMMAND_BIT+25)
#define COMMAND_MISC5	    (COMMAND_BIT+26)
#define COMMAND_MISC6	    (COMMAND_BIT+27)
#define COMMAND_IGNORE   (COMMAND_BIT+28)
#define COMMAND_REPORT	    (COMMAND_BIT+29)
#define COMMAND_G_WEIGHT   (COMMAND_BIT+30)
#define COMMAND_S_WEIGHT   (COMMAND_BIT+31)
#define COMMAND_MISC7	    (COMMAND_BIT+32)
#define COMMAND_GRAD_PREP   (COMMAND_BIT+33)
#define COMMAND_SA	    (COMMAND_BIT+34)
#define COMMAND_SORT	    (COMMAND_BIT+35)
#define COMMAND_CLEAR_AQUEOUS (COMMAND_BIT+36)
#define COMMAND_ANNOTATE    (COMMAND_BIT+37)
#define COMMAND_SURFK	    (COMMAND_BIT+38)
#define COMMAND_OUTER_TRAIN	    (COMMAND_BIT+39)
#define COMMAND_OUTER_DECODE	    (COMMAND_BIT+40)
#define COMMAND_OUTER_ASSESS	    (COMMAND_BIT+41)
#define COMMAND_OUTER_ASSIGN	    (COMMAND_BIT+42)
#define COMMAND_MAKE_FUZZY	    (COMMAND_BIT+43)
#define COMMAND_ENHANCE	    (COMMAND_BIT+44)
#define COMMAND_START	  (COMMAND_BIT+45)
#define COMMAND_STOP	  (COMMAND_BIT+46)
#define COMMAND_RANDOMT	  (COMMAND_BIT+47)
#define COMMAND_GRAD	  (COMMAND_BIT+48)

#define NETWORK_LOLD	  	(NETWORK_BIT+0) 

#define INTERPRETER_LOLD  	(INTERPRETER_BIT+0)
//#define INTERPRETER_SCHEMES	(INTERPRETER_BIT+2) // ten schemes..


struct HMMDirective
{
	int command;
	int network;
	int interpreter;

	int n_data;
	int n_data_space;

	char **comment;
	char **aux_string;
	border_assignment *borders;
	char **raw_observable_string;
	char **raw_class_string;
	double **train_zval;
	int *protein;
	int cur_protein;

	double *weights;

	int **translated_observable_string;
	int *translated_len;

	HMMDirective();
	
	void addData( char *sequence, char *theClass=NULL, double *zvals=NULL, double the_weight=1.0,char *auxData=NULL, char *comment=NULL  );
	void doCommand( int c_id );
	void useTranslator( int trans );
	void useNetwork( int network );
	void translate(void);
};

struct HMMCommands
{
	HMMDirective *allDirectives[MAX_DIRECTIVES];
	int ndirectives;

	void addDirective( HMMDirective *theDirective );
	HMMCommands( FILE *theFile );
};

int commandInterpreter( const char *read );

#endif
