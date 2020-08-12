#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "ctype.h"

#include "util.h"
#include "HMMDirective.h"
extern int globalN_OBSERVABLES;

int isnum( char t )
{
	if( t >= '0' && t <= '9' )
		return 1;
	return 0;
}

HMMDirective::HMMDirective()
{
	command = -1;
	network = -1;
	interpreter = -1;
	
	n_data = 0;
	n_data_space = 10;

	comment = (char **)malloc( sizeof( char *) * n_data_space );
	raw_observable_string = (char **)malloc( sizeof( char *) * n_data_space );
	raw_class_string      = (char **)malloc( sizeof( char *) * n_data_space );
	translated_len      = (int *)malloc( sizeof( int) * n_data_space );
	translated_observable_string      = (int **)malloc( sizeof( int *) * n_data_space );
	aux_string      = (char **)malloc( sizeof( char *) * n_data_space );
	weights = (double *)malloc( sizeof(double) * n_data_space );
	train_zval = (double **)malloc( sizeof(double*) * n_data_space );
	protein = (int *)malloc( sizeof(int) * n_data_space );
	borders = (border_assignment *)malloc( sizeof(border_assignment) * n_data_space );
	cur_protein = 0;
}

void HMMDirective::addData( char *sequence, char *theClass, double *zvals,  double the_weight, char *auxData, char *comment_in )
{
	// sequence may be an aggregate..
		
	int done = 0;
	int cur_spot = 0;

//	if( theClass && theClass[strlen(theClass)-1] == 'B' )
//		theClass[strlen(theClass)-1] = 'X';

	while( cur_spot < strlen(sequence) )
	{
		int tstart = cur_spot;
		int tstop = strlen(sequence);

		if( theClass )
		{
/*			if( theClass[cur_spot] == 'B' ) theClass[cur_spot] = 'X';
			for( int x = tstart; x < tstop; x++ )
				if( theClass[x] == 'B' )
				{
					tstop = x;
					break;
				}
*/
		}
		
			
		if( n_data == n_data_space )
		{
			n_data_space *= 2;
	
			comment = (char **)realloc( comment, sizeof( char *) * n_data_space );
			raw_observable_string = (char **)realloc( raw_observable_string, sizeof( char *) * n_data_space );
			raw_class_string      = (char **)realloc( raw_class_string,      sizeof( char *) * n_data_space );
			aux_string      = (char **)realloc( aux_string,      sizeof( char *) * n_data_space );
			weights =               (double *)realloc( weights, sizeof(double) * n_data_space );
			translated_observable_string = (int **)realloc( translated_observable_string, sizeof(int*) * n_data_space );
			translated_len = (int *)realloc( translated_len, sizeof(int) * n_data_space );
			train_zval = (double **)realloc( train_zval, sizeof(double) * n_data_space );
			protein = (int *)realloc( protein, sizeof(int) * n_data_space );
			borders = (border_assignment *)realloc( borders, sizeof(border_assignment) * n_data_space );
		}
	
		aux_string[n_data] = NULL;
	
		int tlen = tstop - tstart;

		raw_observable_string[n_data] = ( char *) malloc( sizeof(char) * (tlen+1) );
		memcpy( raw_observable_string[n_data], sequence+tstart,tlen*sizeof(char) );
		raw_observable_string[n_data][tlen] = '\0';
		weights[n_data] = the_weight;
	
		if( theClass )
		{
			raw_class_string[n_data] = (char *)malloc( sizeof(char) * (tlen+1) );
			memcpy( raw_class_string[n_data], theClass+tstart, sizeof(char)*(tlen+1) );
			raw_class_string[n_data][tlen] = '\0';
		}

		if( zvals )
		{
			train_zval[n_data] = (double *)malloc( sizeof(double) * tlen );
			memcpy( train_zval[n_data], zvals+tstart, tlen*sizeof(double) );
		}

		if( comment_in )
		{
			comment[n_data] = (char *)malloc( sizeof(char) * (strlen(comment_in)+1) );
			strcpy( comment[n_data], comment_in );
		}
		else
			comment[n_data] = NULL;

		translated_observable_string[n_data] = NULL;
		translated_len[n_data] = 0;
		protein[n_data] = cur_protein; 
		
		n_data++;

		cur_spot = tstop;
	}

	cur_protein++;
}

void HMMDirective::useTranslator( int trans_in )
{
	interpreter = trans_in;
}

void HMMDirective::useNetwork( int network_in )
{
	network = network_in;
}

void HMMDirective::doCommand( int command_in )
{
	command = command_in;
}

int commandInterpreter( const char *read )
{
	int got_big_scheme = 0;

	if( !strncasecmp( read, "PSEUDO", 6 ) )
		return COMMAND_PSEUDO;
	if( !strncasecmp( read, "START", 5 ) )
		return COMMAND_START;
	if( !strncasecmp( read, "SAVE", 4 ) )
		return COMMAND_SAVE;
	if( !strncasecmp( read, "LOAD", 4 ) )
		return COMMAND_LOAD;
	if( !strncasecmp( read, "RANDOMT", 7 ) )
		return COMMAND_RANDOMT;
	if( !strncasecmp( read, "STOP", 4 ) )
		return COMMAND_STOP;
	if( !strncasecmp( read, "TRAIN", 5  ) )
		return COMMAND_TRAIN;
	if( !strncasecmp( read, "GRAD", 4  ) )
		return COMMAND_GRAD;
	if( !strncasecmp( read, "DECODE", 6  ) )
		return COMMAND_DECODE;
	if( !strncasecmp( read, "LOLD", 4 ) )
		return NETWORK_LOLD;
	else if( !strncasecmp( read, "NUM", 3 ) )
	{
		int scheme = atoi( read+3 );

		return INTERPRETER_LOLD + scheme;
	}
	else
	{
		printf("Could not interpret '%s'.\n", read );
		exit(1);
	}

	return -1;

}

int lipid_translator( char tchar )
{
	int nstates = globalN_OBSERVABLES/2;

	if( tchar >= 'A' && tchar <= 'Z' )
		return tchar - 'A';
	if( tchar >= 'a' && tchar <= 'z' )
		return nstates + tchar - 'a';
}

int LoLd_translator( char tchar )
{
	// sum 1:
	// A 1 cholesterol
	// B 1 dopc
	
	// sum 2:
	// C 2 cholesterol
	// D 1 cholesterol 1 DOPC
	// E 2 DOPC

	// sum 3:
	// F 3 cholesterol
	// G ... etc.
	//

	int val = 0;
	
	if( tchar >= 'A' && tchar <= 'Z' )
		return  tchar - 'A';
	else 
		return 26 + tchar - 'a';

}

void HMMDirective::translate( void )
{
	for( int s = 0; s < n_data; s++ )
	{

		if( translated_observable_string[s] ) free( translated_observable_string[s] );
		translated_observable_string[s] = (int *)malloc( sizeof(int) * strlen( raw_observable_string[s] ) );
		translated_len[s] = strlen(raw_observable_string[s]);

		for( int t = 0; t < translated_len[s]; t++ )
		{
			translated_observable_string[s][t] = lipid_translator( raw_observable_string[s][t] ); 
		}	

	}
}


HMMCommands::HMMCommands(FILE *theFile )
{
	ndirectives = 0;
	char *buffer = (char *)malloc( sizeof(char) * 500000 );
	while( !feof(theFile) )
	{
		getLine( theFile, buffer );

		if( feof(theFile) ) break;

		char *t = buffer;

		while( *t == ' ' || *t == '\t' )  t += 1;
		
		if( !strncasecmp( t, "COMMENT", strlen("COMMENT") ) )
			continue;

		char theCommand[256];
		char theInterpreter[256];
		char theNetwork[256];

		sprintf( theInterpreter, "NUM6");
		sprintf( theNetwork, "LoLd");

		int nr = sscanf(buffer, "COMMAND %s", theCommand );//, theInterpreter, theNetwork );
		if( nr != 1 )
		{
	//		printf("Could not read three fields from line '%s'.\n", buffer );
			printf("Could not read command from line '%s'.\n", buffer );
			exit(1); 
		}

		int c_id = commandInterpreter( theCommand );
		int i_id = commandInterpreter( theInterpreter );
		int n_id = commandInterpreter( theNetwork );

		if( !(c_id & COMMAND_BIT) )
		{
			printf("Unknown command '%s'.\n", theCommand );
			exit(1);
		}
		
		if( !(i_id & INTERPRETER_BIT) )
		{
			printf("Unknown interpreter '%s'.\n", theInterpreter );
			exit(1);
		}
		
		if( !(n_id & NETWORK_BIT) )
		{
			printf("Unknown network '%s'.\n", theNetwork );
			exit(1);
		}
		int fp = ftell(theFile);
		getLine(theFile,buffer);


		char *theString;
		char *theClass;

		allDirectives[ndirectives] = new HMMDirective();
		allDirectives[ndirectives]->useTranslator(i_id);
		allDirectives[ndirectives]->useNetwork(n_id);
		allDirectives[ndirectives]->doCommand(c_id);

		int added = 0;
		while( !feof(theFile) && strncasecmp( buffer, "STOP", 4 ) && strncasecmp( buffer, "COMMAND", 7 ) )
		{

			if( feof(theFile) ) break;
			if( c_id == COMMAND_TRAIN || c_id == COMMAND_PSEUDO || c_id == COMMAND_START || c_id == COMMAND_STOP || c_id == COMMAND_DECODE || c_id == COMMAND_GRAD ) 
			{
				theString = (char *)malloc( sizeof(char) * (strlen( buffer )+1) );
				theClass = (char *)malloc( sizeof(char) * (strlen( buffer )+1) );

				int nr = sscanf( buffer, "%s %s #", theString, theClass  );
			
				char *find_comment = buffer;

				while( *find_comment && *find_comment != '#' )
					find_comment += 1;

				char *use_comment = NULL;
				if( *find_comment == '#' )
					use_comment = find_comment + 1;
	
				if( theClass[0] == '#' || nr == 1 )
				{
					for( int x = 0; x < strlen(buffer)+1; x++ )
					{
						theClass[x] = '.';
#ifdef PROTEIN_STATE
						if( theString[x] == 'c' )
							theClass[x] = 'b';
#endif
					}
				}
				allDirectives[ndirectives]->addData( theString, theClass, NULL, 1.0, NULL, use_comment );

	
				free(theString);
				free(theClass);	
			}
			else if( c_id == COMMAND_SAVE || c_id == COMMAND_LOAD )
				allDirectives[ndirectives]->addData( buffer, NULL );
			fp = ftell(theFile);	
			getLine(theFile, buffer );
		}
				
		if( c_id == COMMAND_TRAIN || c_id == COMMAND_DECODE || c_id == COMMAND_GRAD )
			allDirectives[ndirectives]->translate();
	
		ndirectives++;
		
		if( !strncasecmp(buffer, "STOP", 4 ) )
			break;

		if( !strncasecmp( buffer, "COMMAND", 7 ) )
			fseek( theFile, fp, SEEK_SET );

	}
	free(buffer);
}





