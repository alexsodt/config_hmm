
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HMM.h"
#include "util.h"

extern int globalN_OBSERVABLES;

//#define PROTEIN_STATE // may be defined in HMM.h
//#define THREE_STATE
//#define SINGLE_STATE
//#define FOUR_STATE

int getNObservables( const char *fileName )
{
	FILE *theFile = fopen(fileName, "r");
	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", fileName );
		exit(1);
	}		 
	char buffer[4096];
	
	getLine( theFile, buffer );
	
	int nstates;

	sscanf( buffer, "%d", &nstates );

	fclose(theFile);

	return nstates;
}

int MakeGeneralNetwork( HMMSystem &theSystem, int ignore, const char *fileName )
{
	FILE *theFile = fopen(fileName, "r");
	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", fileName );
		exit(1);
	}		 
	char buffer[4096];
	
	getLine( theFile, buffer );
	
	int nstates;
	sscanf( buffer, "%d", &nstates );


	// we add in an "outer" state A for when a lipid can't be put on the surface, for whatever reason. 
//	nstates += 1;

	char label[256];
	sprintf(label, "A");

	HMMState *theStates[nstates];
	HMMState *borderStates[nstates*nstates];	
	char *state_codes = (char *)malloc( sizeof(char) * nstates );

	for( int s1= 0; s1 < nstates; s1++ )
	for( int s2 = 0; s2 < nstates; s2++ )
	{
		borderStates[s1*nstates+s2] = NULL;
	}

	int max_borders = 0;

	int border_state_index = nstates;
	for( int pass = 0; pass < 2; pass++ )
	{
		for( int s = 0; s < nstates; s++ )
		{
			getLine( theFile, buffer );
			char ind;
			int neighbors;
	
			sscanf(buffer, "%c %d", &ind, &neighbors );
			state_codes[s] = ind;
	
			char *border_string = buffer;
			// advance first two fields.
			while( *border_string && (*border_string == ' ' || *border_string == '\t' )) border_string++; 
			while( *border_string && (*border_string != ' ' && *border_string != '\t' )) border_string++; 
			while( *border_string && (*border_string == ' ' || *border_string == '\t' )) border_string++; 
			while( *border_string && (*border_string != ' ' && *border_string != '\t' )) border_string++; 
			while( *border_string && (*border_string == ' ' || *border_string == '\t' )) border_string++; 
	
			if( neighbors > max_borders )
				max_borders = neighbors;
	
			if( pass == 0 )
			{
				sprintf(label, "%c", ind );
				theStates[s] = new HMMState('S', label, 1, s, 0 );	
			}
		
			if( pass == 1 )
			{
				for( int b = 0; b < neighbors; b++ )
				{
					char neighbor = *border_string;

					while( *border_string && (*border_string != ' ' && *border_string != '\t' ) ) border_string++; 
					while( *border_string && (*border_string == ' ' || *border_string == '\t' ) ) border_string++; 

					int s2=-1;

					for( int sx = 0; sx< nstates; sx++ )
					{
						if( theStates[sx]->label[0] == neighbor ) 
							s2 = sx;
					}

					if( s2 == -1 )
					{
						printf("Couldn't find class '%c' in the list.\n", neighbor );
						exit(1);
					}

					if( s == s2 ) continue;

					if( borderStates[s*nstates+s2] == NULL )
					{
						if( neighbor > ind )
							sprintf(label, "Border_%c_%c", ind, neighbor );
						else
							sprintf(label, "Border_%c_%c", neighbor, ind  );

						printf("Making border state %s.\n", label );
						borderStates[s*nstates+s2] = new HMMState('B', label, 1, border_state_index++, 0 );	
						borderStates[s2*nstates+s] = borderStates[s*nstates+s2];	
					}
				}
			
/*	
				sprintf(label, "Border_A_%c", ind );

				borderStates[s*nstates+0] = new HMMState('B', label, 1, border_state_index++, 0 );	
				borderStates[0*nstates+s] = borderStates[s*nstates+0];	
*/
			}
		}	

		if( pass == 1 )
		{
			// add in all transitions.

			for( int s1 = 0; s1 < nstates; s1++ )
			{
				theStates[s1]->addTransition(theStates[s1]);

				for( int s2 = s1; s2 < nstates; s2++ )
				{
					if( borderStates[s1*nstates+s2] )
					{
						theStates[s1]->addTransition(borderStates[s1*nstates+s2]);
						borderStates[s1*nstates+s2]->addTransition(theStates[s1]);
						
						theStates[s2]->addTransition(borderStates[s1*nstates+s2]);
						borderStates[s1*nstates+s2]->addTransition(theStates[s2]);			
					}
				}
			}

		}

		rewind(theFile);
		// number of states
		getLine(theFile, buffer);	
	}	

	fclose(theFile);

	// add transitions between shared states.

	for( int s1 = 0; s1 < nstates; s1++ )
	{	
		char code1 = state_codes[s1];

		for( int s2 = s1+1; s2 < nstates; s2++ )
		{
			char code2 = state_codes[s2];

			if( borderStates[s1*nstates+s2] )
			{
				for( int s3 = 0; s3 < nstates; s3++ )
				{
					if( s3 == s2 || s3 == s1 ) continue;
					if( s3 < s2 || s3 < s1 ) continue;

					if( borderStates[s1*nstates+s3] && borderStates[s2*nstates+s3] )
					{
						borderStates[s1*nstates+s2]->addTransition( borderStates[s1*nstates+s3] );
						borderStates[s1*nstates+s3]->addTransition( borderStates[s1*nstates+s2] );

						borderStates[s1*nstates+s2]->addTransition( borderStates[s2*nstates+s3] );
						borderStates[s2*nstates+s3]->addTransition( borderStates[s1*nstates+s2] );
					}	
				}
			}
		}
	}

	for( int s1 = 0; s1 < nstates; s1++ )
	for( int s2 = s1+1; s2 < nstates; s2++ )
	{
		if( borderStates[s1*nstates+s2] )
			borderStates[s1*nstates+s2]->addTransition( borderStates[s1*nstates+s2] );
	}
	
	for( int s1 = 0; s1 <nstates; s1++ )
		theSystem.addState(theStates[s1]);

	for( int s1 = 0; s1 <nstates; s1++ )
	for( int s2 = s1; s2 <nstates; s2++ )
	{
		if( borderStates[s1*nstates+s2] )
			theSystem.addState(borderStates[s1*nstates+s2]);
	}
	


	theSystem.setup();	

	printf("Generated %d observables.\n", nstates);

	return nstates;
}

