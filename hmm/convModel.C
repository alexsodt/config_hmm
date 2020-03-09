#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HMM.h"
#include "util.h"

#define MAX_FILE_TRANSITIONS    100
#define MAX_FILE_OBSERVABLES    100

struct HMMSystemFileHeader
{
        int nstates;
        int interpreter;
        int model_id;
};

struct HMMStateFileEntry
{
        int file_id;
        char label[256];
        double transmission[MAX_FILE_TRANSITIONS];
        double emission[MAX_FILE_OBSERVABLES];
};

int main( int argc, char **argv )
{
	if( argc != 2 )
	{
		printf("Syntax: convModel model.save\n"); 
	}

	HMMSystemFileHeader theHeader;	

	FILE *theFile = fopen(argv[1],"r");

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	fread( buffer, sizeof(char), 6, theFile );
	rewind(theFile);
	int binary_to_ascii = 1;

	if( !strncasecmp( buffer, "ASCII", 5 ) )
		binary_to_ascii = 0;
	
	if( binary_to_ascii ) 
	{
		fread( &theHeader, sizeof(HMMSystemFileHeader), 1, theFile );
	
		int nstates = theHeader.nstates;
	
		printf("ASCII\n");
		printf("%d\n", theHeader.nstates );
		printf("%d\n", theHeader.interpreter );
		printf("%d\n", theHeader.model_id );

	        for( int s = 0; s < nstates; s++ )
	        {
	                HMMStateFileEntry theEntry;

			fread( &theEntry, sizeof(HMMStateFileEntry), 1, theFile );
	
			printf("%d\n", theEntry.file_id );
			printf("%s\n", theEntry.label );		
			for( int i = 0; i < nstates; i++ )
				printf("%lf ", theEntry.transmission[i] );
			printf("\n");
			for( int i = 0; i < MAX_FILE_OBSERVABLES; i++ )
			{
				double val = theEntry.emission[i];
				if( !(val >0 || val < 1 ) )
					val = 0;
				if( fabs(val) < 1e-250 )
					val = 0;
				printf("%lf ", val );
			}
			printf("\n");
	        }	
	}
	else
	{
		FILE *outFile = fopen("bin.model","wb");

		getLine( theFile, buffer );
		// ASCII
		getLine( theFile, buffer );
		sscanf(buffer, "%d", &(theHeader.nstates) );
		getLine( theFile, buffer );
		sscanf(buffer, "%d", &(theHeader.interpreter) );
		getLine( theFile, buffer );
		sscanf(buffer, "%d", &(theHeader.model_id) );
		
		fwrite( &theHeader, sizeof(HMMSystemFileHeader), 1, outFile );

		int nstates = theHeader.nstates;

	        for( int s = 0; s < nstates; s++ )
	        {
	                HMMStateFileEntry theEntry;
	
			getLine( theFile, buffer );
			sscanf(buffer, "%d", &(theEntry.file_id) );
			getLine( theFile, buffer );
			buffer[255] = '\0';
			strcpy( theEntry.label, buffer );

			getLine(theFile, buffer );
			int nv = readNDoubles( buffer, theEntry.transmission, nstates );
			for( int x = nv; x < MAX_FILE_TRANSITIONS; x++ )
				theEntry.transmission[x] = 0;
			getLine(theFile, buffer );
			nv = readNDoubles( buffer, theEntry.emission, MAX_FILE_OBSERVABLES );
			for( int x = nv; x < MAX_FILE_TRANSITIONS; x++ )
				theEntry.emission[x] = 0;

			fwrite( &theEntry, sizeof(HMMStateFileEntry), 1, outFile );
	        }	
	}
}
