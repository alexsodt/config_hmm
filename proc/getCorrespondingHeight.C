#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

int main( int argc, char **argv )
{
	if( argc < 3 )
	{
		printf("Syntax: getCorrespondingHeight.exe file.heights ref_height\n");
		return -1;
	}

	FILE *theFile = fopen(argv[1],"r");

	if( !theFile )
	{
		printf("Couldn't open file '%s'.\n", argv[1] );
		return -1;
	}

	double ref_height = atof(argv[2]);
	char *buffer = (char *)malloc( sizeof(char) * 1000000 );

	char *best_atom_name = NULL;
	double best_chi2 = 1e10;
	double best_h = 0;
	while( !feof(theFile) )
	{
		getLine(theFile, buffer );

		char *t = (char *)malloc( sizeof(char) * (1+strlen(buffer)) );

		double h;
		int nr = sscanf(buffer, "%s %lf\n", t, &h );	

		if( !best_atom_name || (h-ref_height)*(h-ref_height) < best_chi2 )
		{
			if( best_atom_name ) free(best_atom_name);

			best_atom_name = (char *)malloc( sizeof(char) * (1 + strlen(t) ) );
			strcpy( best_atom_name, t );
			best_chi2 = (h-ref_height)*(h-ref_height);
			best_h = h;
		}
	}
	if( best_atom_name )
		printf("%s %lf\n", best_atom_name, best_h );
}
