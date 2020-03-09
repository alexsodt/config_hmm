#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Fuzzy.h"

#define BUFFA 1


// label surfaces A and B.
void makeSurfaceMask( char *string, double *train_z, char *ABmask )
{
	int len = strlen(string);

	for( int x = 0; x < len; x++ )
		ABmask[x] = '.';

	// positive Z is the 'A' surface.

	char curSurface = 'A';
	char prevTopo   = 'X';

	for( int x = 1; x < len; x++ )
	{
		if( string[x] == '.' || (string[x]=='M'&&string[x-1]!='M') )
		{
			if( string[x-1] != '.' ) // switch?
			{
				if( prevTopo == 'M' ) // this is the case of MMMM....XXXX.... , we need a new surface.
				{
					
					if( curSurface == 'A' ) curSurface = 'B';
					else if( curSurface == 'B' ) curSurface = 'A';
				}					
			}

			ABmask[x] = curSurface;
		}
		else if( (string[x] !='M' && string[x-1]=='M') )
		{
			if( string[x-1] != '.' ) // switch?
			{
				if( prevTopo == 'M' ) // this is the case of MMMM....XXXX.... , we need a new surface.
				{
					if( curSurface == 'A' ) curSurface = 'B';
					else if( curSurface == 'B' ) curSurface = 'A';
				}					
			}

			ABmask[x-1] = curSurface;
		
		}
		else prevTopo = string[x];
	}

	for( int x = 0; x < len-1; x++ )
	{
		if( string[x] == '.' && string[x+1] == 'M' )
			ABmask[x+1] = ABmask[x];
	}
	
	for( int x = len-1; x > 0; x-- )
	{
		if( string[x] == '.' && string[x-1] == 'M' )
			ABmask[x-1] = ABmask[x];
	}
/*
	for( int x = 1; x < len; x++ )
	{
		if( string[x] == '.' )
		{
			// possible surface site.
			if( ABmask[x-1] == 'A' )
				ABmask[x] = 'A';
			else if( ABmask[x-1] == '.' )
			{
				if( string[x-1] == 'M' )
					ABmask[x] = 'A';
			}
		}
	}
	for( int x = len-2; x >= 0; x-- )
	{
		if( string[x] == '.' )
		{
			// possible surface site.
			if( ABmask[x+1] == 'B' )
				ABmask[x] = 'B';
			else if( ABmask[x+1] == '.' )
			{
				if( string[x+1] == 'M' )
					ABmask[x] = 'B';
			}
		}
	}

	for( int x = 0; x < len-1; x++ )
	{
		if( ABmask[x] == '.' && ABmask[x+1] == 'A' )
			ABmask[x] = 'A';
	}
	for( int x = len; x >= 1; x-- )
	{
		if( ABmask[x-1] == 'B' && ABmask[x] == '.' )
			ABmask[x] = 'B';
	}
*/
	// go through and correct.
	

	for( int x = 0; x < len; x++ )
	{
		if( ABmask[x] == 'A' || ABmask[x] == 'B' )
		{		
			char check = ABmask[x];

			int tstart = x;

			while( ABmask[x] == check ) x++;

			int tstop = x;

			int mid = (tstop+tstart)/2;
			
			double z = train_z[mid];

			if( z < -8.0 )
			{
				for( int p = tstart; p < tstop; p++ )
					ABmask[p] = 'B';
			}
			else if( z > 8.0 )
			{
				for( int p = tstart; p < tstop; p++ )
					ABmask[p] = 'A';
			}
			else
			{
				for( int p = tstart; p < tstop; p++ )
					ABmask[p] = '.';
			}
		}
	}



//	printf("string: %s\n", string );
//	printf("mask:   %s\n", ABmask );
}


void makeFuzzy( char *string )
{
	int pt = 0;
	int len = strlen(string);
	int done = 0;

	char cur_code = '?';

	char *str_copy = (char *)malloc( sizeof(char) * (len+1) );

	strcpy( str_copy, string );

	while( !done )
	{
		while( string[pt] == cur_code || string[pt] == '.' )
			pt++;


		if( string[pt] == '\0' ) break;
		// new code.

		int cur_pt = pt;

		while( string[cur_pt] == string[pt] ) cur_pt++;

		int region_length = (cur_pt - pt);
		
		cur_code = string[pt];

		if( cur_code == 'V' || cur_code == 'R' || cur_code == 'E' || cur_code == 'F' ) continue;

		if( (string[pt] == 'M'|| string[pt] == 'N') && region_length > 16 )
		{
			int trim = (region_length-16)/2;
			for( int x = 0; x < trim; x++ )
			{
//				if( pt+x > 3 && pt+x < len-3 )
				if( pt+x >= 0 && pt+x < len )
				{
					string[pt+x] = '.';
					string[cur_pt-x-1] = '.';

				}
			}
		}
		else
		{
		switch( region_length )
		{
			case 1: // do nothing
			case 2: 
				break;
			case 3:
//				if( pt > 3 && pt < len-3 )
				if( pt >= 0 && pt < len )
					string[pt] = '.';
//				if( pt+2 > 3 && pt+2 < len-3 )
				if( pt+2 >= 0 && pt+2 < len )
					string[pt+2] = '.';
				break;
			case 4:
//				if( pt > 3 && pt < len-3 )
				if( pt+3 >= 0 && pt+3 < len )
				string[pt] = '.';
//				if( pt+3 > 3 && pt+3 < len-3 )
				if( pt+3 >= 0 && pt+3 < len )
				string[pt+3] = '.';
				break;
			case 5:
			case 6:
				if( pt >= 0 && pt < len )
				string[pt] = '.';
				if( pt+1 >= 0 && pt+1 < len )
				string[pt+1] = '.';

				if( cur_pt-1 >= 0 && cur_pt-1 < len )
				string[cur_pt-1] = '.';
				if( cur_pt-2 >= 0 && cur_pt-2 < len )
				string[cur_pt-2] = '.';
				break;
			case 7:
			case 8:
				if( pt >= 0 && pt < len )
				string[pt] = '.';
				if( pt+1 >= 0 && pt+1 < len )
				string[pt+1] = '.';
				if( pt+2 >= 0 && pt+2 < len )
				string[pt+2] = '.';

				if( cur_pt-1 >= 0 && cur_pt-1 < len )
				string[cur_pt-1] = '.';
				if( cur_pt-2 >= 0 && cur_pt-2 < len )
				string[cur_pt-2] = '.';
				if( cur_pt-3 >= 0 && cur_pt-3 < len )
				string[cur_pt-3] = '.';
				break;
			default:
				if( pt >= 0 && pt < len )
					string[pt] = '.';
				if( pt+1 >= 0 && pt+1 < len )
					string[pt+1] = '.';
				if( pt+2 >= 0 && pt+2 < len )
					string[pt+2] = '.';
				if( pt+3 >= 0 && pt+3 < len )
					string[pt+3] = '.';

				if( cur_pt-4 >= 0 && cur_pt-4 < len )
					string[cur_pt-4] = '.';
				if( cur_pt-3 >= 0 && cur_pt-3 < len )
					string[cur_pt-3] = '.';
				if( cur_pt-2 >= 0 && cur_pt-2 < len )
					string[cur_pt-2] = '.';
				if( cur_pt-1 >= 0 && cur_pt-1 < len )
					string[cur_pt-1] = '.';
				break;
		}
		}
		
		if( string[pt] == '\0' ) break;
	}	

	int tb = 0;

	while( string[tb] == '.' )
	{	string[tb] = str_copy[tb]; tb++; }
	tb = len-1;
	while( string[tb] == '.' )
	{	string[tb] = str_copy[tb]; tb--; }

	free(str_copy);

}

void Randomize( char *string, char *search )
{

}






