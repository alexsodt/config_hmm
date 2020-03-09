#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int readNDoubles( char *buffer, double *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%lf", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' '))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}



void getLine( FILE *theFile, char *theBuffer )
{
        int i = 0;

        while( !feof(theFile) )
        {
                char tc = fgetc(theFile);
                if( tc != '\n' )//&& i < 19999)
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n')// && i >= 19999 )
		{
		}
		else
                        break;
        }

        theBuffer[i] = '\0';
}

