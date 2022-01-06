// by alex sodt
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int readNInts( char *buffer, int *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%d", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
                        s++;
        }

	return nread;
}


int readNDoubles( char *buffer, double *vals, int nvalues )
{
        int s = 0;

        int slen = strlen(buffer);
	int nread = 0;

        for( int x = 0; x < nvalues; x++ )
        {
                while( s < slen && (buffer[s] == ' ' || buffer[s] == '\t') )// || isalpha(buffer[s])) )
                        s++;

                int nr = sscanf( buffer + s, "%lf", vals+x );

		if( nr != 1 )
			return nread;

		nread++;

                while( s < slen && !(buffer[s] == ' ' || buffer[s] == '\t'))// || isalpha(buffer[s])) )
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

                if( tc != '\n' && i < 999999 )
                {
                        theBuffer[i++] = tc;
                }
                else if( tc != '\n' && i >= 999999 )
		{
		}
		else
                        break;
        }

        theBuffer[i] = '\0';
}

void print5( int val, char *str )
{
        if( val < 10 )
                sprintf(str, "0000%d", val );
        else if( val < 100 )
                sprintf(str, "000%d", val );
        else if( val < 1000 )
                sprintf(str, "00%d", val );
        else if( val < 10000 )
                sprintf(str, "0%d", val );
        else
                sprintf(str, "%d", val );
}

const char *advance_string( const char *t, int nadv )
{
	while( *t && isspace(*t) ) t += 1;
	for( int adv = 0; adv < nadv; adv++ )
	{
		while( *t && !isspace(*t) ) t += 1;
		while( *t && isspace(*t) ) t += 1;
	}

	return t;
}

int decodeString( char *buf, char **out, char **res_out, int nmax )
{
	int done = 0;

	int ns = 0;
	int t = 0;
	int lt = 0;
	int lim = 256;
	char cur[lim+1];
	char res_cur[lim+1];

	while( *(buf+t) )
	{
		while( *(buf+t) == '+' ) t++;

		int got_res = 0;

		while( *(buf+t) && *(buf+t) != '+' && *(buf+t) != '_' )
		{
			if( lt < lim )
			{
				cur[lt] = *(buf+t);	
				lt++;
			}

			t++;
		}
		
		cur[lt] = '\0';

		int res_lt = lt;

		if( *(buf+t) == '_' )
 		{
			t += 1;
			got_res = 1;

			strcpy( res_cur, cur );

			lt = 0;

			while( *(buf+t) && *(buf+t) != '+' )
			{
				if( lt < lim )
				{
					cur[lt] = *(buf+t);	
					lt++;
				}
	
				t++;
			}
		
			cur[lt] = '\0';
		}

 
		if( lt > 0 && ns < nmax )
		{
			out[ns] = (char *)malloc( sizeof(char) * (1 + lt ) );
			strcpy( out[ns], cur );

			if( got_res ) 
			{						
				res_out[ns] = (char *)malloc( sizeof(char) * (1 + res_lt ) );
				strcpy( res_out[ns], res_cur );
			}
			else
			{
				res_out[ns] = (char *)malloc( sizeof(char) * (1 + strlen("any") ) );
				strcpy( res_out[ns], "any" );
			}
			ns++;
		}
			
		lt = 0;

	}

	return ns;
}
