
#include "math.h"
#include "geometry.h"
#include "util.h"

double len3( double *dr )
{
	return sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
}

double dihe( double *r1, double *r2, double *r3, double *r4, double L )
{
	double b1[3] = {r2[0] - r1[0], r2[1] -r1[1], r2[2] - r1[2] };
	double b2[3] = {r3[0] - r2[0], r3[1] - r2[1], r3[2] - r2[2] };
	double b3[3] = {r4[0] - r3[0], r4[1] - r3[1], r4[2] - r3[2] };

	double lb2 = len3(b2);
	double b2xb3[3] = {  b2[1] * b3[2] - b2[2] * b3[1],
			   -(b2[0] * b3[2] - b2[2] * b3[0]),
		             b2[0] * b3[1] - b2[1] * b3[0] };
	double b1xb2[3] = {  b1[1] * b2[2] - b1[2] * b2[1],
			   -(b1[0] * b2[2] - b1[2] * b2[0]),
		             b1[0] * b2[1] - b1[1] * b2[0] };

	double b1db2xb3 = b1[0] * b2xb3[0] + b1[1] * b2xb3[1] + b1[2] * b2xb3[2];
	double xx = b1xb2[0] * b2xb3[0] + b1xb2[1] * b2xb3[1] + b1xb2[2] * b2xb3[2];

	double arg1= lb2 * b1db2xb3;
	double arg2 = xx;

	if( fabs( arg1 ) < 1e-10 )
		arg1 = 1e-10;
	if( fabs( arg2 ) < 1e-10 )
		arg2 = 1e-10;	
 
//	double phi = atan2( lb2 * b1db2xb3, xx );
	double phi = atan2( arg1, arg2 );


	return phi;

/*			
	double dr12[3] = { r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2] };

	while( dr12[0] > L/2 ) dr12[0] -= L;
	while( dr12[0] < -L/2 ) dr12[0] += L;
	while( dr12[1] > L/2 ) dr12[1] -= L;
	while( dr12[1] < -L/2 ) dr12[1] += L;
	
	double l12 = sqrt ( dr12[0]*dr12[0] + dr12[1] * dr12[1] + dr12[2] * dr12[2] );
	
	double dr32[3] = { r3[0] - r2[0], r3[1] - r2[1], r3[2] - r2[2] };

	while( dr32[0] > L/2 ) dr32[0] -= L;
	while( dr32[0] < -L/2 ) dr32[0] += L;
	while( dr32[1] > L/2 ) dr32[1] -= L;
	while( dr32[1] < -L/2 ) dr32[1] += L;
	
	double l32 = sqrt ( dr32[0]*dr32[0] + dr32[1] * dr32[1] + dr32[2] * dr32[2] );
	
	double dr43[3] = { r4[0] - r3[0], r4[1] - r3[1], r4[2] - r3[2] };

	while( dr43[0] > L/2 ) dr43[0] -= L;
	while( dr43[0] < -L/2 ) dr43[0] += L;
	while( dr43[1] > L/2 ) dr43[1] -= L;
	while( dr43[1] < -L/2 ) dr43[1] += L;

	double l43 = sqrt ( dr43[0]*dr43[0] + dr43[1] * dr43[1] + dr43[2] * dr43[2] );

	double cp1[3];

	cp1[0] =   dr12[1] * dr32[2] - dr12[2] * dr32[1];
	cp1[1] = -(dr12[0] * dr32[2] - dr12[2] * dr32[0]);
	cp1[2] =   dr12[0] * dr32[1] - dr12[1] * dr32[0];

	cp1[0] /= (l12*l32);
	cp1[1] /= (l12*l32);
	cp1[2] /= (l12*l32);


	double signcp[3];

	signcp[0] =   dr12[1] * dr32[2] - dr12[2] * dr32[1];
	signcp[1] = -(dr12[0] * dr32[2] - dr12[2] * dr32[1]);
	signcp[2] =   dr12[0] * dr32[1] - dr12[1] * dr32[0];
	
	double cp2[3];

	cp2[0] =   dr32[1] * dr43[2] - dr32[2] * dr43[1];
	cp2[1] = -(dr32[0] * dr43[2] - dr32[2] * dr43[0]);
	cp2[2] =   dr32[0] * dr43[1] - dr32[1] * dr43[0];

	cp2[0] /= -(l12*l32); // (dr23 == -dr32).
	cp2[1] /= -(l12*l32); // (dr23 == -dr32).
	cp2[2] /= -(l12*l32); // (dr23 == -dr32).

	double lcp1 = sqrt( cp1[0]*cp1[0] + cp1[1]*cp1[1] + cp1[2]*cp1[2]);
	double lcp2 = sqrt( cp2[0]*cp2[0] + cp2[1]*cp2[1] + cp2[2]*cp2[2]);

 	double crossz = dr43[0] * signcp[0] + dr43[1] * signcp[1] + dr43[2] * signcp[2];
        double sign = 1;
  
        if( crossz < 0)
		sign = -1;

	double dp = (cp1[0]*cp2[0] + cp1[1] *cp2[1] + cp1[2] * cp2[2])/(lcp1*lcp2);

	if( dp > 1 - 1e-14 )
		dp = 1 - 1e-14;
	if ( dp < -1 + 1e-14 )
		dp = -1 + 1e-14;

	return sign * acos( dp );
*/
}

double bend( double *r1, double *r2, double *r3, double L )
{

	double dr12[3] = { r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2] };
	double dr32[3] = { r3[0] - r2[0], r3[1] - r2[1], r3[2] - r2[2] };

	while( dr12[0] > L/2 ) dr12[0] -= L;
	while( dr12[0] < -L/2 ) dr12[0] += L;
	while( dr12[1] > L/2 ) dr12[1] -= L;
	while( dr12[1] < -L/2 ) dr12[1] += L;
	
	while( dr32[0] > L/2 ) dr32[0] -= L;
	while( dr32[0] < -L/2 ) dr32[0] += L;
	while( dr32[1] > L/2 ) dr32[1] -= L;
	while( dr32[1] < -L/2 ) dr32[1] += L;

	double l12 = sqrt ( dr12[0]*dr12[0] + dr12[1] * dr12[1] + dr12[2] * dr12[2] );
	double l32 = sqrt ( dr32[0]*dr32[0] + dr32[1] * dr32[1] + dr32[2] * dr32[2] );

	double dp = dr12[0] * dr32[0] + dr12[1] * dr32[1] + dr12[2] * dr32[2];
	
	double phi = acos( dp / (l12*l32) );

	return phi;
}
