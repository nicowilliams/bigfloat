/*  Basic complex arithmetic.  Everything assumed to be in Cartesian form
	and made to stay that way.
*/

#include <stdio.h>
#include "bigfloat.h"


/*  set a float to a constant 1  */

void one( FLOAT *x)
{
	null( x);
	x->expnt = 1;
	x->mntsa.e[MS_MNTSA] = 0x40000000;
}

/*  zero out a complex storage location  */

void null_cmplx( COMPLEX *x)
{
	null( & x->real);
	null( & x->imag);
}

/*  add two complex numbers.
	c = a + b 
*/
void add_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX mya, myb;

	copy_cmplx( a, &mya);
	copy_cmplx( b, &myb);
	
	add( &mya.real, &myb.real, &c->real);
	add( &mya.imag, &myb.imag, &c->imag);
}

/*  subtract two complex numbers
	c = a - b
*/
void subtract_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX myb;

	copy_cmplx( b, &myb);
	negate( &myb.real);
	negate( &myb.imag);
	add_cmplx( a, &myb, c);
}

/*  multiply two complex numbers.

	c = (a.real * b.real - a.imag * b.imag) + i(a.imag * b.real + a.real * b.imag)
*/
void multiply_cmplx(  COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX mya, myb;
	FLOAT	temp1, temp2;

	copy_cmplx( a, &mya);
	copy_cmplx( b, &myb);
	multiply( &mya.real, &myb.real, &temp1);
	multiply( &mya.imag, &myb.imag, &temp2);
	subtract( &temp1, &temp2, &c->real);
	multiply( &mya.real, &myb.imag, &temp1);
	multiply( &mya.imag, &myb.real, &temp2);
	add( &temp1, &temp2, &c->imag);
}

/*  divide two complex numbers.
	To keep things in Cartesian form multiply top by
	conjugate of bottom.  Scale result by magnitude of
	bottom.
	
		output is c = ( a * b^*) / |b|
	
	returns 1 if b != 0, 0 if |b| = 0
*/
int divide_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	FLOAT	mag1, mag2;
	COMPLEX	myb;
	
	copy_cmplx( b, &myb);
	multiply( &myb.real, &myb.real, &mag1);
	multiply( &myb.imag, &myb.imag, &mag2);
	add( &mag1, &mag2, &mag1);
	negate( &myb.imag);
	multiply_cmplx( a, &myb, c);
	if( ! divide( &c->real, &mag1, &c->real)) return 0;
	divide( &c->imag, &mag1, &c->imag);
	return 1;
}

/*  compute y = x^k
	where k is a signed integer in range +/-2^31 and x, y are complex.
	returns 1 if ok, 0 if x = 0 and k < 0
*/

int intpwr_cmplx( COMPLEX *x, int k, COMPLEX  *y)
{
	int signflag, n;
	COMPLEX z, t;
	
/*	FLOAT seven, temp;
	
	null(&seven);
	seven.expnt = 3;
	seven.mntsa.e[MS_MNTSA] = 0x70000000;
	square_root( &seven, &seven);
	
/*  initialize Knuth's algorithm A pg 442 semi-numerical algorithms  */

	copy_cmplx( x, &z);
	if ( k < 0 )
	{
		signflag = 1;
		n = -k;
	}
	else
	{
		signflag = 0;
		n = k;
	}
	null_cmplx( &t);
	one( &t.real);
	while (n)
	{
		if ( n & 1 ) multiply_cmplx( &t, &z, &t);
		multiply_cmplx( &z, &z, &z);
		n >>= 1;
	}
	if ( signflag)
	{
		null_cmplx( &z);
		one( &z.real);
		return divide_cmplx( &z, &t, y);
	}
	copy_cmplx( &t, y);
	return 1;
}

/*  compute magnitude of a complex number.  Returns
	FLOAT result.
*/

void magnitude_cmplx( COMPLEX *x, FLOAT *m)
{
	FLOAT x2, y2;
	
	multiply( &x->real, &x->real, &x2);
	multiply( &x->imag, &x->imag, &y2);
	add( &x2, &y2, m);
	square_root( m, m);
}

/*  compute exp(z) for z complex.
	z = x + iy so
	e = exp(z) = exp(x)*(cos(y) + isin(y))
	if x too large, returns 0, otherwise returns e and 1.
	works in place.
*/

int exp_cmplx( COMPLEX *z, COMPLEX *e)
{
	FLOAT	x, y, xp, cy, sy;
	
	copy( &z->real, &x);
	copy( &z->imag, &y);
	if( !exp( &x, &xp) )
	{
		copy( &xp, &e->real);
		null( &e->imag);
		return 0;
	}
	cosine( &y, &cy);
//	printfloat("cos(y)=", &cy);
	sine( &y, &sy);
//	printfloat("sin(y)=", &sy);
	multiply( &xp, &cy, &e->real);
//	printfloat("exp(x)=", &xp);
	multiply( &xp, &sy, &e->imag);
	return 1;
}
