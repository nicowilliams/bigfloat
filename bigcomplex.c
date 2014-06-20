/*  Basic complex arithmetic.  Everything assumed to be in Cartesian form
	and made to stay that way.
*/

#include <stdio.h>
#include "bigfloat.h"


/*  set a float to a constant 1  */

void bf_one( FLOAT *x)
{
	bf_null( x);
	x->expnt = 1;
	x->mntsa.e[MS_MNTSA] = 0x40000000;
}

/*  zero out a complex storage location  */

void bf_null_cmplx( COMPLEX *x)
{
	bf_null( & x->real);
	bf_null( & x->imag);
}

/*  add two complex numbers.
	c = a + b 
*/
void bf_add_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX mya, myb;

	bf_copy_cmplx( a, &mya);
	bf_copy_cmplx( b, &myb);
	
	bf_add( &mya.real, &myb.real, &c->real);
	bf_add( &mya.imag, &myb.imag, &c->imag);
}

/*  subtract two complex numbers
	c = a - b
*/
void bf_subtract_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX myb;

	bf_copy_cmplx( b, &myb);
	bf_negate( &myb.real);
	bf_negate( &myb.imag);
	bf_add_cmplx( a, &myb, c);
}

/*  multiply two complex numbers.

	c = (a.real * b.real - a.imag * b.imag) + i(a.imag * b.real + a.real * b.imag)
*/
void bf_multiply_cmplx(  COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	COMPLEX mya, myb;
	FLOAT	temp1, temp2;

	bf_copy_cmplx( a, &mya);
	bf_copy_cmplx( b, &myb);
	bf_multiply( &mya.real, &myb.real, &temp1);
	bf_multiply( &mya.imag, &myb.imag, &temp2);
	bf_subtract( &temp1, &temp2, &c->real);
	bf_multiply( &mya.real, &myb.imag, &temp1);
	bf_multiply( &mya.imag, &myb.real, &temp2);
	bf_add( &temp1, &temp2, &c->imag);
}

/*  divide two complex numbers.
	To keep things in Cartesian form multiply top by
	conjugate of bottom.  Scale result by magnitude of
	bottom.
	
		output is c = ( a * b^*) / |b|
	
	returns 1 if b != 0, 0 if |b| = 0
*/
int bf_divide_cmplx( COMPLEX *a, COMPLEX *b, COMPLEX *c)
{
	FLOAT	mag1, mag2;
	COMPLEX	myb;
	
	bf_copy_cmplx( b, &myb);
	bf_multiply( &myb.real, &myb.real, &mag1);
	bf_multiply( &myb.imag, &myb.imag, &mag2);
	bf_add( &mag1, &mag2, &mag1);
	bf_negate( &myb.imag);
	bf_multiply_cmplx( a, &myb, c);
	if( ! bf_divide( &c->real, &mag1, &c->real)) return 0;
	bf_divide( &c->imag, &mag1, &c->imag);
	return 1;
}

/*  compute y = x^k
	where k is a signed integer in range +/-2^31 and x, y are complex.
	returns 1 if ok, 0 if x = 0 and k < 0
*/

int bf_intpwr_cmplx( COMPLEX *x, int k, COMPLEX  *y)
{
	int signflag, n;
	COMPLEX z, t;
	
/*	FLOAT seven, temp;
	
	null(&seven);
	seven.expnt = 3;
	seven.mntsa.e[MS_MNTSA] = 0x70000000;
	square_root( &seven, &seven);
	
/*  initialize Knuth's algorithm A pg 442 semi-numerical algorithms  */

	bf_copy_cmplx( x, &z);
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
	bf_null_cmplx( &t);
	bf_one( &t.real);
	while (n)
	{
		if ( n & 1 ) bf_multiply_cmplx( &t, &z, &t);
		bf_multiply_cmplx( &z, &z, &z);
		n >>= 1;
	}
	if ( signflag)
	{
		bf_null_cmplx( &z);
		bf_one( &z.real);
		return bf_divide_cmplx( &z, &t, y);
	}
	bf_copy_cmplx( &t, y);
	return 1;
}

/*  compute magnitude of a complex number.  Returns
	FLOAT result.
*/

void bf_magnitude_cmplx( COMPLEX *x, FLOAT *m)
{
	FLOAT x2, y2;
	
	bf_multiply( &x->real, &x->real, &x2);
	bf_multiply( &x->imag, &x->imag, &y2);
	bf_add( &x2, &y2, m);
	bf_square_root( m, m);
}

/*  compute exp(z) for z complex.
	z = x + iy so
	e = exp(z) = exp(x)*(cos(y) + isin(y))
	if x too large, returns 0, otherwise returns e and 1.
	works in place.
*/

int bf_exp_cmplx( COMPLEX *z, COMPLEX *e)
{
	FLOAT	x, y, xp, cy, sy;
	
	bf_copy( &z->real, &x);
	bf_copy( &z->imag, &y);
	if( !bf_exp( &x, &xp) )
	{
		bf_copy( &xp, &e->real);
		bf_null( &e->imag);
		return 0;
	}
	bf_cosine( &y, &cy);
//	printfloat("cos(y)=", &cy);
	bf_sine( &y, &sy);
//	printfloat("sin(y)=", &sy);
	bf_multiply( &xp, &cy, &e->real);
//	printfloat("exp(x)=", &xp);
	bf_multiply( &xp, &sy, &e->imag);
	return 1;
}
