/********************************************************************************
*																				*
*		The functions in this file use the basic large integer arithmetic sub-	*
*	routines to perform such things as greatest common divisor, modular 		*
*	inversion and modular exponentiation.  Note that many other math packages	*
*	do all these things much better, but this is very easy to follow.			*
*																				*
********************************************************************************/

#include <stdio.h>
#include "bigint.h"

/*  divide a large integer by 2.  A simple shift right operation.  */

void int_div2( x)
BIGINT *x;
{
	INDEX j;
	ELEMENT mask;
	
	INTLOOP(j)
	{
		if (j) mask = ( x->hw[j-1] & 1) ? CARRY : 0;
		else mask = 0;
		x->hw[j] = (x->hw[j] | mask) >> 1;
	}
}

/* compute greatest common divisor using binary method.
	See [DEK, pg 321] Algorithm B for theoretical details.

	Enter with large integers u, v.  
	Returns gcd(u, v) as large integer in w.
*/

void int_gcd(u, v, w)
BIGINT *u, *v, *w;
{
	INDEX  k, i, flag;
	ELEMENT check, carry_bit;
	BIGINT t, U, V;
	
	int_copy( u, &U);
	int_copy( v, &V);
	
/*  find common powers of 2 and eliminate them */

	k = 0;
	
/*  check that both u and v are even */

	while ( !(U.hw[INTMAX] & 1 || V.hw[INTMAX] & 1))   
	{
	/*  increment power of 2 and divide both u and v by 2 */
		k++;					
		int_div2( &U);
		int_div2( &V);
	}
	
/* Now both u and v have been divided by 2^k.  
	If u is odd, set t = v and flag, otherwise t = u and clear flag  */

	if (U.hw[INTMAX] & 1)
	{
		int_copy( &V, &t);
		flag = -1;
	}
	else 
	{
		int_copy( &U, &t);
		flag = 1;
	}
	
	check = 0;
	INTLOOP (i) check |= t.hw[i];
	while (check)
	{
	/* while t is even, divide by 2  */
	
		while ( !(t.hw[INTMAX] & 1)) int_div2( &t);  

/* reset u or v to t depending on sign of flag  */

		if (flag > 0) int_copy( &t, &U);
		else int_copy( &t, &V);
/*		t = u - v;			  core reduction step, gcd remains unchanged  */
		int_sub( &U, &V, &t);
		if (t.hw[0] & MSB_HW)
		{
			flag = -1;
			int_neg( &t);
		}
		else flag = 1;
		check = 0;
		INTLOOP (i) check |= t.hw[i];
	}

	/*  reapply common powers of 2. First do words, then do bits.*/
	
	int_copy( &U, w);
	while ( k > HALFSIZE )
	{
		for (i=0; i<INTMAX; i++) w->hw[i] = w->hw[i+1];
		k -= HALFSIZE;
		w->hw[INTMAX] = 0;
	}
	carry_bit = 0;
	while ( k > 0 )
	{
		INTLOOP (i)
		{
			w->hw[i] = (w->hw[i] << 1) | carry_bit;
			carry_bit = w->hw[i] & CARRY ? 1 : 0;
			w->hw[i] &= LOMASK;
		}
		k--;
	}
}

/*  Binary method for modular exponentiation.  Taken from [DEK,pg 442] Algorithm A.
	Computes z = x^n mod q for x, n and q large integers.
*/

void mod_exp(x, n, q, z)
BIGINT *x, *n, *q, *z;
{
	BIGINT  N, Y, Z, temp, dummy;
	ELEMENT check;
	INDEX   i;
	
/*  initialize variables  */

	int_copy (n, &N);
	int_null( &Y);
	Y.hw[INTMAX] = 1;
	int_copy (x, &Z);

/*  Main loop divides N by 2 each step.  Repeat until N = 0, and return Y as result.  */

	check = 0;
	INTLOOP (i) check |= N.hw[i];
	while (check)
	{

/*  if N is odd, multiply by extra factor of Y */

		if (N.hw[INTMAX] & 1) 
		{
			/*  Y = (Y * Z) % q;  */
			int_mul (&Y, &Z, &temp);
			int_div (&temp, q, &dummy, &Y);
		}
		int_div2( &N);					/* divide N by 2 */
	/*		Z = (Z * Z) % q;		  square Z  */
		int_mul (&Z, &Z, &temp);
		int_div( &temp, q, &dummy, &Z);
		check = 0;
		INTLOOP (i) check |= N.hw[i];
	}
	int_copy (&Y, z);
}

/*  Inversion of numbers in a prime field is similar to solving the linear congruence
	ax = c mod b.  Replace c with 1 and x is the inverse of a mod b.  Taken from [HR,
	pg 309].
	Inputs are large integer a and modulus b.
	Output is x, inverse of a mod b  (ax = 1 mod b).
*/

void mod_inv(a, b, x)
BIGINT *a, *b, *x;
{
	BIGINT  m, n, p0, p1, p2, q, r, temp, dummy;
	ELEMENT check;
	INDEX   sw, i;
	
/*  initialize loop variables  

	sw = 1;
	m = b;
	n = a;
	p0 = 1;
	p1 = m/n;
	q = p1;
	r = m % n;
*/
	sw = 1;
	int_copy( b, &m);
	int_copy( a, &n);
	int_null ( &p0);
	p0.hw[INTMAX] = 1;
	int_div ( &m, &n, &p1, &r);
	int_copy ( &p1, &q);
	
/*  main loop, compute continued fraction  intermediates  */

	check = 0;
	INTLOOP (i) check |= r.hw[i];
	while (check)
	{
		sw = -sw;
		int_copy( &n, &m);
		int_copy( &r, &n);
		int_div( &m, &n, &q, &r);
 /*		p2 = (q * p1 + p0) % b;   core operation of routine  */
 		int_mul( &q, &p1, &temp);
 		int_add( &temp, &p0, &temp);
 		int_div( &temp, b, &dummy, &p2);
 		int_copy( &p1, &p0);
 		int_copy( &p2, &p1);
		check = 0;
		INTLOOP (i) check |= r.hw[i];
	}
	
/*  sw keeps track of sign.  If sw < 0, add modulus to result */

	if (sw < 0) int_sub( b, &p0, x);
	else int_copy( &p0, x);
}

