/*  Basic floating point package with lots of bits.  Use header parameter MNTSA_SIZE
	to change number of bits in mantissa to 32*MNTSA_SIZE.  Works best if even.

  The purpose of this code is to follow the mathematics described in Joe Silverman's
  "Advanced Topics in Arithmetic of Elliptic Curves".  Takes advantage of gcc compiler
  and I hope the underlying hardware.
  
					Author = Mike Rosing
					  date  = Feb. 17, 2000
*/

#include <stdio.h>
#include "bigfloat.h"
#include "multipoly.h"

extern MULTIPOLY	twoxcoef;
extern MULTIPOLY	coscoef;
extern FLOAT		P2;
extern FLOAT		ln2;


/*  copy a floating point value from a to b  */

void copy( FLOAT *a, FLOAT *b)
{
	int i;
	
	OPLOOPd(i) b->mntsa.d[i] = a->mntsa.d[i];
	b->expnt = a->expnt;
}

/* copy a complex value from a to b  */

void copy_cmplx( COMPLEX *a, COMPLEX *b)
{
	copy( &a->real, &b->real);
	copy( &a->imag, &b->imag);
}

void null( FLOAT *a)
{
	int i;
	
	a->expnt = 0;
	OPLOOP(i) a->mntsa.e[i] = 0;
}

/*  negate a value in place  */

void negate( FLOAT *a)
{
	int i;
	unsigned long long int carry;
	
	OPLOOPd(i)  a->mntsa.d[i] ^= ~0;
	a->mntsa.d[0]++;
	if ( a->mntsa.d[0] ) return;
	i=1;
	carry = 1LL;
	while ( carry && ( i < D_SIZE))
	{
		a->mntsa.d[i]++;
		carry = 0LL;
		if ( !a->mntsa.d[i] ) carry = 1LL;
		i++;
	}
}

/*  normalize a floating point number.  Called at the end of add and multiply
	(and divide?),  it shifts results up till sign bit clear and next bit set.
	binary point sits after sign bit, so all fractions are in range 1/4 to 1/2.
	exponent has +/- 2^31 range, on underflow or overflow of exponent
	you ought to flag an error.
*/
void normal(  FLOAT *x)
{
	int	i, j, signflag, upshift, downshift;
	long	xpnt;
	unsigned long mask, uprmask, msb, temp;
	
	signflag = 0;
	xpnt = x->expnt;
	if( x->mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		negate(x);
		signflag = 1;
	}
	for( i=MS_MNTSA; i>=0; i--)  if ( x->mntsa.e[i]) break;
	if ( i<0 )
	{
		x->expnt = 0;  //  result is zero
		return;
	}
	if ( i != MS_MNTSA) 
	{
		xpnt -= 32*(MS_MNTSA - i);
		for ( j=MS_MNTSA; i >= 0; i--)
		{
			x->mntsa.e[j] = x->mntsa.e[i];
			x->mntsa.e[i] = 0;
			j--;
		}
	}
	
/*  find most significant bit.  Same trick as in "degreeof" subroutine  */

	msb = x->mntsa.e[MS_MNTSA];
	downshift = 0;
	mask = ~0;
	for( i=16; i>0; i >>= 1)
	{
		mask ^= mask >> i;
		if ( mask & msb)
		{
			downshift += i;
			msb &= mask;
		}
	}
	if ( downshift == 31 ) // then right shift 1 bit
	{
		for (i=0; i<MS_MNTSA; i++)
		{
			msb = (x->mntsa.e[i+1] & 1) ? SIGN_BIT : 0;
			x->mntsa.e[i] = ( x->mntsa.e[i] >> 1) | msb;
		}
		x->mntsa.e[MS_MNTSA] >>= 1;
		x->expnt = xpnt + 1;
		goto	normlrtn;
	}
	if ( downshift == 30 ) goto normlrtn;
	upshift = 30 - downshift;
	downshift += 2;
	x->mntsa.e[MS_MNTSA] <<= upshift;
	uprmask = ~0 << downshift;
	mask = ~uprmask;
	for( i=MS_MNTSA; i>0; i--)
	{
		x->mntsa.e[i] |= (x->mntsa.e[i-1] & uprmask) >> downshift;
		x->mntsa.e[i-1] = (x->mntsa.e[i-1] & mask) << upshift;
	}

/*  check if input was negative and return  corrected exponent  */

	xpnt -= upshift;
normlrtn:
	x->expnt = xpnt;
	if (signflag) negate (x);
}

/*  compare magnitude of 2 FLOATs.  
	Returns:
		+1 if |a| > |b|
		-1 if |a| < |b|
		  0 if |a| = |b|
*/

int compare( FLOAT *a, FLOAT *b)
{
	FLOAT mya, myb;
	int i;
	
/*  first compare exponents, takes care of most cases  */

	if( a->expnt > b->expnt) return 1;
	if( a->expnt < b->expnt) return -1;

/*  exponents match, check mantissas  */

	copy( a, &mya);
	copy( b, &myb);
	if( mya.mntsa.e[MS_MNTSA] & SIGN_BIT)
		negate( &mya);
	if( myb.mntsa.e[MS_MNTSA] & SIGN_BIT)
		negate( &myb);
	for( i=MS_MNTSA; i >= 0; i--)
	{
		if( mya.mntsa.e[i] > myb.mntsa.e[i]) return 1;
		if( mya.mntsa.e[i] < myb.mntsa.e[i]) return -1;
	}
	return 0;
}

/*  add two floating point numbers a + b = c.
	Makes local copies of data, so any pointers can be the same.
*/
void add( FLOAT *a, FLOAT *b, FLOAT *c)
{
	FLOAT big, small, result;
	int	i, uper, lower, bigsign, resultsign, smallsign;
	unsigned long	msb, himask, lomask;
	unsigned long long carry, temp, carrycheck;

/*  eliminate work if either input is zero  */

	if (iszero (a)) 
	{
		copy (b, c);
		return;
	}
	if (iszero( b))
	{
		copy( a, c);
		return;
	}
	
/* easy to deal with one big and one small, so copy accordingly 
	but equal magnitude and opposite signs causes problems */

	i = compare( a, b);
	if( i == 0)		/*  a = b, very special case  */
	{
		if( a->mntsa.e[MS_MNTSA] ^ b->mntsa.e[MS_MNTSA])
		{
			null( c);
			return;
		}
		copy( a, c);
		c->expnt++;
		return;
	}
	if (i > 0)
	{
		copy( a, &big);
		copy( b, &small);
	}
	else
	{
		copy( b, &big);
		copy( a, &small);
	}
 	if ( (big.expnt - small.expnt) >= (32*MNTSA_SIZE) )
	{
		copy( &big, c);
		return;
	}
	bigsign = big.mntsa.e[MS_MNTSA] & SIGN_BIT ? 1 : 0;
	smallsign = small.mntsa.e[MS_MNTSA] & SIGN_BIT ? 1 : 0;
	
/*  unnormalize small number to align bits, big chunks first  */

	while( (big.expnt - small.expnt) > 31)
	{
		for( i=0; i<MS_MNTSA; i++)
			small.mntsa.e[i] = small.mntsa.e[i+1];
		small.expnt += 32;
		if ( small.mntsa.e[MS_MNTSA] & SIGN_BIT)
			small.mntsa.e[MS_MNTSA] = ~0;
		else
			small.mntsa.e[MS_MNTSA] = 0;
	}

/*  Now move chunks down to finish alignment  */
			
	if( big.expnt > small.expnt)
	{
		uper = big.expnt - small.expnt;
		lower = 32 - uper;
		himask = ~0 << uper;
		lomask = ~himask;
		for( i=0; i<MS_MNTSA; i++)
			small.mntsa.e[i] = (( small.mntsa.e[i] & himask) >> uper) |
							(( small.mntsa.e[i+1] & lomask) << lower)
		;
		msb = ( small.mntsa.e[MS_MNTSA]  & himask) >> uper;
		if ( small.mntsa.e[MS_MNTSA] & SIGN_BIT)
			msb |= ( ~0 << (31 - uper));
		else
			msb &= ( ~0 >> (uper + 1));
		small.mntsa.e[MS_MNTSA] = msb;
	}

/*  mantissas aligned, add everything up.  Propagate carry too.  */

	carry = 0;
	carrycheck = (unsigned long long)SIGN_BIT << 1;
	OPLOOP(i)
	{
		temp = (unsigned long long)big.mntsa.e[i] +
				 (unsigned long long)small.mntsa.e[i] + carry;
		if (temp & carrycheck ) carry = 1;
		else carry = 0;
		result.mntsa.e[i] = temp & 0xffffffff;
	}
	result.expnt = big.expnt;
	resultsign = result.mntsa.e[MS_MNTSA] & SIGN_BIT ? 1 : 0;
	
/*  check for overflow and shift down if needed  */

	if( carry)
	{
		if( (!bigsign && resultsign) || (bigsign && !resultsign) 
			|| (bigsign && smallsign)) 
		{
			himask =  SIGN_BIT;
			for( i = MS_MNTSA; i >= 0; i--)
			{
				msb = result.mntsa.e[i] & 1;
				result.mntsa.e[i] = (result.mntsa.e[i] >> 1) | himask;
				himask = msb ? SIGN_BIT : 0;
			}
		}
		if(bigsign && !resultsign || (bigsign && smallsign)) result.expnt++;

/*  if we have a carry and full null words, then take care of exponent as well */

		i = MS_MNTSA;
		while( i > 0 && (!result.mntsa.e[i])) i--;
		if( (i != MS_MNTSA) &&
				(result.mntsa.e[i] & SIGN_BIT) ) result.expnt++;
	}
	else if( ( !bigsign && resultsign)  ) 
	{
		result.expnt++;
		himask =  0;
		for( i = MS_MNTSA; i >= 0; i--)
		{
			msb = result.mntsa.e[i] & 1;
			result.mntsa.e[i] = (result.mntsa.e[i] >> 1) | himask;
			himask = msb ? SIGN_BIT : 0;
		}
	}
	normal( &result);
	copy( &result, c);
}

/*  because it's useful  c = a - b  */

void subtract( FLOAT *a, FLOAT *b, FLOAT *c)
{
	FLOAT	myb;
	
	copy( b, &myb);
	negate( &myb);
	add( a, &myb, c);
}

/*  round a float to an integer.  add 1/2 and clear fractional bits.  */

void round( FLOAT *a, FLOAT *b)
{
	FLOAT  half;
	int i;
	unsigned long mask;
	
	null( &half);
	half.mntsa.e[MS_MNTSA] = 0x40000000;
	add( a, &half, b);
	i = MS_MNTSA - (b->expnt / 32);
	mask = ~0 << (31 - b->expnt % 32);
	b->mntsa.e[i] &= mask;
	i--;
	while( i >= 0 )
	{
		b->mntsa.e[i] = 0;
		i--;
	}
}
	
/*  multiply two FLOATS to get a third.
	Uses mtheod described by Crenshaw in Embedded Systems Magazine
	March 1997
	Returns c = a * b, 
*/

void multiply( FLOAT *a, FLOAT *b, FLOAT *c)
{
	int	i, j, k, signflag;
	FLOAT	mya, myb;
	unsigned long long int mult, src, dst;

/*  figure out sign of result and use unsigned algorithm  */

	copy( a, &mya);
	copy( b, &myb);
	signflag = 0;
	if( mya.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		signflag = 1;
		negate( &mya);
	}
	if ( myb.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		signflag ^= 1;
		negate( &myb);
	}

/*  compute unnormalized exponent  */

	null( c);
	c->expnt = mya.expnt + myb.expnt + 1;

/*  use longs and multiply up long longs then sum to 
	correct place.  Carry only propagates to next word
	and by order of computation we easily track overflow.
*/
	OPLOOP(i)
	{

/*  compute least significant long.  this will be the high
	order word of the mid point of a full double length
	multiply.  Round down to zero.
*/
		j = MS_MNTSA - i;
		mult = ( (unsigned long long)mya.mntsa.e[i] * 
				(unsigned long long)myb.mntsa.e[j] ) >> 32;
		c->mntsa.d[0] += mult;
	}

/*  compute all full longs and add in correct place  */

	for ( i=1; i<MNTSA_SIZE; i++)  // i is result index
	{
		j = i;
		k = MS_MNTSA;
		while ( j < MNTSA_SIZE )
		{
			src = ((unsigned long long) c->mntsa.e[i]  << 32) | 
					( (unsigned long long)c->mntsa.e[i-1]);
			mult = (unsigned long long)mya.mntsa.e[j] * 
					(unsigned long long)myb.mntsa.e[k];
			dst = src + mult;
			if ( dst < src && (i < MS_MNTSA)) c->mntsa.e[i+1]++;
			j++;
			k--;
			c->mntsa.e[i-1] = dst & 0xffffffffLL;
			c->mntsa.e[i] =  (dst >> 32) & 0xffffffffLL;
		}
	}
	normal( c);
	if (signflag) negate( c);
}

/*  divide FLOATS
	Second cut, use Newton-Raphson method described by
	Oberman and Flynn in CSL-TR-95-675 (Stanford Computer
	systems lab umunhum.stanford.edu/tr/oberman.jul95.tr675.ps.Z)
	computes 1/b = c
	copies data, so any arguments can be the same.
	returns 0 if attempting to divide by zero
	1 otherwise.
*/

int reciprical ( FLOAT *b, FLOAT *c)
{
	FLOAT	myb, x0, x1, two;
	unsigned long long int utop;
	int		i, signflag, digit;
	long		exponent;
		
/*  copy bottom and make positive  */

	signflag = 0;
	copy( b, &myb);
	exponent = myb.expnt;
	if( myb.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		signflag = 1;
		negate( &myb);
	}

/*  check for divide by zero  */

	digit = 0;
	OPLOOP(i)
	{
		if( myb.mntsa.e[i] )
		{
			digit = 1;
			break;
		}
	}
	if ( !digit) return 0;

/*  create constant 2  */

	null( &two);
	two.mntsa.e[MS_MNTSA] = 0X40000000;
	two.expnt = 2;

/*  convert bottom to fraction and find 32 bit first
	guess.  Guess will be normalzied automaticly
	because it started that way.
*/
	null( &x0);
	utop = 0x4000000000000000 / 
		(unsigned long long) myb.mntsa.e[MS_MNTSA];

/*  There is only one possible case for overflow.
	Deal with it.
*/
	if( utop & 0xf00000000) utop = 0xffffffff;
	x0.mntsa.d[MS_MNTSAd] = utop << 31;
	x0.expnt = 1 - myb.expnt;

/*  using guess, compute twice as many bits each step. */

	for( i=0; i<DIVISION_LOOPS; i++)
	{
		multiply( &myb, &x0, &x1);
		subtract( &two, &x1, &x1);
		multiply( &x0, &x1, &x0);
	}

/*  adjust exponent and check sign flag of result */

	copy( &x0, c);
	if (signflag) negate(c);
	return 1;
}

/*  Divide is pretty trivial with reciprical  */

int divide( FLOAT *a, FLOAT *b, FLOAT *c)
{
	FLOAT bottom;
	
	if ( !reciprical( b, &bottom)) return 0;
	multiply( a, &bottom, c);
}

/*  Square root function.  Start with constants.  Taken 
	from Computer Approximations 0293.
	I can't seem to create constants since mntsa is a
	union.  Anybody know how to do that?
*/

static FLOAT p0, p1, p2, q0, q1;

void init_float()
{
	int			degree;
	MULTIPOLY	chebary[ MAXCHEB+1];
	
/*	ascii_to_float("E 0.1767767142", &p0);
	ascii_to_float("E 3.696790108", &p1);
	ascii_to_float("E 3.641977651", &p2);
	ascii_to_float("E 1.287633631", &q0);
	ascii_to_float("E 5.228050594", &q1);
*/	
	p0.expnt = -2;
	p0.mntsa.e[0] = 0xb1f96c1d;
	p0.mntsa.e[1] = 0x4e08339b;
	p0.mntsa.e[2] = 0x4b8e3962;	/*  0.1767767142 */
	p0.mntsa.e[3] = 0x67ab2e18;
	p0.mntsa.e[4] = 0x5f4409fa;
	p0.mntsa.e[5] = 0xa9df1663;
	p0.mntsa.e[6] = 0x5dde9567;
	p0.mntsa.e[7] = 0x5a827a3c;
	
	p1.expnt = 2;
	p1.mntsa.e[0] = 0xa3f42a4b;
	p1.mntsa.e[1] = 0x90ffdbdc;
	p1.mntsa.e[2] = 0x9ec61816;	/*  3.696790108 */
	p1.mntsa.e[3] = 0x6c9416bc;
	p1.mntsa.e[4] = 0x2bc747ec;
	p1.mntsa.e[5] = 0xe49d67c9;
	p1.mntsa.e[6] = 0xc1296f53;
	p1.mntsa.e[7] = 0x764c1ac4;
	
	p2.expnt = 2;
	p2.mntsa.e[0] = 0x66caedd4;
	p2.mntsa.e[1] = 0x7ce7bac2;
	p2.mntsa.e[2] = 0x460c8755;	/*  3.641977651 */
	p2.mntsa.e[3] = 0xdfd03d4d;
	p2.mntsa.e[4] = 0x3ee8fc8b;
	p2.mntsa.e[5] = 0x567e6f7a;
	p2.mntsa.e[6] = 0xf9da54a8;
	p2.mntsa.e[7] = 0x748b14b6;
	
	q0.expnt = 1;
	q0.mntsa.e[0] = 0xb407eb12;
	q0.mntsa.e[1] = 0x65f085f5;
	q0.mntsa.e[2] = 0xd006fa2b;	/*  1.287633631 */
	q0.mntsa.e[3] = 0xddc3d136;
	q0.mntsa.e[4] = 0x4e8ce714;
	q0.mntsa.e[5] = 0xac49b0b4;
	q0.mntsa.e[6] = 0x97fb9afc;
	q0.mntsa.e[7] = 0x526896e3;
	
	q1.expnt = 3;
	q1.mntsa.e[0] = 0xe3263bb7;
	q1.mntsa.e[1] = 0xe6ef73fc;
	q1.mntsa.e[2] = 0x19b79d4b;	/*  5.228050594 */
	q1.mntsa.e[3] = 0x4dbc06bd;
	q1.mntsa.e[4] = 0x73b36731;
	q1.mntsa.e[5] = 0x580a1365;
	q1.mntsa.e[6] = 0x31039445;
	q1.mntsa.e[7] = 0x53a61861;

/*  initialize constants for exp and cosine expansions  */
	
	degree = gen_chebyshev( chebary, MAXCHEB);
	if( degree < MAXCHEB)
	{
		printf("Max degree %d obtained for chebyshev array. \n", degree);
		exit(0);
	}
	degree = calc_2x_coef( chebary, 44, &twoxcoef);
	if( degree < 44)
	{
		printf("Maxdegree %d obtained for 2^x coefficints. \n", degree);
		exit(0);
	}
	degree = calc_cos_coef( chebary, MAXCHEB, &coscoef);
	if( degree < MAXCHEB)
	{
		printf("Maxdegree %d obtained for cos coefficints. \n", degree);
		exit(0);
	}

}

/*  compute real square root of a FLOAT.
	Enter with pointers to source and destination areas.
	They can be the same, and this will work.
	Negative inputs resolved as positive, no errors.
	Uses Heron's algorithm, see Computer Approximations
	pg 90
*/

void square_root( FLOAT *in, FLOAT *out)
{
	FLOAT x, top, bottom, y;
	int i;
	
	char debug[256];
	int j;
	
	copy( in, &x);

/*  check sign and range of input.
	convert to fraction in range 0.25 < x < 1
	and adjust exponent accordingly.
*/
	if( in->mntsa.e[MS_MNTSA] & SIGN_BIT)
		negate( &x);
	if( x.expnt & 1)
	{
		out->expnt = (x.expnt + 1)/2;
		x.expnt = -1;
	}
	else
	{
		out->expnt = x.expnt/2;
		x.expnt = 0;
	}
	
/*  compute 4.7 digits = 15 bits for first estimate.
	Heron's iteration doubles accuracy every step,
	so 5 steps gives 500 bits.  Adjust loop as needed.
*/

	multiply( &p2, &x, &top);
	add( &p1, &top, &top);
	multiply ( &top, &x, &top);
	add ( &p0, &top, &top);
	
	add( &x, &q1, &bottom);
	multiply(  &x, &bottom, &bottom);
	add( &q0, &bottom, &bottom);
	
	divide( &top, &bottom, &y);
	for( i=0; i<5 ; i++)
	{
		divide( &x, &y, &top);
		add( &y, &top, &y);
		y.expnt--;
	}

/*  return properly scaled result  */

	OPLOOP(i) out->mntsa.e[i] = y.mntsa.e[i];
}

/*  convert signed 32 bit integer to a float  */

void int_to_float( int num, FLOAT *x)
{
	null( x);
	x->expnt = 31;
	x->mntsa.e[MS_MNTSA] = num;
	normal( x);
}

