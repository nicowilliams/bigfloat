/****************************************************************
*																*
*		Add on some more complex functions.  These are mostly used to find	*
*	different constants to be used in various expansions of trig and exp		*
*	functions.  Need to combine with polynomials of floats to create all terms	*
*	of a full expansion using Chebyshev polynomials and reducing to standard	*
*	polynomials.													*
*					Author = Mike Rosing							*
*					 date   = May 1, 2000							*
*																*
****************************************************************/

#include "bigfloat.h"
#include "multipoly.h"

void split( FLOAT *x, FLOAT *intprt, FLOAT *frac);

extern RAMDATA ram_block[];

MULTIPOLY	twoxcoef;		/*  2^x expansion series coefficients  */
MULTIPOLY	coscoef;		/*  cos(x)  expansion series coefficients  */
FLOAT		P2;			/*  PI/2  */
FLOAT		ln2;			/*  ln(2)	*/

/*  compute pi to 250 bits  or so.  
	Uses formula arcsin(1/2) = pi/6 = pi/2 - 1 - sum(
		1*3*5*...*(2k-1))/(2^3k (2k+1) k!)
*/
void calcpi( FLOAT *pi)
{
	FLOAT tn, constant;
	int i;
	
	null( pi);
	int_to_float( 1, &tn);
	tn.expnt = -2;
	int_to_float( 3, &constant);
	divide( &tn, &constant, &tn);
	add( pi, &tn, pi);
	for( i=2; i<123; i++)
	{
		int_to_float( 2*i-1, &constant);
		multiply(&constant, &constant, &constant);
		multiply( &constant, &tn, &tn);
		int_to_float( 2*i+1, &constant);
		divide( &tn, &constant, &tn);
		int_to_float( i, &constant);
		divide( &tn, &constant, &tn);
		tn.expnt -= 3;
		add( &tn, pi, pi);
	}
	int_to_float( 1, &constant);
	add( &constant, pi, pi);
	int_to_float( 3, &constant);
	multiply( &constant, pi, pi);
}

/*  Output of above routine is:
tn.expnt = -256
pi = 
E +3.14159265358979323846264338327950288419716939937\
	510582097494459230781640628585258
pi.expnt = 2
pi.mntsa.e[] = {
	 0x1d89cd8c,  0x105df53,  0x4533e63a,  0x94812704,  
	 0xc06e0e68,  0x62633145,  0x10b4611a,  0x6487ed51 
                     }
*/

/*  compute ln(2) to 256 bits.
	ln(2) = 2 sum( 1/( (2k+1) * 3^(2k+1) )
	which converges in about 70+ terms.
*/

void calcln2( FLOAT *ln2)
{
	int	k, epsilon, startxp;
	FLOAT	tk, constant, nine;
	
	null( ln2);
	int_to_float( 3, &constant);
	reciprical( &constant, &tk);	// gives me t0
	startxp = tk.expnt;
	multiply( &tk, &tk, &nine);	// additional 3^2k term
	epsilon = 1;
	k = 0;
	
	while( epsilon > -256)
	{
		add( &tk, ln2, ln2);
		int_to_float( 2*k+1, &constant);
		multiply( &constant, &tk, &tk);
		int_to_float( 2*k+3, &constant);
		divide( &tk, &constant, &tk);
		multiply( &nine, &tk, &tk);
		epsilon = tk.expnt - startxp;
		k++;
	}
	ln2->expnt++;		// final multiply by 2
}

/*  compute y = x^k
	where k is a signed integer in range +/- 2^31 and x, y are FLOAT.
	Copied from complex version.
	Returns 1 if ok, 0 if x = 0 and k<0
*/

int intpwr( FLOAT *x, int k, FLOAT *y)
{
	int signflag, n;
	FLOAT z, t;
	
/*  initialize Knuth's algorithm A pg 442 V2  */

	copy( x, &z);
	if( k<0)
	{
		signflag = -1;
		n = -k;
	}
	else
	{
		signflag = 0;
		n = k;
	}
	int_to_float( 1, &t);
	while( n)
	{
		if( n & 1) multiply( &t, &z, &t);
		multiply( &z, &z, &z);
		n >>= 1;
	}
	if( signflag) return reciprical( &t, y);
	copy( &t, y);
	return 1;
}

/*  This routine computes bessel functions for real arguments
	for nth order to accuracy of 256 bits.  Accuracy is easy to 
	change, assuming storage chages.  Purpose is coefficients
	for exp and trig functions approximated by chebyshev
	polynomials.
	Enter with type = +1 for In(x) and type = -1 for Jn(x),
		FLOAT pointing to x, integer n
	Returns with Float y = Jn(x) or In(x).
	uses |n|
*/

void bessel( int type, int n, FLOAT *x, FLOAT *y)
{
	FLOAT	z2, z4, constant, t1, sum;
	int		startxp, epsilon, j, k;
	
	if( n<0) n = -n;
	copy( x, &z2);
	z2.expnt--;		// divide input by 2
	multiply( &z2, &z2, &z4);		// z^2/4
	if( type < 0) negate( &z4);
	int_to_float( 1, &t1);
	j = n;			// compute 1/n!
	while( j>1)
	{
		int_to_float( j, &constant);
		multiply( &constant, &t1, &t1);
		j--;
	}
	reciprical( &t1, &sum);
	startxp = sum.expnt;
	copy( &sum, &t1);
	
/*  when next term exponent is 256 less than starting
	exponent, we have more bits of accuracy.
	j and k increment by 1 each time to create factorial
	terms 1/k! and 1/(n+j)!
*/

	j = n;
	k = 0;
	epsilon = 1;
	while( epsilon > -250)
	{
		j++;
		k++;
		int_to_float( j, &constant);
		divide( &t1, &constant, &t1);
		multiply( &z4, &t1, &t1);
		int_to_float( k, &constant);
		divide( &t1, &constant, &t1);
		add( &t1, &sum, &sum);
		epsilon = t1.expnt - startxp;
	}
	intpwr( &z2, n, &t1);
	multiply( &t1, &sum, y);
}

/* create table of chebyshev polynomials.
	Enter with maximum degree desired and array large
	enough to hold the MULTIPOLY data pointers.
	Returns array of pointers filled and maximum
	degree successfully created.
	
	T(n,x) = 2xT(n-1, x) - T(n-2, x)
*/

int gen_chebyshev( MULTIPOLY *chebary, int degree)
{
	int 		numgen;
	FLOAT	*fptr;
	MULTIPOLY twox, temp;
	
	if( degree < 0) return 0;
		
	twox.degree = 1;
	if( !get_space( &twox)) return 0;
	fptr = Address(twox);
	null( fptr);
	fptr++;
	int_to_float( 2, fptr);
	
	numgen = 0;
	chebary[0].degree = 0;
	if( !get_space( chebary))  goto chebdie;
	fptr = Address( chebary[0]);
	int_to_float( 1, fptr);
	numgen++;
	if (degree < 1) return 1;
	chebary[1].degree = 1;
	if( !get_space( &chebary[1])) goto chebdie;
	fptr = Address( chebary[1]);
	null( fptr);
	fptr++;
	int_to_float( 1, fptr);
	numgen++;
	
	while( numgen <= degree)
	{
		chebary[numgen].degree = numgen;
		if( !get_space( &chebary[numgen])) break;
		multi_mul( twox, chebary[numgen - 1], &temp);
		multi_sub( temp, chebary[numgen - 2], &chebary[numgen]);
		free_space( &temp);
		numgen++;
	}
chebdie:
	free_space( &twox);
	return numgen;
}

/*  compute coefficients of 2^x using Chebyshev polynomials.
	
	2^x  =  I(o, ln(2)) + sum{ 2*I(n, ln(2))*T(n, x)} n = 1... oo
	
	Enter with pointer to table of cheybshev polynomials, 
	maximum degree of approximation and pointer to where
	you want result.
	Returns with power series in x of above, and max degree
	actually computed.
*/

int calc_2x_coef( MULTIPOLY *cheb,  int maxdegree, 
				MULTIPOLY *twoxcoef)
{
	INDEX		i, j;
	FLOAT		ibesl, *iptr;
	MULTIPOLY	tnterm, sum;
	
	char		test[32];
	
	calcln2( &ln2);
	
	sum.degree = 0;
	if( !get_space( &sum))
	{
		printf(" no space left, calc_2x_coef \n");
		return 0;
	}
	iptr = Address(sum);
	bessel (+1, 0, &ln2, iptr);
	for( i=1; i<=maxdegree; i++)
	{
		bessel( +1, i, &ln2, &ibesl);
		ibesl.expnt++;
		multi_dup( cheb[i], &tnterm);
		for( j = (i&1); j<=i; j += 2)
		{
			iptr = Address( tnterm) + j;
			multiply( &ibesl, iptr, iptr);
		}
		if( !multi_add( tnterm, sum, &sum)) break;
		free_space( &tnterm);
	}
	if( i< maxdegree)  free_space( &tnterm);
	else i = maxdegree;
	multi_dup( sum, twoxcoef);
	free_space( &sum);
	return (i);
}

/*  compute coefficients of cos(x*PI/2) using Chebyshev polynomials.
	
	cos(x*P2)  =   J(0, P2) + 2*sum{(-1)^n*J(2n, P2)*T(2n+1, x)} n = 1... oo
	
	Enter with pointer to table of cheybshev polynomials, 
	maximum degree of approximation and pointer to where
	you want result.
	Returns with power series in x^2 of above, and max degree
	actually computed.  cos( x*PI/2) = sum{ a_j * (x^2)^j}.
	
	NOTE:  All formulas in all the books I found are WRONG.  None have
	the (-1)^n factor, but they act like it's there when constructing
	the polynomial.
*/

int calc_cos_coef( MULTIPOLY *cheb,  int maxdegree, 
				MULTIPOLY *coscoef)
{
	INDEX		i, j, k;
	FLOAT		jbesl, *jptr, *kptr;
	MULTIPOLY	tnterm, sum;
	
	char		test[32];

	calcpi( &P2);
	P2.expnt--;
	sum.degree = 0;
	if( !get_space( &sum))
	{
		printf(" no space left, calc_cos_coef \n");
		return 0;
	}
	jptr = Address(sum);
	bessel(-1, 0, &P2, jptr);
	for( i=1; i<=maxdegree/2; i++)
	{
		k = 2*i;
		bessel( -1, k, &P2, &jbesl);
		jbesl.expnt++;
		if( i&1) negate( &jbesl);
		multi_dup( cheb[k], &tnterm);
		jptr = Address( tnterm);
		for( j = 0; j<=k; j += 2)
		{
			multiply( &jbesl, jptr, jptr);
			jptr += 2;
		}
		if( !multi_add( tnterm, sum, &sum)) break;
		free_space( &tnterm);
	}
	if( i< maxdegree/2)  
	{
		free_space( &tnterm);
		i = 2*i;
	}
	else i = maxdegree;

/*  convert to polynomial in x^2  */

	j = sum.degree/2;
	coscoef->degree = j;
	if( !get_space(coscoef))
	{
		printf(" no space for cosine.\n");
		return 0;
	}
	for( k=0; k<=j; k++)
	{
		kptr = AddressOf( coscoef) + k;
		jptr = Address( sum) + 2*k;
		copy( jptr, kptr);
	}
	free_space( &sum);
	return (i);
}

/*  evaluate simple polynomial.
	Input:  MULTIPOLY coefficients, pointer to x, pointer to y
	Output:  y = F(x)
	y can equal x
*/

void polyeval( MULTIPOLY coef, FLOAT *x, FLOAT *y)
{
	INDEX 	i;
	FLOAT	sum, *cof;
	
	null( &sum);
	for( i=coef.degree; i>0; i--)
	{
		cof = Address( coef) + i;
		add( cof, &sum, &sum);
		multiply( x, &sum, &sum);
	}
	cof = Address( coef);
	add( cof, &sum, y);
}

/*  compute 2^x for x in the range -1 ... 1.
	Most inputs will be in range +/- .5 ... 1 but this routine could handle
	unnormalized inputs.
	Enter with pointer to input and storage for output
	Returns y = 2^x ( works in place, both pointers can be the same)
*/

void twoexp( FLOAT *x, FLOAT *y)
{
	polyeval( twoxcoef, x, y);
}

/*  compute cos( x) for x in range +/- PI/2
	As above, works in place.
	This is a *core* routine, no range checking!
*/

void corecos(FLOAT *x, FLOAT *y)
{
	FLOAT	x2;
	
	divide( x, &P2, &x2);
	multiply( &x2, &x2, &x2);
	polyeval( coscoef, &x2, y);
}

/*  convert a float to a long.  Overflow is max
	possible result.
*/

int	float_to_int( FLOAT *f)
{
	FLOAT	dummy, intprt;
	int		value;
	
	if( f->expnt < 1) return 0;
	if( f->expnt > 31)
	{
		if( f->mntsa.e[MS_MNTSA] & SIGN_BIT)
			return SIGN_BIT;
		return ~SIGN_BIT;
	}
	split( f, &intprt, &dummy);
	if( intprt.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		negate( &intprt);
		value = -(intprt.mntsa.e[MS_MNTSA] >> ( 31 - intprt.expnt));
	}
	else	value = intprt.mntsa.e[MS_MNTSA] >> ( 31 - intprt.expnt);
	return value;
}
	
/*  compute e^x for any x.  |x| > 2^32/ln(2) will overflow
	and return max possible value and 0.
	Otherwise returns y = exp(x) and 1.
	works in place.
*/

int exp( FLOAT *x, FLOAT *y)
{
	FLOAT	z, xp;
	long		xpnt;
	INDEX	i;
	
/*  convert to base 2  */

	divide( x, &ln2, &z);
	
/*  check range is possible to do  */

	if (z.expnt > 32)
	{
		if( x->mntsa.e[MS_MNTSA] & SIGN_BIT)
			null(y);
		else
		{
			OPLOOP(i) y->mntsa.e[MS_MNTSA] = ~0;
			y->mntsa.e[MS_MNTSA] >>= 1;
			y->expnt = ~0 >> 1;
		}
		return 0;
	}
	
/*  we can perform operation, send z mods 2 to core */

	split(&z, &xp, &z);
	xpnt = float_to_int( &xp);
	twoexp( &z, y);
	
/*  next add xpnt to exponent of y  */

	y->expnt +=  xpnt;
	return 1;
}

/*  split a FLOAT into its integer and fractional parts  */

void split( FLOAT *x, FLOAT *intprt, FLOAT *frac)
{
	FLOAT	fracpart;
	INDEX	i, signflag;
	ELEMENT	mask, xpchk;

/*  if number < 1, return just a fraction  */
	
	if( x->expnt <= 0)
	{
		null( intprt);
		copy( x, frac);
		return;
	}

/*  zero out 31 bits of ms word, then one block at a time */

	copy( x, &fracpart);
	signflag = 0;
	if( fracpart.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		negate( &fracpart);
		signflag = 1;
	}
	i = MS_MNTSA;
	xpchk = fracpart.expnt;
	if( xpchk > 31)
	{
		fracpart.mntsa.e[i] = 0;
		i--;
		xpchk -= 31;
	}
	while( (xpchk > 32) && ( i>0 ))
	{
		fracpart.mntsa.e[i] = 0;
		i--;
		xpchk -= 32;
	}

/*  next zero out the remaining integer bits in a left over word  */

	if( i != MS_MNTSA) mask = ~0UL >> xpchk;
	else 	mask = ( ~0UL >> (xpchk+1));
	fracpart.mntsa.e[i] &= mask;
	if( signflag) negate( &fracpart);
	normal( &fracpart);
	subtract( x, &fracpart, intprt);
	copy( &fracpart, frac);
}

/*  compute cosine(x) for any x.
	x values larger than PI*2^200 will be in gross error, so watch out!
	works in place, returns y = cos(x)
*/

void cosine( FLOAT *x, FLOAT *y)
{
	FLOAT	z, PI, dummy, PI3;
	int		cmpr;
	
/*  create 2*PI  */

	copy( &P2, &PI);
	PI.expnt += 2;

/*  check range of input and force modulo 2PI operation  */

	cmpr = compare( x, &PI);
	if( cmpr > 0 )
	{
		divide( x, &PI, &z);
		split( &z, &dummy, &z);
		multiply( &PI, &z, &z);
	}
	else copy( x, &z);
	if( z.mntsa.e[MS_MNTSA] & SIGN_BIT) negate( &z);

/*  z is now in range 0...2PI.  Now convert to range of core cos */

	cmpr = compare( &z, &P2);
	if( cmpr <= 0)
	{
		corecos( &z, y);
		return;
	}
	PI.expnt--;
	add( &PI, &P2, &PI3);	// 3 PI/2
	cmpr = compare( &z, &PI3);
	if( cmpr > 0)
	{
		PI.expnt++;
		subtract( &PI, &z, &z);	// 2PI - x
		corecos( &z, y);
		return;
	}
	subtract( &PI, &z, &z);	// PI - x
	corecos( &z, y);
	negate( y);
}

/*  compute sine(x) for any x.
	same as cosine, jus move arguments around.
	works in place, 
	returns y = sin(x)
*/

void sine( FLOAT *x, FLOAT *y)
{
	FLOAT	z, PI, dummy;
	int		cmpr, signflag;
	
/*  create 2*PI and reduce x modulo 2PI signed  */

	copy( &P2, &PI);
	PI.expnt += 2;
	cmpr = compare( x, &PI);
	if( cmpr > 0)
	{
		divide( x, &PI, &z);
		split( &z, &dummy, &z);
		multiply( &PI, &z, &z);
	}
	else copy( x, &z);
	if( z.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		negate( &z);
		signflag = 1;
	}
	else signflag = 0;

/*  z is no in range 0... 2*PI.
	if bigger than PI, flip sign of result and fold back to 0..PI,
	then subtract PI/2 to put into corecos range.
*/
	PI.expnt--;
	cmpr = compare( &z, &PI);
	if( cmpr > 0)
	{
		signflag ^= 1;
		subtract( &z, &PI, &z);
	}
	subtract( &z, &P2, &z);
	corecos( &z, y);
	if( signflag) negate(y);
}
