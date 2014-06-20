/*  Some basic input output routines to make use of FLOAT variables.  */

#include <stdio.h>
#include "bigfloat.h"

char bf_digitof( FLOAT *frac);

double logbase10of2 = 0.30102999566398119521;

/*  The purpose of this routine is to convert an ascii string to FLOAT.
	Input format makes the exponent important and the mantissa
	secondary, as argued by an old professor of mine long ago.
	
	Format is:  Esxxxxxxxxx\tsxxxxx----.xxxxx----
	
	where initial E can be lower case but must be present.
	s is either + or - or null (not space!)
	\t is space or tab
	exponent can have as many digits as you want but more than 9 is useless
	mantissa can have decimal point anywhere.
	Garbage characters ignored, but if in mantissa assumed to be decimal point.
	
	returns 0 if can't parse number
		     1 if it can.
*/

int bf_ascii_to_float( char *instring, FLOAT *outnum)
{
	int	signflag, exponent;
	char	nextchar;
	FLOAT	digit, ten, scale;
	
/*  look for start of number (E or e)  */

	while ( nextchar = *instring++)
		if ( nextchar == 'e' || nextchar == 'E') break;
	if ( !nextchar ) return 0;

/*  parse exponent.  If no number after exponent return error.
	If blank after the E, exponent assumed 0.
*/
	signflag = 0;
	bf_null( outnum);
	nextchar = *instring++;
	if ( nextchar == '-')
	{
		signflag = -1;
		nextchar = *instring++;
	}
	if ( nextchar == '+' ) nextchar = *instring++;
	exponent = 0;
expnt:
	while ( nextchar >= '0' && nextchar <= '9')
	{
		exponent *= 10;
		exponent += nextchar & 0xf;
		nextchar = *instring++;
	}
	if ( !nextchar) return 0;

/*  ignore garbage characters and keep parsing exponent until death  */

	if ( nextchar != ' ' && nextchar != '\t' )
	{
		nextchar = *instring++;
		goto expnt;
	}
	if ( signflag) exponent = -exponent;

/*  now parse mantissa
	Decimal point can be anywhere.  first check sign, then add in each digit
*/

	while ( nextchar == ' ' || nextchar == '\t')
	{
		nextchar = *instring++;
		if( !nextchar) return 0;
	}
	bf_null( &ten);
	ten.expnt = 4;
	ten.mntsa.e[MS_MNTSA] = 0x50000000;
	signflag = 0;
	if( nextchar == '-' )
	{
		signflag = -1;
		nextchar = *instring++;
	}
	if ( nextchar == '+' ) nextchar = *instring++;
	while ( nextchar >= '0' && nextchar <= '9')
	{
		bf_multiply( outnum, &ten, outnum);
		bf_null ( &digit);
		digit.expnt = 31;
		digit.mntsa.e[MS_MNTSA] = nextchar & 0xf;
		bf_normal ( &digit);
		bf_add( &digit, outnum, outnum);
		nextchar = *instring++;
	}
	if ( !nextchar ) goto scale;

/*  now parse digits after the decimal point.  scale keeps track of how many
	digits we've collected.  You'd have to have a lot of digits to underflow!!
*/
	nextchar = *instring++;
	bf_null( &scale);
	scale.expnt = 1;
	scale.mntsa.e[MS_MNTSA] = 0x40000000;
	bf_divide( &scale, &ten, &scale);
	while (nextchar)
	{
		if( nextchar >= '0' && nextchar <= '9')
		{
			bf_null ( &digit);
			digit.expnt = 31;
			digit.mntsa.e[MS_MNTSA] = nextchar & 0xf;
			bf_normal ( &digit);
			bf_multiply( &digit, &scale, &digit);
			bf_add( &digit, outnum, outnum);
			bf_divide( &scale, &ten, &scale);
		}
		nextchar = *instring++;
	}

/*  combine exponent with mantissa and convert powers of 10 to
	powers of 2.
*/

scale:
	if ( exponent)
	{
		if ( exponent < 0)
		{
			while (exponent)
			{
				bf_divide( outnum, &ten, outnum);
				exponent++;
			}
		}
		else
		{
			while ( exponent )
			{
				bf_multiply( outnum, &ten, outnum);
				exponent--;
			}
		}
	}
	if ( signflag ) bf_negate( outnum);
	return 1;
}

/*  convert  FLOAT to human ascii.  Output format is:

	Esxxxxxxxxx\bsx.xxxxx-----
	
where s is + or -, \b is a space

	Enter with pointer to float and pointer to result space.
	Maximum length of string depends on MNTSA_SIZE in bigfloat.h
	256 bits ~ 77 digits, so 100 bytes is good.  (cut off set at 80)
	Scale room accordingly.  
	Conversion based on Knuth, "Radix Conversion" in Seminumerical
*/

void bf_float_to_ascii( FLOAT *numbr, char *outstring)
{
	long exponent, exp10;
	FLOAT ten, fraction, scale;
	unsigned long mask, lastbit;
	int signflag, i;
	char  *digit, nxtdgt;
	double xpnt;
	
/*  Make this a macro or subroutine?  */

	bf_null( &ten);
	ten.expnt = 4;
	ten.mntsa.e[MS_MNTSA] = 0x50000000;

/*  first character out is always 'E'  */

	digit = outstring;
	*digit++ = 'E';

/*  if number doesn't need an exponent, skip exponent phase  */

	exponent = numbr->expnt;
	bf_copy( numbr, &fraction);
	if( (exponent < -3) || (exponent > 3))
	{
		if ( exponent < -3)
		{
			xpnt = exponent;
			xpnt *= -logbase10of2;
			*digit++ = '-';
			signflag = -1;
		}
		else   /*  positive exponent  */
		{
			xpnt = exponent;
			xpnt *= logbase10of2;
			*digit++ = '+';
			signflag = 0;
		}
		exp10 = xpnt + 0.5;

/*  find msb in exp10 by brute force  */

		lastbit = 0x40000000;
		while  (! (lastbit & exp10) && lastbit) lastbit >>= 1;
		bf_one( &scale);

/*  compute scale = 10^exp10, in binary  */
 
	 	for (mask=1;  mask < lastbit; mask <<= 1 )
	 	{
	 		if (mask & exp10) bf_multiply( &scale, &ten, &scale);
	 		bf_multiply( &ten, &ten, &ten);    // square power of 10
	 	}
	 	bf_multiply( &scale, &ten, &scale);
	 	if ( signflag) bf_multiply( &scale, numbr, &fraction);
	 	else bf_divide( numbr, &scale, &fraction);
	 
 /*  now we have number converted to fraction in range .125+ to 8-
 	We also know the binary representation of the decimal
 	exponent, so let's output the decimal representation.
 */
	 	for( i=0; i<9; i++) digit[i] = '0';
	 	i = 8;
	 	while( i)
	 	{
	 		digit[i] |= exp10 % 10;
	 		exp10 /= 10;
	 		if( !exp10) break;
	 		i--;
	 	}
	 	digit += 9;
	 } /*  end of exponent conversion  */
	 
/*  Make this a macro or subroutine?  */

	bf_null( &ten);
	ten.expnt = 4;
	ten.mntsa.e[MS_MNTSA] = 0x50000000;

/*	 now output fractional part.
	Make positive, then spit out 1 digit of fraction at a time.
*/
	*digit++ = ' ';
	if( fraction.mntsa.e[MS_MNTSA] & SIGN_BIT)
	{
		bf_negate( &fraction);
		*digit++ = '-';
	}
	else *digit++ = '+';
	
/*  spit out first digit then decimal point  */

	*digit++ = bf_digitof( &fraction);
	*digit++ = '.';

/*  multiply result by 10 and output next decimal digit.  
	go till number is gone.
*/
	i=0;
	while( !bf_iszero( &fraction) & i<80)
	{
		bf_multiply( &fraction, &ten, &fraction);
		*digit++ = bf_digitof( &fraction);
		i++;
	}
	*digit++ = 0;
}

/*  this subroutine takes in a FLOAT ( > 0 ) and computes
	the integer part, chopping it off and returning a char.
	Output is normalized back to a fraction.
*/

char bf_digitof( FLOAT *frac)
{
	char num;
	unsigned long mask;
	long xp2 = 31 - frac->expnt;
	
	if ( xp2 > 30 ) return '0';
	mask = ~0 << xp2;
	num = (( mask & frac->mntsa.e[MS_MNTSA]) >> xp2);
	frac->mntsa.e[MS_MNTSA] &= ~mask;
	bf_normal( frac);
	return num | '0';
}

/*  general test routine.  Add more as you need 'em.
	Make inline too.
	
	check if FLOAT is zero.  returns 1 if so, 0 if nonzero
	all normalized numbers represented.  No soft cut off.
	With 32 bit exponent, it's insane to have soft 0.
*/

int bf_iszero( FLOAT *x)
{
	if( x->mntsa.e[MS_MNTSA] ) return 0;
	return 1;
}		

/*  print a text string and float.  For debugging mostly.  */

void bf_printfloat( char *string, FLOAT *numbr)
{
	char	numstrng[128];
	
	printf("%s\n", string);
	bf_float_to_ascii( numbr, numstrng);
	printf("%s\n", numstrng);
}

/*  print complex number to stdout  */

void bf_print_cmplx( char *string, COMPLEX *num)
{
	printf("%s\n", string);
	bf_printfloat("  real part is", &num->real);
	bf_printfloat("  imaginary part is", &num->imag);
	printf("\n");
}
