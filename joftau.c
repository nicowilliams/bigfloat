/*  compute all the coefficients of j-invariant q expansion.
	These are integers, and j is an "elliptic modular function".
	Purpose of this code is to see what it looks like, somehow.
*/

#include <stdio.h>
#include <stdlib.h>
#include "bigfloat.h"
#include "multipoly.h"

extern RAMDATA ram_block[];
extern FLOAT P2;

#define gridsize 512

/*  compute first two terms of j(tau) sequence
	j ~ 1/q + 744

	returns q = exp(2*i*PI*tau) as well as j.
*/
void bf_firstj( COMPLEX *tau, COMPLEX *q, COMPLEX *j)
{
	COMPLEX top, ipi;
	
	bf_null_cmplx( &ipi);
	bf_copy( &P2, &ipi.imag);
	ipi.imag.expnt += 2;
	bf_multiply_cmplx( &ipi, tau, q);
	bf_exp_cmplx( q, q);
	bf_null_cmplx( &top);
	bf_one( &top.real);
	bf_divide_cmplx( &top, q, j);
	bf_int_to_float( 744, &top.real);
	bf_add_cmplx( &top, j, j);
}	

main()
{
	FLOAT o1, dcubed, n, *offset;
	INDEX 	i, j, k, limit;
	MULTIPOLY 	sigma3;
	MULTIPOLY	q24, tau1, tau2;
	MULTIPOLY	joftop, jofbot, joftau;
	FLOAT	*coef, *tsubj, *tnew, *prevc;
	FLOAT	bctop, bcbottom;
	int		shift, maxstore;
	MULTIPOLY	cheb[100];
	
	struct
	{
		ELEMENT x, y;
		COMPLEX start, jt;
	} datablock;
	FILE *svplot;
	COMPLEX tau, jtau, arc[512], q, qn;
	FLOAT theta, dtheta;
	COMPLEX	temp;
	
	limit = 50;
	maxstore = limit+5;
		
	bf_init_ram_space();
	bf_init_float();

/*  create table of sigma_3(n) (sum of cube of all factors
	of n).
*/			
	sigma3.degree = limit;
	if( !bf_get_space( &sigma3)) 
	{
		printf("no space for that much data\n");
		exit(0);
	}
	bf_one( &o1);
	bf_null( &n);
	for( i=1; i<limit; i++)
	{
		bf_add( &o1,  &n, &n);
		bf_multiply( &n, &n, &dcubed);
		bf_multiply( &n, &dcubed, &dcubed);
		for( j=i; j<limit; j+=i)
		{
			offset = Address( sigma3) + j;
			bf_add( &dcubed, offset, offset);
		}
	}
	/*  save to disk here if necessary */
/*	for( i=0; i<limit; i++)
	{
		tsubj = Address( sigma3) + i;
		printf("i=%d\n", i);
		printfloat("sigma3(i) = ", tsubj);
	}

/*  compute Ramanujan's tau function to some ridiculous degree.
	First step is compute coefficients of (1-x)^24 using
	binomial expansion.
*/

	q24.degree = 24;
	if( !bf_get_space( &q24))
	{
		printf("no room for binomial coefficients?\n");
		exit(0);
	}
	coef = Address(q24);
	bf_int_to_float( 1, coef);
	bf_int_to_float( 24, &bctop);
	bf_int_to_float( 1, &bcbottom);
	for( i=1; i<=24; i++)
	{
		prevc = coef; 
		coef = Address(q24) + i;
		bf_multiply( prevc, &bctop, coef);
		bf_divide( coef, &bcbottom, coef);
		bf_negate( coef);
		bf_round( coef, coef);
/*		printf("i= %d\n", i);
		printfloat("  coef =", coef);

/*  decrement top and increment bottom for next term
	in binomial expansion */

		bf_subtract( &bctop, &o1, &bctop);
		bf_add( &o1, &bcbottom, &bcbottom);
	}

/*  now compute Ramanujan's tau.
	Keep two versions, a source and a destination.
	work back and forth multiplying source by
	binomial coefficients and sum to power shifted
	location.  Seeking coefficients for each power 
	of q:  q*product( 1- q^n)^24  n = 1 to infinity.
*/

	tau1.degree = maxstore;
	tau2.degree = maxstore;
	if( !bf_get_space( &tau1))
	{
		printf("can't allocate first tau block.\n");
		exit(0);
	}
	if( !bf_get_space( &tau2))
	{
		printf("can't allocate second tau block.\n");
		exit(0);
	}
	coef = Address( q24);
	tnew = Address( tau1);
	bf_multi_copy( 25, coef, tnew);
	k = 2;
	while( k<maxstore )  /*  loop over products  */
	{
/*  multiply by 1, first coefficient of each product term  */

		prevc = Address( tau1);
		tnew = Address( tau2);
		bf_multi_copy( maxstore, prevc, tnew);

/*  for each coefficient in 24 term product of next term,
	multiply by every term in previous product and sum
	to shifted location.  */
	
		for( i=1; i <= 24; i++)
		{
			coef = Address( q24) + i;
			shift = k*i;
			if( shift > maxstore) continue;
			for( j=0; j<= k*24; j++)
			{
				if( j + shift > maxstore) continue;
				tsubj = Address( tau1) + j;
				tnew = Address( tau2) + j + shift;
				bf_multiply( tsubj, coef, &bctop);
				bf_add( &bctop, tnew, tnew);
				bf_round( tnew, tnew);
			}
		}
		k++;

/*  now flip things over and go back to tau 1  */

		if( k>maxstore) break;
		prevc = Address( tau2);
		tnew = Address( tau1);
		bf_multi_copy( maxstore, prevc, tnew);
		for( i=1; i<=24; i++)
		{
			coef = Address(q24) + i;
			shift = k*i;
			if( shift > maxstore) continue;
			for( j=0; j<= k*24; j++)
			{
				if( j + shift > maxstore) continue;
				tsubj = Address( tau2) + j;
				tnew = Address( tau1) + j + shift;
				bf_multiply( tsubj, coef, &bctop);
				bf_add( &bctop, tnew, tnew);
				bf_round( tnew, tnew);
			}
		}
		k++;
	}

/*  compute top polynomial of joftau 
	 (1 + 240*sum(sigma3(n)*q^n))^3 */

	bf_multi_dup( sigma3, &joftop);
	coef = Address( joftop);
	bf_int_to_float( 1, coef);
	bf_int_to_float( 240, &bctop);
	for( i=1; i<limit; i++)
	{
		coef = Address( joftop) + i;
		bf_multiply( &bctop, coef, coef);
	}
	bf_power_mul( joftop, joftop, &joftau);  
	bf_power_mul( joftop, joftau, &joftop);

/*  finally compute joftau coefficients */

	bf_power_div( joftop, tau1, &joftau);
/*	for( i=0; i<limit; i++)
	{
		tsubj = Address( joftau) + i;
		printf("i=%d\n", i);
		printfloat("joftau(i) = ", tsubj);
	}
*/
/*  region F is defined as | Re(tau) | < 1/2 and || tau || > 1.
	For each point tau in F, find j(tau).
	save binary data to disk.  Format is (x, y) and (start,
	end) where x and y are integer indexes, and start, end
	are complex points.  Initial start is an arc along
	tau = 1.
*/
	svplot = fopen("joftau.complex", "wb");
	if( !svplot)
	{
		printf( "can't create output file\n");
		exit(0);
	}

/*  create an arc along bottom of F  */

	bf_copy( &P2, &theta);
	theta.expnt += 2;  // create 2PI/3
	bf_int_to_float( 3, &n);
	bf_divide( &theta, &n, &theta);
	bf_copy( &theta, &dtheta);
	dtheta.expnt--;		// create PI/3/511
	bf_int_to_float( gridsize-1, &n);
	bf_divide( &dtheta, &n, &dtheta);
	for( i=0; i<gridsize; i++)
	{
		bf_cosine( &theta, &arc[i].real);
		bf_sine( &theta, &arc[i].imag);
		bf_subtract( &theta, &dtheta, &theta);
	}

/*  Move arc up F, and find j(tau) for each point  */

	bf_int_to_float( 10, &bctop);
	bf_int_to_float( gridsize, &n);
	bf_divide( &bctop, &n, &bctop);  // upper limit of F
	for( i=0; i<gridsize; i++)
	{
		printf("i= %d\n", i);
		bf_int_to_float(i, &n);
		bf_multiply( &bctop, &n, &tau.imag);
		bf_null( &tau.real);
		datablock.y = i;
		for( j=0; j<gridsize; j++)
		{
			datablock.x = j;

			bf_add_cmplx( &tau, &arc[j], &datablock.start);
//			print_cmplx("data block start", &datablock.start);
			
/*  compute j(tau) for this point.  Note power of q = index - 1 */

			bf_firstj( &datablock.start, &q, &jtau);
//			print_cmplx("q = exp(2 i PI tau)", &q);
//			print_cmplx("first terms", &jtau);
			bf_copy_cmplx( &q, &qn);
			for( k=2; k<limit; k++)
			{
				offset = Address( joftau) + k;
				bf_multiply( offset, &qn.real, &temp.real);
				bf_multiply( offset, &qn.imag, &temp.imag);
				bf_add_cmplx( &temp, &jtau, &jtau);
				bf_multiply_cmplx( &q, &qn, &qn);
			}
			
/*  save data point to disk  */

//			print_cmplx("j(tau) = ", &jtau);
			bf_copy_cmplx( &jtau, &datablock.jt);
			if( fwrite(  &datablock, sizeof( datablock), 1, svplot) <= 0)
				printf("can't write to disk\n");
		}
	}
	fclose( svplot);
	printf("all done!\n\n");
}			


