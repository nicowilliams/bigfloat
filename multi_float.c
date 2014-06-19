/********************************************************************************
*																				*
*		The purpose of this file is basic manipulation of power series.   These are infinite		*
*    series as opposed to multivariate polynomials of finite degree.  I use the same memory			*
*    allocation and all degrees start at zero, so offset to a block corresponds to coefficient of x to	*
*    that power.																		*
*    Algorithms from Knuth page 506 (section 4.7)											*
*																				*
*									Author = Mike Rosing							*
*									  date  =  April 11, 2000							*
*																				*
********************************************************************************/

#include "bigfloat.h"
#include "multipoly.h"

extern RAMDATA ram_block[];

/*  add two power series.  C = A + B
	Creates new output space.  takes care of memory management,
	but ASSUMES A and B were previously allocated.  
	Returns 1 on success, 0 on failure to create C.
*/

int power_add( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C)
{
	ELEMENT		i, j, big, small;
	MULTIPOLY	result;
	FLOAT		*Result, *Aptr, *Bptr;
	
/*  find out which polynomial is larger and allocate that much space  */

	if (A.degree > B.degree)
	{
		small = B.degree;
		big = A.degree;
	}
	else
	{
		small = A.degree;
		big = B.degree;
	}
	result.degree = big;
	if (!get_space( &result)) return 0;

/*  now add the shorter length amounts together  */

	Result = Address(result);
	Aptr = Address( A);
	Bptr = Address( B);
	for( i=0; i<= small; i++)
	{
		add( Aptr, Bptr, Result);
		Result++;
		Aptr++;
		Bptr++;
	}

	if( A.degree == big)
		multi_copy( big - small, Aptr, Result);
	else
		multi_copy( big - small, Bptr, Result);
				
/*  take care of memory management  */

	if ( A.memdex == C->memdex ) free_space( &A);
	if ( B.memdex == C->memdex ) free_space( &B);
	C->degree = result.degree;
	C->memdex = result.memdex;
	return 1;
}

/*  Multiply two power series. Uses Knuth's construction
	from page 506. of Semi Numerical Algorithms (sect. 4.7)
	computes C = A*B 
	returns 0 if C can't be allocated.
*/

int power_mul( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C)
{
	ELEMENT		i, k;
	MULTIPOLY	result;
	FLOAT		*Aptr, *Bptr, *Result;
	FLOAT		temp;
	
/*  create space for result */

	if (A.degree > B.degree)
		result.degree = A.degree;
	else
		result.degree = B.degree;
	if ( !get_space( &result) ) return 0;
	for( i=0; i<=result.degree; i++)
	{
		Result = Address(result) + i;
		null( Result);
		for( k=0; k<=i; k++)
		{
			if( k <= A.degree) Aptr = Address( A) + k;
			else continue;
			if( i-k < B.degree) Bptr = Address( B) + i - k;
			else continue;
			multiply( Aptr, Bptr, &temp);
			add( &temp, Result, Result);
		}
	}
	
/*  take care of memory management  */

	if (A.memdex == C->memdex) free_space( &A);
	if (B.memdex == C->memdex) free_space( &B);
	C->memdex = result.memdex;
	C->degree = result.degree;
	return 1;
}

/*  Subroutine to test if a FLOAT value is 0.
	returns 0 if all bits in field are 0, first non-zero ELEMENT
	otherwise.
*/

ELEMENT zero_check( FLOAT *z)
{
	INDEX	i;
	
	if( z->expnt) return z->expnt;
	for( i= MS_MNTSA; i >= 0; i++)
		if( z->mntsa.e[i] ) return z->mntsa.e[i];
	return 0;
}
	
/*  Power series division.  Resulting degree is minimum
	degree of two source polynomials.
	
		C = A/B
	returns 1 if ok, 0 if no space for C
*/

int power_div( MULTIPOLY A, MULTIPOLY B,  MULTIPOLY *C)
{
	ELEMENT		n, k;
	FLOAT		temp, *Result, *Aptr, *Bptr;
	MULTIPOLY	result;

/*  check to see if attempting to divide by 0  */

	Bptr = Address( B);
	if (!zero_check(Bptr)) return (0);

/*  create space for result  */

	if( A.degree < B.degree) result.degree = A.degree;
	else result.degree = B.degree;
	if( !get_space(&result)) return 0;

/*  do first term  */

	Aptr = Address(A);
	Bptr = Address( B);
	Result = Address( result);
	divide( Aptr, Bptr, Result);

/*  do rest of terms  */

	for( n=1; n<= result.degree; n++)
	{
		Result = Address( result) + n;
		for( k=0; k<n; k++)
		{
			Aptr = Address( result) +k;
			Bptr = Address( B) + n - k;
			multiply( Aptr, Bptr, &temp);
			add( &temp, Result, Result);
		}
		Aptr = Address( A) + n;
		subtract( Aptr, Result, Result);
		Bptr = Address( B);
		divide( Result, Bptr, Result);
	}
		
/*  take care of memory management  */

	if (A.memdex == C->memdex) free_space( &A);
	if (B.memdex == C->memdex) free_space( &B);
	C->memdex = result.memdex;
	C->degree = result.degree;
	return 1;
}

/*  NOTE: the following routines work with *polynomials*, which are similar
	to but not the same as *power series* above.  The main difference is
	that in a power series we don't care about higher order terms and let
	them fall off the end because they are too small.  With polynomials
	we keep every term and assume things get larger.
*/

/*  add two multivariate polynomials.  C = A + B
	Creates new output space.  takes care of memory management,
	but ASSUMES A and B were previously allocated.  
	Returns 1 on success, 0 on failure to create C.
*/

int multi_add( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C)
{
	ELEMENT		i, j;
	MULTIPOLY	shortpoly, result, longpoly;
	FLOAT		*Result, *Aptr, *Bptr;
	
/*  find out which polynomial is larger and allocate that much space  */

	if (A.degree > B.degree)
	{
		shortpoly.degree = B.degree;
		shortpoly.memdex = B.memdex;
		longpoly.degree = A.degree;
		longpoly.memdex = A.memdex;
	}
	else
	{
		shortpoly.degree = A.degree;
		shortpoly.memdex = A.memdex;
		longpoly.degree = B.degree;
		longpoly.memdex = B.memdex;
	}
	result.degree = longpoly.degree;
	if (!get_space( &result)) return 0;

/*  now add the shorter length amounts together  */

	Result = Address(result);
	Aptr = Address( A);
	Bptr = Address( B);
	for( i=0; i<= shortpoly.degree; i++)
	{
		add(  Aptr, Bptr, Result);
		Result++;
		Aptr++;
		Bptr++;
	}

/*  and copy the rest  */

	Result = Address(result) + shortpoly.degree + 1;
	Aptr = Address(longpoly) + shortpoly.degree + 1;
	multi_copy( longpoly.degree - shortpoly.degree, Aptr, Result);
	Result = Address(result);
	while( result.degree)
		if( !zero_check( &Result[result.degree])) result.degree--;
		else break;
		
/*  take care of memory management  */

	if ( A.memdex == C->memdex ) free_space( &A);
	if ( B.memdex == C->memdex ) free_space( &B);
	C->degree = result.degree;
	C->memdex = result.memdex;
	return 1;
}

/*  subtract two multivariate polynomials.  C = A - B
	Negate all components of B and then call add.  
	Simplest way to deal with it.
*/

int multi_sub( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C)
{
	INDEX		i;
	FLOAT		*ptr;
	MULTIPOLY	temp;
	
	multi_dup( B, &temp);
	 
	for( i=0; i<=temp.degree; i++)
	{
		ptr = Address(temp) + i;
		negate( ptr);
	}
	multi_add( A, temp, C);
	free_space( &temp);	
}

/*  Multiply two multivariate polynomials. Uses Knuth's construction
	from page 399 of Semi Numerical Algorithms (sect. 4.6)
	computes C = A*B 
	returns 0 if C can't be allocated.
*/

int multi_mul( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C)
{
	ELEMENT		i, j, k;
	MULTIPOLY	shortmulti, longmulti, result;
	FLOAT		*Short, *Long, *Result;
	FLOAT		temp;
	
/*  create space for result using sum of degrees of source polynomials */

	result.degree = A.degree + B.degree;
	if ( !get_space( &result) ) return 0;
	if (A.degree > B.degree)
	{
		longmulti.memdex = A.memdex;
		longmulti.degree = A.degree;
		shortmulti.memdex = B.memdex;
		shortmulti.degree = B.degree;
	}
	else
	{
		longmulti.memdex = B.memdex;
		longmulti.degree = B.degree;
		shortmulti.memdex = A.memdex;
		shortmulti.degree = A.degree;
	}
	for( k=0; k<shortmulti.degree; k++)
	{
		Result = Address(result) + k;
		null( Result);
		for( i=0; i<=k; i++)
		{
			j = k - i;
			Short = Address( shortmulti) + i;
			Long = Address( longmulti) + j;
			multiply(  Short, Long,  &temp);
			Result = Address( result) + k;
			add( Result, &temp, Result);
		}
	}
	for( k=shortmulti.degree; k<longmulti.degree; k++)
	{
		Result = Address(result) + k;
		null( Result);
		for( i=0; i<=shortmulti.degree; i++)
		{
			j = k - i;
			Short = Address(shortmulti) + i;
			Long = Address( longmulti) + j;
			multiply( Short,  Long, &temp);
			Result = Address( result) + k;
			add( Result, &temp, Result);
		}
	}
	for( k=longmulti.degree; k<=result.degree; k++)
	{
		Result = Address( result) + k;
		null( Result);
		for( i = k-longmulti.degree; i <= shortmulti.degree;  i++)
		{
			j = k - i;
			Short = Address( shortmulti) + i;
			Long = Address( longmulti) + j;
			multiply( Short, Long, &temp);
			Result = Address( result) + k;
			add( Result, &temp, Result);
		}
	}
	Result = Address( result);
	while( result.degree)
		if( !zero_check( &Result[result.degree])) result.degree--;
		else break;
	
/*  take care of memory management  */

	if (A.memdex == C->memdex) free_space( &A);
	if (B.memdex == C->memdex) free_space( &B);
	C->memdex = result.memdex;
	C->degree = result.degree;
	return 1;
}
