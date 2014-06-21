
#include "bigfloat.h"

typedef struct
{
	ELEMENT		degree;		/*  degree of polynomial  */
	ELEMENT		memdex;		/*  index into ram_block array  */
} MULTIPOLY;

typedef struct
{
	ELEMENT	flag;		/*  marks block as free or used  for crunch time */
	ELEMENT	up;		/*  next index to ram_block in address increasing order */
	ELEMENT	down;	/*  previous index to ram_block in adress decreasing order */
	ELEMENT	size;		/*  number of FIELD2N blocks allocated for this ram_block */
	FLOAT	*start;	/*  address in pool_mem of this data block  */
}  RAMDATA;

#define  Address(v)		ram_block[(v).memdex].start
#define AddressOf(v)	ram_block[(v)->memdex].start

/*  MAXCHEB determines number of Chebyshev polynomials to use and
	amounts to maximum degree of exp and cos expansions.  
	Make even for cosine.
*/

#define	MAXCHEB		54

/*
int mbf_multi_div( MULTIPOLY Top, MULTIPOLY Bottom, 
				MULTIPOLY *Quotient, MULTIPOLY *Remainder);
int	mbf_multi_gcd( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *gcd);
int mbf_gen_division_polynomial( MULTIPOLY *f, int length, CURVE curv);
int mbf_gen_xmodf( MULTIPOLY f, MULTIPOLY *xmod);
int mbf_xqmodf( MULTIPOLY f, MULTIPOLY *xmodf, MULTIPOLY *xq);
*/

/* Function prototypes produced by cproto(1) */

/* multi_float.c */
int mbf_power_add(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
int mbf_power_mul(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
ELEMENT mbf_zero_check(FLOAT *);
int mbf_power_div(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
int mbf_multi_add(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
int mbf_multi_sub(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
int mbf_multi_mul(MULTIPOLY, MULTIPOLY, MULTIPOLY *);
/* multi_mem.c */
void mbf_init_ram_space(void);
int mbf_get_space(MULTIPOLY *);
void mbf_crunch_ram(void);
void mbf_free_space(MULTIPOLY *);
void mbf_multi_copy(ELEMENT, FLOAT *, FLOAT *);
int mbf_multi_dup(MULTIPOLY, MULTIPOLY *);
/* multipoly.h */
