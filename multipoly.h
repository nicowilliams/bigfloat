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

int get_space(MULTIPOLY *newpoly);
void multi_copy( ELEMENT length, FLOAT *source, FLOAT *destination);
void crunch_ram();
int multi_dup( MULTIPOLY from,  MULTIPOLY *to);
ELEMENT zero_check( FLOAT *z);

/*int multi_add( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C);
int multi_mul( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *C);
ELEMENT zero_check( FIELD2N z);
int multi_div( MULTIPOLY Top, MULTIPOLY Bottom, 
				MULTIPOLY *Quotient, MULTIPOLY *Remainder);
int	multi_gcd( MULTIPOLY A, MULTIPOLY B, MULTIPOLY *gcd);
int gen_division_polynomial( MULTIPOLY *f, int length, CURVE curv);
int gen_xmodf( MULTIPOLY f, MULTIPOLY *xmod);
int xqmodf( MULTIPOLY f, MULTIPOLY *xmodf, MULTIPOLY *xq);
*/