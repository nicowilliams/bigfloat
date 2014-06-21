/*** bigint.h ***/

#include "field2n.h"

/*  The following are used by the multiprecision integer package.
	This really is very crude.  See J. W. Crenshaw "Programmers
	Toolbox", Embedded Systems Programming Dec. 1996 for why you
	want to do this in assembler.
*/

#define	HALFSIZE	(WORDSIZE/2)
#define	HIMASK		(-1L<<HALFSIZE)
#define LOMASK		(~HIMASK)
#define CARRY		(1L<<HALFSIZE)
#define MSB_HW		(CARRY>>1)
#define	INTMAX		(4*MAXLONG-1)
#define MAXSTRING	(MAXLONG*WORDSIZE/3)

#define	INTLOOP(i)	for(i=INTMAX;i>=0;i--)

typedef struct {
	ELEMENT		hw[4*MAXLONG];
}  BIGINT;

/* Function prototypes produced by cproto(1) */

/* bigint.c */
void bi_int_null(BIGINT *);
void bi_int_copy(BIGINT *, BIGINT *);
void bi_field_to_int(FIELD2N *, BIGINT *);
void bi_int_to_field(BIGINT *, FIELD2N *);
void bi_int_neg(BIGINT *);
void bi_int_add(BIGINT *, BIGINT *, BIGINT *);
void bi_int_sub(BIGINT *, BIGINT *, BIGINT *);
void bi_int_mul(BIGINT *, BIGINT *, BIGINT *);
void bi_int_div(BIGINT *, BIGINT *, BIGINT *, BIGINT *);
void bi_ascii_to_bigint(char *, BIGINT *);
void bi_bigint_to_ascii(BIGINT *, char *);
/* int_functions.c */
void bi_int_div2(BIGINT *);
void bi_int_gcd(BIGINT *, BIGINT *, BIGINT *);
void bi_mod_exp(BIGINT *, BIGINT *, BIGINT *, BIGINT *);
void bi_mod_inv(BIGINT *, BIGINT *, BIGINT *);
