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

void bi_int_null();
void bi_int_copy();
void bi_field_to_int();
void bi_int_to_field();
void bi_int_neg();
void bi_int_add();
void bi_int_sub();
void bi_int_mul();
void bi_int_div();
void bi_ascii_to_bigint();
void bi_bigint_to_ascii();
void bi_int_gcd();
void bi_mod_exp();
void bi_mod_inv();
void bi_int_div2();
