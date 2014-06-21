/*  The idea behind bigfloat is to get enough precision to plot interesting
	regions of "fractal" or Hausdorf spaces.  Change precision with the
	MNTSA_SIZE define, this is the number of 32 bit longwords used to
	store the mantissa.  Should be even.

			Author = Mike Rosing
			  date  = feb. 14, 2000
*/

#define	MNTSA_SIZE	8
#define	D_SIZE		(MNTSA_SIZE/2)
#define	SIGN_BIT		0x80000000
#define	MSB			0x40000000
#define	MS_MNTSA	(MNTSA_SIZE-1)
#define	MS_MNTSAd	(D_SIZE-1)
#define	DIVISION_LOOPS	3
/*  division loops = log_2(MNTSA_SIZE), see division() */

typedef struct
{
	long		expnt;
	union
	{	/*  data stored in little endian order  */
		unsigned long long d[D_SIZE];
		unsigned long e[MNTSA_SIZE];
	} mntsa;
} FLOAT;

typedef struct
{
	FLOAT	real;
	FLOAT	imag;
} COMPLEX;

typedef	short int	INDEX;
typedef	unsigned long	ELEMENT;

#define	OPLOOPd(i)		for(i=0; i<D_SIZE ; i++)

#define	OPLOOP(i)			for(i=0; i<MNTSA_SIZE; i++)

/* Function prototypes produced by cproto(1) */

/* bigcomplex.c */
void bf_one(FLOAT *);
void bf_null_cmplx(COMPLEX *);
void bf_add_cmplx(COMPLEX *, COMPLEX *, COMPLEX *);
void bf_subtract_cmplx(COMPLEX *, COMPLEX *, COMPLEX *);
void bf_multiply_cmplx(COMPLEX *, COMPLEX *, COMPLEX *);
int bf_divide_cmplx(COMPLEX *, COMPLEX *, COMPLEX *);
int bf_intpwr_cmplx(COMPLEX *, int, COMPLEX *);
void bf_magnitude_cmplx(COMPLEX *, FLOAT *);
int bf_exp_cmplx(COMPLEX *, COMPLEX *);
/* bigfloat.c */
void bf_copy(FLOAT *, FLOAT *);
void bf_copy_cmplx(COMPLEX *, COMPLEX *);
void bf_null(FLOAT *);
void bf_negate(FLOAT *);
void bf_normal(FLOAT *);
int bf_compare(FLOAT *, FLOAT *);
void bf_add(FLOAT *, FLOAT *, FLOAT *);
void bf_subtract(FLOAT *, FLOAT *, FLOAT *);
void bf_round(FLOAT *, FLOAT *);
void bf_multiply(FLOAT *, FLOAT *, FLOAT *);
int bf_reciprical(FLOAT *, FLOAT *);
int bf_divide(FLOAT *, FLOAT *, FLOAT *);
void bf_init_float(void);
void bf_square_root(FLOAT *, FLOAT *);
void bf_int_to_float(int, FLOAT *);
/* bigfloat.h */
/* bigfunc.c */
void bf_calcpi(FLOAT *);
void bf_calcln2(FLOAT *);
int bf_intpwr(FLOAT *, int, FLOAT *);
void bf_bessel(int, int, FLOAT *, FLOAT *);
int bf_gen_chebyshev(MULTIPOLY *, int);
int bf_calc_2x_coef(MULTIPOLY *, int, MULTIPOLY *);
int bf_calc_cos_coef(MULTIPOLY *, int, MULTIPOLY *);
void bf_polyeval(MULTIPOLY, FLOAT *, FLOAT *);
void bf_twoexp(FLOAT *, FLOAT *);
void bf_corecos(FLOAT *, FLOAT *);
int bf_float_to_int(FLOAT *);
int bf_exp(FLOAT *, FLOAT *);
void bf_split(FLOAT *, FLOAT *, FLOAT *);
void bf_cosine(FLOAT *, FLOAT *);
void bf_sine(FLOAT *, FLOAT *);
/* bigio.c */
int bf_ascii_to_float(char *, FLOAT *);
void bf_float_to_ascii(FLOAT *, char *);
char bf_digitof(FLOAT *);
int bf_iszero(FLOAT *);
void bf_printfloat(char *, FLOAT *);
void bf_print_cmplx(char *, COMPLEX *);
