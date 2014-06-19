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
