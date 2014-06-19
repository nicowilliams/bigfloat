/********************************************************
*														*
*		Build up my own memory control system.  Problem is with	*
*	bugs in chunk_alloc, and I'd rather have the ram be process		*
*	controlled rather than on the heap anyway.  I'll let the experts	*
*	fix egcs and Linux.										*
*														*
*					Author = Mike Rosing					*
*					  date  = Sept. 5, 1999					*
*     Changed to work with floats instead of FIELD2N 4/7/00			*
*														*
********************************************************/

#include <stdio.h>
#include "bigfloat.h"
#include "multipoly.h"

#define block_limit   100000

RAMDATA  ram_block[block_limit];

#define POOL_SIZE	1024*1024
FLOAT	 pool_mem[POOL_SIZE];

#define IndexListLength	block_limit

ELEMENT ramIndexList[ IndexListLength];
ELEMENT StartRamIndex, StopRamIndex, IndexListInUse;

/*  chunk_alloc dies so rather than debug OS problems I've created my own
	memory management.  Designed to handle multivariate polynomial math,
	this memory manager  uses a globally allocated pool, a linked list of
	pointers into that poool, and an "index list" used to track unused pointers.
	Space is allocated from the pool for each MULTIPOLY requested and an
	inex into the ram_block array is returned.  This index does not change
	for the life of the variable.  This allows me to crunch ram and physically
	move data without the mathematics being aware of it.
*/

/*  initialize ram allocation system  */

void init_ram_space()
{
	StartRamIndex = 0;
	StopRamIndex = 0;
	IndexListInUse = 0;
	ramIndexList[0] = 1;
	ram_block[0].up = 0;
	ram_block[0].down = 0;
	ram_block[0].flag = 0;		/*  this is permenent  */
	ram_block[0].size = POOL_SIZE;
	ram_block[0].start = pool_mem;
}

/*  ram_block[0] is special.  It always points to the last free ram block
	available and its up pointer ill always point to the first memory
	block and its down pointer will reference the next to last block.
	The purpose of the up and down indecies is to keep ram allocation
	order indepenedent of ram_block index.
*/

/*  get a chunk of ram from the pool for a MULTIPOLY. Enter with degree
	value preset.  Returns with index into ram_block which contains a 
	pointer to ram and value 1 if space available.  Returns garbage index
	and value 0 if not enough space available for the MULTIPOLY.
	If no space immediately available off end of pool, whole thing is 
	crunched once.  If still no space, you get a zero.
	Purpose of ramIndexList is to reuse ram_block areas removed during
	a crunch.
*/

int get_space( MULTIPOLY *newpoly)
{
	ELEMENT  need, ramindex;
	
	need = (newpoly->degree + 1)*(sizeof(FLOAT)/sizeof(ELEMENT));
	if( ram_block[0].size < need)
	{
		crunch_ram();
		if ( ram_block[0].size < need)
		{
			newpoly->memdex = -1;
			return 0;
		}
	}
	if( IndexListInUse)
	{
		ramindex = ramIndexList[ StartRamIndex];
		StartRamIndex++;
		if( StartRamIndex >= StopRamIndex) IndexListInUse = 0;
	}
	else
	{
		ramindex = ramIndexList[StopRamIndex];
		ramIndexList[StopRamIndex]++;
	}
	
	if( ramindex >= block_limit) return 0;  /*  we ran of index space!  */
	newpoly->memdex = ramindex;
	ram_block[ramindex].start = ram_block[0].start;
	ram_block[0].start += need;
	ram_block[0].size -= need;
	ram_block[ramindex].size = need;
	ram_block[ramindex].up = 0;
	ram_block[ramindex].down  = ram_block[0].down;
	ram_block[ram_block[0].down].up = ramindex;
	ram_block[0].down = ramindex;
	ram_block[ramindex].flag = 1;
	return 1;
}

/*  crunch memory:  Squeeze all the free space out of pool_mem.  I do this
	in a crude way but the purpose is to isolate the math code from memory 
	management problems. 
	ram_block[] contains a pointer to a block of FIELD2N (called .start), the
	number of FIELD2N's (.size) and .up, .down number containing the pointer
	of the next and previous memory blocks respectively.
	The ramIndexList[] tells me which indecies into ram_block have been
	removed after two free spaces have been concatenated.  It's a list of
	NULL pointers in ram_block.
*/

void  crunch_ram()
{
	ELEMENT	index, up, down, up2, size;
	FLOAT	*from, *to;
	
/*  if this is not first time here, reset ramIndexList pointers.
	use cyclic buffer to avoid this step.
*/
	if( IndexListInUse)
	{
		index = 0;
		while( StartRamIndex <= StopRamIndex)
		{
			ramIndexList[index] = ramIndexList[StartRamIndex];
			StartRamIndex++;
			index++;
		}
		StartRamIndex = 0;
		StopRamIndex = index - 1;
	}

/*  next concatenate any free space regions  */

	index = 1;
	while( index)
	{
		if( ram_block[index].flag)  index = ram_block[index].up;  /*  look for a free block*/
		else
		{
			up = ram_block[index].up;
			if( !ram_block[up].flag)  /* two free blocks in a row at this point */
			{
				ramIndexList[StopRamIndex+1] = ramIndexList[StopRamIndex];
				ramIndexList[StopRamIndex] = index;
				IndexListInUse = 1;
				StopRamIndex++;
				ram_block[up].size += ram_block[index].size;
				ram_block[up].start = ram_block[index].start;
				down = ram_block[index].down;
				ram_block[down].up = up;
				ram_block[up].down = down;
				index = up;
			}
			else
			{
/*  move a used block down in address over a free block and change pointers around */

				down = ram_block[index].down;
				up2 = ram_block[up].up;
				from = ram_block[up].start;
				to = ram_block[index].start;
				size = ram_block[up].size;
				multi_copy( size, from, to);
				ram_block[up].start = to;
				ram_block[index].start += size;
				ram_block[down].up = up;
				ram_block[up].down = down;
				ram_block[index].up = up2;
				ram_block[index].down = up;
				ram_block[up2].down = index;
			}
		}/*  ram_block[index].flag set/clear */
	}/*  index traversal  */
}/* crunch_ram  */

/*  last and simplest routine.
	Mark a block of ram as free space.  
	Nothing happens until crunch time.
*/
void free_space( MULTIPOLY *x)
{
	ram_block[x->memdex].flag = 0;
}

/*  copy a multivariate polynomial to another place.
	This is a physical copy, both regions are ASSUMED to actually be able
	to hold the data !!
Note: exceptionally useful when multiplying by x^k.
*/

void multi_copy( ELEMENT length, FLOAT *source, FLOAT *destination)
{
	INDEX  i;
	
	if (!length) return;
	for( i=0; i<length; i++)
	{
		copy( source, destination);
		source++;
		destination++;
	}
}

/*  The following seems like a better interface to the higher level math.
	Leave what works alone for now, but this may be able to eliminate
	the above routine completely.
*/

int multi_dup( MULTIPOLY from, MULTIPOLY *to)
{
	INDEX	i;
	FLOAT	*f, *t;
	
	to->degree = from.degree;
	if( !get_space( to)) return 0;
	f = Address( from);
	t = AddressOf( to);
	for( i=0; i<=from.degree; i++)
	{
		copy( f, t);
		f++;
		t++;
	}
	return 1;
}
/*main()
{
	MULTIPOLY	test1, test2, test3;
	INDEX		i;
	
	printf("start\n");

	init_ram_space();
	test1.degree = 25;
	get_space( &test1);
	
	printf("Start address is %x, start free space is %x\n", pool_mem, 
			ram_block[0].start);
	printf("size of space is %d, and it is marked with flag = %d\n",
			ram_block[test1.memdex].size, ram_block[test1.memdex].flag);
	printf("and the pointer to the ram block for test1 is %x\n",
			ram_block[test1.memdex].start);
	test2.degree = 63;
	get_space(&test2);
	test3.degree = 127;
	get_space(&test3);
	free_space(&test2);
	crunch_ram();
	printf("index  flag  up   down   size   address\n");
	for(i=0; i<5; i++)
		printf("  %d     %d    %d     %d      %d         %x\n", i, ram_block[i].flag,
			ram_block[i].up, ram_block[i].down, 
			ram_block[i].size, ram_block[i].start);
	printf("finish\n");
}
	*/