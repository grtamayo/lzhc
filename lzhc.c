/*
	---- An LZ77 Coding with Hashing ----
	
	Filename:               lzhc.C
	Decompression program:  lzhcx.c
	Date:                   (August 19-23, 2006) 2/2/2022

	Description:

		LZ77/LZSS with *hashing*.
		
		(1)  The table indices are determined by the hash value of the first
		*three* characters, spreading the list nodes into a "wider" table.
		Consequently, the size of the list in the table entries is minimized,
		allowing fewer list node search operations.

		(2)  An advantage of "single-character" hashing (using 256 linked-list
		heads) is that the first character is always determined, so the search
		can quickly start at the 2nd character. An n-character hash (of size
		>= 2) does not have this property and the search always begins at the
		first character. Interestingly, hashing is still much faster.

	May 17, 2008:
		:   lzlist.c --> lzhash.c
		:   allowed up to 16K sliding window
		:   used a faster LZ77 searching technique

	To compile:
		bcc32 -O2 -w lzhhc
		bcc32c -O2 -w lzhhc
		g++ -O3 lzhhc.c -o lzhhc
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "utypes.h"
#include "gtbitio2.c"
#include "lzhash.c"
#include "mtf.c"
#include "ucodes2.c"

/* the decompressor's must also equal these values. */
#define NUM_POS_BITS     17                 /* 12..21 tested working */
#define WIN_BUFSIZE      (1<<NUM_POS_BITS)
#define WIN_MASK         (WIN_BUFSIZE-1)
#define HASH_SHIFT       (NUM_POS_BITS-8)

/* accepted minimum size for a string is 3 */
#define NUM_LEN_BITS     (NUM_POS_BITS-1)
#define PAT_BUFSIZE      (1<<NUM_LEN_BITS)
#define MIN_LEN            4
#define EOF_CODE         256
#define MTF_CODES_N      257

#define HASH_BYTES_N     3
/* 3-byte hash */
#define hash_w(pos) \
	(((win_buf[(pos) & WIN_MASK] << HASH_SHIFT) ^ (win_buf[((pos)+1) & WIN_MASK]<<1) \
		^ (win_buf[((pos)+2) & WIN_MASK]<<4)) & WIN_MASK)
#define hash_p(pos) \
	(((pattern[(pos)] << HASH_SHIFT) ^ (pattern[((pos)+1) % PAT_BUFSIZE]<<1) \
		^ (pattern[((pos)+2) % PAT_BUFSIZE]<<4)) % PAT_BUFSIZE)

typedef struct {
	char algorithm[4];
} file_stamp;

typedef struct {
	unsigned int pos, len;
} dpos_t;

dpos_t dpos;
unsigned char win_buf[ WIN_BUFSIZE ];      /* the "sliding" window buffer. */
unsigned char pattern[ PAT_BUFSIZE ];      /* the "look-ahead" buffer. */
int win_cnt = 0, pat_cnt = 0, buf_cnt = 0;  /* some counters. */

void copyright( void );
dpos_t search ( uchar w[], uchar p[] );
void put_codes( dpos_t *dpos );
void send_EOF_code( void );

int main( int argc, char *argv[] )
{
	float ratio = 0.0;
	int i;
	file_stamp fstamp;

	if ( argc != 3 ) {
		fprintf(stderr, "\n Usage: lzhc infile outfile");
		copyright();
		return 0;
	}
	
	clock_t start_time = clock();
	
	if ( (gIN = fopen( argv[1], "rb" )) == NULL ) {
		fprintf(stderr, "\nError opening input file.");
		return 0;
	}

	if ( (pOUT = fopen( argv[2], "wb" )) == NULL ) {
		fprintf(stderr, "\nError opening output file.");
		return 0;
	}
	init_put_buffer();

	fprintf(stderr, "\n----[ A Lempel-Ziv (LZ77) Coding with Hashing Implementation ]----\n");
	fprintf(stderr, "\nWindow Buffer size used  = %15lu bytes", (ulong) WIN_BUFSIZE );
	fprintf(stderr, "\nLook-Ahead Buffer size   = %15lu bytes", (ulong) PAT_BUFSIZE );

	fprintf(stderr, "\n\nName of input file : %s", argv[1] );
	
	/* Write the FILE STAMP. */
	rewind( pOUT );
	strcpy( fstamp.algorithm, "LZ7" );
	fwrite( &fstamp, sizeof(file_stamp), 1, pOUT );
	nbytes_out = sizeof(file_stamp);
	
	/* start Compressing to output file. */
	fprintf(stderr, "\n\nCompressing...");

	/* initialize the table of pointers. */
	if ( !alloc_lzhash(WIN_BUFSIZE) ) goto halt_prog;

	/* set the sliding-window to all zero (0) values. */
	memset( win_buf, 0, WIN_BUFSIZE );

	/* initialize the search list. */
	for ( i = 0; i < WIN_BUFSIZE; i++ ) {
		lznext[i] = LZ_NULL;
		lzprev[i] = LZ_NULL;
		insert_lznode( hash_w(i), i );
	}

	/* make sure to rewind the input file */
	rewind(gIN);

 	/* fill the pattern buffer. */
	buf_cnt = fread( pattern, 1, PAT_BUFSIZE, gIN );

	/* initialize the input buffer. */
	init_get_buffer();
	nbytes_read = buf_cnt;
	
	alloc_mtf(MTF_CODES_N);
	
	/* compress */
	while ( buf_cnt > 0 ) {	/* look-ahead buffer not empty? */
		dpos = search( win_buf, pattern );

		/* encode prefix bits. */
		if ( dpos.len > MIN_LEN ) {  /* more than MIN_LEN match? */
			put_ONE();           /* yes, send a 1 bit. */
		}
		else if ( dpos.len == MIN_LEN ) { /* exactly MIN_LEN matching characters? */
			put_ZERO();          /* yes, send a 0 bit. */
			put_ONE();           /* and a 1 bit. */
		}
		else {	/* less than MIN_LEN matching characters. */
			put_ZERO();          /* send a 0 bit. */
			put_ZERO();          /* one more 0 bit to indicate a no match. */
		}

		/* encode window position or len codes. */
		put_codes( &dpos );
	}
	send_EOF_code();        /* Send EOF code */
	flush_put_buffer();
	fprintf(stderr, "complete.");
	
	/* get infile's size and get compression ratio. */
	nbytes_read = get_nbytes_read();
	
	fprintf(stderr, "\n\nName of output file: %s", argv[2] );
	fprintf(stderr, "\nLength of input file     = %15llu bytes", nbytes_read );
	fprintf(stderr, "\nLength of output file    = %15llu bytes", nbytes_out );
	
	ratio = (((float) nbytes_read - (float) nbytes_out) /
		(float) nbytes_read ) * (float) 100;
	fprintf(stderr, "\nCompression ratio:         %15.2f %% in %3.2f secs.\n", ratio, 
		(double) (clock()-start_time) / CLOCKS_PER_SEC );
	
	halt_prog:

	free_put_buffer();
	free_get_buffer();
	free_lzhash();
	if ( gIN ) fclose( gIN );
	if ( pOUT ) fclose( pOUT );

	return 0;
}

void copyright( void )
{
	fprintf(stderr, "\n\n Written by: Gerald Tamayo, 2006\n");
}

/*
This function searches the sliding window buffer for the largest
"string" stored in the pattern buffer.

The function uses an "array of pointers" to singly-linked
lists, which contain the various occurrences or "positions" of a
particular character in the sliding-window.
*/
dpos_t search( uchar w[], uchar p[] )
{
	register int i, j, k;
	dpos_t dpos = { 0, 0 };

	/* point to start of lzhash[ index ] */
	i = lzhash[ hash_p(pat_cnt) & WIN_MASK ];

	if ( buf_cnt > 1 ) while ( i != LZ_NULL ) {
		/* ---- FAST LZ77 SEARCH ---- */
		/*
		First, match the "context" string (the current longest string or
		the partial match) plus 1 "suffix" byte (ie., the first byte tested
		for a match) from right to left...

		The context length (i.e., the current match, dpos.len)) is a
		"skip count" as in Boyer-Moore search; thus our approach here
		does not need to prepare a "skip table" for the symbols.

		dpos.len points to the first suffix symbol; if it is a mismatch,
		the search can end immediately.  We verify the *context* string
		first since the *suffix* string may match completely but the context
		string may not, and hence would be a mismatch for the whole string.
		*/
		if ( (j=pat_cnt+dpos.len) >= PAT_BUFSIZE ) j -= PAT_BUFSIZE;
		k = dpos.len;
		do {
			if ( p[j] != w[ (i+k) & WIN_MASK ] ) {
				goto skip_search;  /* allows fast search. */
			}
			if ( j-- == 0 ) j = PAT_BUFSIZE-1;
		} while ( (--k) >= 0 );

		/* then match the rest of the "suffix" string from left to right. */
		if ( (j=pat_cnt+dpos.len+1) >= PAT_BUFSIZE ) j -= PAT_BUFSIZE;
		k = dpos.len+1;
		if ( k < buf_cnt ) do {
			if ( p[j] != w[ (i+k) & WIN_MASK ] ) {
				break;
			}
			/* rotate if necessary. */
			if ( (++j) == PAT_BUFSIZE ) j = 0;
		} while ( (++k) < buf_cnt );

		if ( k > dpos.len ) {
			/* if greater than previous length, record it. */
			dpos.pos = i;
			dpos.len = k;

			/* maximum match, end the search. */
			if ( k == buf_cnt ) break;
		}

		skip_search:

		/* point to next occurrence of this hash index. */
		i = lznext[i];
	}

	return dpos;
}

/*
Transmits a length/position pair of codes according
to the match length received.

When we receive a match length of 0, we quickly set the length
code to 1 (we have to "slide" through the window buffer at least
one character at a time).

Due to the algorithm, we only encode the match length if it is
greater than MIN_LEN. Next, a byte or a "position code" is
transmitted.

Then this function properly performs the "sliding" part by
copying the matched characters to the window buffer; note that
the linked list is also updated.

Finally, it "gets" characters from the input file according
to the number of matching characters.
*/
void put_codes( dpos_t *dpos )
{
	int i, j, k;

	if ( dpos->len > MIN_LEN ) {	/* encode len code only if > MIN_LEN. */
		k = dpos->len - (MIN_LEN+1);  /* map (MIN_LEN+1) to 0, (MIN_LEN+2) to 1, ... */
		#define MFOLD 2         /* m = 2 works for this type of encoding */
		i = k >> MFOLD;  /* "fold" the suffix string. (10/21/2008) */
		while ( i-- ) {         /*   encode only a part of the unary codes. */
			put_ONE();
		}
		put_nbits( (k % (1 << MFOLD)) << 1, MFOLD+1 );
	}
	
	/* encode position for match len >= MIN_LEN. */
	if ( dpos->len >= MIN_LEN ) {
		k = dpos->pos;
		put_nbits( k, NUM_POS_BITS );
	}
	else {
		dpos->len = 1;
		/* emit just the byte. */
		k = (unsigned char) pattern[pat_cnt];
		/* Implemented VL coding for better compression. (1/12/2010) */
		put_vlcode(mtf(k), 3);
	}
	
	/* ---- if its a match, then "slide" the buffer. ---- */
	/*
	the string at the "left" of win_buf[win_cnt] contains the
	byte at win_buf[win_cnt] so its hash value may change too.
	*/
	if ( (j=win_cnt-(HASH_BYTES_N-1)) < 0 ) {
		/* record the left-most string index (k). */
		j = k = WIN_BUFSIZE+j;
	}
	else k = j;

	/* remove the strings (i.e., positions) from the hash list. */
	for ( i = 0; i < (dpos->len+(HASH_BYTES_N-1)); i++, j++ ) {
		delete_lznode( hash_w(j), j & WIN_MASK );
	}

	j = win_cnt;
	for ( i = 0; i < dpos->len; i++, j++ ) {
		/* write the character to the window buffer. */
		*(win_buf +((j) & (WIN_MASK)) ) =
			*(pattern + ((pat_cnt+i) % PAT_BUFSIZE));
	}

	/* with the new characters, rehash at this position. */
	j = k;
	for ( i = 0; i < (dpos->len+(HASH_BYTES_N-1)); i++, j++ ) {
		insert_lznode( hash_w(j), j & WIN_MASK );
	}

	/* get dpos.len bytes */
	for ( i = 0; i < (dpos->len); i++ ){
		if( (k=gfgetc()) != EOF ) {
			*(pattern + ((pat_cnt+i) % PAT_BUFSIZE)) =
				(uchar) k;
		}
		else { break; }
	}

	/* update counters. */
	buf_cnt -= (dpos->len-i);
	win_cnt = (win_cnt+dpos->len) & WIN_MASK;
	if ( (pat_cnt = (pat_cnt+dpos->len)) >= PAT_BUFSIZE ) {
		pat_cnt -= PAT_BUFSIZE;
	}
}

void send_EOF_code( void )
{	
	put_ZERO(); put_ZERO();  /* send prefix bits */
	/* Implemented VL coding for better compression. (1/12/2010) */
	put_vlcode(mtf(EOF_CODE), 3);
}

