/*
	Filename:		lzhcx.c	(this is the eXtractor), 2006/2022
	Encoder :		lzhc.c

	NOTES:

	-- Decompression in LZ77 is faster since you just have to
	extract the bytes from the window buffer using the pos and
	len variables.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "utypes.h"
#include "gtbitio2.c"
#include "ucodes2.c"
#include "mtf.c"

#define NUM_POS_BITS     17                 /* 12..21 tested working */
#define WIN_BUFSIZE      (1<<NUM_POS_BITS)
#define WIN_MASK         (WIN_BUFSIZE-1)

#define NUM_LEN_BITS     (NUM_POS_BITS-1)
#define PAT_BUFSIZE      (1<<NUM_LEN_BITS)
#define MIN_LEN            4
#define EOF_CODE         256
#define MTF_CODES_N      257

typedef struct {
	char algorithm[4];
} file_stamp;

typedef struct {
	uint pos, len;
} dpos_t;

dpos_t dpos;
uchar win_buf[ WIN_BUFSIZE ];
uchar pattern[ PAT_BUFSIZE ];
uint win_cnt = 0;

void copyright( void );

int main( int argc, char *argv[] )
{	
	int i, k;
	FILE *out;
	file_stamp fstamp;

	if ( argc != 3 ) {
		fprintf(stderr, "\n Usage: lzhcx infile outfile");
		copyright();
		return 0;
	}
	
	clock_t start_time = clock();
	
	if ( (gIN = fopen( argv[1], "rb" )) == NULL ) {
		fprintf(stderr, "\nError opening input file.");
		return 0;
	}
	fread( &fstamp, sizeof(file_stamp), 1, gIN );
	init_get_buffer();

	if ( (out = fopen( argv[2], "wb" )) == NULL ) {
		fprintf(stderr, "\nError opening output file.");
		goto halt_prog;
	}

	fprintf(stderr, "\n Name of input  file : %s", argv[1] );
	fprintf(stderr, "\n Name of output file : %s", argv[2] );

	fprintf(stderr, "\n\n Decompressing...");

	/* initialize sliding-window. */
	memset( win_buf, 0, WIN_BUFSIZE );
	
	alloc_mtf(MTF_CODES_N);
	
	/* start decompression. */
	while ( 1 ) {
		if ( get_bit() == 1 ){	/* match len > MIN_LEN? */
			for ( k = 0; get_bit(); k++ ) ;
			#define MFOLD    2
			k <<= MFOLD;
			k += get_nbits(MFOLD);
			
			/* map true values: 0 to (MIN_LEN+1), 1 to (MIN_LEN+2), and so on. */
			k += (MIN_LEN+1);
			dpos.len = k;

			/* get window position. */
			dpos.pos = get_nbits( NUM_POS_BITS );

			/* if its a match, then "slide" the window buffer. */
			i = dpos.len;
			while ( i-- ) {
				/* copy byte. */
				pattern[i] = win_buf[ (dpos.pos+i) & WIN_MASK ];
			}
			i = dpos.len;
			while ( i-- ) {
				win_buf[ (win_cnt+i) & WIN_MASK ] = pattern[ i ];
			}
			
			/* output string. */
			fwrite( pattern, dpos.len, 1, out );
			win_cnt = (win_cnt + dpos.len) & WIN_MASK;
		}
		else {
			switch ( get_bit() ){
			case 0:
			
			/* get VL-coded byte and output it. */
			k = get_mtf_c(get_vlcode(3));
			if ( k == EOF_CODE ) goto exit_loop;   /* EOF code! */
			fputc( k, out );
			
			/* slide the window buffer. */
			win_buf[ win_cnt ] = (unsigned char) k;
			if ( (++win_cnt) == WIN_BUFSIZE ) win_cnt = 0;
			break;

			case 1:

			/* get window position. */
			dpos.pos = get_nbits( NUM_POS_BITS );
			dpos.len = MIN_LEN;    /* we know len == MIN_LEN. */
			
			/* if its a match, then "slide" the window buffer. */
			i = dpos.len;
			while ( i-- ) {
				/* copy byte. */
				pattern[i] = win_buf[ (dpos.pos+i) & WIN_MASK ];
			}
			i = dpos.len;
			while ( i-- ) {
				win_buf[ (win_cnt+i) & WIN_MASK ] = pattern[ i ];
			}
			
			/* output string. */
			fwrite( pattern, dpos.len, 1, out );
			win_cnt = (win_cnt + dpos.len) & WIN_MASK;
			
			break;
			}	/* end switch */
		}
	}
	
	exit_loop:
	
	fprintf(stderr, "done, in %3.2f secs\n", 
        (double)(clock()-start_time) / CLOCKS_PER_SEC);

	halt_prog:
	
	free_get_buffer();
	free_mtf_table();
	if ( gIN ) fclose( gIN );
	if ( out ) fclose( out );

	return 0;
}

void copyright( void )
{
	fprintf(stderr, "\n\n Written by: Gerald Tamayo, 2006/2022\n");
}
