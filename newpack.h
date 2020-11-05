////////////////////////////////////////////////////////////////////////////
//                           **** NEWPACK ****                            //
//                        New Lossless Compressor                         //
//                   Copyright (c) 2020 David Bryant.                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

#ifndef NEWPACK_H_
#define NEWPACK_H_

int newpack_compress (FILE *infile, FILE *outfile, FILE *verbose, int block_size, int interleave, int flags);
int newpack_decompress (FILE *infile, FILE *outfile, FILE *verbose, int use_model_only);

#define INTERLEAVE_SEARCH       0
#define MAX_INTERLEAVE          16

#define HISTORY_DEPTH_MASK      0x7
#define EXHAUSTIVE_FLAG         0x8
#define TRY_RLE_FLAG            0x10
#define FORCE_RLE_FLAG          0x20
#define TRY_NUM_FLAG            0x40
#define FORCE_NUM_FLAG          0x80
#define VERY_VERBOSE_FLAG       0x1000

#define NUM_FLAGS               (TRY_NUM_FLAG | FORCE_NUM_FLAG)
#define RLE_FLAGS               (TRY_RLE_FLAG | FORCE_RLE_FLAG)

#define MAX_HISTORY_BITS        30      // determines ram footprint of bins_mask and history lookup table
#define MAX_USED_BINS           1048576 // determines ram footprint of histogram and probability tables
#define MAX_PROBABILITY         0x3FFF

#endif
