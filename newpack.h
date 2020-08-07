////////////////////////////////////////////////////////////////////////////
//                           **** NEWPACK ****                            //
//                        New Lossless Compressor                         //
//                   Copyright (c) 2020 David Bryant.                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

#ifndef NEWPACK_H_
#define NEWPACK_H_

int newpack_compress (FILE *infile, FILE *outfile, FILE *verbose, int block_size, int history_bytes, int very_verbose);
int newpack_decompress (FILE *infile, FILE *outfile, FILE *verbose, int use_model_only);

#endif
