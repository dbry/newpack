////////////////////////////////////////////////////////////////////////////
//                           **** NEWPACK ****                            //
//                        New Lossless Compressor                         //
//                   Copyright (c) 2020 David Bryant.                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// main.c

// This module implements the command-line program to experiment with the
// new compression algorithms.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>

#include "newpack.h"

#define HISTORY_DEPTH_DEFAULT 5
#define MAX_BLOCK_SIZE 100000000
#define DEF_BLOCK_SIZE 20000000
#define MIN_BLOCK_SIZE 100000

static const char *usage = "\n"
" NEWPACK  Experimental General-Purpose Lossless Compressor  Version 0.0.2\n"
" Copyright (c) 2020 David Bryant. All Rights Reserved.\n\n"
" Usage:   NEWPACK [-d] [-options] infile outfile\n"
"           specify '-' for stdin or stdout\n\n"
" Options: -d     = decompression (default is compression)\n"
"          -a     = use all specialty modes (equivalent to -irn)\n"
"          -[1-7] = probability model depth (default = 5 for 20 MB block)\n"
"          -0     = no history employed in model (symbol frequency only)\n"
"          -e     = exhaustive search (very slow for very little gain)\n"
"          -i     = automatically detect interleave and use if better\n"
"          -i<n>  = force specified interleave stride (1-16, default = 1)\n"
"          -b<n>  = specify block size (1 - 100 MB or 100000+ bytes,\n"
"                   default = 20 MB, does not override history depth\n"
"          -l     = use long blocks (100 MB; history depth set to -6)\n"
"          -s     = use short blocks (4 MB; history depth set to -4)\n"
"          -t     = use tiny blocks (1 MB; history depth set to -3)\n"
"          -n     = try numerical data type preprocessing and use if better\n"
"          -nn    = force numerical data type preprocessing (i.e., deltas)\n"
"          -r     = try simple run-length encoding and use if better\n"
"          -rr    = force simple run-length encoding (really just for testing)\n"
"          -m     = decode random output based solely on model\n"
"          -vv    = very verbose messaging (include internal details)\n"
"          -v     = verbose messaging\n\n"
" Warning: EXPERIMENTAL - DON'T EVEN THINK OF ARCHIVING WITH THIS!!\n";

int main (argc, argv) int argc; char **argv;
{
    int history_depth = HISTORY_DEPTH_DEFAULT, exhaustive_mode = 0, numerical_mode = 0, rle_mode = 0;
    int error_count = 0, decompress_mode = 0, interleave = 1, verbose = 0, use_model_only = 0;
    int history_depth_override = 0, block_size_override = 0;
    char *infilename = NULL, *outfilename = NULL;
    int block_size = DEF_BLOCK_SIZE;
    FILE *infile, *outfile;

    while (--argc) {
#ifdef _WIN32
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {
                    case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
                        history_depth = **argv & HISTORY_DEPTH_MASK;
                        history_depth_override = 1;
                        break;

                    case 'A': case 'a':
                        interleave = INTERLEAVE_SEARCH;
                        rle_mode = numerical_mode = 1;
                        break;

                    case 'M': case 'm':
                        use_model_only = 1;

                    case 'D': case 'd':
                        decompress_mode = 1;
                        break;

                    case 'E': case 'e':
                        exhaustive_mode = 1;
                        break;

                    case 'R': case 'r':
                        rle_mode++;
                        break;

                    case 'N': case 'n':
                        numerical_mode++;
                        break;

                    case 'B': case 'b':
                        block_size_override = block_size = strtol (++*argv, argv, 10);
                        --*argv;

                        if (block_size < MIN_BLOCK_SIZE)
                            block_size *= 1000000;

                        if (block_size < MIN_BLOCK_SIZE || block_size > MAX_BLOCK_SIZE) {
                            fprintf (stderr, "block size must be 1 - 100 MB or 100000+ bytes\n");
                            ++error_count;
                        }

                        break;

                    case 'I': case 'i':
                        interleave = strtol (++*argv, argv, 10);
                        --*argv;

                        if (interleave < 0 || interleave > 16) {
                            fprintf (stderr, "interleave must be 1 - 16\n");
                            ++error_count;
                        }

                        break;

                    case 'L': case 'l':
                        if (!block_size_override)
                            block_size = MAX_BLOCK_SIZE;

                        if (!history_depth_override)
                            history_depth = 6;

                        break;

                    case 'S': case 's':
                        if (!block_size_override)
                            block_size = 4000000;

                        if (!history_depth_override)
                            history_depth = 4;

                        break;

                    case 'T': case 't':
                        if (!block_size_override)
                            block_size = 1000000;

                        if (!history_depth_override)
                            history_depth = 3;

                        break;

                    case 'V': case 'v':
                        verbose++;
                        break;

                    default:
                        fprintf (stderr, "illegal option: %c !\n", **argv);
                        ++error_count;
                }
        else if (!infilename)
            infilename = *argv;
        else if (!outfilename)
            outfilename = *argv;
        else {
            fprintf (stderr, "extra option: %s !\n", *argv);
            error_count++;
        }
    }

    if (!outfilename || error_count) {
        printf ("%s", usage);
        return 0;
    }

    if (strcmp (infilename, "-")) {
        infile = fopen (infilename, "rb");

        if (!infile) {
            fprintf (stderr, "can't open %s for input\n", infilename);
            return 1;
        }
    }
    else {
#ifdef _WIN32
        _setmode (_fileno (stdin), O_BINARY);
#endif
        infile = stdin;
    }

    if (strcmp (outfilename, "-")) {
        outfile = fopen (outfilename, "wb");

        if (!outfile) {
            fprintf (stderr, "can't open %s for output\n", outfilename);
            return 1;
        }
    }
    else {
#ifdef _WIN32
        _setmode (_fileno (stdout), O_BINARY);
#endif
        outfile = stdout;
    }

    if (!error_count) {
        time_t start_time = time (NULL);

        if (decompress_mode) {
            if (verbose) fprintf (stderr, "decompressing %s to %s ...\n", infilename, outfilename);
            newpack_decompress (infile, outfile, verbose ? stderr : NULL, use_model_only);
        }
        else {
            int flags = history_depth;
            if (verbose) fprintf (stderr, "compressing %s to %s ...\n", infilename, outfilename);

            if (verbose > 1)
                flags |= VERY_VERBOSE_FLAG;

            if (exhaustive_mode)
                flags |= EXHAUSTIVE_FLAG;

            if (numerical_mode)
                flags |= numerical_mode > 1 ? FORCE_NUM_FLAG : TRY_NUM_FLAG;

            if (rle_mode)
                flags |= rle_mode > 1 ? FORCE_RLE_FLAG : TRY_RLE_FLAG;

            newpack_compress (infile, outfile, verbose ? stderr : NULL, block_size, interleave, flags);
        }

        fclose (infile);
        fclose (outfile);

        if (verbose) fprintf (stderr, "compression time: %d seconds\n", (int)(time (NULL) - start_time));
    }

    return 0;
}

