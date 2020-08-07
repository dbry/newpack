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

#define BLOCK_SIZE 16*1024*1024
#define HISTORY_BYTES -1

static const char *usage = "\n"
" NEWPACK  Experimental Lossless Compressor  Version 0.0.1\n"
" Copyright (c) 2020 David Bryant. All Rights Reserved.\n\n"
" Usage:   NEWPACK [-options] infile outfile\n"
"           specify '-' for stdin or stdout\n\n"
" Options: -d     = decompress\n"
"          -h<n>  = history bytes (compress only, else find best)\n"
"          -b<n>  = block size (compress only, default = 16 MB)\n"
"          -m     = decode random output based solely on model\n"
"                   (works for the first block of file only)\n"
"          -vv    = very verbose (include internal details)\n"
"          -v     = verbose\n\n"
" Warning: EXPERIMENTAL - DON'T EVEN THINK OF ARCHIVING WITH THIS!!\n";

static int decompress (FILE *infile, FILE *outfile, FILE *verbose, int use_model_only)
{
    char preamble [12];

    if (fread (preamble, 1, 12, infile) != 12 || strncmp (preamble, "newpack0.0", 10) ||
        (preamble [10] != 'F') || (preamble [11] != '1')) {
            fprintf (stderr, "not a valid newpack 0.0 file!\n");
            return 1;
    }

    if (preamble [10] == 'F')
        return newpack_decompress (infile, outfile, verbose, use_model_only);
    else {
        fprintf (stderr, "not a valid newpack 0.0 file!\n");
        return 1;
    }
}

int main (argc, argv) int argc; char **argv;
{
    int error_count = 0, decompress_mode = 0, history_bytes = HISTORY_BYTES, verbose = 0, use_model_only = 0;
    char *infilename = NULL, *outfilename = NULL;
    FILE *infile, *outfile;
    int block_size = BLOCK_SIZE;

    while (--argc) {
#ifdef _WIN32
        if ((**++argv == '-' || **argv == '/') && (*argv)[1])
#else
        if ((**++argv == '-') && (*argv)[1])
#endif
            while (*++*argv)
                switch (**argv) {
                    case 'M': case 'm':
                        use_model_only = 1;

                    case 'D': case 'd':
                        decompress_mode = 1;
                        break;

                    case 'B': case 'b':
                        block_size = strtol (++*argv, argv, 10);
                        --*argv;

                        if (block_size < 1000 || block_size > 16*1024*1024) {
                            fprintf (stderr, "block size must be 1KB - 16MB\n");
                            ++error_count;
                        }

                        break;

                    case 'H': case 'h':
                        history_bytes = strtol (++*argv, argv, 10);
                        --*argv;

                        if (history_bytes < 0 || history_bytes > 8) {
                            fprintf (stderr, "history bytes must be 0 - 8\n");
                            ++error_count;
                        }

                        break;

                    case 'V': case 'v':
                        verbose++;
                        break;

                    default:
                        fprintf (stderr, "illegal option: %c !", **argv);
                        ++error_count;
                }
        else if (!infilename)
            infilename = *argv;
        else if (!outfilename)
            outfilename = *argv;
        else {
            fprintf (stderr, "extra option: %s !", *argv);
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
            fprintf (stderr, "can't open %s for input", infilename);
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
            fprintf (stderr, "can't open %s for output", outfilename);
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
            decompress (infile, outfile, verbose ? stderr : NULL, use_model_only);
        }
        else {
            if (verbose) fprintf (stderr, "compressing %s to %s ...\n", infilename, outfilename);
            newpack_compress (infile, outfile, verbose ? stderr : NULL, block_size, history_bytes, verbose > 1);
        }

        fclose (infile);
        fclose (outfile);

        if (verbose) fprintf (stderr, "processing time: %d seconds\n", (int)(time (NULL) - start_time));
    }
}

