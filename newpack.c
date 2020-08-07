////////////////////////////////////////////////////////////////////////////
//                           **** NEWPACK ****                            //
//                        New Lossless Compressor                         //
//                   Copyright (c) 2020 David Bryant.                     //
//                          All Rights Reserved.                          //
//      Distributed under the BSD Software License (see license.txt)      //
////////////////////////////////////////////////////////////////////////////

// newpack.c

// This module implements an experimental general-purpose lossless data
// compressor, the genesis of which was the "fast" DSD compression mode
// I developed for WavPack. I discovered early on that several regular
// lossless compressors (e.g. bzip2) do a surprisingly decent job on DSD
// (1-bit PCM) audio files, and after I developed the "fast" DSD mode of
// WavPack I discovered that it could often do a decent job on other types
// of data (like text).
//
// This implementation pre-scans each buffer of data to be compressed and
// creates a probability model for each byte which is based on some number
// of previous bytes in the stream. This model is compressed with simple
// RLE schemes and sent to the output file, and then the actual data is
// encoded using the model a range coder and sent last in the stream.
//
// The model can incorporate from 0 to 8 previous bytes in the prediction,
// however memory contraints limit the number that are actually practical.
// For text files around 2-3 previous bytes are generally optimum, and for
// binary files 1-2 bytes work best. Only if the file uses a very small
// symbol set can 4 or more bytes be used successfully in the model.
//
// Because the model can become very large, better compression is achieved
// with large blocks of data (assuming the data is somewhat homogeneous).
// This version allows block sizes from 1 Kbyte to 16 Mbytes.

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <math.h>

#define MAX_BYTES_PER_BIN   1280    // maximum bytes for the value lookup array (per bin)
                                    //  such that the total storage per bin = 2K (also
                                    //  counting probabilities and summed_probabilities)

#define MAX_HISTORY_BINS    2097152 // determines ram footprint of histogram
#define MAX_USED_BINS       65536   // determines ram footprint of probability tables (2K each on decode)
#define MAX_PROBABILITY     0xA0    // set to 0xff to disable RLE encoding for probabilities table
#define NO_COMPRESSION      0xFF    // send this for HISTORY_BITS to disable compression for block

// #define CAPTURE_BINS_MASK
// #define CAPTURE_PROB_TABLES

static int compress_buffer (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, char history_bytes, unsigned char *code_map_ptr);

/*------------------------------------------------------------------------------------------------------------------------*/

int newpack_compress (FILE *infile, FILE *outfile, FILE *verbose, int block_size, int history_bytes, int very_verbose)
{
    long long total_bytes_read = 0, total_bytes_written = 12;
    unsigned char *inbuffer = malloc (block_size);
    FILE *details = very_verbose ? verbose : NULL;

    fwrite ("newpack0.0", 1, 10, outfile);
    fputc ('F', outfile);
    fputc ('1', outfile);

    while (1) {
        int bytes_read = fread (inbuffer, 1, block_size, infile), bytes_written, bc;
        unsigned char code_map [32], *bp = inbuffer;
        int best_history_bytes = history_bytes;

        if (!bytes_read)
            break;

        total_bytes_read += bc = bytes_read;
        memset (code_map, 0, sizeof (code_map));

        while (bc--) {
            code_map [*bp >> 3] |= 1 << (*bp & 7);
            bp++;
        }

#if defined CAPTURE_BINS_MASK || defined CAPTURE_PROB_TABLES
        int outbytes = bytes_read + (bytes_read << 1) + 65536;
#else
        int outbytes = bytes_read + (bytes_read >> 3) + 65536;
#endif
        char *outbuffer = malloc (outbytes);

        if (history_bytes == -1) {
            char *trial_buffer = malloc (outbytes);
            int trial_bytes;

            for (bytes_written = trial_bytes = 0; trial_bytes <= 8; ++trial_bytes) {
                int trial_result = compress_buffer (inbuffer, bytes_read, trial_buffer, outbytes, details, trial_bytes, code_map);

                if (!trial_result)
                    continue;

                if (!bytes_written || trial_result < bytes_written) {
                    char *temp = outbuffer;

                    best_history_bytes = trial_bytes;
                    bytes_written = trial_result;
                    outbuffer = trial_buffer;
                    trial_buffer = temp;

                    if (trial_result + 1000 < outbytes)
                        outbytes = trial_result + 1000;
                }
            }

            free (trial_buffer);
        }
        else
            while (1) { 
                bytes_written = compress_buffer (inbuffer, bytes_read, outbuffer, outbytes, details, best_history_bytes, code_map);

                if (bytes_written)
                    break;
                else if (best_history_bytes)
                    best_history_bytes--;
                else {
                    fprintf (stderr, "fatal encoding error, compress_buffer() always returns zero!\n");
                    break;
                }
            }

        total_bytes_written += fwrite (outbuffer, 1, bytes_written, outfile);
        free (outbuffer);
    }

    free (inbuffer);

    if (verbose) fprintf (verbose, "compression complete: %lld bytes --> %lld bytes (%.2f%%)\n",
        total_bytes_read, total_bytes_written, total_bytes_written * 100.0 / total_bytes_read);

    return 0;
}

#if (MAX_PROBABILITY < 0xff)

// 0x00                     = terminate (pad with zeros)
// 0x01 - MAX_PROBABILITY   = non-zero probability values
// MAX_PROBABILITY+1 - 0xFF = zero count + MAX_PROBABILITY

static int rle_encode_bytes (unsigned char *src, int bcount, char *outfile)
{
    int max_rle_zeros = 0xff - MAX_PROBABILITY;
    int outbytes = 0, zcount = 0;

    while (bcount--) {
        if (*src) {
            while (zcount) {
                if (outfile) *outfile++ = MAX_PROBABILITY + (zcount > max_rle_zeros ? max_rle_zeros : zcount);
                zcount -= (zcount > max_rle_zeros ? max_rle_zeros : zcount);
                outbytes++;
            }

            if (outfile) *outfile++ = *src++; else src++;
            outbytes++;
        }
        else {
            zcount++;
            src++;
        }
    }

    while (zcount) {
        if (outfile) *outfile++ = MAX_PROBABILITY + (zcount > max_rle_zeros ? max_rle_zeros : zcount);
        zcount -= (zcount > max_rle_zeros ? max_rle_zeros : zcount);
        outbytes++;
    }

    if (outfile) *outfile++ = 0;
    outbytes++;

    return outbytes;
}

#endif

// 0x00        = terminate (pad with zeros)
// 0x01 - 0xFE = 0 to 253 zeros followed by one
// 0xFF        = 254 zeros (and no one)

static int rle_encode_bits (unsigned char *src, int bcount, char *outfile)
{
    int outbytes = 0, zcount = 0, index = 0;

    for (index = 0; index < bcount << 3; ++index)
        if (src [index >> 3] & (1 << (index & 7))) {
            while (zcount >= 254) {
                if (outfile) *outfile++ = 0xFF;
                zcount -= 254;
                outbytes++;
            }

            if (outfile) *outfile++ = zcount + 1;
            outbytes++;
            zcount = 0;
        }
        else
            zcount++;

    if (outfile) *outfile++ = 0;
    outbytes++;

    return outbytes;
}

static void calculate_probabilities (int hist [256], unsigned char probs [256], unsigned short prob_sums [256], int num_codes)
{
    int sum_values = 0, max_hits = 0, i;
    double scaler;

    for (i = 0; i < num_codes; ++i)
        if (hist [i] > max_hits) max_hits = hist [i];

    if (max_hits == 0) {    // this shouldn't happen!
        memset (probs, 0, sizeof (*probs) * num_codes);
        memset (prob_sums, 0, sizeof (*prob_sums) * num_codes);
        return;
    }

    if (max_hits > MAX_PROBABILITY)
        scaler = (double) MAX_PROBABILITY / max_hits;
    else
        scaler = 1.0;

    for (i = 0; i < num_codes; ++i) {
        int value;

        if (hist [i]) {
            value = floor (hist [i] * scaler + 0.5);

            if (value > MAX_PROBABILITY)
                value = MAX_PROBABILITY;
            else if (!value)
                value = 1;
        }
        else
            value = 0;

        prob_sums [i] = sum_values += value;
        probs [i] = value;
    }

#if 0   // this code reduces probability values when they are completely redundant (i.e., common divisor), but
        // this doesn't really happen often enough to make it worthwhile

    if (min_value > 1) {
        for (i = 0; i < num_codes; ++i)
            if (probs [i] % min_value)
                break;

        if (i == num_codes) {
            for (i = 0; i < num_codes; ++i) {
                prob_sums [i] /= min_value;
                probs [i] /= min_value;
            }

            // fprintf (stderr, "fixed min_value = %d, divisor = %d, probs_sum = %d\n", min_value, divisor, prob_sums [num_codes-1]);
        }
    }
#endif
}

static int compress_buffer (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, char history_bytes, unsigned char *code_map_ptr)
{
    int bc = num_samples, num_codes = 0, history_bins, bins_mask_size, used_bins, p0, p1, i;
    unsigned int low = 0, high = 0xffffffff, mult;
    unsigned char code_map [32], code_xlate [256];
    unsigned char *bp = buffer, *bins_mask;
    unsigned short *summed_probabilities;
    unsigned char *probabilities;
    int total_summed_probabilities = 0;
    int (*histogram) [256];
    char *outp = outbuffer;

    if (outbytes < 100)
        return 0;

    memcpy (code_map, code_map_ptr, sizeof (code_map));

    for (num_codes = i = 0; i < sizeof (code_map); i++)
        num_codes += __builtin_popcount (code_map [i]);

    if (num_codes < 250) {
        int incode, outcode;

        for (incode = outcode = 0; incode < 256; ++incode)
            if (code_map [incode >> 3] & (1 << (incode & 7)))
                code_xlate [incode] = outcode++;
    }
    else {
        num_codes = 256;
        for (i = 0; i < 256; ++i)
            code_xlate [i] = i;
    }

    for (history_bins = 1, i = history_bytes; i; i--)
        if ((history_bins *= num_codes) > MAX_HISTORY_BINS)
            return 0;

    bins_mask_size = (history_bins + 7) / 8;
    bins_mask = calloc (1, bins_mask_size);
    histogram = calloc (1, sizeof (*histogram) * history_bins);
    p0 = p1 = 0;

    int history_mod = history_bins / num_codes;

    while (bc--) {
        bins_mask [p0 >> 3] |= 1 << (p0 & 7);
        histogram [p0] [code_xlate [*bp]]++;

        if (num_codes != 256) {
            if (history_bytes == 0)
                p0 = *bp++ * 0;
            else if (history_bytes == 1)
                p0 = code_xlate [*bp++];
            else
                p0 = (p0 % history_mod) * num_codes + code_xlate [*bp++];
        }
        else
            p0 = ((p0 << 8) | code_xlate [*bp++]) & (history_bins - 1);
    }

    for (used_bins = i = 0; i < bins_mask_size; ++i)
        used_bins += __builtin_popcount (bins_mask [i]);

    if (used_bins > MAX_USED_BINS) {
        free (bins_mask);
        free (histogram);
        return 0;
    }

    probabilities = malloc (sizeof (*probabilities) * used_bins * num_codes);
    summed_probabilities = malloc (sizeof (*summed_probabilities) * used_bins * num_codes);
    int *hist_xlate = malloc (sizeof (int) * history_bins), hist_index;

    for (hist_index = p0 = 0; p0 < history_bins; p0++)
        if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
            calculate_probabilities (histogram [p0], probabilities + hist_index * num_codes, summed_probabilities + hist_index * num_codes, num_codes);
            total_summed_probabilities += summed_probabilities [hist_index * num_codes + num_codes - 1];
            hist_xlate [p0] = hist_index++;
        }

    free (histogram);

    // This code detects the case where the required value lookup tables grow silly big and cuts them back down. This would
    // normally only happen with large blocks or poorly compressible data. The target is to guarantee that the total memory
    // required for all three decode tables will be 2K bytes per history bin.

    while (total_summed_probabilities > used_bins * MAX_BYTES_PER_BIN) {
        int max_sum = 0, sum_values = 0, largest_bin = 0;

        for (p0 = 0; p0 < used_bins; ++p0)
            if (summed_probabilities [p0 * num_codes + num_codes - 1] > max_sum) {
                max_sum = summed_probabilities [p0 * num_codes + num_codes - 1];
                largest_bin = p0;
            }

        total_summed_probabilities -= max_sum;
        p0 = largest_bin;

        for (p1 = 0; p1 < num_codes; ++p1)
            summed_probabilities [p0 * num_codes + p1] = sum_values += probabilities [p0 * num_codes + p1] = (probabilities [p0 * num_codes + p1] + 1) >> 1;

        total_summed_probabilities += summed_probabilities [p0 * num_codes + num_codes - 1];
    }

    memcpy (outp, &num_samples, sizeof (num_samples)); outp += sizeof (num_samples); outbytes -= sizeof (num_samples);
    memcpy (outp, &history_bytes, sizeof (history_bytes)); outp += sizeof (history_bytes); outbytes -= sizeof (history_bytes);

    unsigned char num_codes_byte = num_codes, max_probability = MAX_PROBABILITY;
    memcpy (outp, &num_codes_byte, sizeof (num_codes_byte)); outp += sizeof (num_codes_byte); outbytes -= sizeof (num_codes_byte);

    if (num_codes != 256) {
        memcpy (outp, code_map, sizeof (code_map)); outp += sizeof (code_map); outbytes -= sizeof (code_map);
    }

    memcpy (outp, &max_probability, sizeof (max_probability)); outp += sizeof (max_probability); outbytes -= sizeof (max_probability);

    int encoded_bins_mask_size = rle_encode_bits (bins_mask, bins_mask_size, NULL);

#ifdef CAPTURE_BINS_MASK
    if (0) {
#else
    if (encoded_bins_mask_size < bins_mask_size) {
#endif
        if (outbytes < encoded_bins_mask_size + 100) {
            free (bins_mask); free (probabilities); free (summed_probabilities); free (hist_xlate);
            return 0;
        }

        *outp++ = 1; outbytes--;
        outp += rle_encode_bits (bins_mask, bins_mask_size, outp); outbytes -= encoded_bins_mask_size;

        if (verbose) fprintf (verbose, "%d history bytes, %d of %d bins used, bin mask size = %d, compressed = %d\n",
            history_bytes, used_bins, history_bins, bins_mask_size, encoded_bins_mask_size);
    }
    else {
        if (outbytes < bins_mask_size + 100) {
            free (bins_mask); free (probabilities); free (summed_probabilities); free (hist_xlate);
            return 0;
        }

        *outp++ = 0; outbytes--;
        memcpy (outp, bins_mask, bins_mask_size); outp += bins_mask_size; outbytes -= bins_mask_size;

        if (verbose) fprintf (verbose, "%d history bytes, %d of %d bins used, bin mask size = %d\n",
            history_bytes, used_bins, history_bins, bins_mask_size);

#ifdef CAPTURE_BINS_MASK
        FILE *capture = fopen ("capture.bin", "wb");

        if (capture) {
            fwrite (bins_mask, 1, bins_mask_size, capture);
            fclose (capture);
        }
#endif
    }

    free (bins_mask);

#if (MAX_PROBABILITY < 0xFF) && !defined (CAPTURE_PROB_TABLES)
    int encoded_prob_size = rle_encode_bytes (probabilities, used_bins * num_codes, NULL);

    if (outbytes < encoded_prob_size + 100) {
        free (probabilities); free (summed_probabilities); free (hist_xlate);
        return 0;
    }

    rle_encode_bytes (probabilities, used_bins * num_codes, outp); outp += encoded_prob_size; outbytes -= encoded_prob_size;
    if (verbose) fprintf (verbose, "%d of 256 codes used, probability tables size = %d, compressed = %d\n",
        num_codes, used_bins * num_codes, encoded_prob_size);
#else
    if (outbytes < used_bins * num_codes + 100) {
        free (probabilities); free (summed_probabilities); free (hist_xlate);
        return 0;
    }

    memcpy (outp, probabilities, used_bins * num_codes); outp += used_bins * num_codes; outbytes -= used_bins * num_codes;
    if (verbose) fprintf (verbose, "%d of 256 codes used, probability tables size = %d, not compressed\n",
        num_codes, used_bins * num_codes);

#ifdef CAPTURE_PROB_TABLES
        FILE *capture = fopen ("capture.bin", "wb");

        if (capture) {
            fwrite (probabilities, 1, used_bins * num_codes, capture);
            fclose (capture);
        }
#endif

#endif

    char *encoding_start = outp;

    bp = buffer;
    bc = num_samples;
    p0 = p1 = 0;

    while (bc--) {
        int bin = hist_xlate [p0], bindex = bin * num_codes;
        int code = code_xlate [*bp++];

        mult = (high - low) / summed_probabilities [bindex + num_codes - 1];

        if (!mult) {
            high = low;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24; outbytes--;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            mult = (high - low) / summed_probabilities [bindex + num_codes - 1];
        }

        if (code)
            low += summed_probabilities [bindex + code - 1] * mult;

        high = low + probabilities [bindex + code] * mult - 1;

        while ((high >> 24) == (low >> 24)) {
            *outp++ = high >> 24; outbytes--;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        if (num_codes != 256) {
            if (history_bytes == 0)
                p0 = 0;
            else if (history_bytes == 1)
                p0 = code;
            else
                p0 = (p0 % history_mod) * num_codes + code;
        }
        else
            p0 = ((p0 << 8) | code) & (history_bins - 1);

        if (outbytes < 100) {
            free (hist_xlate); free (summed_probabilities); free (probabilities);
            return 0;
        }
    }

    high = low;

    while ((high >> 24) == (low >> 24)) {
        *outp++ = high >> 24; outbytes--;
        high = (high << 8) | 0xff;
        low <<= 8;
    }

    if (verbose) fprintf (verbose, "encoded data bytes = %d, total bytes written = %d\n",
        (int) (outp - encoding_start), (int) (outp - outbuffer));

    free (hist_xlate);
    free (summed_probabilities);
    free (probabilities);
    return outp - outbuffer;
}

/*------------------------------------------------------------------------------------------------------------------------*/

#define BUFFER_SIZE (1024*1024)

int newpack_decompress (FILE *infile, FILE *outfile, FILE *verbose, int use_model_only)
{
    unsigned char *input_buffer = malloc (BUFFER_SIZE), *inp = input_buffer, *inpx = input_buffer;
    unsigned char *probabilities, **value_lookup, *lookup_buffer, history_bytes, max_probability;
    long long total_bytes_read = 12, total_bytes_written = 0;
    unsigned char code_map [32], code_xlate [256], *vp;
    int num_samples, history_bins, num_codes;
    unsigned short *summed_probabilities;

    while (1) {
        unsigned int low = 0, high = 0xffffffff, mult, value, sum_values;
        char *output_buffer, *outp;
        int p0, p1, i;

        for (i = 0; i < 4; ++i) {
            if (inp == inpx) {
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                if (inp == inpx)
                    break;
            }

            ((char *) &num_samples) [i] = *inp++;
        }

        if (i != 4)
            break;

        if (inp == inpx) {
            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (inp == inpx)
                break;
        }

        history_bytes = *inp++;

        if (history_bytes == NO_COMPRESSION) {
            outp = output_buffer = malloc (num_samples);

            while (num_samples--) {
                if (inp == inpx) {
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    if (inp == inpx)
                        break;
                }

                *outp++ = *inp++;
            }

            total_bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
            free (output_buffer);
            continue;
        }

        if (inp == inpx) {
            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (inp == inpx)
                break;
        }

        unsigned char num_codes_byte;
        num_codes_byte = *inp++;
        num_codes = num_codes_byte ? num_codes_byte : 256;

        if (num_codes != 256) {
            int incode, outcode;

            for (i = 0; i < sizeof (code_map); ++i) {
                if (inp == inpx) {
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    if (inp == inpx)
                        break;
                }

                code_map [i] = *inp++;
            }

            if (i != sizeof (code_map))
                break;

            for (incode = outcode = 0; incode < 256; ++incode)
                if (code_map [incode >> 3] & (1 << (incode & 7)))
                    code_xlate [outcode++] = incode;
        }
        else
            for (i = 0; i < 256; ++i)
                code_xlate [i] = i;

        for (history_bins = 1, i = history_bytes; i; i--)
            if ((history_bins *= num_codes) > MAX_HISTORY_BINS) {
                fprintf (stderr, "fatal decoding error, history bins = %d\n", history_bins);
                return 1;
            }

        if (inp == inpx) {
            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (inp == inpx)
                break;
        }

        max_probability = *inp++;

        int bin_mask_size = (history_bins + 7) / 8, used_bins = 0, bit_index;
        unsigned char *bins_mask = malloc (bin_mask_size), *ptr, byte;

        if (inp == inpx) {
            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            if (inp == inpx)
                break;
        }

        byte = *inp++;

        if (byte) {
            memset (bins_mask, 0, bin_mask_size);
            bit_index = -1;

            while (1) {
                if (inp == inpx) {
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    if (inp == inpx)
                        break;
                }

                byte = *inp++;

                if (!byte)
                    break;

                if (byte < 0xFF) {
                    bit_index += byte;

                    if (bit_index >= bin_mask_size << 3)
                        break;

                    bins_mask [bit_index >> 3] |= 1 << (bit_index & 7);
                }
                else
                    bit_index += 0xFE;
            }
        }
        else {
            for (i = 0; i < bin_mask_size; ++i) {

                if (inp == inpx) {
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                    if (inp == inpx)
                        break;
                }

                bins_mask [i] = *inp++;
            }

            if (i != bin_mask_size)
                break;
        }

        for (i = 0; i < bin_mask_size; ++i)
            used_bins += __builtin_popcount (bins_mask [i]);

        if (used_bins > MAX_USED_BINS) {
            fprintf (stderr, "fatal decoding error, used bins = %d\n", used_bins);
            return 1;
        }

        lookup_buffer = malloc (used_bins * MAX_BYTES_PER_BIN);
        value_lookup = malloc (sizeof (*value_lookup) * used_bins);
        summed_probabilities = malloc (sizeof (*summed_probabilities) * used_bins * num_codes);
        ptr = probabilities = malloc (sizeof (*probabilities) * used_bins * num_codes);

        if (max_probability < 0xff) {
            while (1) {
                if (inp == inpx)
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                if (*inp > max_probability) {
                    int zcount = *inp++ - max_probability;

                    while (zcount--)
                        *ptr++ = 0;
                }
                else if (*inp)
                    *ptr++ = *inp++;
                else {
                    inp++;
                    break;
                }
            }

            if (ptr != probabilities + used_bins * num_codes) {
                fprintf (stderr, "fatal decoding error, read %d bytes, expected %d\n", (int)(ptr - probabilities), used_bins * num_codes);
                return 1;
            }
        }
        else {
            for (i = 0; i < used_bins * num_codes; ++i) {
                if (inp == inpx)
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                *ptr++ = *inp++;
            }
        }

        outp = output_buffer = malloc (num_samples);
        int *hist_xlate = malloc (sizeof (int) * history_bins), hist_index = 0;
        memset (hist_xlate, -1, sizeof (int) * history_bins);
        unsigned char *lb_ptr = lookup_buffer;
        int total_summed_probabilities = 0;

        for (p0 = 0; p0 < history_bins; ++p0)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
                for (sum_values = i = 0; i < num_codes; ++i)
                    summed_probabilities [hist_index * num_codes + i] = sum_values += probabilities [hist_index * num_codes + i];

                if ((total_summed_probabilities += sum_values) > used_bins * MAX_BYTES_PER_BIN) {
                    fprintf (stderr, "summed probabilities exceeded limit!\n");
                    return 1;
                }

                value_lookup [hist_index] = lb_ptr;

                for (i = 0; i < num_codes; i++) {
                    int c = probabilities [hist_index * num_codes + i];

                    while (c--)
                        *lb_ptr++ = i;
                }

                hist_xlate [p0] = hist_index++;
            }

        for (i = 4; i--;) {
            if (inp == inpx)
                total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

            value = (value << 8) | *inp++;
        }

        free (bins_mask);
        int history_mod = history_bins / num_codes;
        p0 = p1 = 0;

        while (num_samples--) {
            int bin = hist_xlate [p0], bindex = bin * num_codes;

            if (use_model_only) {
                static unsigned long long kernel = 0x3141592653589793;

                if (bin == -1)
                    bin = hist_xlate [p0 = 0];

                kernel = ((kernel << 4) - kernel) ^ 1;
                kernel = ((kernel << 4) - kernel) ^ 1;
                kernel = ((kernel << 4) - kernel) ^ 1;
                i = value_lookup [bin] [(unsigned int)(kernel >> 32) % summed_probabilities [bindex + num_codes - 1]];
            }
            else {
                if (bin == -1) {
                    fprintf (stderr, "attempt to access unused history bin %d!\n", p0);
                    return 1;
                }

                mult = (high - low) / summed_probabilities [bindex + num_codes - 1];

                if (!mult) {
                    for (i = 4; i--;) {
                        if (inp == inpx)
                            total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                        value = (value << 8) | *inp++;
                    }

                    low = 0;
                    high = 0xffffffff;
                    mult = high / summed_probabilities [bindex + num_codes - 1];

                    if (!mult) {
                        fprintf (stderr, "fatal decoding error, mult = 0!\n");
                        return 1;
                    }
                }

                i = (value - low) / mult;

                if (i >= summed_probabilities [bindex + num_codes - 1]) {
                    fprintf (stderr, "fatal decoding error, index too big!\n");
                    return 1;
                }

                if ((i = value_lookup [bin] [i]))
                    low += summed_probabilities [bindex + i - 1] * mult;

                high = low + probabilities [bindex + i] * mult - 1;
            }

            *outp++ = code_xlate [i];

            if (num_codes != 256) {
                if (history_bytes == 0)
                    p0 = 0;
                else if (history_bytes == 1)
                    p0 = i;
                else
                    p0 = (p0 % history_mod) * num_codes + i;
            }
            else
                p0 = ((p0 << 8) | i) & (history_bins - 1);

            while ((high >> 24) == (low >> 24)) {
                if (inp == inpx)
                    total_bytes_read += (inpx = input_buffer + fread (inp = input_buffer, 1, BUFFER_SIZE, infile)) - inp;

                value = (value << 8) | *inp++;
                high = (high << 8) | 0xff;
                low <<= 8;
            }
        }

        total_bytes_written += fwrite (output_buffer, 1, outp - output_buffer, outfile);
        free (summed_probabilities);
        free (probabilities);
        free (lookup_buffer);
        free (value_lookup);
        free (output_buffer);
        free (hist_xlate);

        if (use_model_only)
            break;
    }

    free (input_buffer);

    if (verbose) fprintf (verbose, "decompression complete: %lld bytes --> %lld bytes (%.2f%%)\n",
        total_bytes_read, total_bytes_written, total_bytes_read * 100.0 / total_bytes_written);

    return 0;
}
