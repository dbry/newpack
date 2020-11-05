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
// creates a probability model for each byte which is based on an up to
// 24-bit hash of immediately previous bytes in the stream. This model
// consists of a bitmask representing that hash values that are actually
// encountered and a concatenation of the probability tables for each
// encountered hash value. These two structures are recusively compressed
// and sent to the output file. Then the actual data is encoded using the
// model and a range coder and sent last in the stream.
//
// To create the hash, each byte is converted to a new representation with
// from 1 to 8 bits, and then these are concatenated into a 32-bit word
// that is finally truncated to the desired size (1-24 bits). The codes are
// reordered from most to least common to minimize collisions in the
// hash.
//
// Because the model can become very large, better compression is achieved
// with large blocks of data (assuming the data is somewhat homogeneous).

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <math.h>

#include "newpack.h"

#if MAX_PROBABILITY > 0xFF
typedef unsigned short prob_t;
#else
typedef unsigned char prob_t;
#endif

static unsigned short log2linear [256], linear2log [MAX_PROBABILITY+1], tables_generated;

static void generate_log_tables (void);
static unsigned int checksum (void *data, int bcount);
static int compress_buffer_best (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, int flags);
static int compress_buffer_size (unsigned char *buffer, int num_samples, FILE *verbose, int hash_bits, int *history_bits, int *histogram, int flags);
static int compress_buffer (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, int hash_bits, int *history_bits, int *histogram, int flags);
static int decompress_interleave (unsigned char **input, int *input_size, unsigned char *output, int output_size, int use_model_only);
static int decompress_buffer (unsigned char **input, int *input_size, unsigned char *output, int output_size, int use_model_only);
static int rle_encode (void *input, int input_size, void *output);
static int rle_decode (unsigned char **input, int *input_size, void *output, int output_size);
static void interleave (void *buffer, int num_bytes, int stride, void *source);
static void deinterleave (void *buffer, int num_bytes, int stride, void *source);
static unsigned int hashed_code (unsigned int input_code, unsigned int num_hash_values);
static int peak_interleave (unsigned char *buffer, int num_bytes, int max_interleave, FILE *verbose);
static void convert2deltas (void *buffer, int byte_count);
static void convert2sums (void *buffer, int byte_count);


/* ----------------------------- compression --------------------------------*/

int newpack_compress (FILE *infile, FILE *outfile, FILE *verbose, int block_size, int interleave, int flags)
{
    int outbuffer_size = block_size + (block_size >> 3) + 65536, i;
    char *outbuffer = malloc (outbuffer_size), *outbuffer2 = NULL;
    FILE *details = (flags & VERY_VERBOSE_FLAG) ? verbose : NULL;
    long long total_bytes_read = 0, total_bytes_written = 0;
    unsigned char *inbuffer = malloc (block_size);
    int interleave_hits [MAX_INTERLEAVE] = { 0 };
    int total_blocks = 0;

    if (interleave == INTERLEAVE_SEARCH)
        outbuffer2 = malloc (outbuffer_size);

    generate_log_tables ();

    if (interleave < 0 || interleave > MAX_INTERLEAVE)
        interleave = 1;

    if (interleave > 1)
        while (block_size % interleave)
            block_size--;

    while (1) {
        int bytes_read = fread (inbuffer, 1, block_size, infile);
        int interleave_to_use = interleave, interleave_size, interleave_outsize;
        int outbytes = 0, straight_outbytes = 0, numeric_mode_count = 0;
        int flags_to_use = flags;

        if (!bytes_read)
            break;

        total_blocks++;

        if (details)
            fprintf (details, "\nblock %d containing %d bytes:\n", total_blocks, bytes_read);

        if (interleave_to_use == INTERLEAVE_SEARCH)
            interleave_to_use = peak_interleave (inbuffer, bytes_read, MAX_INTERLEAVE, details);

        total_bytes_read += bytes_read;

        if (interleave_to_use > 1 && bytes_read > interleave_to_use * 256) {
            if (interleave == INTERLEAVE_SEARCH) {

                if (flags & FORCE_NUM_FLAG) {
                    convert2deltas (inbuffer, bytes_read);
                    straight_outbytes = compress_buffer_best (inbuffer, bytes_read, outbuffer2, bytes_read + (bytes_read >> 3) + 65536, details, flags);
                    convert2sums (inbuffer, bytes_read);
                }
                else {
                    straight_outbytes = compress_buffer_best (inbuffer, bytes_read, outbuffer2, bytes_read + (bytes_read >> 3) + 65536, details, flags & ~TRY_NUM_FLAG);

                    if (flags & TRY_NUM_FLAG) {
                        int numerical_outbytes;

                        convert2deltas (inbuffer, bytes_read);
                        numerical_outbytes = compress_buffer_best (inbuffer, bytes_read, outbuffer, bytes_read + (bytes_read >> 3) + 65536, details,
                            (flags & HISTORY_DEPTH_MASK) ? flags - 1 : flags);
                        convert2sums (inbuffer, bytes_read);

                        if (numerical_outbytes < straight_outbytes) {
                            char *temp = outbuffer2;
                            outbuffer2 = outbuffer;
                            outbuffer = temp;
                            straight_outbytes = numerical_outbytes;
                            numeric_mode_count++;
                        }
                    }
                }
            }

            deinterleave (inbuffer, bytes_read, interleave_to_use, NULL);

            if (flags_to_use & HISTORY_DEPTH_MASK)
                flags_to_use--;

            if (interleave_to_use > 10 && (flags_to_use & HISTORY_DEPTH_MASK))
                flags_to_use--;
        }
        else
            interleave_to_use = 1;

        interleave_size = bytes_read / interleave_to_use;
        interleave_outsize = interleave_size + (interleave_size >> 3) + 65536;

        char *outbuffer3 = NULL;

        if (flags & TRY_NUM_FLAG)
            outbuffer3 = malloc (interleave_outsize);

        for (i = 0; i < interleave_to_use; ++i) {
            int bytes_to_process = (i < interleave_to_use - 1) ? interleave_size : bytes_read - interleave_size * i;
            int straight_outsize, numerical_outsize;

            if (flags & FORCE_NUM_FLAG) {
                convert2deltas (inbuffer + interleave_size * i, bytes_to_process);
                outbytes += straight_outsize = compress_buffer_best (inbuffer + interleave_size * i, bytes_to_process,
                    outbuffer + outbytes, outbuffer_size - outbytes, details, flags_to_use);
            }
            else {
                straight_outsize = compress_buffer_best (inbuffer + interleave_size * i, bytes_to_process,
                    outbuffer + outbytes, outbuffer_size - outbytes, details, flags_to_use & ~TRY_NUM_FLAG);

                if (flags & TRY_NUM_FLAG) {
                    convert2deltas (inbuffer + interleave_size * i, bytes_to_process);
                    numerical_outsize = compress_buffer_best (inbuffer + interleave_size * i, bytes_to_process,
                        outbuffer3, interleave_outsize, details, (flags & HISTORY_DEPTH_MASK) ? flags - 1 : flags);

                    if (numerical_outsize < straight_outsize) {
                        memcpy (outbuffer + outbytes, outbuffer3, numerical_outsize);
                        outbytes += numerical_outsize;
                        numeric_mode_count++;
                    }
                    else
                        outbytes += straight_outsize;
                }
                else
                    outbytes += straight_outsize;
            }
        }

        free (outbuffer3);

        if (details) {
            if (straight_outbytes)
                fprintf (details, "interleaved(%d) = %d (depth=%d), straight = %d (depth=%d), ratio = %.2f%%\n",
                    interleave_to_use, outbytes, flags_to_use & HISTORY_DEPTH_MASK, straight_outbytes,
                    flags & HISTORY_DEPTH_MASK, outbytes * 100.0 / straight_outbytes);
            else if (interleave_to_use != 1)
                fprintf (details, "interleaved(%d) = %d (depth=%d)\n",
                    interleave_to_use, outbytes, flags_to_use & HISTORY_DEPTH_MASK);

            if (flags & TRY_NUM_FLAG)
                fprintf (details, "numeric mode used %d of %d times\n", numeric_mode_count,
                    interleave_to_use + (interleave_to_use > 1 && interleave == INTERLEAVE_SEARCH));
        }

        total_bytes_written += fwrite ("newpack0.2", 1, 10, outfile);
        total_bytes_written += fwrite (&bytes_read, 1, sizeof (bytes_read), outfile);

        if (straight_outbytes && straight_outbytes < outbytes) {
            straight_outbytes++;
            total_bytes_written += fwrite (&straight_outbytes, 1, sizeof (straight_outbytes), outfile);
            fputc (1, outfile);
            total_bytes_written++;
            total_bytes_written += fwrite (outbuffer2, 1, straight_outbytes - 1, outfile);
            interleave_hits [0]++;
        }
        else {
            outbytes++;
            total_bytes_written += fwrite (&outbytes, 1, sizeof (outbytes), outfile);
            fputc (interleave_to_use, outfile);
            total_bytes_written++;
            total_bytes_written += fwrite (outbuffer, 1, outbytes - 1, outfile);
            interleave_hits [interleave_to_use - 1]++;
        }
    }

    if (details) fprintf (details, "\n");

    if (details && interleave == INTERLEAVE_SEARCH) {
        if (interleave_hits [0]) fprintf (details, "interleave not applied in %d of %d blocks\n", interleave_hits [0], total_blocks);
        for (i = 1; i < MAX_INTERLEAVE; ++i)
            if (interleave_hits [i])
                fprintf (details, "interleave stride of %d applied in %d of %d blocks\n", i+1, interleave_hits [i], total_blocks);
    }

    if (verbose) fprintf (verbose, "compression complete: %lld bytes --> %lld bytes (%.2f%%), %d blocks used\n",
        total_bytes_read, total_bytes_written, total_bytes_written * 100.0 / total_bytes_read, total_blocks);

    free (inbuffer);
    free (outbuffer);
    free (outbuffer2);

    return 0;
}

static int compress_buffer_best (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, int flags)
{
    int bytes_written, history_bits = 0, hash_bits, rle_encode_bytes = 0;
    int num_codes = 0, bc = num_samples;
    unsigned char *bp = buffer;
    char *outp = outbuffer;
    int histogram [256];

    outp++;         // reserve first byte for flags
    outbytes--;

    if (flags & RLE_FLAGS) {
        rle_encode_bytes = rle_encode (buffer, num_samples, NULL);

        if (rle_encode_bytes > outbytes)
            rle_encode_bytes = 0;

        if (flags & FORCE_RLE_FLAG) {
            if (rle_encode_bytes) {
                if (outbuffer) {
                    rle_encode (buffer, num_samples, outp);
                    outp [-1] = flags;
                }

                return rle_encode_bytes + 1;
            }
            else
                return 0;
        }
    }

    memset (histogram, 0, sizeof (histogram));

    while (bc--)
        if (!histogram [*bp++]++)
            num_codes++;

    if (outbuffer && !(flags & HISTORY_DEPTH_MASK))
        bytes_written = compress_buffer (buffer, num_samples, outp, outbytes, NULL, 0, &history_bits, histogram, flags & ~NUM_FLAGS);
    else
        bytes_written = compress_buffer_size (buffer, num_samples, NULL, 0, &history_bits, histogram, flags & ~NUM_FLAGS);

    if (!(flags & HISTORY_DEPTH_MASK)) {
        if (rle_encode_bytes && rle_encode_bytes <= bytes_written) {
            if (outbuffer) {
                rle_encode (buffer, num_samples, outp);
                outp [-1] = flags;
            }

            if (verbose) fprintf (verbose, "rle encoding used, %d bytes <= %d bytes\n", rle_encode_bytes, bytes_written);
            return rle_encode_bytes + 1;
        }

        if (verbose && rle_encode_bytes)
            fprintf (verbose, "rle encoding not used, %d bytes > %d bytes\n", rle_encode_bytes, bytes_written);

        if (outbuffer)
            outp [-1] = flags & ~RLE_FLAGS;

        return bytes_written + 1;
    }

    if (!bytes_written) {
        fprintf (stderr, "fatal error: compress_buffer() returned zero when it shouldn't!\n");
        return 0;
    }

    int least_bytes, best_hash_bits, best_history_bits, estimated_bytes, search_dir;
    int max_hash_bits = num_codes == 1 ? 0 : 32 - __builtin_clz (num_codes - 1);

    history_bits = MAX_HISTORY_BITS;
    estimated_bytes = compress_buffer_size (buffer, num_samples, verbose, 1, &history_bits, histogram, flags & ~NUM_FLAGS);
    best_history_bits = history_bits;
    best_hash_bits = search_dir = 1;
    least_bytes = estimated_bytes;

    if (verbose) fprintf (verbose, "hash bits = %d, history_bits = %d, size = %d\n", 1, history_bits, estimated_bytes);

    history_bits = MAX_HISTORY_BITS;
    estimated_bytes = compress_buffer_size (buffer, num_samples, verbose, max_hash_bits, &history_bits, histogram, flags & ~NUM_FLAGS);

    if (verbose) fprintf (verbose, "hash bits = %d, history_bits = %d, size = %d\n", max_hash_bits, history_bits, estimated_bytes);

    if (estimated_bytes < least_bytes) {
        least_bytes = estimated_bytes;
        best_history_bits = history_bits;
        best_hash_bits = max_hash_bits;
        search_dir = -1;
    }

    for (hash_bits = best_hash_bits + search_dir; hash_bits < max_hash_bits && hash_bits > 1; hash_bits += search_dir) {
        history_bits = MAX_HISTORY_BITS;
        estimated_bytes = compress_buffer_size (buffer, num_samples, verbose, hash_bits, &history_bits, histogram, flags & ~NUM_FLAGS);

        if (verbose) fprintf (verbose, "hash bits = %d, history_bits = %d, size = %d\n", hash_bits, history_bits, estimated_bytes);

        if (estimated_bytes > least_bytes + (least_bytes / 10) + 1000)
            break;

        if (estimated_bytes < least_bytes) {
            least_bytes = estimated_bytes;
            best_history_bits = history_bits;
            best_hash_bits = hash_bits;
        }
    }

    if (bytes_written <= least_bytes) {
        best_history_bits = best_hash_bits = 0;
        least_bytes = bytes_written;
    }

    if (verbose) fprintf (verbose, "best: hash bits = %d, history_bits = %d, size = %d\n", best_hash_bits, best_history_bits, least_bytes);

    if (rle_encode_bytes && rle_encode_bytes <= least_bytes) {
        if (outbuffer) {
            rle_encode (buffer, num_samples, outp);
            outp [-1] = flags;
        }

        if (verbose) fprintf (verbose, "rle encoding used, %d bytes <= %d bytes\n", rle_encode_bytes, least_bytes);
        return rle_encode_bytes + 1;
    }

    if (verbose && rle_encode_bytes)
        fprintf (verbose, "rle encoding not used, %d bytes > %d bytes\n", rle_encode_bytes, least_bytes);

    if (outbuffer) {
        least_bytes = compress_buffer (buffer, num_samples, outp, outbytes, verbose, best_hash_bits, &best_history_bits, histogram, flags & ~NUM_FLAGS);
        outp [-1] = flags & ~RLE_FLAGS;
    }

    return least_bytes + 1;
}

static void calculate_probabilities (int *hist, prob_t *probs, unsigned short *prob_sums, int num_codes)
{
    int code_hits = 0, min_value = 0, max_value = 0, sum_values = 0, i;

    for (i = 0; i < num_codes; ++i)
        if (hist [i]) {
            if (code_hits++) {
                if (hist [i] < min_value) min_value = hist [i];
                if (hist [i] > max_value) max_value = hist [i];
                sum_values += hist [i];
            }
            else
                sum_values = min_value = max_value = hist [i];
        }

    if (!code_hits) {        // this shouldn't happen!
        memset (probs, 0, sizeof (prob_t) * num_codes);
        memset (prob_sums, 0, sizeof (*prob_sums) * num_codes);
        return;
    }

    if (min_value == max_value) {
        for (sum_values = i = 0; i < num_codes; ++i)
            prob_sums [i] = sum_values += probs [i] = !!hist [i];

        return;
    }

    int target_min = 8 * min_value / max_value + 1;
    double scaler = 1.0;

    if (sum_values > MAX_PROBABILITY)
        scaler = (double) MAX_PROBABILITY / sum_values;

    if (min_value * scaler > target_min)
        scaler = (double) target_min / min_value;

    if (scaler == 1.0)
        for (sum_values = i = 0; i < num_codes; ++i)
            prob_sums [i] = sum_values += probs [i] = log2linear [linear2log [hist [i]]];
    else
        for (sum_values = i = 0; i < num_codes; ++i) {
            if (hist [i]) {
                int prob = floor (hist [i] * scaler + 0.5);
                if (!prob) prob = 1;
                prob_sums [i] = sum_values += probs [i] = log2linear [linear2log [prob]];
            }
            else
                prob_sums [i] = sum_values += probs [i] = 0;
        }

    return;
}

static double calculate_probabilities_size (int *hist, prob_t *probs, int num_codes)
{
    int code_hits = 0, min_value = 0, max_value = 0, sum_values = 0, total_hits, i;
    double encoding_bits = 0.0;

    for (i = 0; i < num_codes; ++i)
        if (hist [i]) {
            if (code_hits++) {
                if (hist [i] < min_value) min_value = hist [i];
                if (hist [i] > max_value) max_value = hist [i];
                sum_values += hist [i];
            }
            else
                sum_values = min_value = max_value = hist [i];
        }

    total_hits = sum_values;

    if (!code_hits) {        // this shouldn't happen!
        memset (probs, 0, sizeof (prob_t) * num_codes);
        return encoding_bits;
    }

    if (min_value == max_value) {
        for (i = 0; i < num_codes; ++i)
            probs [i] = !!hist [i];

        return log (code_hits / 1.0) / log (2) * total_hits;
    }

    int target_min = 8 * min_value / max_value + 1;
    double scaler = 1.0;

    if (sum_values > MAX_PROBABILITY)
        scaler = (double) MAX_PROBABILITY / sum_values;

    if (min_value * scaler > target_min)
        scaler = (double) target_min / min_value;

    if (scaler == 1.0)
        for (sum_values = i = 0; i < num_codes; ++i)
            sum_values += probs [i] = log2linear [linear2log [hist [i]]];
    else
        for (sum_values = i = 0; i < num_codes; ++i) {
            if (hist [i]) {
                int prob = floor (hist [i] * scaler + 0.5);
                if (!prob) prob = 1;
                sum_values += probs [i] = log2linear [linear2log [prob]];
            }
            else
                sum_values += probs [i] = 0;
        }

    for (i = 0; i < num_codes; ++i)
        if (hist [i])
            encoding_bits += log ((double) sum_values / probs [i]) / log (2) * hist [i];

    return encoding_bits;
}

static int compress_buffer (unsigned char *buffer, int num_samples, char *outbuffer, int outbytes, FILE *verbose, int hash_bits, int *history_bits, int *src_histogram, int flags)
{
    int bc = num_samples, num_codes = 0, history_bins, history_mask, bins_mask_size, used_bins, p0, p1, i;
    unsigned char hash_xlate [256], code_xlate [256], codes [256], *bp = buffer, *bins_mask;
    unsigned int total_summed_probabilities = 0, low = 0, high = 0xffffffff, mult;
    int max_history_bits = 0, max_used_bins = 1;
    unsigned short *summed_probabilities;
    unsigned short *hist_xlate_16 = NULL;
    unsigned int *hist_xlate_32 = NULL;
    char *outp = outbuffer;
    prob_t *probabilities;

    if (flags & HISTORY_DEPTH_MASK) {
        max_used_bins = 1 << ((flags & HISTORY_DEPTH_MASK) * 2 + 6);
        max_history_bits = (flags & HISTORY_DEPTH_MASK) * 3 + 9;
    }

    if (*history_bits > max_history_bits)
        *history_bits = max_history_bits;

    memcpy (outp, &num_samples, sizeof (num_samples));
    outp += sizeof (num_samples);
    outbytes -= sizeof (num_samples);

    unsigned char hash_bits_byte = hash_bits, history_bits_byte = *history_bits;

    memcpy (outp, &hash_bits_byte, sizeof (hash_bits_byte));
    outp += sizeof (hash_bits_byte);
    outbytes -= sizeof (hash_bits_byte);

    memcpy (outp, &history_bits_byte, sizeof (history_bits_byte));
    outp += sizeof (history_bits_byte);
    outbytes -= sizeof (history_bits_byte);

    history_mask = (history_bins = 1 << *history_bits) - 1;

    if (*history_bits) {
        memset (code_xlate, -1, sizeof (code_xlate));

        while (num_codes < 256) {
            int max_hits = 0, next_code = 0;

            for (i = 0; i < 256; ++i)
                if (src_histogram [i] > max_hits && code_xlate [i] == 0xff) {
                    max_hits = src_histogram [i];
                    next_code = i;
                }

            if (max_hits) {
                code_xlate [next_code] = num_codes;
                codes [num_codes++] = next_code;
            }
            else
                break;
        }

        if (outbytes + num_codes < 100)
            return 0;

        for (i = 0; i < 256; ++i)
            hash_xlate [i] = hashed_code (i, 1 << hash_bits);

        bins_mask_size = (history_bins + 7) / 8;
        bins_mask = calloc (1, bins_mask_size);
        p0 = p1 = 0;

        while (bc--) {
            bins_mask [p0 >> 3] |= 1 << (p0 & 7);
            p0 = ((p0 << hash_bits) | hash_xlate [code_xlate [*bp++]]) & history_mask;
        }

        bins_mask [p0 >> 3] |= 1 << (p0 & 7);   // make sure final code is in bins_mask even though no ptable entry required

        for (used_bins = i = 0; i < bins_mask_size; ++i)
            used_bins += __builtin_popcount (bins_mask [i]);

        if (used_bins > max_used_bins) {
            fprintf (stderr, "compress_buffer(): maximum used bin exceeded, %d\n", used_bins);
            return 0;
        }

        if (used_bins > 65536)
            hist_xlate_32 = malloc (sizeof (*hist_xlate_32) * history_bins);
        else
            hist_xlate_16 = malloc (sizeof (*hist_xlate_16) * history_bins);

        int hist_index;

        for (hist_index = p0 = 0; p0 < history_bins; p0++)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
                if (hist_xlate_32)
                    hist_xlate_32 [p0] = hist_index++;
                else
                    hist_xlate_16 [p0] = hist_index++;
            }

        int *histogram = calloc (1, num_codes * used_bins * sizeof (*histogram));

        bp = buffer;
        bc = num_samples;
        p0 = p1 = 0;

        while (bc--) {
            histogram [(hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0]) * num_codes + code_xlate [*bp]]++;
            p0 = ((p0 << hash_bits) | hash_xlate [code_xlate [*bp++]]) & history_mask;
        }

        probabilities = malloc (sizeof (prob_t) * used_bins * num_codes);
        summed_probabilities = malloc (sizeof (*summed_probabilities) * used_bins * num_codes);

        for (p0 = 0; p0 < used_bins; p0++) {
            calculate_probabilities (histogram + p0 * num_codes, probabilities + p0 * num_codes, summed_probabilities + p0 * num_codes, num_codes);
            total_summed_probabilities += summed_probabilities [p0 * num_codes + num_codes - 1];
        }

        free (histogram);

        unsigned char num_codes_byte = num_codes;

        memcpy (outp, &num_codes_byte, sizeof (num_codes_byte));
        outp += sizeof (num_codes_byte);
        outbytes -= sizeof (num_codes_byte);

        memcpy (outp, codes, num_codes);
        outp += num_codes;
        outbytes -= num_codes;

        if (verbose) fprintf (verbose, "%d codes: %02x most common, %02x least common\n", num_codes, codes [0], codes [num_codes - 1]);

        int compressed_bit_mask_size = compress_buffer_best (bins_mask, bins_mask_size, outp, outbytes, NULL,
            (flags & HISTORY_DEPTH_MASK) > 1 ? TRY_RLE_FLAG | 1 : TRY_RLE_FLAG);

        if (!compressed_bit_mask_size) {
            free (bins_mask); free (summed_probabilities); free (probabilities);
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            return 0;
        }

        outbytes -= compressed_bit_mask_size;
        outp += compressed_bit_mask_size;

        if (verbose) fprintf (verbose, "%d hash bits, %d of %d bins used, bin mask size = %d, compressed = %d\n",
            hash_bits, used_bins, history_bins, bins_mask_size, compressed_bit_mask_size);

        int ptable_encode_size = used_bins * num_codes * sizeof (unsigned char);
        unsigned char *pared_ptable = malloc (ptable_encode_size), *dst = pared_ptable;
        prob_t *src = probabilities;
        int next_code;

        for (p0 = 0; p0 < history_bins; p0++)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7)))
                for (next_code = 0; next_code < num_codes; ++next_code) {
                    p1 = ((p0 << hash_bits) | hash_xlate [next_code]) & history_mask;
                    if (bins_mask [p1 >> 3] & (1 << (p1 & 7)))
                        *dst++ = linear2log [*src++];
                    else
                        src++;
                }

        ptable_encode_size = (int)((dst - pared_ptable) * sizeof (unsigned char));
        free (bins_mask);

        int compressed_ptable_size = compress_buffer_best (pared_ptable, ptable_encode_size, outp, outbytes,
            NULL, (flags & HISTORY_DEPTH_MASK) > 1 ? TRY_RLE_FLAG | 1 : TRY_RLE_FLAG);

        if (!compressed_ptable_size) {
            free (pared_ptable); free (summed_probabilities); free (probabilities);
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            return 0;
        }

        outbytes -= compressed_ptable_size;
        outp += compressed_ptable_size;

        if (verbose) fprintf (verbose, "%d of 256 codes used, probability tables size = %d, compressed = %d\n",
            num_codes, ptable_encode_size, compressed_ptable_size);

        free (pared_ptable);
    }
    else {
        unsigned char log_ptable [256];
        num_codes = 256;

        for (i = 0; i < 256; ++i)
            hash_xlate [i] = code_xlate [i] = i;

        probabilities = malloc (sizeof (prob_t) * 256);
        summed_probabilities = malloc (sizeof (*summed_probabilities) * 256);
        calculate_probabilities (src_histogram, probabilities, summed_probabilities, 256);

        for (i = 0; i < 256; ++i)
            log_ptable [i] = linear2log [probabilities [i]];

        int compressed_ptable_size = compress_buffer_best (log_ptable, sizeof (log_ptable), outp, outbytes, NULL, FORCE_RLE_FLAG);

        if (!compressed_ptable_size) {
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            free (summed_probabilities); free (probabilities);
            return 0;
        }

        outbytes -= compressed_ptable_size;
        outp += compressed_ptable_size;
    }

    char *encoding_start = outp;

    bp = buffer;
    bc = num_samples;
    p0 = p1 = 0;

    while (bc--) {
        int bindex = history_mask ? (hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0]) * num_codes : 0;
        int code = code_xlate [*bp++];

        mult = (high - low) / summed_probabilities [bindex + num_codes - 1];

        if (!mult) {
            high = low;

            while ((high >> 24) == (low >> 24)) {
                *outp++ = high >> 24;
                outbytes--;
                high = (high << 8) | 0xff;
                low <<= 8;
            }

            mult = (high - low) / summed_probabilities [bindex + num_codes - 1];
        }

        if (code)
            low += summed_probabilities [bindex + code - 1] * mult;

        high = low + probabilities [bindex + code] * mult - 1;

        while ((high >> 24) == (low >> 24)) {
            *outp++ = high >> 24;
            outbytes--;
            high = (high << 8) | 0xff;
            low <<= 8;
        }

        if (history_mask)
            p0 = ((p0 << hash_bits) | hash_xlate [code]) & history_mask;

        if (outbytes < 100) {
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            free (summed_probabilities); free (probabilities);
            return 0;
        }
    }

    high = low;

    while ((high >> 24) == (low >> 24)) {
        *outp++ = high >> 24;
        outbytes--;
        high = (high << 8) | 0xff;
        low <<= 8;
    }

    if (verbose) fprintf (verbose, "encoded data bytes = %d, total bytes written = %d (%.2f%%)\n",
        (int) (outp - encoding_start), (int) (outp - outbuffer), (outp - outbuffer) * 100.0 / num_samples);

    free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
    free (summed_probabilities);
    free (probabilities);
    return outp - outbuffer;
}

static int compress_buffer_size (unsigned char *buffer, int num_samples, FILE *verbose, int hash_bits, int *history_bits, int *src_histogram, int flags)
{
    int bc = num_samples, num_codes = 0, smallest_size = 0, history_bins, history_mask, bins_mask_size, used_bins, p0, p1, i;
    unsigned char hash_xlate [256], code_xlate [256], /* codes [256], */ *bp = buffer, *bins_mask;
    int max_history_bits = 0, max_used_bins = 1;
    unsigned short *hist_xlate_16 = NULL;
    unsigned int *hist_xlate_32 = NULL;

    if (flags & HISTORY_DEPTH_MASK) {
        max_used_bins = 1 << ((flags & HISTORY_DEPTH_MASK) * 2 + 6);
        max_history_bits = (flags & HISTORY_DEPTH_MASK) * 3 + 9;
    }

    if (*history_bits > max_history_bits)
        *history_bits = max_history_bits;

    if (!*history_bits) {
        int outbytes = sizeof (num_samples) + 2, encoding_bytes, ptable_bytes;
        unsigned char log_ptable [256];
        prob_t probabilities [256];

        encoding_bytes = ceil ((32.0 + calculate_probabilities_size (src_histogram, probabilities, 256)) / 8.0);

        for (i = 0; i < 256; ++i)
            log_ptable [i] = linear2log [probabilities [i]];

        ptable_bytes = compress_buffer_best (log_ptable, sizeof (log_ptable), NULL, sizeof (log_ptable) * 2, NULL, FORCE_RLE_FLAG);

        return outbytes + ptable_bytes + encoding_bytes;
    }

    memset (code_xlate, -1, sizeof (code_xlate));

    while (num_codes < 256) {
        int max_hits = 0, next_code = 0;

        for (i = 0; i < 256; ++i)
            if (src_histogram [i] > max_hits && code_xlate [i] == 0xff) {
                max_hits = src_histogram [i];
                next_code = i;
            }

        if (max_hits)
            code_xlate [next_code] = num_codes++;
        else
            break;
    }

    for (i = 0; i < 256; ++i)
        hash_xlate [i] = hashed_code (i, 1 << hash_bits);

    history_mask = (history_bins = 1 << *history_bits) - 1;
    bins_mask_size = (history_bins + 7) / 8;
    bins_mask = calloc (1, bins_mask_size);

    while (1) {
        bp = buffer;
        bc = num_samples;
        used_bins = p0 = p1 = 0;

        while (bc--) {
            if (!(bins_mask [p0 >> 3] & 1 << (p0 & 7))) {
                bins_mask [p0 >> 3] |= 1 << (p0 & 7);

                if (++used_bins > max_used_bins)
                    break;
            }

            p0 = ((p0 << hash_bits) | hash_xlate [code_xlate [*bp++]]) & history_mask;
        }

        if (used_bins <= max_used_bins)
            break;

        history_mask = (history_bins /= 2) - 1;
        bins_mask_size = (history_bins + 7) / 8;
        memset (bins_mask, 0, bins_mask_size);
        --*history_bits;
    }

    if (used_bins > 65536)
        hist_xlate_32 = malloc (sizeof (*hist_xlate_32) * history_bins);
    else
        hist_xlate_16 = malloc (sizeof (*hist_xlate_16) * history_bins);

    int hist_index;

    for (hist_index = p0 = 0; p0 < history_bins; p0++)
        if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
            if (hist_xlate_32)
                hist_xlate_32 [p0] = hist_index++;
            else
                hist_xlate_16 [p0] = hist_index++;
        }

    int *histogram = calloc (1, num_codes * used_bins * sizeof (*histogram));

    bp = buffer;
    bc = num_samples;
    p0 = p1 = 0;

    while (bc--) {
        histogram [(hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0]) * num_codes + code_xlate [*bp]]++;
        p0 = ((p0 << hash_bits) | hash_xlate [code_xlate [*bp++]]) & history_mask;
    }

    int initial_history_bits = *history_bits;
    int sizes_by_history_bits [32] = { 0 };

    while (1) {
        prob_t *probabilities = malloc (sizeof (prob_t) * used_bins * num_codes);
        int outbytes = 0;

        double encoding_bits = 32.0;

        for (p1 = p0 = 0; p0 < history_bins; p0++)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
                int bin = hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0];
                encoding_bits += calculate_probabilities_size (histogram + bin * num_codes, probabilities + p1 * num_codes, num_codes);
                p1++;
            }

        outbytes += sizeof (num_samples) + 2;   // num_samples + hash_bits + history_bits;
        outbytes += 1 + num_codes;              // num_codes + that many codes

        int compressed_bit_mask_size = compress_buffer_best (bins_mask, bins_mask_size, NULL, bins_mask_size * 2 + 1024, NULL,
            TRY_RLE_FLAG | ((flags & EXHAUSTIVE_FLAG) && (flags & HISTORY_DEPTH_MASK) > 1));

        if (!compressed_bit_mask_size) {
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            free (bins_mask); free (probabilities);
            return 0;
        }

        outbytes += compressed_bit_mask_size;

        int ptable_encode_size = used_bins * num_codes * sizeof (unsigned char);
        unsigned char *pared_ptable = malloc (ptable_encode_size), *dst = pared_ptable;
        prob_t *src = probabilities;
        int next_code;

        for (p0 = 0; p0 < history_bins; p0++)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7)))
                for (next_code = 0; next_code < num_codes; ++next_code) {
                    p1 = ((p0 << hash_bits) | hash_xlate [next_code]) & history_mask;
                    if (bins_mask [p1 >> 3] & (1 << (p1 & 7)))
                        *dst++ = linear2log [*src++];
                    else
                        src++;
                }

        ptable_encode_size = (int)((dst - pared_ptable) * sizeof (unsigned char));

        int compressed_ptable_size = compress_buffer_best (pared_ptable, ptable_encode_size, NULL, ptable_encode_size * 2 + 1024, NULL,
            TRY_RLE_FLAG | ((flags & EXHAUSTIVE_FLAG) && (flags & HISTORY_DEPTH_MASK) > 1));

        if (!compressed_ptable_size) {
            free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
            free (pared_ptable); free (probabilities);
            return 0;
        }

        outbytes += compressed_ptable_size;
        free (probabilities);
        free (pared_ptable);

        outbytes += (int) ceil (encoding_bits / 8.0);
        sizes_by_history_bits [*history_bits] = outbytes;

        if (*history_bits <= 1 || *history_bits <= hash_bits)
            break;

        if (!smallest_size || outbytes < smallest_size)
            smallest_size = outbytes;

        for (p0 = 0, p1 = history_bins / 2; p0 < history_bins / 2; p0++, p1++) {
            int lower = bins_mask [p0 >> 3] & (1 << (p0 & 7));
            int upper = bins_mask [p1 >> 3] & (1 << (p1 & 7));

            if (lower || upper) {
                if (lower && upper) {
                    int *lower_hist = histogram + (hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0]) * num_codes;
                    int *upper_hist = histogram + (hist_xlate_32 ? hist_xlate_32 [p1] : hist_xlate_16 [p1]) * num_codes;
                    int entries = num_codes;

                    while (entries--)
                        *lower_hist++ += *upper_hist++;

                    used_bins--;
                }
                else if (upper) {
                    bins_mask [p0 >> 3] |= 1 << (p0 & 7);
                    if (hist_xlate_32)
                        hist_xlate_32 [p0] = hist_xlate_32 [p1];
                    else
                        hist_xlate_16 [p0] = hist_xlate_16 [p1];
                }
            }
        }

        --*history_bits;
        history_mask = (history_bins /= 2) - 1;
        bins_mask_size = (history_bins + 7) / 8;
    }

    free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
    free (histogram);
    free (bins_mask);

    for (i = initial_history_bits; i; i--)
        if (sizes_by_history_bits [i] &&
            sizes_by_history_bits [i] < sizes_by_history_bits [*history_bits])
                *history_bits = i;

    return sizes_by_history_bits [*history_bits];
}

/* ----------------------------- decompression --------------------------------*/

int newpack_decompress (FILE *infile, FILE *outfile, FILE *verbose, int use_model_only)
{
    int total_bytes_read = 0, total_bytes_written = 0;

    generate_log_tables ();

    while (1) {
        int output_size, input_size, bytes_read;
        char preamble [10];

        bytes_read = fread (preamble, 1, sizeof (preamble), infile);

        if (!bytes_read)
            break;

        if (bytes_read != sizeof (preamble) || strncmp (preamble, "newpack0.2", sizeof (preamble)) ||
            fread (&output_size, 1, sizeof (output_size), infile) != sizeof (output_size) ||
            fread (&input_size, 1, sizeof (input_size), infile) != sizeof (input_size)) {
                fprintf (stderr, "not a valid newpack 0.2 file!\n");
                return 1;
        }

        total_bytes_read += sizeof (preamble) + sizeof (input_size) + sizeof (output_size);
        unsigned char *enc_buffer = malloc (input_size), *enc_ptr = enc_buffer;
        unsigned char *dec_buffer = malloc (output_size);

        if (fread (enc_buffer, 1, input_size, infile) != input_size) {
            fprintf (stderr, "can't read %d byte-block from file!\n", input_size);
            return 1;
        }

        total_bytes_read += input_size;

        int decoded_bytes = decompress_interleave (&enc_ptr, &input_size, dec_buffer, output_size, use_model_only);

        if (decoded_bytes != output_size) {
            fprintf (stderr, "decode error, expected %d decoded bytes, got %d\n", output_size, decoded_bytes);
            return 1;
        }

        total_bytes_written += fwrite (dec_buffer, 1, decoded_bytes, outfile);
        free (enc_buffer);
        free (dec_buffer);
    }

    if (verbose && total_bytes_written)
        fprintf (verbose, "bytes read / written = %d / %d (%.2f%%)\n", total_bytes_read, total_bytes_written,
            total_bytes_read * 100.0 / total_bytes_written);

    return 0;
}

static int decompress_buffer (unsigned char **input, int *input_size, unsigned char *output, int output_size, int use_model_only)
{
    unsigned char **value_lookup, *lookup_buffer, history_bits, hash_bits;
    int flags, num_samples, history_bins, num_codes, used_bins, ptable_entries;
    unsigned char code_xlate [256], hash_xlate [256];
    unsigned short *summed_probabilities;
    prob_t *probabilities;

    unsigned int low = 0, high = 0xffffffff, value = 0, mult, sum_values;
    unsigned char *outp = output;
    int p0, p1, i;

    flags = *(*input)++;
    --*input_size;

    if (flags & RLE_FLAGS) {
        int outbytes = rle_decode (input, input_size, output, output_size);

        if (flags & NUM_FLAGS)
            convert2sums (output, outbytes);

        return outbytes;
    }

    for (i = 0; i < 4; ++i)
        ((char *) &num_samples) [i] = *(*input)++;

    *input_size -= 4;

    hash_bits = *(*input)++;
    --*input_size;
    history_bits = *(*input)++;
    --*input_size;

    if (hash_bits > 8 || hash_bits > history_bits || history_bits > MAX_HISTORY_BITS) {
        fprintf (stderr, "fatal decoding error, hash bits = %d, history bits = %d\n", hash_bits, history_bits);
        return -1;
    }

    if (history_bits) {
        unsigned char num_codes_byte;
        num_codes_byte = *(*input)++;
        num_codes = num_codes_byte ? num_codes_byte : 256;
        --*input_size;

        for (i = 0; i < num_codes; ++i) {
            code_xlate [i] = *(*input)++;
            --*input_size;
        }

        for (i = 0; i < 256; ++i)
            hash_xlate [i] = hashed_code (i, 1 << hash_bits);
    }
    else {
        for (i = 0; i < 256; ++i)
            hash_xlate [i] = code_xlate [i] = i;

        num_codes = 256;
    }

    history_bins = 1 << history_bits;

    int bins_mask_size = ((1 << history_bits) + 7) / 8;
    unsigned char *bins_mask = malloc (bins_mask_size);
    int bins_mask_decode_bytes;

    if (history_bits)
        bins_mask_decode_bytes = decompress_buffer (input, input_size, bins_mask, bins_mask_size, 0);
    else
        bins_mask_decode_bytes = bins_mask [0] = 1;

    if (bins_mask_decode_bytes != bins_mask_size) {
        fprintf (stderr, "fatal decoding error, bin mask size wrong, wanted %d, got %d\n",
            bins_mask_size, bins_mask_decode_bytes);
        return -1;
    }

    for (used_bins = i = 0; i < bins_mask_size; ++i)
        used_bins += __builtin_popcount (bins_mask [i]);

    if (used_bins > MAX_USED_BINS) {
        fprintf (stderr, "fatal decoding error, used bins = %d\n", used_bins);
        return -1;
    }

    ptable_entries = used_bins * num_codes;
    summed_probabilities = malloc (sizeof (*summed_probabilities) * ptable_entries);
    probabilities = malloc (sizeof (prob_t) * ptable_entries);

    int prob_decode_bytes;

    prob_decode_bytes = decompress_buffer (input, input_size, (unsigned char *) probabilities, ptable_entries * sizeof (prob_t), 0);

    if (prob_decode_bytes < ptable_entries * sizeof (prob_t)) {
        prob_t *unpared_ptable = malloc (ptable_entries * sizeof (prob_t)), *dst = unpared_ptable;
        unsigned char *src = (unsigned char *) probabilities;
        int next_code, unpared_size;

        for (p0 = 0; p0 < history_bins; p0++)
            if (bins_mask [p0 >> 3] & (1 << (p0 & 7)))
                for (next_code = 0; next_code < num_codes; ++next_code) {
                    p1 = ((p0 << hash_bits) | hash_xlate [next_code]) & (history_bins - 1);
                    if (bins_mask [p1 >> 3] & (1 << (p1 & 7)))
                        *dst++ = log2linear [*src++];
                    else
                        *dst++ = 0;
                }

        unpared_size = (int)((dst - unpared_ptable) * sizeof (prob_t));

        if (prob_decode_bytes == (int)((src - (unsigned char *) probabilities) * sizeof (unsigned char)) &&
            unpared_size == ptable_entries * sizeof (prob_t)) {
                memcpy (probabilities, unpared_ptable, ptable_entries * sizeof (prob_t));
                prob_decode_bytes = unpared_size;
            }

        free (unpared_ptable);
    }

    if (prob_decode_bytes != ptable_entries * sizeof (prob_t)) {
        fprintf (stderr, "fatal decoding error, probabilities size wrong, wanted %d, got %d\n",
            (int)(ptable_entries * sizeof (prob_t)), prob_decode_bytes);
        return -1;
    }

    int total_summed_probabilities = 0;

    for (p0 = 0; p0 < ptable_entries; ++p0)
        total_summed_probabilities += probabilities [p0];

    lookup_buffer = malloc (total_summed_probabilities);
    value_lookup = malloc (sizeof (*value_lookup) * used_bins);

    unsigned short *hist_xlate_16 = NULL;
    unsigned int *hist_xlate_32 = NULL;
    int hist_index = 0;

    if (used_bins > 65536)
        hist_xlate_32 = malloc (sizeof (*hist_xlate_32) * history_bins);
    else
        hist_xlate_16 = malloc (sizeof (*hist_xlate_16) * history_bins);

    unsigned char *lb_ptr = lookup_buffer;

    for (p0 = 0; p0 < history_bins; ++p0)
        if (bins_mask [p0 >> 3] & (1 << (p0 & 7))) {
            for (sum_values = i = 0; i < num_codes; ++i)
                summed_probabilities [hist_index * num_codes + i] = sum_values += probabilities [hist_index * num_codes + i];

            value_lookup [hist_index] = lb_ptr;

            for (i = 0; i < num_codes; i++) {
                int c = probabilities [hist_index * num_codes + i];

                while (c--)
                    *lb_ptr++ = i;
            }

            if (hist_xlate_32)
                hist_xlate_32 [p0] = hist_index++;
            else
                hist_xlate_16 [p0] = hist_index++;
        }

    for (i = 4; i--;) {
        value = (value << 8) | *(*input)++;
        --*input_size;
    }

    free (bins_mask);
    p0 = p1 = 0;

    int scount = num_samples;

    while (scount--) {
        int bin = hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0], bindex = bin * num_codes;

        mult = (high - low) / summed_probabilities [bindex + num_codes - 1];

        if (!mult) {
            for (i = 4; i--;) {
                value = (value << 8) | *(*input)++;
                --*input_size;
            }

            low = 0;
            high = 0xffffffff;
            mult = high / summed_probabilities [bindex + num_codes - 1];

            if (!mult) {
                fprintf (stderr, "fatal decoding error, mult = 0!\n");
                return -1;
            }
        }

        i = (value - low) / mult;

        if (i >= summed_probabilities [bindex + num_codes - 1]) {
            fprintf (stderr, "fatal decoding error, index too big!\n");
            return -1;
        }

        if ((i = value_lookup [bin] [i]))
            low += summed_probabilities [bindex + i - 1] * mult;

        high = low + probabilities [bindex + i] * mult - 1;
        *outp++ = code_xlate [i];
        p0 = ((p0 << hash_bits) | hash_xlate [i]) & (history_bins - 1);

        while ((high >> 24) == (low >> 24)) {
            value = (value << 8) | *(*input)++;
            high = (high << 8) | 0xff;
            --*input_size;
            low <<= 8;
        }
    }

    if (use_model_only) {
        static unsigned long long kernel = 0x3141592653589793;

        outp = output;
        p0 = p1 = 0;

        while (num_samples--) {
            int bin = hist_xlate_32 ? hist_xlate_32 [p0] : hist_xlate_16 [p0], bindex = bin * num_codes;

            if (bin >=used_bins || !summed_probabilities [bindex + num_codes - 1]) {
                if (hist_xlate_32)
                    bindex = (bin = hist_xlate_32 [p0 = 0]) * num_codes;
                else
                    bindex = (bin = hist_xlate_16 [p0 = 0]) * num_codes;
            }

            kernel = ((kernel << 4) - kernel) ^ 1;
            kernel = ((kernel << 4) - kernel) ^ 1;
            kernel = ((kernel << 4) - kernel) ^ 1;
            i = value_lookup [bin] [(unsigned int)(kernel >> 32) % summed_probabilities [bindex + num_codes - 1]];
            *outp++ = code_xlate [i];
            p0 = ((p0 << hash_bits) | hash_xlate [i]) & (history_bins - 1);
        }
    }

    free (hist_xlate_32 ? (void *) hist_xlate_32 : (void *) hist_xlate_16);
    free (summed_probabilities);
    free (probabilities);
    free (lookup_buffer);
    free (value_lookup);

    if (flags & NUM_FLAGS)
        convert2sums (output, (int)(outp - output));

    return outp - output;
}

static int decompress_interleave (unsigned char **input, int *input_size, unsigned char *output, int output_size, int use_model_only)
{
    int total_bytes_written = 0, interleave_value, i;

    interleave_value = *(*input)++;
    --*input_size;

    if (interleave_value < 1 || interleave_value > 16) {
        fprintf (stderr, "fatal decoding error, interleave = %d!\n", interleave_value);
        return -1;
    }

    for (i = 0; i < interleave_value; ++i)
        total_bytes_written += decompress_buffer (input, input_size, output + total_bytes_written, output_size - total_bytes_written, use_model_only);

    if (interleave_value > 1)
        interleave (output, total_bytes_written, interleave_value, NULL);

    return total_bytes_written;
}

/* ---------------------------- hashing ---------------------------- */

static unsigned int hashed_code (unsigned int input_code, unsigned int num_hash_values)
{
    return ((0 - ((input_code / num_hash_values) & 1)) ^ input_code) % num_hash_values;
}

/* ---------------------------- logrithms ---------------------------- */

#define KNEE 35

static void generate_log_tables (void)
{
    double value = MAX_PROBABILITY;
    int straight = 0, i;

    for (i = 255; i >= 0; --i) {
        int integer = (int) floor (value + 0.5);

        if (integer > i && !straight)
            log2linear [i] = integer;
        else {
            log2linear [i] = i;
            straight = 1;
        }

        // printf ("%d: %d (%.3f)\n", i, log2linear [i], value);
        linear2log [log2linear [i]] = i;

        value -= value / KNEE;
    }

    int low_log;

    for (i = 1; i <= MAX_PROBABILITY; ++i)
        if (linear2log [i])
            low_log = linear2log [i];
        else {
            int low_log_error = log2linear [low_log] - i;
            int high_log_error = log2linear [low_log + 1] - i;

            if (abs (low_log_error) <= abs (high_log_error))
                linear2log [i] = low_log;
            else
                linear2log [i] = low_log + 1;
        }

    double max_error = 0.0;
    int worst_value = 0;

    for (i = 0; i <= MAX_PROBABILITY; ++i)
        if (log2linear [linear2log [i]] != i) {
            double error = abs (log2linear [linear2log [i]] - i) / i;

            if (error > max_error) {
                max_error = error;
                worst_value = i;
            }
        }

    // fprintf (stderr, "max error = %d --> %d --> %d, %.2f%%\n",
    //     worst_value, linear2log [worst_value], log2linear [linear2log [worst_value]], max_error * 100.0);

    // fprintf (stderr, "ave error = %.2f%%\n", 50.0 / KNEE);

    tables_generated = 1;
}

/* ---------------------------- checksum ---------------------------- */

static unsigned int checksum (void *data, int bcount)
{
    unsigned char *dptr = data;
    unsigned int sum = -1;

    while (bcount--)
        sum = (sum * 3) + *dptr++;

    return sum;
}

/* -------------------------- Run-Length Encoder ---------------------------- */

int rle_encode (void *input, int input_size, void *output)
{
    unsigned char *inptr = input, *outptr = output;
    int dups, lits, bc = input_size, index = 0;

    while (bc) {
        if (bc > 1 && inptr [0] == inptr [1]) {
            for (dups = 1; dups < 127 && bc > dups + 1 && inptr [0] == inptr [dups + 1]; dups++);

            if (output) {
                outptr [index++] = 0x80 + dups;
                outptr [index++] = *inptr;
            }
            else
                index += 2;

            inptr += dups + 1;
            bc -= dups + 1;
        }
        else {
            for (lits = 1; lits < 127 && lits < bc; ++lits)
                if (lits >= 3 && inptr [lits] == inptr [lits - 1] && inptr [lits] == inptr [lits - 2]) {
                    lits -= 2;
                    break;
                }

            if (output)
                outptr [index++] = 0x00 + lits;
            else
                index++;

            bc -= lits;

            if (output)
                while (lits--)
                    outptr [index++] = *inptr++;
            else {
                index += lits;
                inptr += lits;
            }
        }
    }

    if (output)
        outptr [index++] = 0;
    else
        index++;

    return index;
}

int rle_decode (unsigned char **input, int *input_size, void *output, int output_size)
{
    unsigned char *outlim = (unsigned char *) output + output_size;
    unsigned char *outp = output;

    while (*input_size && **input) {
        if (**input & 0x80) {
            int count = (*(*input)++ & 0x7f) + 1;
            int value = *(*input)++;

            *input_size -= 2;

            if (outp + count > outlim) {
                fprintf (stderr, "rle_decode(): output buffer too small!\n");
                return -1;
            }

            while (count--)
                *outp++ = value;
        }
        else {
            int count = *(*input)++ & 0x7f;

            if (count > --*input_size) {
                fprintf (stderr, "rle_decode(): bytes exhausted!\n");
                return -1;
            }

            if (outp + count > outlim) {
                fprintf (stderr, "rle_decode(): output buffer too small!\n");
                return -1;
            }

            while (count--) {
                *outp++ = *(*input)++;
                --*input_size;
            }
        }
    }

    if (*input_size < 1) {
        fprintf (stderr, "rle_decode(): missing termination!\n");
        return -1;
    }

    ++*input;
    --*input_size;

    return outp - (unsigned char *) output;
}

/* ----------------------------- interleaving ------------------------------- */

static int peak_interleave (unsigned char *buffer, int num_bytes, int max_interleave, FILE *verbose)
{
    int interleave_factors [max_interleave+1], ave_factor = 0, max_factor, min_factor, peak, j;
    unsigned char *ptr;

    memset (interleave_factors, 0, sizeof (interleave_factors));

    for (ptr = buffer; ptr < buffer + num_bytes; ++ptr) {
        int max_len = (buffer + num_bytes - ptr - 1) / 3, n;

        for (n = 1; n <= max_interleave && n < max_len; ++n)
            if (*ptr == ptr [n] && *ptr == ptr [n+n] && *ptr == ptr [n+n+n])
                interleave_factors [n]++;
    }

    for (peak = j = 1; j <= max_interleave; ++j) {
        ave_factor += interleave_factors [j] = floor (interleave_factors [j] / sqrt (sqrt (j)) + 0.5);
        if (j == 1)
            min_factor = max_factor = interleave_factors [j];
        else {
            if (interleave_factors [j] > max_factor) {
                max_factor = interleave_factors [j];
                peak = j;
            }

            if (interleave_factors [j] < min_factor)
                min_factor = interleave_factors [j];
        }
    }

    if (!max_factor || peak == 1)
        return 1;

    ave_factor = (ave_factor + max_interleave - 1) / max_interleave;

    if (max_factor < ave_factor * 2 || num_bytes / max_factor > 1000 ||
        (peak == max_interleave && max_factor < interleave_factors [peak-1] * 2) ||
        (peak < max_interleave && max_factor < interleave_factors [peak-1] + interleave_factors [peak+1]))
            return 1;

    if (verbose) {
        char str [max_interleave*32];
        str [0] = 0;

        for (j = 1; j <= max_interleave; ++j) {
            if (interleave_factors [j] == max_factor)
                sprintf (str+strlen(str), " +%d+", interleave_factors [j]);
            else if (interleave_factors [j] == min_factor)
                sprintf (str+strlen(str), " -%d-", interleave_factors [j]);
            else
                sprintf (str+strlen(str), " %d", interleave_factors [j]);
        }

        fprintf (verbose, "peak_i = %d, factors =%s, ave = %d, ratio = %d\n", peak, str, ave_factor, num_bytes / max_factor);
    }

    return peak;
}

static void interleave (void *buffer, int num_bytes, int stride, void *source)
{
    unsigned char *bufptr = buffer, *temp, *tptr;
    int index = 0, bc = num_bytes;

    if (!source || source == buffer) {
        tptr = temp = malloc (num_bytes);
        memcpy (temp, bufptr, num_bytes);
    }
    else
        tptr = source;

    while (bc--) {
        bufptr [index] = *tptr++;
        if ((index += stride) >= num_bytes)
            index = (index % stride) + 1;
    }

    if (!source || source == buffer)
        free (temp);
}

static void deinterleave (void *buffer, int num_bytes, int stride, void *source)
{
    unsigned char *bufptr = buffer, *temp;
    int index = 0, bc = num_bytes;

    if (!source || source == buffer) {
        temp = malloc (num_bytes);
        memcpy (temp, buffer, num_bytes);
    }
    else temp = source;

    while (bc--) {
        *bufptr++ = temp [index];
        if ((index += stride) >= num_bytes)
            index = (index % stride) + 1;
    }

    if (!source || source == buffer)
        free (temp);
}

/* ----------------------------- numeric delta ------------------------------- */

static void convert2deltas (void *buffer, int byte_count)
{
    if (byte_count > 1) {
        unsigned char *bptr = buffer, last = 0;
        int bc = byte_count;

        while (bc--)
            last += *bptr++ -= last;
    }
}

static void convert2sums (void *buffer, int byte_count)
{
    if (byte_count > 1) {
        unsigned char *bptr = buffer, sum = 0;
        int bc = byte_count;

        while (bc--)
            sum = *bptr++ += sum;
    }
}
