/* The MIT License

    Copyright (C) 2016 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include <htslib/sam.h>
#include "version.h"

void usage(FILE *fp) {
    fprintf(fp, "\n");
    fprintf(fp, "Program: capmq\n");
    fprintf(fp, "Version: %s (using htslib %s)\n", CAPMQ_VERSION, hts_version());
    fprintf(fp, "About:   cap mapping quality (MAPQ) to the specified value\n");
    fprintf(fp, "Usage:   capmq [options] in-file out-file\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -C INT            Cap mapping quality scores at INT\n");
    fprintf(fp, "  -s                Store original MAPQ in om:i aux tag\n");
    fprintf(fp, "  -r                Restore original MAPQ from om:i aux tag\n");
    fprintf(fp, "  -I fmt(,opt...)   Input format and format-options [auto].\n");
    fprintf(fp, "  -O fmt(,opt...)   Output format and format-options [SAM].\n");
    fprintf(fp, "\n");
    fprintf(fp,
"Standard htslib format options apply. So to create a CRAM file with lossy\n\
template names enabled and a larger number of sequences per slice, try:\n\
\n\
    capmq -O cram,lossy_names,seqs_per_slice=100000\n\
\n");
}

int main(int argc, char **argv) {
    samFile *in, *out = NULL;
    htsFormat in_fmt = {0};
    htsFormat out_fmt = {0};
    bam_hdr_t *header;
    bam1_t *b = NULL;
    int ret, opt, q, capQ = 0, storeQ = 0, restoreQ = 0;
    uint8_t *om = NULL;

    while ((opt = getopt(argc, argv, "I:O:C:srh")) != -1) {
    switch (opt) {
        case 'I':
            hts_parse_format(&in_fmt, optarg);
            break;

        case 'O':
            hts_parse_format(&out_fmt, optarg);
            break;

        case 'C':
            capQ = atoi(optarg);
            break;

        case 's':
            storeQ = 1;
            break;

        case 'r':
            restoreQ = 1;
            break;

        case 'h':
            usage(stdout);
            return 1;

        default:
            usage(stderr);
            return 1;
        }
    }
    
    if (!capQ && !restoreQ) {
        usage(stderr);
        return 1;
    }
    
    char *fnin = optind < argc ? argv[optind++] : "-";
    if (!(in = sam_open_format(fnin, "r", &in_fmt))) {
        perror(argv[optind]);
        return 1;
    }

    char mode[5] = "w";
    char *fnout = optind < argc ? argv[optind++] : "-";
    sam_open_mode(mode+1, fnout, NULL);

    if (!(out = sam_open_format(fnout, mode, &out_fmt))) {
        perror("(stdout)");
        return 1;
    }

    if (!(header = sam_hdr_read(in))) {
        fprintf(stderr, "Failed to read file header\n");
        return 1;
    }
    if (out && sam_hdr_write(out, header) != 0) {
        fprintf(stderr, "Failed to write file header\n");
        return 1;
    }

    b = bam_init1();
    if (!b) {
        fprintf(stderr, "Failed to allocate bam struct\n");
        return 1;
    }

    while ((ret = sam_read1(in, header, b)) >= 0) {
        if (b->core.tid >= 0) {
            if (restoreQ) {
                om = bam_aux_get(b, "om");
                if (om) {
                    b->core.qual = bam_aux2i(om);
                    bam_aux_del(b, om);
                    om = 0;
                }
            } else {
                q = b->core.qual;
                if (b->core.qual > capQ) {
                    if (storeQ) bam_aux_append(b, "om", 'i', 4, (uint8_t*)&q);
                    b->core.qual = capQ;
                }
            }
        }
        if (sam_write1(out, header, b) < 0) {
            fprintf(stderr, "Failed to write to output file\n");
            return 1;
        }
    }
    if (ret < -1) {
        fprintf(stderr, "Error reading input.\n");
        return 1;
    }
    bam_destroy1(b);

    bam_hdr_destroy(header);

    if (sam_close(in) != 0) {
        fprintf(stderr, "Error while closing input fd\n");
        return 1;
    }

    if (out && sam_close(out) != 0) {
        fprintf(stderr, "Error while closing output fd\n");
        return 1;
    }

    return 0;
}
