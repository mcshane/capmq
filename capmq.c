/* The MIT License

    Copyright (C) 2016-2017 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>
            Jennifer Liddle <js10@sanger.ac.uk>

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
#include <stdbool.h>
#include <getopt.h>
#include <math.h>
#include <errno.h>

#include <htslib/sam.h>
#include <cram/sam_header.h>
#include "version.h"

// Read Group / Capping Quality pair (for -g or -G options)
typedef struct {
    char *rg;   // read group
    uint8_t  capQ;  // map quality
} rgv_t;

// Array of RG/Value pairs (for -g or -G options)
typedef struct {
    int end;    // number of entries
    int max;    // maximum number of entries
    rgv_t **rgv; // list of RG/Value pairs
} rgva_t;

// initialise the rgv array
static rgva_t *rgva_init(int max)
{
    rgva_t *rgva = calloc(1, sizeof(rgva_t));
    rgva->end = 0;
    rgva->max = max;
    rgva->rgv = calloc(max, sizeof(rgv_t *));
    return rgva;
}

static int rgv_cmp(const void *p1, const void *p2)
{
    rgv_t *rgv1 = *(rgv_t **)p1;
    rgv_t *rgv2 = *(rgv_t **)p2;
    return strcmp(rgv1->rg, rgv2->rg);
}

static int rgv_keycmp(const void *p1, const void *p2)
{
    char *key = (char *)p1;
    rgv_t *rgv2 = *(rgv_t **)p2;
    return strcmp(key, rgv2->rg);
}

static void rgva_sort(rgva_t *rgva)
{
    qsort(rgva->rgv, rgva->end, sizeof(rgv_t*), rgv_cmp);
}

static rgv_t *rgva_find(rgva_t *rgva, char *rg)
{
    void *result = bsearch(rg, rgva->rgv, rgva->end, sizeof(rgv_t *), rgv_keycmp);
    return result ? *(rgv_t **)result : NULL;
}


// push entry onto rgv array
static void rgva_push(rgva_t *rgva, rgv_t *rgv)
{
    if (rgva->end == rgva->max) {
        // expand the array
        rgva->max *= 2;
        rgva->rgv = realloc(rgva->rgv, rgva->max * sizeof(rgv_t *));
    }
    rgva->rgv[rgva->end] = rgv;
    rgva->end++;
}

// Global options
typedef struct {
    bool verbose;
    uint8_t capQ;
    bool storeQ;
    bool restoreQ;
    bool freemix;
    uint8_t minQ;
    samFile *in;
    samFile *out;
    char *argv_list;
    rgva_t *rgva;
} opts_t;

static void free_opts(opts_t *opts)
{
    int n;
    if (!opts) return;
    if (opts->in) sam_close(opts->in);
    if (opts->out) sam_close(opts->out);
    free(opts->argv_list);
    if (opts->rgva) {
        for (n=0; n < opts->rgva->end; n++) {
            rgv_t *rgv = opts->rgva->rgv[n];
            free(rgv->rg);
            free(rgv);
        }
    }
}

// Convert freemix value to quality value, or die trying
static inline uint8_t f2q(double f)
{
    double qd;
    if (f == 0) {
      return 0;
    } else if (f < 0) {
      fprintf(stderr, "cannot interpret negative freemix value %f\n", f);
    }
    qd = -10.0 * log10(f);
    if (qd < 0.0) {
      return 0;
    }
    if (qd > 255.0) {
      return 255;
    }
    return (uint8_t) qd;
}

/*
 * Get uint8_t from string, or die trying
 */
static uint8_t uint8_from_str(char *str)
{
    uint8_t retval;
    unsigned long int val;
    errno = 0;
    val = strtoul(str, NULL, 10);
    if ((errno == ERANGE && (val == ULONG_MAX || val == 0)) 
        || (errno != 0 && val == 0)) {
        fprintf(stderr, "failed to parse string `%s' as unsigned long\n", str);
        perror("strtoul");
        exit(1);
    }
    if (val > 255) {
        fprintf(stderr, "ERROR: value of string `%s' was out of range (must be between 0 and 255 inclusive)\n", str);
        exit(1);
    } else {
        retval = (uint8_t)val;
    }
    return retval;
}

/*
 * Get valid capq from string, or die trying
 */
static uint8_t capq_from_str(char *str, bool freemix)
{
    uint8_t capQ;

    if (freemix) {
      double val;
      errno = 0;
      val = strtod(str, NULL);
      if (errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALL)) {
          fprintf(stderr, "overflow parsing freemix argument `%s' as double\n", str);
          perror("strtod");
          exit(1);
      } else if (errno == ERANGE && val == 0) {
          fprintf(stderr, "WARNING: underflow parsing freemix argument `%s' as double (will set capQ to 255)\n", str);
          capQ = 255;
      } else if (errno != 0 && val == 0) {
          fprintf(stderr, "failed to parse freemix argument `%s' as double\n", str);
          perror("strtod");
          exit(1);
      } else {
          capQ = f2q(val);
      }
    } else {
      capQ = uint8_from_str(str);
    }
    return capQ;
}

/*
 * Parse an RG:val pair
 * Format is RG:val
 *
 * where RG  is the Read Group
 *       val is the capQ value (or freemix value)
 */
static void parse_rgv(rgva_t *rgva, char *arg, bool freemix)
{
    char *argstr = strdup(arg);
    char *s = strrchr(argstr,':');
    if (s) {
        rgv_t *rgv = calloc(1, sizeof(rgv_t));
        *s=0;
        rgv->rg = strdup(argstr);
        rgv->capQ = capq_from_str(s+1, freemix);
        rgva_push(rgva,rgv);
    }
    free(argstr);
}

/*
 * Parse RG, val pairs from a tab delimited file
 */
static void parse_gfile(char *fname, opts_t *opts)
{
    char *buf = NULL;
    size_t n;
    FILE *fh = fopen(fname,"r");
    if (!fh) {
        fprintf(stderr,"ERROR: Can't open file %s: %s\n", fname, strerror(errno));
        exit(1);
    }

    while (getline(&buf, &n, fh) > 0) {
        if (buf[strlen(buf)-1] == '\n') buf[strlen(buf)-1]=0;   // remove trailing lf
        if (*buf && *buf!='#') {    // ignore blank lines and comments
            char *s = strchr(buf,'\t');
            if (s) {
                *s=':'; parse_rgv(opts->rgva, buf, opts->freemix);
            }
        }
        free(buf); buf=NULL;
    }
    fclose(fh);
}

/*
 * convert SAM_hdr to bam_hdr
 */
static void sam_hdr_unparse2(SAM_hdr *sh, bam_hdr_t *h)
{
    free(h->text);
    sam_hdr_rebuild(sh);
    h->text = strdup(sam_hdr_str(sh));
    h->l_text = sam_hdr_length(sh);
    sam_hdr_free(sh);
}

/*
 * Display usage information
 */
static void usage(FILE *fp) {
    fprintf(fp, "\n");
    fprintf(fp, "Program: capmq\n");
    fprintf(fp, "Version: %s (using htslib %s)\n", CAPMQ_VERSION, hts_version());
    fprintf(fp, "About:   cap mapping quality (MAPQ) to the specified value\n");
    fprintf(fp, "Usage:   capmq [options] in-file out-file\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -C max              Cap MAPQ at max (default: 255)\n");
    fprintf(fp, "  -S                  Do not store original MAPQ in om:i aux tag\n");
    fprintf(fp, "  -r                  Restore original MAPQ from om:i aux tag\n");
    fprintf(fp, "  -v                  verbose\n");
    fprintf(fp, "  -g RG:max           Cap MAPQ for read group IDs.\n");
    fprintf(fp, "                      This can be specified more than once, and if specified\n");
    fprintf(fp, "                      will overide the -C paramater for those read groups.\n");
    fprintf(fp, "  -G filename         As for -g, but group ID/max value pairs are read from\n");
    fprintf(fp, "                      a tab delimited file.\n");
    fprintf(fp, "  -f                  The values to -C, -g or in the file specified with -G\n");
    fprintf(fp, "                      are NOT maximum MAPQ scores, but estimated fraction of\n");
    fprintf(fp, "                      contamination (e) from which to calculate the maximum\n");
    fprintf(fp, "                      MAPQ as int(10*log10(1/e)).\n");
    fprintf(fp, "  -m min              Minimum MAPQ. Do not set the calculated quality\n");
    fprintf(fp, "                      to less than this value. Only used with -f\n");
    fprintf(fp, "                      (default: 0)\n");
    fprintf(fp, "  -I fmt(,opt...)     Input format and format-options [auto].\n");
    fprintf(fp, "  -O fmt(,opt...)     Output format and format-options [SAM].\n");
    fprintf(fp, "\n");
    fprintf(fp,
"Standard htslib format options apply. So to create a CRAM file with lossy\n\
template names enabled and a larger number of sequences per slice, try:\n\
\n\
    capmq -O cram,lossy_names,seqs_per_slice=100000\n\
\n");
}

/*
 * Parse and validate command line arguments
 */
static opts_t *parse_args(int argc, char **argv)
{
    htsFormat in_fmt = {0};
    htsFormat out_fmt = {0};
    int opt;

    opts_t* opts = calloc(sizeof(opts_t), 1);
    if (!opts) { perror("cannot allocate option parsing memory"); return NULL; }

    opts->rgva = rgva_init(255);
    opts->argv_list = stringify_argv(argc, argv);
    opts->storeQ = true;
    opts->minQ = 0;
    opts->capQ = 255;

    // a bit hacky, but I need to know if -f is in effect before parsing -g or -G or -C
    if (strstr(opts->argv_list,"-f")) opts->freemix = true;

    while ((opt = getopt(argc, argv, "m:g:G:I:O:C:sSrhvf")) != -1) {
    switch (opt) {
        case 'I': hts_parse_format(&in_fmt, optarg);
                  break;

        case 'O': hts_parse_format(&out_fmt, optarg);
                  break;

        case 'C': opts->capQ = capq_from_str(optarg, opts->freemix);
                  break;

        case 's': opts->storeQ = true;
                  break;

        case 'S': opts->storeQ = false;
                  break;

        case 'r': opts->restoreQ = true;
                  break;

        case 'f': opts->freemix = true;
                  break;

        case 'v': opts->verbose = true;
                  break;

        case 'g': parse_rgv(opts->rgva, optarg, opts->freemix);
                  break;

        case 'G': parse_gfile(optarg,opts);
                  break;

        case 'm': opts->minQ = uint8_from_str(optarg);
                  break;

        case 'h': usage(stdout);
                  return 0;

        default:
            fprintf(stderr,"Unknown option: '%c'\n", opt);
            usage(stderr);
            return NULL;
        }
    }

    if (opts->capQ == 255 && !opts->restoreQ && !opts->rgva->end) {
        fprintf(stderr, "Nothing to do!\n");
        return NULL;
    }

    char *fnin = optind < argc ? argv[optind++] : "-";
    if (!(opts->in = sam_open_format(fnin, "r", &in_fmt))) {
        perror(argv[optind]);
        return NULL;
    }

    char mode[5] = "w";
    char *fnout = optind < argc ? argv[optind++] : "-";
    sam_open_mode(mode+1, fnout, NULL);

    if (!(opts->out = sam_open_format(fnout, mode, &out_fmt))) {
        perror("(stdout)");
        return NULL;
    }

    rgva_sort(opts->rgva);

    if (opts->freemix) {
      if (opts->capQ < opts->minQ) {
          if (opts->verbose) {
              fprintf(stderr, 
                      "Default mapping quality cap calculated from freemix (%d) "
                      "was lower than the minimum specifed by `-m' (%d), using "
                      "the latter as the default mapping quality cap.\n", 
                      opts->capQ, opts->minQ);
          }
          opts->capQ = opts->minQ;
      }
      if (opts->rgva) {
          int n;
          for (n=0; n < opts->rgva->end; n++) {
              rgv_t *rgv = opts->rgva->rgv[n];
              if (rgv->capQ < opts->minQ) {
                  if (opts->verbose) {
                      fprintf(stderr, 
                              "Mapping quality cap calculated from freemix (%d) "
                              "for read group `%s' was lower than the minimum "
                              "specifed by `-m' (%d), using the latter as the "
                              "mapping quality cap for this read group.\n",
                              rgv->capQ, rgv->rg, opts->minQ);
                  }
                  rgv->capQ = opts->minQ;
              }
          }
      }
    }

    return opts;
}

/*
 * Process the file
 */
int capq(opts_t *opts)
{
    bam_hdr_t *header;
    bam1_t *b = NULL;
    int ret;
    uint8_t *om = NULL;
    uint8_t *rg = NULL;
    rgv_t *rgv = NULL;
    int n;

    if (opts->verbose) {
        fprintf(stderr, "Capping mapping qualities of %s to a maximum of %d by default\n", opts->in->fn, opts->capQ);
        if (opts->rgva) {
            for (n=0; n < opts->rgva->end; n++) {
                rgv_t *rgv = opts->rgva->rgv[n];
                fprintf(stderr, "Capping mapping qualities to a maximum of %d for read group %s\n", rgv->capQ, rgv->rg); 
            }
        }
    }

    // read header
    if (!(header = sam_hdr_read(opts->in))) {
        fprintf(stderr, "Failed to read file header\n");
        return 1;
    }

    // Add @PG line to header
    SAM_hdr *sh = sam_hdr_parse_(header->text,header->l_text);
    sam_hdr_add_PG(sh, "capmq",
                   "VN", CAPMQ_VERSION,
                   "CL", opts->argv_list,
                   "DS", "cap map quality values",
                   NULL, NULL);
    sam_hdr_unparse2(sh,header);

    // write new header
    if (opts->out && sam_hdr_write(opts->out, header) != 0) {
        fprintf(stderr, "Failed to write file header\n");
        return 1;
    }

    b = bam_init1();
    if (!b) {
        fprintf(stderr, "Failed to allocate bam struct\n");
        return 1;
    }

    // Loop over each read in the BAM file
    while ((ret = sam_read1(opts->in, header, b)) >= 0) {
        if (b->core.tid >= 0) {
            om = bam_aux_get(b, "om");

            // the restore option overrides everything else
            if (opts->restoreQ) {
                if (om) {
                    b->core.qual = bam_aux2i(om);   // restore quality
                    bam_aux_del(b, om);             // delete om tag
                }
            } else {
                uint8_t capQ = opts->capQ;
                rg = bam_aux_get(b, "RG");
                if (rg) {
                    // handle -g / -G options
                    char *rgtag = bam_aux2Z(rg);
                    rgv = rgva_find(opts->rgva, rgtag);
                    if (rgv) capQ = rgv->capQ;
                }
                if (b->core.qual > capQ) {
                    if (opts->storeQ && !om) {
                        // handle -s option
                        int q = b->core.qual;
                        bam_aux_append(b, "om", 'i', 4, (uint8_t*)&q);
                    }
                    b->core.qual = capQ;
                }
            }
        }
        if (sam_write1(opts->out, header, b) < 0) {
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

    return 0;
}

/*
 * parse arguments and do things with them
 */
int main(int argc, char *argv[])
{
    int ret = 1;
    opts_t* opts = parse_args(argc, argv);
    if (opts) {
        ret = capq(opts);
    }
    free_opts(opts);
    return ret;
}

