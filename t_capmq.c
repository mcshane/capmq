/* The MIT License

    Copyright (C) 2017 Genome Research Ltd.

    Author: Jennifer Liddle <js10@sanger.ac.uk>

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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

bool verbose = true;

int run_test(char *cmd, char *expected)
{
    char *lineptr = NULL;
    char res_str[1024];
    size_t n;
    int nlines=0;
    int nPG=0;
    int om[4]={0,0,0,0};
    int q[4];

    if (verbose) fprintf(stderr, "Testing command: %s\n", cmd);

    FILE *f = popen(cmd, "r");
    if (!f) {
        fprintf(stderr,"Can't run command: %s\n", cmd);
        return 1;
    }

    while (getline(&lineptr, &n, f) != -1) {
        if (strncmp(lineptr,"@PG",3) == 0) nPG++;;
        if (*lineptr != '@') {  // read line, not header line
            char *s;
            if (nlines < 4) {
                s = strstr(lineptr, "\t4");     // look for quality score
                q[nlines] = s ? atoi(s+1) : -1;
                s = strstr(lineptr, "om:i:");
                om[nlines] = s ? atoi(s+5) : -1;
            }
            nlines++;
        }
    }
    free(lineptr);

    pclose(f);
    sprintf(res_str,"%d %d om[%d,%d,%d,%d] q[%d,%d,%d,%d]",
                     nlines, nPG, om[0], om[1], om[2], om[3], q[0], q[1], q[2], q[3]);

    if (strcmp(res_str,expected)) {
        fprintf(stderr,"Expected: %s\n",expected);
        fprintf(stderr,"Actual:   %s\n",res_str);
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int pass=0, fail=0;

    // this should do nothing except add @PG, because no MQ is over 100
    if (run_test("./capmq -C100 test1.sam","4 1 om[-1,-1,-1,-1] q[45,46,47,48]")) fail++; else pass++;

    // cap all values, and create om tags
    if (run_test("./capmq -C40 -s test1.sam","4 1 om[45,46,47,48] q[40,40,40,40]")) fail++; else pass++;

    // cap and restore. End result should be unchanged 
    if (run_test("./capmq -C40 -s test1.sam | capmq -r","4 2 om[-1,-1,-1,-1] q[45,46,47,48]")) fail++; else pass++;

    // cap and cap and restore. End result should still be unchanged 
    if (run_test("./capmq -C41 -s test1.sam | capmq -C 5 -s | capmq -r","4 3 om[-1,-1,-1,-1] q[45,46,47,48]")) fail++; else pass++;

    // read groups, no default cap
    if (run_test("./capmq -gb:41 -ga:40 -gx:42 -s test1.sam","4 1 om[45,46,47,48] q[40,40,41,41]")) fail++; else pass++;

    // read groups, with default cap
    if (run_test("./capmq -C40 -ga:41 -s test1.sam","4 1 om[45,46,47,48] q[41,41,40,40]")) fail++; else pass++;

    // cap value using freemix
    if (run_test("./capmq -C0.00005 -f test1.sam","4 1 om[-1,-1,-1,-1] q[43,43,43,43]")) fail++; else pass++;

    // read groups using freemix
    if (run_test("./capmq -gb:0.00005 -f -ga:0.0001 test1.sam","4 1 om[-1,-1,-1,-1] q[40,40,43,43]")) fail++; else pass++;

    // read groups using freemix and -m
    if (run_test("./capmq -m41 -gb:0.00005 -f -ga:0.0001 test1.sam","4 1 om[-1,-1,-1,-1] q[41,41,43,43]")) fail++; else pass++;

    // read groups from file
    if (run_test("./capmq -m41 -f -G test1.txt test1.sam","4 1 om[-1,-1,-1,-1] q[41,41,43,43]")) fail++; else pass++;

    printf("Passed %d tests\n", pass);
    if (fail) printf("FAILED %d tests\n", fail);

    return fail;
}

