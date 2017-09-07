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

typedef int (*content_test_func)(FILE *f, char *expected);

int sam_content_test(FILE *f, char *expected) {
    char res_str[1024];
    size_t n;
    int nPG=0;
    int om[6]={0,0,0,0,0,0};
    int q[6];
    int nlines = 0;
    char *lineptr = NULL;
    while (getline(&lineptr, &n, f) != -1) {
        if (strncmp(lineptr,"@PG",3) == 0) nPG++;;
        if (*lineptr != '@') {  // read line, not header line
            char *s;
            char *col;
            int colnum;
            if (nlines < 6) {
                col = strtok(lineptr, "\t");
                for (colnum=0; colnum < 4; colnum++) {
                    col = strtok(NULL, "\t");
                }
                q[nlines] = col!=NULL ? atoi(col) : -1;
                for (colnum=0; colnum < 7; colnum++) {
                    strtok(NULL, "\t");
                }
                while ((col = strtok(NULL, "\t")) != NULL) {
                  s = strstr(col, "om:i:");
                  om[nlines] = s!=NULL ? atoi(s+5) : -1;
                }
            }
            nlines++;
        }
    }
    free(lineptr);

    sprintf(res_str,"%d %d om[%d,%d,%d,%d,%d,%d] q[%d,%d,%d,%d,%d,%d]",
                    nlines, nPG, om[0], om[1], om[2], om[3], om[4], om[5], q[0], q[1], q[2], q[3], q[4], q[5]);

    if (strcmp(res_str,expected)) {
        fprintf(stderr,"Expected: %s\n",expected);
        fprintf(stderr,"Actual:   %s\n",res_str);
        return 1;
    }

    return 0;
}

int content_contains_test(FILE *f, char *expected) {
    size_t n;
    int nlines = 0;
    char *lineptr = NULL;

    while (getline(&lineptr, &n, f) != -1) {
        nlines++;
        char *s;
        s = strstr(lineptr, expected);
        if (s != NULL) {
            return 0;
        }
    }
    free(lineptr);

    fprintf(stderr, "Expected string not found in any output line: %s\n", expected);
    return 1;
}

int run_test(char *cmd, char *expected, int expected_status, content_test_func test)
{
    if (verbose) fprintf(stderr, "Testing command: %s\n", cmd);

    FILE *f = popen(cmd, "r");
    if (!f) {
        fprintf(stderr,"Can't run command: %s\n", cmd);
        return 1;
    }

    int result = test(f, expected);

    int status = pclose(f);
    if (WIFEXITED(status)) {
      int exit_status = WEXITSTATUS(status);
      if (exit_status != expected_status) {
          fprintf(stderr, "Expected exit status: %d\n", expected_status);
          fprintf(stderr, "Actual exit status: %d\n", exit_status);
          return 1;
      }
    } else {
      fprintf(stderr, "Program did not exit normally\n");
      return 1;
    }

    return result;
}

int main(int argc, char *argv[])
{
    int pass=0, fail=0;

    // this should do nothing except add @PG, because no MQ is over 100
    if (run_test("./capmq -C100 test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[45,46,47,48,4,4]",0,&sam_content_test)) fail++; else pass++;

    // cap all values, and create om tags
    if (run_test("./capmq -C40 test1.sam","6 1 om[45,46,47,48,-1,-1] q[40,40,40,40,4,4]",0,&sam_content_test)) fail++; else pass++;

    // this should cap all mapping qualities at 0
    if (run_test("./capmq -C0 test1.sam","6 1 om[45,46,47,48,4,4] q[0,0,0,0,0,0]",0,&sam_content_test)) fail++; else pass++;

    // cap and restore. End result should be unchanged
    if (run_test("./capmq -C40 test1.sam | ./capmq -r","6 2 om[-1,-1,-1,-1,-1,-1] q[45,46,47,48,4,4]",0,&sam_content_test)) fail++; else pass++;

    // cap and cap and restore. End result should still be unchanged
    if (run_test("./capmq -C41 test1.sam | ./capmq -C 5 | ./capmq -r","6 3 om[-1,-1,-1,-1,-1,-1] q[45,46,47,48,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups, no default cap
    if (run_test("./capmq -gb:41 -ga:40 -gx:42 test1.sam","6 1 om[45,46,47,48,-1,-1] q[40,40,41,41,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups, with default cap
    if (run_test("./capmq -C40 -ga:41 test1.sam","6 1 om[45,46,47,48,-1,-1] q[41,41,40,40,4,4]",0,&sam_content_test)) fail++; else pass++;

    // cap value using freemix
    if (run_test("./capmq -S -C0.00005 -f test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[43,43,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups using freemix
    if (run_test("./capmq -S -gb:0.00005 -f -ga:0.0001 test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[40,40,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups using freemix and -m
    if (run_test("./capmq -m41 -S -gb:0.00005 -f -ga:0.0001 test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[41,41,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups from file
    if (run_test("./capmq -S -f -G test-b-a.txt test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[40,40,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups from file - override `a' with minimum
    if (run_test("./capmq -m41 -S -f -G test-b-a.txt test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[41,41,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // read groups from file - RG `a' not matched
    if (run_test("./capmq -m41 -S -f -G test-b-ai.txt test1.sam","6 1 om[-1,-1,-1,-1,-1,-1] q[45,46,43,43,4,4]",0,&sam_content_test)) fail++; else pass++;

    // this should do nothing and say so
    if (run_test("./capmq test1.sam 2>&1","Nothing to do",1,&content_contains_test)) fail++; else pass++;

    printf("Passed %d tests\n", pass);
    if (fail) printf("FAILED %d tests\n", fail);

    return fail;
}

