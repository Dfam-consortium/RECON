#include "seqlist.h"
#include "msps.h"
#include "minunit.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> // for getopt

// Function prototypes
int GetOutputInterval(char *, int, int);

// Example MinUnit unit test
MU_TEST(test_check) {
    mu_check(1 == 1);
}

// Example MinUnit test suite
MU_TEST_SUITE(test_suite) {
    MU_RUN_TEST(test_check);
}

// Print help message for command-line usage
void print_usage(const char *prog_name) {
    printf("Usage: %s [-s seq_list] [-m msp_file] [-n num_output] [-t] [-h]\n", prog_name);
    printf("Options:\n");
    printf("  -s  Path to sequence name list file (required)\n");
    printf("  -m  Path to MSP file (required)\n");
    printf("  -n  Number of output files (optional, default is 1)\n");
    printf("  -t  Run unit tests\n");
    printf("  -h  Show this help message\n");
}

int main(int argc, char *argv[]) {
    FILE *msp_file = NULL, *seq_list = NULL; /* input files */
    FILE **output, *err, *msp_no; /* output files */
    int noof = 1, i, opt; /* noof == Number Of Output Files */
    char line[150], output_name[50];
    MSP_t cur;
    int img_ct = -1;
    int run_tests = 0; // Flag to indicate running unit tests

    // Parse command-line arguments using getopt
    while ((opt = getopt(argc, argv, "s:m:n:th")) != -1) {
        switch (opt) {
            case 's':
                seq_list = fopen(optarg, "r");
                if (!seq_list) {
                    fprintf(stderr, "Error: Cannot open sequence name list file '%s'.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                msp_file = fopen(optarg, "r");
                if (!msp_file) {
                    fprintf(stderr, "Error: Cannot open MSP file '%s'.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'n':
                noof = atoi(optarg);
                if (noof < 1) {
                    fprintf(stderr, "Error: Number of output files must be greater than 0.\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 't':
                run_tests = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // If -t option is set, run unit tests
    if (run_tests) {
        MU_RUN_SUITE(test_suite);
        MU_REPORT();
        return MU_EXIT_CODE;
    }

    // Ensure mandatory options are provided
    if (!seq_list || !msp_file) {
        fprintf(stderr, "Error: Both sequence name list (-s) and MSP file (-m) are required.\n");
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // Open error and summary files
    err = fopen("images/errors", "w");
    if (!err) {
        fprintf(stderr, "Error: Cannot open 'images/errors' for output.\n");
        exit(EXIT_FAILURE);
    }
    msp_no = fopen("summary/ori_msp_no", "w");
    if (!msp_no) {
        fprintf(stderr, "Error: Cannot open 'summary/ori_msp_no' for output.\n");
        exit(EXIT_FAILURE);
    }

    // Allocate memory for output files and open them
    output = (FILE **)malloc(noof * sizeof(FILE *));
    for (i = 0; i < noof; i++) {
        sprintf(output_name, "images/spread%d", i + 1);
        output[i] = fopen(output_name, "w");
        if (!output[i]) {
            fprintf(stderr, "Error: Cannot open '%s' for output.\n", output_name);
            exit(EXIT_FAILURE);
        }
    }

    // Load sequence names into global variable (GetSeqNames implementation in seqlist.h)
    GetSeqNames(seq_list);

    // Process MSP file and spread data across output files
    while (fgets(line, sizeof(line), msp_file)) {
        if (scan_msp(&cur, line)) {
            fprintf(err, "Wrong format: %s\n", line);
            continue;
        }

        img_ct++;
        i = GetOutputInterval(cur.query.frag.seq_name, seq_no, noof);
        if (i < 0) {
            fprintf(err, "Sequence name %s not found in the sequence list.\n", cur.query.frag.seq_name);
            exit(EXIT_FAILURE);
        }
        fprintf(output[i], "%d %d %s %d %d\n", img_ct, cur.score, cur.query.frag.seq_name, cur.query.frag.lb, cur.query.frag.rb);

        img_ct++;
        i = GetOutputInterval(cur.sbjct.frag.seq_name, seq_no, noof);
        if (i < 0) {
            fprintf(err, "Sequence name %s not found in the sequence list.\n", cur.sbjct.frag.seq_name);
            exit(EXIT_FAILURE);
        }
        fprintf(output[i], "%d %d %s %d %d\n", img_ct, cur.score, cur.sbjct.frag.seq_name, cur.sbjct.frag.lb, cur.sbjct.frag.rb);
    }

    // Write final image count to the summary file
    fprintf(msp_no, "%d\n", img_ct / 2 + 1);

    // Clean up
    fclose(err);
    fclose(msp_no);
    for (i = 0; i < noof; i++) {
        fclose(output[i]);
    }
    free(output);

    exit(EXIT_SUCCESS);
}

int GetOutputInterval(char *seq_name, int seq_no, int noof) {
    int interval = GetSeqIndex(0, seq_no - 1, seq_name) / (seq_no / noof);
    if (interval < 0) {
        printf("Sequence name %s not found in the sequence list.\n", seq_name);
        return -1;
    }
    if (interval > noof - 1) interval = noof - 1; /* Ensure the interval doesn't exceed output files */
    return interval;
}
