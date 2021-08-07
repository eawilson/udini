#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <getopt.h>
#include <ctype.h>
#include <errno.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "zlib.h"


const int NAME = 0;
const int SEQ = 1;
const int PLUS = 2;
const int QUAL = 3;
const unsigned GZBUFFER_LEN = 1024 * 1024;
const size_t LINE_LEN = 1024;


bool endswith(const char *text, const char *suffix);


int main (int argc, char **argv) {
    /*
     * 
     */
    const char *output_filename = "-", *stats_filename = "stats.json", *umi_sequences = NULL;
    char **input_filenames = NULL, *umi_buffer = NULL, **valid_umis = NULL, *lines[2][4] = {NULL}, *umis[2] = {NULL}, ch = '\0';
    bool interleaved = false, whitespace = true, validated_umis = true, invalid_short_read = false, invalid_n_read = false;
    int umi_length = 0, umi_stem_length = 0, min_read_length = 50;
    int lineno = 0, readno = 0, consecutive_ns = 0, max_consecutive_ns = 2, i = 0;
    int best = 0, nextbest = 0, edit_distance = 0, j = 0;
    size_t input_filenames_len = 0, umi_buffer_len = 0, valid_umis_len = 0, current_input_filename = 0;
    size_t total_reads = 0, invalid_umi_reads[2] = {0}, invalid_short_reads = 0, invalid_n_reads = 0, length = 0;
    size_t line_lens[2][4] = {0}, seq_len = 0, qual_len = 0;
    gzFile output_fastq = NULL, input_fastqs[2] = {NULL};
    FILE *stats_file = NULL;

    // variable needed by strtol
    char *endptr = NULL;
    long val = 0;
    // variables needed by getopt_long
    int option_index = 0, c = 0;
    static struct option long_options[] = {{"output", required_argument, 0, 'o'},
                                           {"umi", required_argument, 0, 'u'},
                                           {"umi-length", required_argument, 0, 'l'},
                                           {"umi-stem-length", required_argument, 0, 'k'},
                                           {"umi-sequences", required_argument, 0, 'q'},
                                           {"min-read-length", required_argument, 0, 'm'},
                                           {"max-consecutive-ns", required_argument, 0, 'n'},
                                           {"stats", required_argument, 0, 's'},
                                           {"interleaved", no_argument, 0, 'i'},
                                           {"help", no_argument, 0, 'h'},
                                           {0, 0, 0, 0}};

    // Parse optional arguments
    while (c != -1) {
        c = getopt_long(argc, argv, "o:u:l:k:q:m:n:s:ih", long_options, &option_index);

        switch (c) {
            case 'o':
                if ((!endswith(optarg, ".fastq")) && !endswith(optarg, ".fastq.gz") && (strcmp(optarg, "-") != 0)) {
                    fprintf(stderr, "Error: Output file must be of type fastq or fastq.gz\n");
                    exit(EXIT_FAILURE);
                    }
                output_filename = optarg;
                break;
                
            case 'u':
                if (strcmp(optarg, "thruplex") == 0) {
                    umi_length = 6;
                    umi_stem_length = 11;
                    }
                else if (strcmp(optarg, "thruplex_hv") == 0) {
                    umi_length = 7;
                    umi_stem_length = 1;
                    umi_sequences = "AAGCTGA ACAACGA ACTCGTA ATTGCTC CGAGTAC CGCTAAT "
                                    "CTAGTAG GACATCG GTCTCTG TACCTCA TCTGGTA TGAACGG "
                                    "ACGACTC ATCTGGA CAATAGC CCTAGGT CGTCTCA GAGTCTC "
                                    "GGCAATG TCCACTA TCTCCAT TGTCAAC TGTGTCT TTGTAGT";
                    }
                else if (strcmp(optarg, "prism") == 0) {
                    umi_length = 8;
                    umi_stem_length = 0;
                    umi_sequences = "GAGACGAT GCACAACT TTCCAAGG GCGTCATT CGCATGAT "
                                    "GAAGGAAG ACGGAACA ACTGAGGT CGGCTAAT TGAAGACG "
                                    "GCTATCCT GTTACGCA TGGACTCT AGCGTGTT ATCCAGAG "
                                    "GATCGAGT CTTAGGAC TTGCGAAG GTGCCATA CTGTTGAC "
                                    "TCGCTGTT GATGTGTG TTCGTTGG ACGTTCAG AAGCACTG "
                                    "TTGCAGAC GTCGAAGA CAATGTGG ACCACGAT ACGACTTG "
                                    "GATTACCG ACTAGGAG";
                    }
                else {
                    fprintf(stderr, "Error: Unsupported umi type: %s\n", optarg);
                    exit(EXIT_FAILURE);
                    }
                break;
                           
            case 'l':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg) || (val > INT_MAX) || (val < 0)) {
                    fprintf(stderr, "Error: Invalid --umi-length\n");
                    exit(EXIT_FAILURE);
                    }
                umi_length = (int)val;
                break;
           
            case 'k':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg) || (val > INT_MAX) || (val < 0)) {
                    fprintf(stderr, "Error: Invalid --umi-length\n");
                    exit(EXIT_FAILURE);
                    }
                umi_stem_length = (int)val;
                break;

            case 'q':
                for (i = 0; optarg[i] != '\0'; ++i) {
                    if (optarg[i] != 'A' && optarg[i] != 'C' && optarg[i] != 'G' && optarg[i] != 'T' && optarg[i] != ' ') {
                        fprintf(stderr, "Error: Invalid --umi-sequences\n");
                        exit(EXIT_FAILURE);
                       }
                    }
                umi_sequences = optarg;
                break;
           
            case 'm':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg) || (val > INT_MAX) || (val < 0)) {
                    fprintf(stderr, "Error: Invalid --min-read-length\n");
                    exit(EXIT_FAILURE);
                    }
                min_read_length = (int)val;
                break;
           
            case 'n':
                errno = 0;
                val = strtol(optarg, &endptr, 10);
                if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) || (endptr == optarg) || (val > INT_MAX) || (val < 0)) {
                    fprintf(stderr, "Error: Invalid --max-consecutive-ns\n");
                    exit(EXIT_FAILURE);
                    }
                max_consecutive_ns = (int)val;
                break;
               
            case 's':
                if (!endswith(optarg, ".json")) {
                    fprintf(stderr, "Error: Stats file must be of type json\n");
                    exit(EXIT_FAILURE);
                    }
                stats_filename = optarg;
                break;
            
            case 'i':
                interleaved = true;
                break;
                
            case 'h':
                // help message
                fprintf(stderr, "Program: udini\n");
                fprintf(stderr, "Version: 2.0.0\n");
                fprintf(stderr, "Usage:   udini [options] input_fastq(s)\n");
                fprintf(stderr, "Options: -o, --output STR              Name of output fastq/fastq.gz or '-' for \n");
                fprintf(stderr, "                                       stdout (default '-')\n");
                fprintf(stderr, "         -u, --umi STR                 Known UMI type (one of 'thruplex', \n");
                fprintf(stderr, "                                       'thruplex_hv' or 'prism')\n");
                fprintf(stderr, "         -l, --umi-length INT          UMI length\n");
                fprintf(stderr, "         -k, --umi-stem-length INT     UMI stem length\n");
                fprintf(stderr, "         -q, --umi-sequences STR       Whitespace separated list of valid UMI \n");
                fprintf(stderr, "                                       sequences\n");
                fprintf(stderr, "         -m, --min-read-length INT     Filter all read pairs with either read \n");
                fprintf(stderr, "                                       shorter than min length (default 50)\n");
                fprintf(stderr, "         -n, --max-consecutive-ns INT  Filter all read pairs with either read \n");
                fprintf(stderr, "                                       containing a greater number of \n");
                fprintf(stderr, "                                       consecutive Ns (default 2)\n");
                fprintf(stderr, "         -s, --stats STR               Name of output statistics file \n");
                fprintf(stderr, "                                       (default 'stats.json')\n");
                fprintf(stderr, "         -i, --interleaved             Each input fastq contains alternating \n");
                fprintf(stderr, "                                       R1/R2 reads\n");
                fprintf(stderr, "         -h, --help                    Display this message\n");
                exit(EXIT_SUCCESS);
                
            case '?':
                // unknown option, getopt_long already printed an error message.
                exit(EXIT_FAILURE);
            }
        }
    
    input_filenames_len = argc - optind;
    input_filenames = argv + optind;
    if (input_filenames_len == 0) {
        fprintf(stderr, "Error: No input file supplied\n");
        exit(EXIT_FAILURE);
        }
    if ((input_filenames_len % 2) && !interleaved) {
        fprintf(stderr, "Error: Must be an even number of alternating R1/R2 input fastqs unless --interleaved is selected\n");
        exit(EXIT_FAILURE);
        }
    for (i = 0; i <input_filenames_len; ++i) {
        if (strcmp(input_filenames[i], "-") == 0) {
            if (input_filenames_len > 1) {
                fprintf(stderr, "Error: Stdin cannot be mixed with other inputs\n");
                exit(EXIT_FAILURE);
                }
            }
        else if (!endswith(input_filenames[i], ".fastq") && !endswith(input_filenames[i], ".fastq.gz")) {
            fprintf(stderr, "Error: Input file must be of type fastq or fastq.gz\n");
            exit(EXIT_FAILURE);
            }
        }
    
    
    if (umi_sequences) {
        umi_buffer_len = strlen(umi_sequences);
        if ((umi_buffer = malloc(umi_buffer_len + 1)) == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for umi sequnces buffer\n");
            exit(EXIT_FAILURE);
            }
        memcpy(umi_buffer, umi_sequences, umi_buffer_len + 1);
        
        // We will need less pointers than this but given the low number any additonal work to calculate the exact number is excessive
        if ((valid_umis = calloc(umi_buffer_len, sizeof(char *))) == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for valid umis\n");
            exit(EXIT_FAILURE);
            }
        
        for (i = 0; i < umi_buffer_len; ++i) {
            if (isspace(umi_buffer[i])) {
                umi_buffer[i] = '\0';
                whitespace = true;
                }
            else if (whitespace) {
                valid_umis[valid_umis_len++] = umi_buffer + i;
                whitespace = false;
                }
            }
        
        for (i = 0; i < valid_umis_len; ++i) {
            if (!umi_length) {
                umi_length = strlen(valid_umis[i]);
                }
            else if (strlen(valid_umis[i]) != umi_length) {
                fprintf(stderr, "Error: Not all UMI sequences are %i nt long\n", umi_length);
                exit(EXIT_FAILURE);
                }
            }
        }
    
    
    // Make sure we have a true min line length to prevent buffer overflows when removing umis
    if (min_read_length < 2 * (umi_length + umi_stem_length) + 1) {
        min_read_length = 2 * (umi_length + umi_stem_length) + 1;
        }
    
    
    if (strcmp(output_filename, "-") == 0) {
        if ((output_fastq = gzdopen(dup(fileno(stdout)), "wbT")) == NULL) {
            fprintf(stderr, "Error: Unable to open stdout for writing\n");
            exit(EXIT_FAILURE);
            }
        }
    else if (endswith(output_filename, ".fastq")) {
        if ((output_fastq = gzopen(output_filename, "wbT")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s for writing\n", input_filenames[current_input_filename]);
            exit(EXIT_FAILURE);
            }
        }
    else { // fastq.gz
        if ((output_fastq = gzopen(output_filename, "wb")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s for writing\n", input_filenames[current_input_filename]);
            exit(EXIT_FAILURE);
            }
        if (gzbuffer(output_fastq, GZBUFFER_LEN) == -1) {
            fprintf(stderr, "Error: Unable to set gzbuffer size\n");
            exit(EXIT_FAILURE);
            }
        }
    
    
    for (i = 0; i < 8; ++i) {
        line_lens[i / 4][i % 4] = LINE_LEN;
        if ((lines[i / 4][i % 4] = malloc(line_lens[i / 4][i % 4] + 1)) == NULL) {
            fprintf(stderr, "Error: Unable to allocate memory for line buffer\n");
            exit(EXIT_FAILURE);
            }
        }
    
    
    for (current_input_filename = 0; current_input_filename < input_filenames_len; ++current_input_filename) {
        
        for (readno = 0;; ++readno) {
            if (strcmp(input_filenames[current_input_filename], "-") == 0) {
                if ((input_fastqs[readno] = gzdopen(dup(fileno(stdin)), "rb")) == NULL) {
                    fprintf(stderr, "Error: Unable to open stdin for reading\n");
                    exit(EXIT_FAILURE);
                    }
                }
            else {
                if ((input_fastqs[readno] = gzopen(input_filenames[current_input_filename], "rb")) == NULL) {
                    fprintf(stderr, "Error: Unable to open %s for reading\n", input_filenames[current_input_filename]);
                    exit(EXIT_FAILURE);
                    }
                }
            if (gzbuffer(input_fastqs[readno], GZBUFFER_LEN) == -1) {
                fprintf(stderr, "Error: Unable to set gzbuffer size\n");
                exit(EXIT_FAILURE);
                }
            
            if (readno == 1) {
                break;
                }
        
            if (interleaved) {
                input_fastqs[1] = input_fastqs[0];
                break;
                }
            else {
                ++current_input_filename;
                }
            }
        
        
        for (;; ++total_reads) {
            for (i = 0; i < 8; ++i) {
                readno = i / 4;
                lineno = i % 4;

                length = 0;
                while (true) {
                    if (gzgets(input_fastqs[readno], lines[readno][lineno] + length, line_lens[readno][lineno] - length) == NULL) {
                        if (gzeof(input_fastqs[readno])) {
                            if (i == 0 && gzgetc(input_fastqs[1]) == -1 && gzeof(input_fastqs[1])) {
                                i = 8;
                                break;
                                }
                            fprintf(stderr, "Error: Truncated fastq\n");
                            }
                        else {
                            fprintf(stderr, "Error: Unable to read fastq\n");
                            }
                        exit(EXIT_FAILURE);
                        
                        }
                    length += strlen(lines[readno][lineno] + length);
                    if (lines[readno][lineno][length - 1] == '\n') {
                        lines[readno][lineno][--length] = '\0';
                        if (lines[readno][lineno][length - 1] == '\r') {
                            lines[readno][lineno][--length] = '\0';
                            }
                        break;
                        }
                    
                    line_lens[readno][lineno] += LINE_LEN;
                    if ((lines[readno][lineno] = realloc(lines[readno][lineno], line_lens[readno][lineno] + 1)) == NULL) {
                        fprintf(stderr, "Error: Unable to reallocate memory for line buffer.\n");
                        exit(EXIT_FAILURE);
                        }
                    }
                }
            
            if (i == 9) {
                for (readno = 0; readno < 2; ++readno) {
                    if (gzclose(input_fastqs[readno]) != Z_OK) {
                        fprintf(stderr, "Error: Unable to successfully close fastq\n");
                        exit(EXIT_FAILURE);
                        }
                    if (interleaved) {
                        break;
                        }
                    }
                break;
                }
            
            
            invalid_n_read = false;
            invalid_short_read = false;
            for (readno = 0; readno < 2; ++readno) {
                seq_len = strlen(lines[readno][SEQ]);
                qual_len = strlen(lines[readno][QUAL]);
                if (seq_len != qual_len) {
                    fprintf(stderr, "Error: Sequence and quality differ in length\n");
                    exit(EXIT_FAILURE);
                    }
                
                if (seq_len < min_read_length) {
                    invalid_short_read = true;
                    }
                
                if (max_consecutive_ns) {
                    consecutive_ns = 0;
                    for (i = 0; i < seq_len; ++i) {
                        if (lines[readno][SEQ][i] == 'N') {
                            if (++consecutive_ns > max_consecutive_ns) {
                                invalid_n_read = true;
                                break;
                                }
                            }
                        else {
                            consecutive_ns = 0;
                            }
                        }
                    }        
                }
                
            if (invalid_n_read || invalid_short_read) {
                invalid_n_reads += invalid_n_read;
                invalid_short_reads += invalid_short_read;
                continue;
                }
            
            
            // Remove any whitespace and fastq comment after name
            for (readno = 0; readno < 2; ++readno) {
                for (i = 0; lines[readno][NAME][i] != '\0'; ++i) {
                    if (isspace(lines[readno][NAME][i])) {
                        lines[readno][NAME][i] = '\0';
                        break;
                        }
                    }
                }
            
            if (strcmp(lines[0][NAME], lines[1][NAME]) != 0) {
                fprintf(stderr, "Error: Mismatched paired reads, names don't match\n");
                exit(EXIT_FAILURE);
                }
            
            
            if (umi_length) {
                validated_umis = true;
                for (readno = 0; readno < 2; ++readno) {
                    if (umi_sequences) {
                        best = umi_length;
                        nextbest = umi_length;
                        for (j = 0; j < valid_umis_len; ++j) {
                            edit_distance = 0;
                            for (i = 0; i < umi_length; ++i) {
                                if (lines[readno][SEQ][i] != valid_umis[j][i]) {
                                    if (++edit_distance == nextbest) {
                                        break;
                                        }
                                    }
                                }
                            
                            if (edit_distance < best) {
                                nextbest = best;
                                best = edit_distance;
                                umis[readno] = valid_umis[j];
                                if (edit_distance == 0) {
                                    break;
                                    }
                                }
                            else if (edit_distance < nextbest) {
                                nextbest = edit_distance;
                                }
                            }
                        
                        if (best > 1 || (best == 1 && nextbest < 3)) {
                            ++invalid_umi_reads[readno];
                            validated_umis = false;
                            }
                        }
                    else {
                        umis[readno] = lines[readno][SEQ];
                        }
                    }
                
                if (validated_umis) {
                    for (readno = 0; readno < 2; ++readno) {
                        if (gzprintf(output_fastq, "%s RX:Z:%.*s-%.*s\tQX:Z:%.*s %.*s\n%s\n+\n%s\n",
                                                   lines[readno][NAME],
                                                   umi_length,
                                                   umis[0],
                                                   umi_length,
                                                   umis[1],
                                                   umi_length,
                                                   lines[0][QUAL],
                                                   umi_length,
                                                   lines[1][QUAL],
                                                   lines[readno][SEQ] + umi_length + umi_stem_length,
                                                   lines[readno][QUAL] + umi_length + umi_stem_length) < 0) {
                            fprintf(stderr, "Error: Unable to write to output fastq\n");
                            exit(EXIT_FAILURE);
                            }
                        }
                    }
                }
            else {
                for (readno = 0; readno < 2; ++readno) {
                    if (gzprintf(output_fastq, "%s\n%s\n+\n%s\n",
                                               lines[readno][NAME],
                                               lines[readno][SEQ] + umi_stem_length,
                                               lines[readno][QUAL] + umi_stem_length) < 0) {
                        fprintf(stderr, "Error: Unable to write to output fastq\n");
                        exit(EXIT_FAILURE);
                        }
                    }      
                }
            }
        }        
    
    
    if (gzclose(output_fastq) != Z_OK) {
        fprintf(stderr, "Error: Unable to successfully close output fastq\n");
        exit(EXIT_FAILURE);
        }
    

    if ((stats_file = fopen(stats_filename, "r+")) == NULL) {
        if ((stats_file = fopen(stats_filename, "w")) == NULL) {
            fprintf(stderr, "Error: Unable to open %s for writing\n", stats_filename);
            exit(EXIT_FAILURE);
            }
        fprintf(stats_file, "{\n");
        }
    else {
        fseek(stats_file, 0, SEEK_END);
        if (ftell(stats_file) == 0) {
            fprintf(stats_file, "{}\n");
            }
        
        fseek(stats_file, -1, SEEK_CUR);
        while (true) {
            ch = fgetc(stats_file);
            if (ch != '}' && !isspace(ch)) {
                if (ch != '{') {
                    fprintf(stats_file, ",");
                    }
                break;
                }
            if (fseek(stats_file, -2, SEEK_CUR) == -1) {
                fprintf(stderr, "Error: Malformed stats file\n");
                exit(EXIT_FAILURE);
                }
            }
        }

    fprintf(stats_file, "\n    \"total_reads\": %zu,\n", total_reads);
    fprintf(stats_file, "    \"invalid_short_reads\": %zu,\n", invalid_short_reads);
    fprintf(stats_file, "    \"invalid_n_reads\": %zu", invalid_n_reads);
    if (umi_length) {
        fprintf(stats_file, ",\n    \"invalid_umi_reads_r1\": %zu,\n", invalid_umi_reads[0]);
        fprintf(stats_file, "    \"invalid_umi_reads_r2\": %zu", invalid_umi_reads[1]);
        }
    fprintf(stats_file, "\n}\n");

    if (fclose(stats_file) == EOF) {
        fprintf(stderr, "Error: Unable to close stats file\n");
        exit(EXIT_FAILURE);
        }
    
    free(umi_buffer);
    for (i = 0; i <8; ++i) {
        free(lines[i / 4][i % 4]);
        }
    }



bool endswith(const char *text, const char *suffix) {
    int offset = strlen(text) - strlen(suffix);
    return (offset >= 0 && strcmp(text + offset, suffix) == 0);
    }


