#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "util_lib.h"

// which corresponds to either 0 (done) or 1 (empty) buffer
void enqueue(ThreadData *td, unsigned int buffer, unsigned int which) {

    QUEUE *item = (QUEUE *) malloc(sizeof(QUEUE));
    item->td = td;
    item->buffer = buffer;
    item->next = 0;

    QUEUE **head = 0;
    QUEUE **tail = 0;

    if (which) {
        head = &q_empty_head;
        tail = &q_empty_tail;
    } else {
        head = &q_done_head;
        tail = &q_done_tail;
    }

//    printf("INSIDE ENQ, head %d tail %d\n", *head, *tail);
    if (!*head) {
        *head = item;
        *tail = item;
    } else {
        (*tail)->next = item;
        *tail = item;
    }
}

void printq(unsigned int which) {
    QUEUE *head;

    if (which) head = q_empty_head;
    else head = q_done_head;

//    printf("QUEUE %d head %d\n", which, head);
    while (head) {
        printf("QUEUE %d has id %d buffer %d\n", which, head->td->id, head->buffer);
        head = head->next;
    }

}

// which corresponds to either 0 (done) or 1 (empty) buffer
QUEUE *deq(unsigned int which) {

    QUEUE **head = 0;
    if (which) head = &q_empty_head;
    else head = &q_done_head;

    if (!*head) return 0;

    QUEUE *temp = *head;
    *head = (*head)->next;
    return temp;
}

// which corresponds to either 0 (done) or 1 (empty) buffer
void cutnpaste_q(QUEUE **dest, unsigned int which) {

    QUEUE **head = 0;

    if (which) {
        head = &q_empty_head;
    } else {
        head = &q_done_head;
    }

    if (!*head) {
        *dest = 0;
        return;
    }

    *dest = *head;
    *head = 0;
}

/* This is the viterbi score matrix for this sequence. */
double **dmatrix(int num_col) {
    int i;
    double **m = calloc(NUM_STATE, sizeof(double *));

    if (!m) {
        fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");
        exit(EXIT_FAILURE);
    }

    for (i=0; i < NUM_STATE; i++) {
        m[i] = calloc(num_col, sizeof(double));
        if (!m[i]) {
            fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
            exit(EXIT_FAILURE);
        }
    }
    return m;
}


/* NUM_STATE = NUM_STATE. This is the viterbi score matrix for this sequence. */
int **imatrix(int num_col) {
    int i;
    int **m = calloc(NUM_STATE, sizeof(int *));

    if (!m) {
        fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in imatrix()");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < NUM_STATE; i++) {
        m[i] = calloc(num_col, sizeof(int));
        if (!m[i]) {
            fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i ," in imatrix()");
            exit(EXIT_FAILURE);
        }
    }
    return m;
}

int *ivector(int nh) {
    int *v = calloc(nh, sizeof(int));

    if (!v) {
        fprintf(stderr, "%s\n", "ERROR: Allocation failure in ivector()");
        exit(EXIT_FAILURE);
    }

    return v;
}

void free_dmatrix(double **m) {
    int i;

    for (i=NUM_STATE-1; i>=0; i--) {
        free(m[i]);
        m[i] = NULL;
    }
    free(m);
    m = NULL;
}


void free_imatrix(int **m) {
    int i;

    for (i=NUM_STATE-1; i>=0; i--) {
        free(m[i]);
        m[i] = NULL;
    }
    free(m);
    m = NULL;
}


HMM_StateTransition hmm_state_transition_parse (const char *tr) {
    if (strcmp(tr, "MM") == 0)
        return TR_MM;
    if (strcmp(tr, "MI") == 0)
        return TR_MI;
    if (strcmp(tr, "MD") == 0)
        return TR_MD;
    if (strcmp(tr, "II") == 0)
        return TR_II;
    if (strcmp(tr, "IM") == 0)
        return TR_IM;
    if (strcmp(tr, "DD") == 0)
        return TR_DD;
    if (strcmp(tr, "DM") == 0)
        return TR_DM;
    if (strcmp(tr, "GE") == 0)
        return TR_GE;
    if (strcmp(tr, "GG") == 0)
        return TR_GG;
    if (strcmp(tr, "ER") == 0)
        return TR_ER;
    if (strcmp(tr, "RS") == 0)
        return TR_RS;
    if (strcmp(tr, "RR") == 0)
        return TR_RR;
    if (strcmp(tr, "ES") == 0)
        return TR_ES;
    if (strcmp(tr, "ES1") == 0)
        return TR_ES;

    return -1;
}

Nucleotide nucleotide_parse (char nt) {
    if (nt == 'A' || nt == 'a')
        return NUCL_A;
    if (nt == 'C' || nt == 'c')
        return NUCL_C;
    if (nt == 'G' || nt == 'g')
        return NUCL_G;
    if (nt == 'T' || nt == 't')
        return NUCL_T;

    return NUCL_INVALID;
}

Nucleotide nucleotide_complement (Nucleotide nt) {
    if (nt == NUCL_INVALID)
        return NUCL_INVALID;

    return 3 - nt;
}

int trinucleotide (Nucleotide a, Nucleotide b, Nucleotide c) {
    if (c == NUCL_INVALID)
        return 0;
    if (b == NUCL_INVALID)
        return c;
    if (a == NUCL_INVALID)
        return (b << 2) + c;

    return (a << 4) + (b << 2) + c;
}

static int trinucleotide_pep (Nucleotide a, Nucleotide b, Nucleotide c) {
    if (a == NUCL_INVALID || b == NUCL_INVALID || c == NUCL_INVALID)
        return 64;

    return (a << 4) + (b << 2) + c;
}

// Calculates the reverse complement of the bases in `dna` and saves it in `reverse_complement`
void get_rc_dna(const Nucleotide dna[], int dna_len, char *reverse_complement) {
    static char NUCL_CHARACTERS_RC[] = "TGCAN";
    int i;

    for (i = 0; i < dna_len; i++)
        reverse_complement[dna_len-i-1] = NUCL_CHARACTERS_RC[dna[i]];
}

void get_protein(const Nucleotide dna[], int dna_len, char *protein, Strand strand) {
    int i;

    /* If our DNA sequence is not a multiple of 3, throw away the last nucleotides */
    dna_len -= dna_len % 3;

    if (strand == FORWARD_STRAND) {
        for (i = 0; i < dna_len; i += 3)
            protein[i/3] = translation_table[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
    } else {
        for (i = 0; i < dna_len; i += 3)
            protein[(dna_len-i)/3-1] = translation_table_rc[trinucleotide_pep(dna[i], dna[i+1], dna[i+2])];
    }

    /* Don't forget the string terminator */
    protein[(dna_len/3)] = '\0';
}

void print_usage() {

    puts("USAGE: FGS++ -s [seq_file_name] -m [max_mem_use] -o [output_file_name] -w [1 or 0] -t [train_file_name] -p [thread_num] -e [1 or 0] -d [1 or 0] ");
    puts("EXAMPLE USAGE: FGS++ -s example/NC_000913-454.fna -o output -w 0 -t 454_5 -p 16 ");
    puts("MINIMAL USAGE: FGS++ -s [seq_file_name] -o [output_file_name] -w [1 or 0] -t [train_file_name] ");
    puts("INFO: FragGeneScan++ will only output the amino acid files by default. To obtain the meta information set -e 1 and for the DNA files set -d 1\n");
    puts("    Mandatory parameters:");
    puts("       -s [seq_file_name]:    sequence file name including the full path");
    puts("                              for standard input specify as stdin");
    puts("       -o [output_file_name]: output file name including the full path");
    puts("                              for standard output specify as stdout");
    puts("       -w [1 or 0]:           1 if the sequence file has complete genomic sequences");
    puts("       		                0 if the sequence file has short sequence reads");
    puts("       -t [train_file_name]:  file name that contains model parameters; this file should be in the -r directory");
    puts("                           Note that four files containing model parameters already exist in the \"train\" directory");
    puts("                           [complete] for complete genomic sequences or short sequence reads without sequencing error");
    puts("                           [sanger_5] for Sanger sequencing reads with about 0.5% error rate");
    puts("                           [sanger_10] for Sanger sequencing reads with about 1% error rate");
    puts("                           [454_5] for 454 pyrosequencing reads with about 0.5% error rate");
    puts("                           [454_10] for 454 pyrosequencing reads with about 1% error rate");
    puts("                           [454_30] for 454 pyrosequencing reads with about 3% error rate");
    puts("                           [illumina_1] for Illumina sequencing reads with about 0.1% error rate");
    puts("                           [illumina_5] for Illumina sequencing reads with about 0.5% error rate");
    puts("                           [illumina_10] for Illumina sequencing reads with about 1% error rate\n");
    puts("    Optional flags");
    puts("       -r [train_file_dir]:    Full path of the directory containing the training model files.");
    puts("       -p [thread_num]         The number of threads used by FragGeneScan++; default is 1 thread.");
    puts("       -e [1 or 0]             Output metadata for sequences.");
    puts("       -d [1 or 0]             Output DNA file.");
    puts("       -m [max_mem_usage]      Maximum amount of memory to be used by the application, in megabytes, default 1024 for 1GB");
    puts("       -x [translation_table]  Which translation table to use (default: 11)");
}

void stopMemset(char *ptr, int length) {
    int i;
    for (i=0; i<length; i++) {
        if (ptr[i] == '\0') {
            return;
        }
        ptr[i] = '\0';
    }
}
