/**
 * @file hmm.h
 *
 * @brief The main header file.
 */

#ifndef __HMM_H
#define __HMM_H

#include <ctype.h>
#include <ctype.h>
#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fasta.h"

/** The queue for worker threads that are waiting for the writer thread */
#define DONE_Q 0
/** The queue for worker threads that are waiting for input */
#define EMPTY_Q 1

/* Some precomputed mathematical constants */
#define LOG_53 -0.63487827243
#define LOG_16 -1.83258146375
#define LOG_30 -1.20397280433
#define LOG_25 -1.38629436112
#define LOG_83 -0.18632957819
#define LOG_10 -2.30258509299
#define LOG_07 -2.65926003693
#define LOG_95 -0.05129329438

/** A nucleotide */
typedef enum {
    /** Adenine (A) */
    NUCL_A,
    /** Cytosine (C) */
    NUCL_C,
    /** Guanine (G) */
    NUCL_G,
    /** Thymine (T) */
    NUCL_T,

    NUCL_INVALID
} Nucleotide;

/** The amount of possible states */
#define NUM_STATE 29

/** The possible states of the HMM */
typedef enum {
    NOSTATE = -1,
    S_STATE,
    E_STATE,
    R_STATE,
    S_STATE_1,
    E_STATE_1,
    M1_STATE,
    M2_STATE,
    M3_STATE,
    M4_STATE,
    M5_STATE,
    M6_STATE,
    M1_STATE_1,
    M2_STATE_1,
    M3_STATE_1,
    M4_STATE_1,
    M5_STATE_1,
    M6_STATE_1,
    I1_STATE,
    I2_STATE,
    I3_STATE,
    I4_STATE,
    I5_STATE,
    I6_STATE,
    I1_STATE_1,
    I2_STATE_1,
    I3_STATE_1,
    I4_STATE_1,
    I5_STATE_1,
    I6_STATE_1,
} HMM_State;

/** The amount of possible state transitions */
#define NUM_TRANSITIONS 14

/**
 * The transition types (used in e.g. HMM->tr)
 * The first letter is the from state, the second letter the to state.
 *
 * The first 7 transitions occur in gene regions, which has states
 * M[atch], I[nsertion] and D[eletion]
 *
 * For example, TR_MD is a transition from a match state to a deletion state.
 */
typedef enum {
    TR_MM,
    TR_MI,
    TR_MD,
    TR_II,
    TR_IM,
    TR_DD,
    TR_DM,
    TR_GE,
    TR_GG,
    TR_ER,
    TR_RS,
    TR_RR,
    TR_ES,
    TR_ES1,
} HMM_StateTransition;

char *translation_table[65];
char *translation_table_rc[65];

// semaphores
#ifdef __APPLE__
typedef sem_t* SEM_T;
#elif __linux
typedef sem_t SEM_T;
#define sem_wait(x) sem_wait(&x)
#define sem_post(x) sem_post(&x)
#endif

/**
 * Should always be used when accessing/modifying a queue.
 */
SEM_T sema_Q;
/**
 * Used by the reader thread to make sure a worker thread doesn't try to handle input
 * before it's fully read.
 */
SEM_T sema_R;

SEM_T sema_r;
SEM_T sema_w;


typedef struct {

    /**
     * pi[1..N] pi[i] is the initial state distribution.
     */
    double  pi[NUM_STATE];

    /**
     * The transition probabilities (see also ::StateTransition)
     */
    double tr[NUM_TRANSITIONS];

    /**
     * The transition probability from a lowest-level state  to a lowest-level state
     */
    double e_M_1[6][16][4];
    double e_M[6][16][4];

    double tr_R_R[4][4];
    double tr_I_I[4][4];
    double tr_M_I[4][4];

    double tr_S[61][64];
    double tr_E[61][64];
    double tr_S_1[61][64];
    double tr_E_1[61][64];

    double S_dist[6];  /*sigma, mu,alpha, sigma_r, mu_r, alpha_r */
    double E_dist[6];
    double S1_dist[6];
    double E1_dist[6];
} HMM;

typedef struct {

    double trans[44][6][16][4];
    double rtrans[44][6][16][4];
    double noncoding[44][4][4];
    double start[44][61][64];
    double stop[44][61][64];
    double start1[44][61][64];
    double stop1[44][61][64];

    double S_dist[44][6];
    double E_dist[44][6];
    double S1_dist[44][6];
    double E1_dist[44][6];

} TRAIN;

/**
 * The data that can be used by each worker thread separately.
 */
typedef struct {
    /** The ID for this thread, a unique number. */
    unsigned int id;

    /** The hidden Markov model that will be used in this thread. */
    HMM *hmm;

    bool wholegenome;
    bool format;

    unsigned int *output_num_sequences;
    unsigned int *input_num_sequences;

    /** The FASTA headers. */
    char*** record_headers;
    /** The FASTA nucleotide sequences. */
    char*** record_sequences;
    /** The lengths of the FASTA nucleotide sequences. */
    int** record_sequences_lens;

    /** The outputs for printing to the meta file (with the putative gene coordinates. */
    char*** output_buffer;
    /** The outputs for the .faa-file (with all putative genes translated). */
    char*** aa_buffer;
    /** The outputs for the .ffn-file (with the putative genes). */
    char*** dna_buffer;

    /** The predicted insertions */
    int* insert;
    /** The predicted deletions */
    int* c_delete;

    /* Temporary buffers */
    char* dna;
    char* dna1;
    char* dna_f;
    char* dna_f1;
    char* protein;
    char* temp_str;

    SEM_T sema_r;
    SEM_T sema_w;
} ThreadData;

ThreadData *thread_datas;

HMM hmm;
TRAIN train;

/**
 * Creates the data for a worker thread
 */
void thread_data_init(ThreadData* td, unsigned int id);
int read_seq_into_buffer(FastaFile* fp, ThreadData *td, unsigned int buf, bool initial_input);

void get_prob_from_cg(HMM *hmm, TRAIN *train, char *O, int len_seq);
void get_train_from_file(char *filename, HMM *hmm_ptr, char *mfilename, char *mfilename1, char *nfilename,
                         char *sfilename,char *pfilename,char *s1filename,char *p1filename, char *dfilename, TRAIN *train_ptr);
void viterbi(HMM *hmm_ptr, char *O, char* output_buffer, char* aa_buffer, char *dna_buffer,
             char *sequence_head, bool whole_genome, bool format, int len_seq,
             char* dna_ptr, char* dna1_ptr, char* dna_f_ptr, char* dna_f1_ptr, char* protein_ptr,
             int* insert_ptr, int* c_delete_ptr, char* temp_str_ptr);

void free_hmm(HMM *hmm);

/**
 * Calculates the reverse-complement of given DNA
 *
 * @param dna The string we want to take the reverse-complement of.
 * @param[out] rc_dna The output buffer for the reverse-complement.
 */
void get_rc_dna(char *dna, int dna_len, char *rc_dna);

void get_rc_dna_indel(char* dna_f, int dna_len, char* dna_f1);
void get_corrected_dna(char *dna, char *dna_f);
void print_usage();

/**
 * @brief Prints the output of a gene to the different buffers.
 *
 * Writes a predicted gene and corresponding information to the different buffers.
 *
 * Putative gene information
 * -------------------------
 * @param codon_start What strand the gene lies on (1, -1 or 0 if unknown).
 * @param start_t The start position of the gene.
 * @param end_t The ending position of the gene.
 * @param frame The frame number.
 * @param dna The complete DNA input string.
 * @param dna_seq The complete DNA nucleotide sequence.
 * @param dna_len The length of `dna`and `dna_seq`.
 * @param dna_f The (supposedly?) formatted nucleotide sequence of the gene (on the sense strand).
 * @param insertions The list of the insertions
 * @param insertions_len The length of `insertions`
 * @param deletions The list of the deletions
 * @param deletions_len The length of `deletions`
 *
 * The output buffers
 * ------------------
 * @param output_buffer The buffer for the .out-file, consisting of putative genes coordinates.
 * @param aa_buffer The buffer for the .faa-file, the FASTA file that lists the translated proteins
 *                  of the putative genes.
 * @param dna_buffer The buffer for the .fnn-file, the FASTA file that lists the putative genes.
 *
 * Other output information
 * -------------------------
 * @param sequence_head_short The sequence header for the FASTA files (without '>').
 * @param format Whether to use a formatted sequence (see `dna_f` and `rc_dna_f`)
 * @param multiple (> 0) if multiple genes were found
 *
 * Temporary output buffers (already allocated)
 * ------------------------
 * @param rc_dna buffer which you can fill with the reverse complement of the dna.
 * @param rc_dna_f buffer which you can fill with the formatted(?) reverse complement of the dna.
 * @param protein buffer which you can fill with the translated gene.
 * @param temp_str_ptr General-purpose string buffer
 */
void print_outputs(int codon_start, int start_t, int end_t, int frame,
                   char *output_buffer, char *aa_buffer, char *dna_buffer,
                   char *sequence_head_short,
                   char *dna, int dna_len, Nucleotide dna_seq[], char *rc_dna, char* dna_f, char *rc_dna_f, char *protein,
                   int *insertions, int *deletions, int insertions_len, int deletions_len,
                   bool format, char *temp_str_ptr, unsigned int multiple);

// helper functions to cleanup the main function
void setTrainDirectory(char* train_path);
void conductWork();

/**
 * Initializes the writer thread and the worker threads.
 */
void initializeThreads();

void destroySemaphores();
void initializeSemaphores();
void setupProgram(int argc, char** argv);
#endif
