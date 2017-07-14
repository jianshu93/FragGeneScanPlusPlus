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
#include <math.h>
#include <pthread.h>
#include <semaphore.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "fasta.h" // STRINGLEN 

#define DONE_Q 0
#define EMPTY_Q 1
#define max_dbl 10000000000.0

#define LOG_53 -0.63487827243
#define LOG_16 -1.83258146375
#define LOG_30 -1.20397280433
#define LOG_25 -1.38629436112
#define LOG_83 -0.18632957819
#define LOG_10 -2.30258509299
#define LOG_07 -2.65926003693
#define LOG_95 -0.05129329438

#define NUM_STATE 29

#define NOSTATE -1
#define S_STATE 0
#define E_STATE 1
#define R_STATE 2
#define S_STATE_1 3
#define E_STATE_1 4
#define M1_STATE 5
#define M2_STATE 6
#define M3_STATE 7
#define M4_STATE 8
#define M5_STATE 9
#define M6_STATE 10
#define M1_STATE_1 11
#define M2_STATE_1 12
#define M3_STATE_1 13
#define M4_STATE_1 14
#define M5_STATE_1 15
#define M6_STATE_1 16
#define I1_STATE 17
#define I2_STATE 18
#define I3_STATE 19
#define I4_STATE 20
#define I5_STATE 21
#define I6_STATE 22
#define I1_STATE_1 23
#define I2_STATE_1 24
#define I3_STATE_1 25
#define I4_STATE_1 26
#define I5_STATE_1 27
#define I6_STATE_1 28

#define TR_MM 0
#define TR_MI 1
#define TR_MD 2
#define TR_II 3
#define TR_IM 4
#define TR_DD 5
#define TR_DM 6
#define TR_GE 7
#define TR_GG 8
#define TR_ER 9
#define TR_RS 10
#define TR_RR 11
#define TR_ES 12
#define TR_ES1 13

char hmm_file[STRINGLEN];
char aa_file[STRINGLEN];
char seq_file[STRINGLEN];
char out_file[STRINGLEN];
char dna_file[STRINGLEN];
char train_file[STRINGLEN];
char mstate_file[STRINGLEN];
char rstate_file[STRINGLEN];
char nstate_file[STRINGLEN];
char sstate_file[STRINGLEN];
char pstate_file[STRINGLEN];
char s1state_file[STRINGLEN];     /* stop codon of gene in - stand */
char p1state_file[STRINGLEN];
char dstate_file[STRINGLEN];
char train_dir[STRINGLEN];

// semaphores
#ifdef __APPLE__
typedef sem_t* SEM_T;
#elif __linux
typedef sem_t SEM_T;
#define sem_wait(x) sem_wait(&x)
#define sem_post(x) sem_post(&x)
#endif

SEM_T sema_Q;
SEM_T sema_R;
SEM_T sema_r;
SEM_T sema_w;


typedef struct {

    double  pi[29];    /* pi[1..N] pi[i] is the initial state distribution. */

    double tr[14];                 /* transition probability from a (delete/insert/match) state to a state */

    double e_M_1[6][16][4];      /* transition probability from a lowest-level state  to a  lowest-level state*/
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
 * The data that can be used by each thread separately (without locking).
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

void thread_data_init(ThreadData* td);
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
 * Translates the given DNA sequence.
 *
 * @param dna The sequence we want to translate.
 * @param[out] protein The output buffer for the protein.
 * @param strand What strand to translate.
 */
void get_protein(char *dna, char *protein, int strand);

/**
 * Calculates the reverse-complement of given DNA
 *
 * @param dna The string we want to take the reverse-complement of.
 * @param[out] rc_dna The output buffer for the reverse-complement.
 */
void get_rc_dna(char *dna, char *rc_dna);

void get_rc_dna_indel(char* dna_f, char* dna_f1);
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
 * @param dna The nucleotide sequence of the gene (on the sense strand).
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
                   char* output_buffer, char* aa_buffer, char* dna_buffer,
                   char* sequence_head_short,
                   char* dna, char* rc_dna, char* dna_f, char* rc_dna_f, char* protein,
                   int* insertions, int* deletions, int insertions_len, int deletions_len,
                   bool format, char* temp_str_ptr, unsigned int multiple);

// helper functions to cleanup the main function
void setTrainDirectory(char* train_path);
void conductWork();
void initializeThreads();
void destroySemaphores();
void initializeSemaphores();
void setupProgram(int argc, char** argv);
#endif
