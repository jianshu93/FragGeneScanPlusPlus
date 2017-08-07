/**
 * @file util_lib.h
 */

#ifndef __UTIL_LIB_H
#define __UTIL_LIB_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "hmm.h"

double **dmatrix(int num_col);
int **imatrix(int num_col);
int *ivector(int nh);

void free_dmatrix(double **m);
void free_imatrix(int **m);

/**
 * Parses the given string into a HMM_StateTransition.
 */
HMM_StateTransition hmm_state_transition_parse (const char *nt);

/**
 * Parses the given character as a nucleotide.
 */
Nucleotide nucleotide_parse (char nt);

/**
 * Returns the reverse complement of the given nucleotide.
 */
Nucleotide nucleotide_complement (Nucleotide nt);

/**
 * Returns the trinucleotide from combining the given nucleotides.
 */
int trinucleotide (Nucleotide a, Nucleotide b, Nucleotide c);

/**
 * Translates the given DNA sequence.
 *
 * @param dna The sequence we want to translate.
 * @param dna The length of the sequence.
 * @param[out] protein The output buffer for the protein.
 * @param strand What strand to translate.
 */
void get_protein(const Nucleotide dna[], int dna_len, char *protein, Strand strand);

/**
 * Print out the usage info for the executable.
 */
void print_usage();

/**
 * A queue for worker thread buffers.
 */
typedef struct BufferQueue {
    /** The next thread in the queue */
    struct BufferQueue *next;
    /** The worker thread */
    ThreadData *td;
    /** The specific buffer */
    unsigned int buffer;
} QUEUE;

/** The queue of empty buffers, waiting for input by the reader thread. */
QUEUE *q_empty_head;
QUEUE *q_empty_tail;
/** The queue of finished buffers, waiting for the writer thread to handle output. */
QUEUE *q_done_head;
QUEUE *q_done_tail;

/**
 * Prints out the contents of the given buffer.
 */
void printq(unsigned int which);

/**
 * Enqueues the given worker thread's buffer to the given queue
 */
void enq(ThreadData *td, unsigned int buffer, unsigned int which);

/**
 * Dequeues the first element of the given queue.
 */
QUEUE *deq(unsigned int which);

/**
 * Loads the requested queue into the first arguments' address.
 */
void cutnpaste_q(QUEUE **dest, unsigned int which);

/**
 * A lazy version of memset().
 */
void stopMemset(char *ptr, int length);

#endif
