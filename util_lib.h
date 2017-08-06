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

HMM_StateTransition hmm_state_transition_parse (char *nt);

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
void get_protein(Nucleotide dna[], int dna_len, char *protein, Strand strand);
void print_usage();

typedef struct q {
    struct q *next;
    ThreadData *td;
    unsigned int buffer;
} QUEUE;

QUEUE *q_empty_head;
QUEUE *q_empty_tail;
QUEUE *q_done_head;
QUEUE *q_done_tail;


void printq(unsigned int which);
void enq(ThreadData *td, unsigned int buffer, unsigned int which);
QUEUE *deq(unsigned int which);

void cutnpaste_q(QUEUE **dest, unsigned int which);

void stopMemset(char *ptr, int length);

#endif
