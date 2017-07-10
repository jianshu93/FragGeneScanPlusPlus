/* fasta.h
 * Declarations for simple FASTA i/o library
 * SRE, Sun Sep  8 05:37:38 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.h,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */

#ifndef FASTA_H
#define FASTA_H

#include <stdio.h>

/* 1,000,000 chars => 1 MB */
#define STRINGLEN 1000000

typedef struct fastafile_s {
    FILE *fp;
    char  buffer[STRINGLEN];
} FastaFile;

FastaFile *fasta_file_new(char *seqfile);

/**
 * @brief Reads the next record from the FASTA file.
 *
 * Reads a record FASTA file into a header buffer and sequence buffer.
 *
 * @note `out_seq` and `out_header` need to be freed by the caller.
 *
 * @param out_seq An unallocated buffer to contain the sequence. Freed by the caller.
 * @param out_header An unallocated buffer to contain the header. Freed by the caller.
 * @param out_seq_len Will be filled with the length of the sequence.
 *
 * @returns 1 if more sequences can be read in the file, 0 otherwise.
 */
int fasta_file_read_record(FastaFile *fp, char **out_seq, char **out_name, int *out_L);

void fasta_file_free(FastaFile *ffp);

#endif
