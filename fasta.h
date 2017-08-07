/**
 * @file fasta.h
 *
 * @brief Declarations for simple FASTA I/O operations
 *
 * Purpose:  A very rudimentary FASTA file reading API. Designed
 *           for simplicity and clarity, not for robustness.
 *
 *           A code example:
 * \code
 *     FastaFile *ffp;
 *     char *seq, *name;
 *     int seqlen;
 *
 *     ffp = fasta_file_new(seqfile);
 *     while (fasta_file_read_record(ffp, &seq, &name, &seqlen) {
 *         // do stuff with sequence;
 *         free(name);
 *         free(seq);
 *     }
 *     fasta_file_free(ffp);
 * \endcode
 *
 * Commentary:
 *           The basic problem with reading FASTA files is that there is
 *           no end-of-record indicator. When you're reading sequence n,
 *           you don't know you're done until you've read the header line
 *           for sequence n+1, which you won't parse 'til later (when
 *           you're reading in the sequence n+1). One common trick for
 *           this is to implement a one-line "lookahead" buffer that you
 *           can peek at, before parsing later.
 *
 *           This buffer is kept in a small structure (a FastaFile), rather
 *           than in a static char[] in the function. This allows
 *           us to have multiple FASTA files open at once. The static approach
 *           would only allow us to have one file open at a time. ANSI C
 *           predates the widespread use of parallel programming. It was
 *           not overly concerned about the drawbacks of statics. Today,
 *           though, you should keep in mind that you may someday want to
 *           turn your program into a multithreaded, parallel program, and
 *           all functions in parallelized code must be "reentrant": able to
 *           be called a second time - with different arguments,
 *           and while the code in the first function call is still executing! -
 *           without overwriting or corrupting any static storage in the
 *           function. Statics have fewer uses now (for example, to
 *           test that some initialization code for a function is run once
 *           and only once.)
 *
 * Limitations:
 *           There is no error handling, for clarity's sake. Also,
 *           the parser is brittle. Improper FASTA files (for instance,
 *           blank lines between records) will cause unexpected
 *           behavior. Real file parsers are more complex.
 *           In real life, they have to deal with absolutely anything the user might
 *           pass as a "FASTA file"; and either parse it correctly,
 *           or detect that it's an invalid format and fail cleanly.
 *
 *           Lines are read in from the file using ANSI C's fgets(). fgets()
 *           requires a maximum buffer length (here, STRINGLEN, which is
 *           defined as 512 in bio5495.h). Some FASTA files have very long
 *           description lines, however; notably the NCBI NR database. Static
 *           limitations on things like line or sequence lengths should be
 *           avoided. An example of a replacement for fgets() that dynamically
 *           allocates its buffer size and allows any line length is
 *           SQUID's sre_fgets().
 *
 *           We use ANSI C's strtok() to parse the sequence name out of the line.
 *           strtok() is deprecated in modern programs because it is not threadsafe.
 *           (See comments above.) An example of a threadsafe version is
 *           SQUID's sre_strtok().
 *
 * Returns:
 */

#ifndef FASTA_H
#define FASTA_H

#include <stdbool.h>
#include <stdio.h>

/**
 * The default length for a string buffer.
 * By default, set to 1,000,000 chars => 1 MB
 */
#define STRINGLEN 1000000

/** A FASTA-formatted file */
typedef struct {
    /** The file pointer to the FASTA file */
    FILE *fp;
    /** The buffer that is used to read data into */
    char buffer[STRINGLEN];
} FastaFile;

/**
 * Creates a new FastaFile.
 *
 * @param seqfile The path of the FASTA file.
 * @returns a FastaFile pointer, or NULL on failure (e.g. if the file doesn't exist, or isn't readable).
 */
FastaFile *fasta_file_new(const char *seqfile);

/**
 * Reads the next record of the FASTA file into a header and sequence buffer.
 *
 * @param fp The FASTA file to read from.
 * @param[out] out_seq An unallocated buffer to contain the sequence. Freed by the caller.
 * @param[out] out_header An unallocated buffer to contain the header. Freed by the caller.
 * @param[out] out_seq_len Will be filled with the length of the sequence.
 *
 * @returns true if more sequences can be read in the file, false otherwise.
 */
bool fasta_file_read_record(FastaFile *fp, char **out_seq, char **out_header, int *out_seq_len);

/**
 * Frees the resource for the FastaFile.
 */
void fasta_file_free(FastaFile *ffp);

#endif
