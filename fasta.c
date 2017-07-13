/* Simple API for FASTA file reading
 * for Bio5495/BME537 Computational Molecular Biology
 * SRE, Sun Sep  8 05:35:11 2002 [AA2721, transatlantic]
 * CVS $Id: fasta.c,v 1.1 2003/10/05 18:43:39 eddy Exp $
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "fasta.h"

FastaFile *fasta_file_new(char *seqfile) {

    FastaFile *ffp;
    ffp = malloc(sizeof(FastaFile));

    if (strcmp(seqfile, "stdin") == 0) {
        ffp->fp = stdin;
    } else {
        ffp->fp = fopen(seqfile, "r");              /* Assume seqfile exists & readable!   */
        if (ffp->fp == NULL) {
            free(ffp);
            return NULL;
        }
        if ((fgets(ffp->buffer, STRINGLEN, ffp->fp)) == NULL) {
            free(ffp);
            return NULL;
        }
    }

    return ffp;
}

int fasta_file_read_record(FastaFile *ffp, char **out_seq, char **out_header, int *out_seq_len) {
    char *s, *header, *seq;
    int n, nalloc;

    /* Peek at the lookahead buffer; check if it's a valid FASTA header. */
    if (ffp->buffer[0] != '>')
        return 0;

    /* Parse out the header */
    s  = strtok(ffp->buffer+1, "\n");
    //header = malloc(sizeof(char) * (strlen(s)+1));
    header = malloc(sizeof(char) * 1024);
    strcpy(header, s);

    /* Everything else 'til the next descline is the sequence.
     * Note the idiom for dynamic reallocation of seq as we
     * read more characters, so we don't have to assume a maximum
     * sequence length.
     */
    seq = calloc(1024, sizeof(char));
    n = 0;
    nalloc = 128;
    while (fgets(ffp->buffer, STRINGLEN, ffp->fp)) {
        if (ffp->buffer[0] == '>')
            break;	/* We've reached the next header */

        for (s = ffp->buffer; *s != '\0'; s++) {
            if (!isalpha(*s))
                continue;  /* accept any alphabetic character */
            seq[n] = *s;                  /* store the character, bump length n */
            n++;
            if (nalloc == n) {        /* are we out of room in seq? if so, expand */
                /* (remember, need space for the final '\0')*/
                nalloc += 128;
                seq = realloc(seq, sizeof(char) * nalloc);
            }
        }
    }
    seq[n] = '\0';

    *out_header = header;
    *out_seq = seq;
    *out_seq_len = n;
    return 1;
}

void fasta_file_free(FastaFile *ffp) {
    fclose(ffp->fp);
    free(ffp);
}
