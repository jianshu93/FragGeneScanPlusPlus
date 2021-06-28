/**
 * @file run_hmm.h
 */

#ifndef __RUN_HMM_H
#define __RUN_HMM_H

#include <stdio.h>
#include "hmm.h"

void writeDNA();
void writeMeta();
void writeAminoAcids(FILE *aa_outfile_fp, ThreadData *td, unsigned int buffer);

void parseArguments(int argc, char **argv);
void checkFiles();
void setMemoryLimits();
void checkOutputFiles();


void *writerThread(void *args);
void readerThread();
void *workerThread(void *_thread_datas);

#endif
