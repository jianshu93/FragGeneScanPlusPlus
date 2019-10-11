#include "hmm.h"
#include "run_hmm.h"
#include "translation_tables.h"
#include "util_lib.h"

/* User-specified options */
bool wholegenome = false;
bool output_dna = false;
bool output_meta = false;
bool verbose = false;
int chunk_size = 0;
unsigned int threadnum = 1;
char aa_file[STRINGLEN];
char seq_file[STRINGLEN];
char out_file[STRINGLEN];
char dna_file[STRINGLEN];
char train_dir[STRINGLEN];

/* Files located in the train directory */
char hmm_file[STRINGLEN];
char train_file[STRINGLEN];
char mstate_file[STRINGLEN];
char rstate_file[STRINGLEN];
char nstate_file[STRINGLEN];
char sstate_file[STRINGLEN];
char pstate_file[STRINGLEN];
char s1state_file[STRINGLEN];     /* stop codon of gene in - stand */
char p1state_file[STRINGLEN];
char dstate_file[STRINGLEN];

unsigned int MAX_SEQS_PER_BUFFER;

/* READER THREAD */
/** The input sequence file */
FastaFile *fp = NULL;
/** Whether we're done reading the whole input */
bool num_reads_flag = false;
/** The nr of sequences we've read */
long read_counter = 0;
off_t stopped_at_fpos; // tracks how far we've read in the input

/* WRITER THREAD */
/** The writer thread, the thread that waits for finished worker threads and handles the output */
pthread_t writer_thread;
FILE *outfile_fp;
FILE *dna_outfile_fp;
/** The nr of sequences the writer thread has outputted */
long writer_counter = 0;

void parseArguments(int argc, char **argv) {
    /* read command line argument */
    //!! This argument reading should all be encapsulated in a single function, this will make reading the code much easier, right now we have to always move around it.
    if (argc <= 8) {
        fprintf(stderr, "ERROR: You missed some parameters for input\n");
        print_usage();
        exit(EXIT_FAILURE);
    }

    int c;
    int translation_table_id = 11;

    while ((c = getopt(argc, argv, "dec:o:p:r:s:t:vw:x:")) != -1) {
        switch (c) {
        case 'd':
            output_dna = true;
            break;
        case 'e':
            output_meta = true;
            break;
        case 'c':
            chunk_size = atoi(optarg);
            break;
        case 'o':
            strcpy(out_file, optarg);
            strcpy(aa_file, out_file);
            strcat(aa_file, ".faa");
            strcpy(dna_file, out_file);
            strcat(dna_file, ".ffn");
            break;
        case 'p':
            threadnum = atoi(optarg);
            if (threadnum < 1) {
                fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'r':
            setTrainDirectory(optarg);
            break;
        case 's':
            strcpy(seq_file, optarg);

            if (strcmp(seq_file, "stdin") == 0) {
                break;
            } else if (access(seq_file, F_OK)==-1) {
                fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 't':
            strcpy(train_file, optarg);
            strcpy(hmm_file, train_dir);
            strcat(hmm_file, train_file);
            if (access(hmm_file, F_OK)==-1) {
                fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'v':
            verbose = true;
            break;
        case 'w':
            wholegenome = atoi(optarg);
            if (wholegenome != 0 && wholegenome != 1) {
                fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
                print_usage();
                exit(EXIT_FAILURE);
            }
            break;
        case 'x':
            translation_table_id = atoi(optarg);
            break;
        default:
            break;
        }

        if (!parse_translation_tables(translation_table_id, &translation_table, &translation_table_rc)) {
            fprintf(stderr, "ERROR: No translation table with number %s\n", optarg);
            exit(EXIT_FAILURE);
        }
    }
}

void checkFiles() {
    /* check whether the specified files exist */
    if (access(mstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: Forward prob. file [%s] does not exist\n", mstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(rstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: Backward prob. file [%s] does not exist\n", rstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(nstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: noncoding prob. file [%s] does not exist\n", nstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(sstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: start prob. file [%s] does not exist\n", sstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(pstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: stop prob. file [%s] does not exist\n", pstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(s1state_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: start1 prob. file [%s] does not exist\n", s1state_file);
        exit(EXIT_FAILURE);
    }
    if (access(p1state_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: stop1 prob. file [%s] does not exist\n", p1state_file);
        exit(EXIT_FAILURE);
    }
    if (access(dstate_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: pwm dist. file [%s] does not exist\n", dstate_file);
        exit(EXIT_FAILURE);
    }
    if (access(hmm_file, F_OK)==-1) {
        fprintf(stderr, "ERROR: hmm file [%s] does not exist\n", hmm_file);
        exit(EXIT_FAILURE);
    }
}

void setMemoryLimits() {
    /* check for mem limit, allocate buffer */
    if (chunk_size <= 0) {
        printf("Minimum chunk size specified invalid, defaulting to 1\n");
        chunk_size = 1;
    }

    // 5 stands for the number of buffers we are currently using per thread
    MAX_SEQS_PER_BUFFER = chunk_size;
}

void checkOutputFiles() {

    // remove them, if they already exist
    remove(aa_file);
    if (output_meta) remove(out_file);
    if (output_dna) remove(dna_file);
}

void setTrainDirectory(char *train_path) {
    strcpy(train_dir, train_path);
    strcat(train_dir, "/");
    strcpy(mstate_file, train_dir);
    strcat(mstate_file, "gene");
    strcpy(rstate_file, train_dir);
    strcat(rstate_file, "rgene");
    strcpy(nstate_file, train_dir);
    strcat(nstate_file, "noncoding");
    strcpy(sstate_file, train_dir);
    strcat(sstate_file, "start");
    strcpy(pstate_file, train_dir);
    strcat(pstate_file, "stop");
    strcpy(s1state_file, train_dir);
    strcat(s1state_file, "stop1");
    strcpy(p1state_file, train_dir);
    strcat(p1state_file, "start1");
    strcpy(dstate_file, train_dir);
    strcat(dstate_file, "pwm");
}

void initializeSemaphores() {

#ifdef __APPLE__
    sem_unlink("/sema_Q");
    sem_unlink("/sema_R");
    sem_unlink("/sema_F");
    sem_unlink("/sema_r");
    sem_unlink("/sema_w");

    if ((sema_Q = sem_open("/sema_Q", O_CREAT, 0644, 1)) == SEM_FAILED ||
            (sema_R = sem_open("/sema_R", O_CREAT, 0644, 1)) == SEM_FAILED ||
            (sema_F = sem_open("/sema_F", O_CREAT, 0644, 1)) == SEM_FAILED ||
            (sema_r = sem_open("/sema_r", O_CREAT, 0644, 1)) == SEM_FAILED ||
            (sema_w = sem_open("/sema_w", O_CREAT, 0644, 1)) == SEM_FAILED) {
        perror("ERROR: sem_open");
        exit(EXIT_FAILURE);
    }

#elif __linux
    sem_init(&sema_Q, 0, 1);
    sem_init(&sema_R, 0, 1);
    sem_init(&sema_F, 0, 1);
    sem_init(&sema_r, 0, 0);
    sem_init(&sema_w, 0, 0);
#endif
}

void destroySemaphores() {

#ifdef __APPLE__
    sem_unlink("/sema_Q");
    sem_unlink("/sema_R");
    sem_unlink("/sema_F");
    sem_unlink("/sema_r");
    sem_unlink("/sema_w");

    char name[40];
    int j;
    for (j = 0; j<threadnum; j++) {
        sprintf(name, "/sema_r%d", j);
        sem_unlink(name);

        sprintf(name, "/sema_w%d", j);
        sem_unlink(name);
    }

#elif __linux
    sem_destroy(&sema_Q);
    sem_destroy(&sema_R);
    sem_destroy(&sema_F);
    sem_destroy(&sema_r);
    sem_destroy(&sema_w);
#endif

}

void initializeThreads() {
    unsigned int i, j;

    // allocate memory for each thread only once!
    log_debug("Allocating memory for all threads...\n");

    pthread_t *thread = calloc(threadnum, sizeof(pthread_t));
    thread_datas = calloc(threadnum, sizeof(ThreadData));
    for (i = 0; i < threadnum; i++)
        thread_data_init(thread_datas + i, i);

    log_debug("Allocated memory for all threads!\n");


    log_debug("Starting the writer thread...\n");
    pthread_create(&writer_thread, 0, writerThread, 0);


    log_debug("Opening the sequence file...\n");
    fp = fasta_file_new(seq_file);
    if (!fp) {
        printf("ERROR! Could not open seqence file %s for reading...!\n", seq_file);
        exit(EXIT_FAILURE);
    }


    log_debug("Giving workers initial inputs...\n");
    for (j = 0; j < threadnum; j++) {
        for (i = 0; i < 2; i++) {
            if ((stopped_at_fpos = read_seq_into_buffer(fp, thread_datas + j, i, true)) != 0) {
                sem_post(thread_datas[j].sema_r);
            }
        }
    }

    log_debug("Starting worker threads...\n");
    for (j = 0; j < threadnum; j++)
        pthread_create(&thread[j], 0, workerThread, (void *)(thread_datas+j));
}

void readerThread() {
    // master loop - while we haven't exhausted reading the file yet
    while (stopped_at_fpos!=0) {
        sem_wait(sema_r);

        /* Fetch the queue of threads that are waiting for input */
        sem_wait(sema_Q);
        QUEUE *temp;
        cutnpaste_q(&temp, EMPTY_Q);
        sem_post(sema_Q);

        /* Iterate over the queue */
        for (; temp != NULL; temp = temp->next) {
            sem_wait(sema_R);
            stopped_at_fpos = read_seq_into_buffer(fp,  temp->td, temp->buffer, false);

            /* We've fully read the input sequence file */
            if (stopped_at_fpos == 0) {
                sem_wait(sema_F);
                num_reads_flag = true;
                sem_post(sema_F);
            }

            sem_post(sema_R);

            /* Tell the worker thread we're done reading input */
            sem_post(temp->td->sema_r);
        }
    }

    log_debug("Finished handing out all the work...\n");
    fasta_file_free(fp);

    sem_post(sema_w); /* ensure it doesn't block */
    num_reads_flag = true;
    sem_post(sema_F);
}

int main (int argc, char **argv) {
    setTrainDirectory("train");

    parseArguments(argc, argv);

    checkFiles();

    setMemoryLimits();

    checkOutputFiles();

    initializeSemaphores();

    log_debug("Max number of sequences per thread : %d\n", MAX_SEQS_PER_BUFFER);

    /* read all initial model */
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

    /* prepare all of the worker threads as well as the writer thread */
    initializeThreads();

    /* The master thread becomes the reader threadn which reads the rest of the input sequence file */
    readerThread();

    /* Now wait for the writer thread to finish before exiting */
    pthread_join(writer_thread, NULL);

    /* destroy the semaphores */
    destroySemaphores();

    printf("Run finished with %d threads.\n", threadnum);

    return EXIT_SUCCESS;
}

int read_seq_into_buffer(FastaFile *ffp, ThreadData *thread_data, unsigned int buf, bool initial_input) {
    char *seq, *name;
    int seq_len;
    unsigned int count = 0, i;

    if (! initial_input) {

        for (i = 0; i < MAX_SEQS_PER_BUFFER; i++) {
            free(thread_data->record_headers[buf][i]);
            free(thread_data->record_sequences[buf][i]);
        }
    }

    while ((count < MAX_SEQS_PER_BUFFER) && fasta_file_read_record(ffp, &seq, &name, &seq_len)) {
        thread_data->record_headers[buf][count] = name;
        thread_data->record_sequences[buf][count] = seq;
        thread_data->record_sequences_lens[buf][count] = seq_len;
        read_counter++;
        count++;
    }

    thread_data->input_num_sequences[buf] = count;

    return count;
}

void thread_data_init(ThreadData *td, unsigned int id) {
    unsigned int i, j;
    // Initialize thread data structure

    td->hmm = calloc(1, sizeof(HMM));
    memcpy(td->hmm, &hmm, sizeof(HMM));

    td->wholegenome = wholegenome;

    td->id = id;

#ifdef __APPLE__
    char name[40];
    sprintf(name, "/sema_r%d", td->id);
    sem_unlink(name);

    if ((td->sema_r = sem_open(name, O_CREAT, 0644, 0)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

    sprintf(name, "/sema_w%d", td->id);
    sem_unlink(name);
    if (( td->sema_w = sem_open(name, O_CREAT, 0644, 2)) == SEM_FAILED ) {
        perror("sem_open");
        exit(EXIT_FAILURE);
    }

#elif __linux
    sem_init(&td->sema_r, 0, 0);
    sem_init(&td->sema_w, 0, 2);
#endif

    // TODO : refactor to as many single large malloc calls as possible
    td->output_num_sequences = calloc(2, sizeof(unsigned int));
    td->input_num_sequences = calloc(2, sizeof(unsigned int));

    td->record_headers = malloc(sizeof(char **) * 2);
    td->record_sequences = malloc(sizeof(char **) * 2);
    td->record_sequences_lens = malloc(sizeof(int *) * 2);

    td->output_buffer = malloc(sizeof(char **) * 2);
    td->aa_buffer = malloc(sizeof(char **) * 2);
    td->dna_buffer = malloc(sizeof(char **) * 2);

    td->dna = calloc(STRINGLEN, sizeof(char));
    td->dna1 = calloc(STRINGLEN, sizeof(char));
    td->protein = calloc(STRINGLEN, sizeof(char));
    td->temp_str = calloc(STRINGLEN, sizeof(char));

    td->insert = malloc(sizeof(int) * STRINGLEN);
    td->c_delete = malloc(sizeof(int) * STRINGLEN);

    for (i = 0; i < 2; i++) {
        td->record_headers[i] = malloc(sizeof(char *) * MAX_SEQS_PER_BUFFER);
        td->record_sequences[i] = malloc(sizeof(char *) * MAX_SEQS_PER_BUFFER);
        td->record_sequences_lens[i] = malloc(sizeof(int) * MAX_SEQS_PER_BUFFER);
        td->output_buffer[i] = malloc(sizeof(char *) * MAX_SEQS_PER_BUFFER);
        td->aa_buffer[i] = malloc(sizeof(char *) * MAX_SEQS_PER_BUFFER);
        td->dna_buffer[i] = malloc(sizeof(char *) * MAX_SEQS_PER_BUFFER);

        for (j = 0; j < MAX_SEQS_PER_BUFFER; j++) {
            td->aa_buffer[i][j] = calloc(STRINGLEN, sizeof(char));
            td->dna_buffer[i][j] = calloc(STRINGLEN, sizeof(char));
            td->output_buffer[i][j] = calloc(STRINGLEN, sizeof(char));
        }
    }
}

void writeOutputFiles(FILE *aa_outfile_fp, ThreadData *td, unsigned int buffer) {
    if (output_meta)
        writeMeta();
    if (output_dna)
        writeDNA();

    writeAminoAcids(aa_outfile_fp, td, buffer);
    log_debug("Wrote results for thread %d, buffer %d.\n", td->id, buffer );
}

void writeDNA() {
    dna_outfile_fp = fopen(dna_file, "a");
    if (!dna_outfile_fp) {
        printf("ERROR: Could not open dna output file %s for writing!\n", dna_file);
        exit(EXIT_FAILURE);
    }
}

void writeMeta() {
    outfile_fp = fopen(out_file, "a");
    if (!outfile_fp) {
        printf("ERROR: Could not open meta output file %s for writing!\n", out_file);
        exit(EXIT_FAILURE);
    }
}

void writeAminoAcids(FILE *aa_outfile_fp, ThreadData *td, unsigned int buffer) {
    unsigned int j;

    for (j = 0; j < td->output_num_sequences[buffer]; j++) {
        writer_counter++;
        char *ptrc;
        if (td->aa_buffer[buffer][j][0]!=0) {
            ptrc=td->aa_buffer[buffer][j];
            while (*ptrc!='\0') {
                if (*ptrc=='\t') *ptrc='>';
                ptrc++;
            }
            fprintf(aa_outfile_fp, ">%s", td->aa_buffer[buffer][j]);
        }

        //!! Why are we clearing dna and output buff?
        stopMemset(td->output_buffer[buffer][j], STRINGLEN);
        stopMemset(td->aa_buffer[buffer][j], STRINGLEN);
        stopMemset(td->dna_buffer[buffer][j], STRINGLEN);
    }
}

FILE *openFilePointers() {
    FILE *aa_outfile_fp;

    aa_outfile_fp = (strcmp(out_file, "stdout") == 0)? stdout : fopen(aa_file, "a");
    if (!aa_outfile_fp) {
        printf("ERROR: Could not open aa output file %s for writing!\n", aa_file);
        exit(EXIT_FAILURE);
    }

    return aa_outfile_fp;
}

void closeFilePointers( FILE **aa_outfile_fp, FILE **outfile_fp, FILE **dna_outfile_fp ) {

    fclose(*aa_outfile_fp);
    if (output_meta)
        fclose(*outfile_fp);
    if (output_dna)
        fclose(*dna_outfile_fp);

    *aa_outfile_fp = NULL;
    *outfile_fp = NULL;
    *dna_outfile_fp = NULL;
}

void *writerThread(void *args) {
    FILE *aa_outfile_fp = openFilePointers();

    sem_wait(sema_F);
    while (true) {
        QUEUE *temp;

        sem_wait(sema_w);
        sem_post(sema_F);

        /* Fetch the queue of worker threads that are done */
        sem_wait(sema_Q);
        cutnpaste_q(&temp, DONE_Q);
        sem_post(sema_Q);

        /* Iterate over the queue */
        for (; temp != NULL; temp = temp->next) {
            sem_wait(sema_R);

            ThreadData *td = temp->td;
            writeOutputFiles(aa_outfile_fp, td, temp->buffer);

            sem_post(sema_R);

            /* Tell the worker thread we've written the output,
             * so it can continue using its buffers */
            sem_post(td->sema_w);
        }

        /* Check if we're done, i.e. the whole input was read and processed*/
        sem_wait(sema_F);
        if (num_reads_flag && writer_counter == read_counter)
            break;
    }
    sem_post(sema_F);

    /* We're done writing, so close the output files */
    closeFilePointers(&aa_outfile_fp, &outfile_fp, &dna_outfile_fp );

    return NULL;
}

void runViterbiOnBuffers(ThreadData *td, unsigned int b) {
    unsigned int i;

    for (i = 0; i < td->input_num_sequences[b]; i++) {
        get_prob_from_cg(td->hmm, &train, td->record_sequences[b][i], td->record_sequences_lens[b][i]);

        if (td->record_sequences[b][i] && td->record_headers[b][i] ) {

            viterbi(td->hmm, td->record_sequences[b][i], td->output_buffer[b][i],
                    td->aa_buffer[b][i], td->dna_buffer[b][i],
                    td->record_headers[b][i], td->wholegenome, td->record_sequences_lens[b][i],
                    td->dna, td->dna1, td->protein, td->insert, td->c_delete, td->temp_str);
        }
    }

}

void *workerThread(void *_thread_datas) {

    ThreadData *td = (ThreadData *)_thread_datas;
    unsigned int b = 0;

    while (true) {
        /* Wait until the reader thread allows us to start */
        sem_wait(td->sema_r);

        /* Wait until the writer thread is done writing output
         * from the previous buffers before continuing */
        sem_wait(td->sema_w);

        runViterbiOnBuffers(td, b);
        td->output_num_sequences[b] = td->input_num_sequences[b];

        log_debug("Thread %d buffer %d done work on %d sequences!\n",
                  td->id, b, td->input_num_sequences[b]);

        /* Critical section: queues */
        sem_wait(sema_Q);
        enqueue(td, b, EMPTY_Q);
        enqueue(td, b, DONE_Q);
        sem_post(sema_Q);

        /* Block again until more was read */
        sem_post(sema_r);
        sem_post(sema_w);

        b = (b + 1) % 2;
    }

    return (void *) 0;
}
