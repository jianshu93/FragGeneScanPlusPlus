CC=gcc
CFLAGS=-g -m64 -w -lm -lpthread
LDFLAGS=-lm -lpthread

FGSpp: run_hmm.o hmm_lib.o util_lib.o fasta.o
	$(CC) $(CFLAGS) -o $@ $^

run_hmm.o: hmm.h run_hmm.h translation_tables.h util_lib.h
fasta.o: fasta.h
hmm_lib.o: util_lib.h
util_lib.o: util_lib.h

clean:
	rm -f *.o FGSpp

.PHONY: clean
