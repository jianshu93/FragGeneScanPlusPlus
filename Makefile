CC=gcc
CFLAGS=-O3 -m64 -w -lm -lpthread
LDFLAGS=-lm -lpthread

SRC = fasta.c hmm_lib.c run_hmm.c util_lib.c
OBJ = ${SRC:.c=.o}

all: FGSpp

FGSpp: ${OBJ}
	$(CC) -o $@ ${OBJ} ${LDFLAGS}

.c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o FGSpp

.PHONY: clean all
