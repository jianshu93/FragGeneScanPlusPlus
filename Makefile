CC=	gcc 
CFLAG= -O3 -march=native -Wall
FLAGS= -lm -lpthread
HEADER=	util_lib.h fasta.h run_hmm.h
SRCS=	util_lib.c hmm_lib.c run_hmm.c fasta.c 
OBJ=	util_lib.o hmm_lib.o run_hmm.o  fasta.o 
EXEC=  FGS++
ASTYLE_FLAGS= --mode=c --style=google -s4 -n -H -k3

all:  $(OBJ) $(HEADER)
	$(CC)  $(CFLAG) -o $(EXEC) $(OBJ) $(FLAGS)

%.o:%.c
	$(CC) $(CFLAG) -c $< 

clean:
	rm -rf $(OBJ) $(EXEC)

astyle:
	astyle $(ASTYLE_FLAGS) $(SRCS) $(HEADER)

doc: $(HEADER)
	doxygen
