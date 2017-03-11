#
# linux & Mac
#LIBS = -llapack -lblas -lm
LIBS = -lm


CFLAGS = -O2
CC = gcc -c
LD = gcc
#CFLAGS = -fast
#CC = icc -c
#LD = icc


HEADERS = treeGlobal.h  treecode.h
OBJ = readin.o treeGlobal.o treecode.o

all: coulomb

%.o: %.c $(HEADERS) Makefile
	$(CC) $(CFLAGS) -o $@ $<

coulomb: coulomb.o $(OBJ)
	$(LD) -o coulomb coulomb.o $(OBJ) $(LIBS) $(LIBSLA)

paratest: paratest.o $(OBJ)
	$(LD) -o paratest paratest.o $(OBJ) $(LIBS) $(LIBSLA)

clean:
	\rm -f *.o
