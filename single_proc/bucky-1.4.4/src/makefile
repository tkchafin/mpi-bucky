#
# To compile, type "make" or make "all"
# To remove files, type "make clean"
#
OBJS = TaxonSet.o Alias.o mbsumtree.o TGM.o Quartets.o
TARGET = all

CC = g++
CFLAGS = -O3 -I.

.SUFFIXES: .C .o 

all: bucky mbsum

bucky: bucky.C bucky.h $(OBJS)
	$(CC) $(CFLAGS) -o bucky bucky.C $(OBJS)

mbsum: mbsum.C mbsumtree.o
	$(CC) $(CFLAGS) -o mbsum mbsum.C mbsumtree.o

.C.o:
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(OBJS) bucky mbsum
