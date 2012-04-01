CC = g++


OBJS = MethodA.o MethodB.o MethodC.o MethodD.o Alignment.o AlignmentReader.o bitfile.o utils.o

all: readzip

readzip: $(OBJS) readzip.o
	$(CC) -o readzip readzip.o $(OBJS)



clean:
	rm -f core *.o *~ readzip
