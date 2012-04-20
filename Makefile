CC = g++
CCFLAGS = -Os


OBJS = MethodA.o MethodB.o MethodC.o MethodD.o Alignment.o AlignmentReader.o bitfile.o utils.o

all: readzip

readzip: $(OBJS) readzip.o
	$(CC) $(CCFLAGS) -o readzip readzip.o $(OBJS)
MethodA.o:
	$(CC) $(CCFLAGS) -c MethodA.cpp 
MethodB.o:
	$(CC) $(CCFLAGS) -c MethodB.cpp 
MethodC.o:
	$(CC) $(CCFLAGS) -c MethodC.cpp 
MethodD.o:
	$(CC) $(CCFLAGS) -c MethodD.cpp 
utils.o:
	$(CC) $(CCFLAGS) -c utils.cpp 
AlignmentReader.o:
	$(CC) $(CCFLAGS) -c AlignmentReader.cpp 
Alignment.o:
	$(CC) $(CCFLAGS) -c Alignment.cpp 
bitfile.o:
	$(CC) $(CCFLAGS) -c bitfile.cpp 

clean:
	rm -f core *.o *~ readzip
