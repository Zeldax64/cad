CC = mpicc
CFLAGS = -Wall -O3 -fopenmp -march=native -funroll-loops
PROG = mpi-mm.elf
SRCDIR = src
PROCS = 4
THREADS = 4

sources = $(wildcard $(SRCDIR)/*.c)
headers = $(wildcard $(SRCDIR)/*.h)
objects = $(patsubst %.c, %.o, $(sources))

run: all
	mpiexec -n $(PROCS) ./$(PROG) $(THREADS)
	#./../copy-files $(PROG)
	#cp $(PROG) ./../../../$(PROG)
	#./exec.sh $(THREADS)

all: $(PROG)

$(PROG): $(objects)
	$(CC) -o $(PROG) $(objects) $(CFLAGS)
	objdump -D $(PROG) > mpi-mm.dump

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(CFLAGS)

.PHONY: clean valgrind

valgrind: $(PROG)
	valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=valgrind-out.txt \
         mpiexec -n $(PROCS) ./$(PROG) $(THREADS)

clean:
	rm -rf $(PROG) $(objects)
