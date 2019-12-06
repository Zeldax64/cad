CC = mpicc
CFLAGS = -Wall -O3 -fopenmp
PROG = mpi-mm.elf
SRCDIR = src
PROCS = 1
THREADS = 1

sources = $(wildcard $(SRCDIR)/*.c)
headers = $(wildcard $(SRCDIR)/*.h)
objects = $(patsubst %.c, %.o, $(sources))

run: $(PROG)
	mpiexec -n $(PROCS) ./$(PROG) $(THREADS)

$(PROG): $(objects)
	$(CC) -o $(PROG) $(objects) $(CFLAGS) 	

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
