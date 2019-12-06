CC = mpicc
CFLAGS = -Wall -O3 -g -fopenmp
PROG = mpi-mm.elf
SRCDIR = src
PROCS = 2#$$(nproc)
THREADS = 2

sources = $(wildcard $(SRCDIR)/*.c)
headers = $(wildcard $(SRCDIR)/*.h)
objects = $(patsubst %.c, %.o, $(sources))

run: $(PROG)
	#@echo "Program Start:"
	mpiexec -n $(PROCS) ./$(PROG) $(THREADS)

$(PROG): $(objects)
	#echo "Compiling..."
	$(CC) -o $(PROG) $(objects) $(CFLAGS) 	

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(CFLAGS) 

.PHONY: clean

clean:
	rm -rf $(PROG) $(objects)
