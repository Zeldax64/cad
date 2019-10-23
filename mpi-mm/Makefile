CC = mpicc
CFLAGS = -Wall -O0 -g
PROG = mpi-mm.elf
SRCDIR = src

sources = $(wildcard $(SRCDIR)/*.c)
headers = $(wildcard $(SRCDIR)/*.h)
objects = $(patsubst %.c, %.o, $(sources))

run: $(PROG)
	#@echo "Program Start:"
	#./$(PROG)
	mpiexec -n $$(nproc) ./$(PROG)

$(PROG): $(objects)
	#echo "Compiling..."
	$(CC) -o $(PROG) $(objects) $(CFLAGS) 	

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(CFLAGS) 

.PHONY: clean

clean:
	rm -rf $(PROG) $(objects)
