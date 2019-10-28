#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

// Number of OpenMP threads used
#define THREADS 4
// Size of the array to be sorted
#define NUM_KEYS 20
// Maximum value+1 found in the array. 
#define MAX_KEY 5

void rank(int* key_array, const int SIZE);
void rand_key(int* key_array, int size);
void print_keys(int* keys, int size);
void parallel_buckets(int size, int *key_array, int *buckets);
void unbucket(int *keys, int *buckets);

int main () {
	int size = NUM_KEYS;
	int *keys = NULL;
	int *buckets = NULL;
	int *reduced_buckets = NULL;
	int world_rank, world_size, elements_per_proc;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if(world_rank == 0) {
		keys = (int*)malloc(size*sizeof(int));
		rand_key(keys, size);
		printf("Input: ");		
		print_keys(keys, size);
		reduced_buckets = (int*)malloc(MAX_KEY*sizeof(int)); 
	}

	buckets = (int*)malloc(MAX_KEY*sizeof(int)); 
	
	elements_per_proc = NUM_KEYS/world_size;
	int *sub_keys = (int*)malloc(elements_per_proc*sizeof(int));

	// Divide input data among MPI processes.
	MPI_Scatter(keys, 
				elements_per_proc,
				MPI_INT,
				sub_keys,
				elements_per_proc,
				MPI_INT,
				0,
				MPI_COMM_WORLD);

	// Compute buckets.
	parallel_buckets(elements_per_proc, sub_keys, buckets);

	// Return buckets to process 0.
	MPI_Reduce(
		buckets,
		reduced_buckets,
		MAX_KEY,
		MPI_INT,
		MPI_SUM,
		0,
		MPI_COMM_WORLD
		);

	// Process 0 finishes the sorting.
	if(world_rank == 0) {
		unbucket(keys, reduced_buckets);
		printf("Unbucket: ");
		print_keys(keys, size);
		
		free(keys);
	}
	
	free(sub_keys);
	MPI_Finalize();
}

// OpenMP bucket sort.
void rank(int* key_array, const int SIZE) {
	int i, k;
	int *key_buff_ptr;

	int** key_buff1;

	key_buff1 = (int**)malloc(THREADS*sizeof(int*));

	for(int i = 0; i < THREADS; ++i) {
		key_buff1[i] = (int*) calloc(MAX_KEY, sizeof(int));
	}
	key_buff_ptr = key_buff1[0];

	#pragma omp parallel private(i,k)
	{
		int *work_buff;
		int myid = 0, num_procs = 1;

		myid = omp_get_thread_num();
		num_procs = omp_get_num_threads();

		work_buff = key_buff1[myid];

		for( i=0; i<MAX_KEY; i++ ) {
			work_buff[i] = 0;
		}

		#pragma omp for nowait
		for( i=0; i<SIZE; i++ ) {
			work_buff[key_array[i]]++;
		}

		for( i=0; i<MAX_KEY-1; i++ ) 
			work_buff[i+1] += work_buff[i];

		#pragma omp barrier
		for( myid=1; myid<num_procs; myid++ ) {
			#pragma omp for //nowait
			for( i=0; i<MAX_KEY; i++ )
				key_buff_ptr[i] += key_buff1[myid][i];
		}

		#pragma omp for
		for( k=0; k<MAX_KEY; k++ ) {
			i = (k==0)? 0 : key_buff_ptr[k-1];
			while ( i<key_buff_ptr[k] )
				key_array[i++] = k;
		}
	} /*omp parallel*/

	// Free data.
	for(int i = 0; i < THREADS; ++i) {
		free(key_buff1[i]);
	}
	free(key_buff1);
}

// Generate a random array of a given size.
void rand_key(int* key_array, int size) {
	for(int i = 0; i < size; ++i) {
		key_array[i] = rand() % MAX_KEY;
	}
}

// Print array of a given size.
void print_keys(int* keys, int size) {
	for(int i = 0; i < size; ++i) {
		printf("%d ", keys[i]);
	}
	printf("\n");
}

// Put array elements into buckets and return them.
void parallel_buckets(int size, int *key_array, int *buckets) {
	int i;
	int *key_buff_ptr;

	int** key_buff1;

	key_buff1 = (int**)malloc(THREADS*sizeof(int*));
	for(int i = 0; i < THREADS; ++i) {
		key_buff1[i] = (int*)calloc(MAX_KEY, sizeof(int));
	}
	key_buff_ptr = key_buff1[0];

	#pragma omp parallel private(i)
	{
		int *work_buff;
		int myid = 0, num_procs = 1;

		myid = omp_get_thread_num();
		num_procs = omp_get_num_threads();

		work_buff = key_buff1[myid];

		for( i=0; i<MAX_KEY; i++ ) {
			work_buff[i] = 0;
		}

		#pragma omp for nowait
		for( i=0; i<size; i++ ) {
			work_buff[key_array[i]]++;
		}

		for( i=0; i<MAX_KEY-1; i++ ) 
			work_buff[i+1] += work_buff[i];

		#pragma omp barrier
		for( myid=1; myid<num_procs; myid++ ) {
			#pragma omp for
			for( i=0; i<MAX_KEY; i++ )
				key_buff_ptr[i] += key_buff1[myid][i];
		}
	} /*omp parallel*/

	// Copy return.
	for(i = 0; i < MAX_KEY; ++i) {
		buckets[i] = key_buff_ptr[i];
	}

	// Free data.
	for(int i = 0; i < THREADS; ++i) {
		free(key_buff1[i]);
	}
	free(key_buff1);
}

// Unbucket array elements and finish sorting.
void unbucket(int *keys, int *buckets) {
	int i, k;

	#pragma omp parallel for private(i, k)
	for(k = 0; k < MAX_KEY; k++) {
		i = (k==0)? 0 : buckets[k-1];
		while (i<buckets[k])
			keys[i++] = k;
	}
}