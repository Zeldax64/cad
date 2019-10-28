#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define THREADS 4
#define NUM_KEYS 80
#define MAX_KEY 30

void rank(int* key_array);
void rand_key(int* key_array, int size);
void print_keys(int* keys, int size);

int main () {
	int size = NUM_KEYS;
	
	int *keys;
	keys = (int*) malloc(size*sizeof(int));
	rand_key(keys, THREADS*NUM_KEYS);

	print_keys(keys, size); // Print iput array.
	rank(keys);
	print_keys(keys,size); // Print output array.

	free(keys);
}

void rank(int* key_array) {
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
		for( i=0; i<NUM_KEYS; i++ ) {
			work_buff[key_array[i]]++;
		}

		for( i=0; i<MAX_KEY-1; i++ ) 
			work_buff[i+1] += work_buff[i];

		#pragma omp barrier
		for( myid=1; myid<num_procs; myid++ ) {
			#pragma omp for nowait
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
}

// Helper functions.
void rand_key(int* key_array, int size) {
	for(int i = 0; i < size; ++i) {
		key_array[i] = rand() % MAX_KEY;
	}
}

void print_keys(int* keys, int size) {
	for(int i = 0; i < size; ++i) {
		printf("%d ", keys[i]);
	}
	printf("\n");
}
