#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "linear.h"
#include "blockmul.h"

#define N_SIZE 10

// Reference multiply
void multiply(uint32_t m, uint32_t n, uint32_t p, float **a, float **b, float **c);

int main() {
//	// Initialize the MPI environment
//	MPI_Init(NULL, NULL);
//
//	// Get the number of processes
//	int world_size;
//	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//	// Get the rank of the process
//	int world_rank;
//	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//
//	// Get the name of the processor
//	char processor_name[MPI_MAX_PROCESSOR_NAME];
//	int name_len;
//	MPI_Get_processor_name(processor_name, &name_len);
//
//	// Print off a hello world message
//	printf("Hello world from processor %s, rank %d out of %d processors\n",
//	processor_name, world_rank, world_size);
//
//	// Finalize the MPI environment.
//	MPI_Finalize();
//
//	return 0;


	uint32_t n = N_SIZE;
	/*
	float **A, **B, **C, **C_ref;

	A = rand_mat_f(n);
	B = rand_mat_f(n);

	C_ref = (float**) malloc(n * sizeof(float*));
	for(uint32_t i = 0; i < n; ++i) {
		C_ref[i] = (float*) calloc(n, sizeof(float*));
	}


	C = (float**) malloc(n * sizeof(float*));
	for(uint32_t i = 0; i < n; ++i) {
		C[i] = (float*) calloc(n, sizeof(float*));
	}

	multiply(n, n, n, A, B, C_ref);

	transpose_mat(n, B);

	//mul_mat(n, A, B, C);
	mul_mat_blk(n, 1, A, B, C);
	int exit = cmp_mat(n, C, C_ref);

	if(exit != 0) {
		printf("Matrizes diferentes!\n");
	}

	free_mat(n, A);
	free_mat(n, B);
	free_mat(n, C_ref);

	return exit;
	*/


	fmat_t* A = init_fmat(5, 1);
	fmat_t* B = init_fmat(1, 5);
	fmat_t* C;

	rand_fmat(A);
	rand_fmat(B);
//
//	print_fmat(A);
//	print_fmat(B);
//
//
//	free_fmat(A);
//	free_fmat(B);
	C = mpi_blk_mul(A, B, 3, 3);
	print_fmat(C);
	return 0;
}

void multiply(uint32_t m, uint32_t n, uint32_t p, float **a, float **b, float **c) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

