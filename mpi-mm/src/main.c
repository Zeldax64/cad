#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "linear.h"
#include "blockmul.h"

#define N_SIZE 10

// Reference multiply
void multiply(uint32_t m, uint32_t n, uint32_t p, float **a, float **b, float **c);
void master();
void slave();


int main() {
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if(world_rank == 0) {
		master();
	}
	else {
		slave();
	}

	// Finalize the MPI environment.
	MPI_Finalize();

	return 0;
}

int verify(fmat_t *A, fmat_t *B, fmat_t *C) {
	fmat_t *C_ref = init_fmat(A->lines, B->cols);
	
	multiply(A->lines, A->cols, B->cols, A->mat, B->mat, C_ref->mat);
	
	if(!cmp_fmat(C, C_ref)) {
		printf("Matrices aren't equal!\n");
		
		printf("Res\n");
		print_fmat(C);
		printf("------\n");
		printf("Reference\n");
		print_fmat(C_ref);
		free_fmat(C_ref);

		return 1;
	}

	free_fmat(C_ref);
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

void master() {
	int world_size;
	fmat_t *A = init_fmat(N_SIZE, N_SIZE);
	fmat_t *B = init_fmat(N_SIZE, N_SIZE);

	fmat_t *C = NULL;

	rand_fmat(A);
	rand_fmat(B);

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	C = mpi_blk_mul(A, B, 10, 10);

	verify(A, B, C);

	printf("Master exit\n");

	free_fmat(A);
	free_fmat(B);
	free_fmat(C);
}

void slave() {
	int my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if(my_rank == 1) {
		while(1) {
			if(has_task()) {
				recv_task();
			}
			else {
				break;
			}
		}
	}
}
