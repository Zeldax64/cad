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


//	uint32_t n = N_SIZE;
//
//	fmat_t* A = init_fmat(5, 1);
//	fmat_t* B = init_fmat(1, 5);
//	fmat_t* C, *C_ref;
//
//	rand_fmat(A);
//	rand_fmat(B);
//
//	C = mpi_blk_mul(A, B, 3, 3);
//
//	C_ref = init_fmat(5, 5);
//	multiply(A->lines, A->cols, B->cols, A->mat, B->mat, C_ref->mat);
//	
//	if(!cmp_fmat(C, C_ref)) {
//		printf("Matrices aren't equal!\n");
//	}
//
//	free_fmat(A);
//	free_fmat(B);
//	free_fmat(C);
//	free_fmat(C_ref);

//	return 0;
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
	fmat_t *fmat = init_fmat(3, 5);

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	for(int i = 1; i < world_size; ++i) {
		MPI_Send_fmat(fmat, i);
	}
}

void slave() {
	fmat_t* fmat;

	fmat = MPI_Recv_fmat(0);

	print_fmat(fmat);
}
