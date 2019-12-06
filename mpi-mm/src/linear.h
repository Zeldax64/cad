#pragma once

#include <stdio.h>

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>

#include <mpi.h>

#include "global.h"

typedef struct farr {
	uint32_t size;

	float *arr;
} farr_t;

typedef struct fmat {
	uint32_t lines;
	uint32_t cols;
	bool transposed;

	float **mat;

} fmat_t;

// Public
fmat_t* init_fmat(uint32_t m, uint32_t n);
void free_fmat(fmat_t* mat);
fmat_t* mul_fmat(fmat_t* a, fmat_t* b);
void rand_fmat(fmat_t* mat);
void transpose_fmat(fmat_t* mat);
void print_fmat(fmat_t* mat);
fmat_t* slice_fmat_lines(fmat_t* fmat, uint32_t lines);
bool cmp_fmat(fmat_t* a, fmat_t* b);

// Float matrix utilities
float* zero_arr_f(uint32_t n);
float* rand_arr_f(uint32_t n);
float** zero_mat_f(uint32_t m, uint32_t n);
float** rand_mat_f(uint32_t m, uint32_t n);
void free_mat(uint32_t n, float** mat);
float mul_arr(uint32_t n, float *a, float *b);

// Matrix operations
void mul_mat(uint32_t n, float **a, float **b, float **c);
void mul_mat_blk(uint32_t n, uint32_t blk_size, float **a, float **b, float **c);
fmat_t* blk_mul(float **a, float **b, uint32_t blk_x, uint32_t blk_y, uint32_t arr_size);
void set_blk(uint32_t i, uint32_t j, uint32_t blk_x, uint32_t blk_y, float **blk, float **c); 

// MPI functions
void MPI_Send_fmat(fmat_t* fmat, int dest);
fmat_t* MPI_Recv_fmat(int src);
