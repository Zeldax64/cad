#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

typedef struct farray {
	uint32_t size;

	float *arr;
} farr_t;

typedef struct fmatrix {
	uint32_t width;
	uint32_t height;

	float **mat;

} fmat_t;

float* rand_arr_f(uint32_t n);
float** rand_mat_f(uint32_t n);
void free_mat(uint32_t n, float** mat);
void transpose_mat(uint32_t n, float **mat);
void print_mat(uint32_t n, float** mat);
float mul_arr(uint32_t n, float *a, float *b);
void mul_mat(uint32_t n, float **a, float **b, float **c);
int cmp_mat(uint32_t n, float **a, float **b);
void mul_mat_blk(uint32_t n, uint32_t blk_size, float **a, float **b, float **c);
float** blk_mul(uint32_t n, uint32_t blk_x, uint32_t blk_y, float **a, float **b);
void set_blk(uint32_t i, uint32_t j, uint32_t blk_x, uint32_t blk_y, float **blk, float **c); 
