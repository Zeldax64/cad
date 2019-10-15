#include "linear.h"

float* rand_arr_f(uint32_t n) {
	float* arr;
	
	arr = (float *) malloc(n * sizeof(float));
	
	for(uint32_t i = 0; i < n; ++i) {
		arr[i] = rand();
	}

	return arr;
}

float** rand_mat_f(uint32_t n) {
	float** mat;

	mat = (float**) malloc(n * sizeof(float*));

	for(uint32_t i = 0; i < n; ++i) {
		mat[i] = rand_arr_f(n);
	}

	return mat;
}

void free_mat(uint32_t n, float** mat) {
	for(uint32_t i = 0; i < n; i++) {
		free(mat[i]);
	}
	free(mat);
}

void print_mat(uint32_t n, float** mat) {
	for(uint32_t i = 0; i < n; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			printf("%f ", mat[i][j]);
		}
		printf("\n");
	}
}

void mul_mat(uint32_t n, float **a, float **b, float **c) {
	for(uint32_t i = 0; i < n; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			c[i][j] = mul_arr(n, a[i], b[j]);
		}
	}
}

float mul_arr(uint32_t n, float *a, float *b) {
	float res = 0.0f;
	for(uint32_t i = 0; i < n; ++i) {
		res += a[i]*b[i];
	}
	return res;
}

void transpose_mat(uint32_t n, float **mat) {
	float aux;
	
	for(uint32_t i = 0; i < n; ++i) {
		for(uint32_t j = 0; j < i; ++j) {
			aux = mat[i][j];
			mat[i][j] = mat[j][i];
			mat[j][i] = aux;
		}
	}
}


int cmp_mat(uint32_t n, float **a, float **b) {
	for(uint32_t i = 0; i < n; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			if(a[i][j] != b[i][j]) {
				return 1;
			}
		}
	}

	return 0;
}


float** blk_mul(uint32_t n, uint32_t blk_x, uint32_t blk_y, float **a, float **b) {
	float **block;

	block = (float**) malloc(blk_x * sizeof(float*));
	for(uint32_t i = 0; i < blk_x; ++i) {
		block[i] = (float*) calloc(blk_y, sizeof(float*));
	}

	for(uint32_t i = 0; i < blk_x; ++i) {
		for(uint32_t j = 0; j < blk_y; ++j) {
			block[i][j] = mul_arr(n, a[i], b[j]);
		}
	}

	return block;
}

void set_blk(uint32_t i, uint32_t j, uint32_t blk_x, uint32_t blk_y, float **blk, float **c) {
	for(uint32_t k = 0; k < blk_x; k++) {
		for(uint32_t l = 0; l < blk_y; l++) {
			c[i+k][j+l] = blk[k][l];
		}
	}
}

void mul_mat_blk(uint32_t n, uint32_t blk_size, float **a, float **b, float **c) {
	uint32_t blk_x, blk_y;
	float **block;

	for(uint32_t i = 0; i < n; i += blk_size) {
		blk_x = (n-i) > blk_size ? blk_size : n-i; 
		for(uint32_t j = 0; j < n; j += blk_size) {
			blk_y = (n-j) > blk_size ? blk_size : n-j; 
			
			block = blk_mul(n, blk_x, blk_y, &a[i], &b[j]);
			
			set_blk(i, j, blk_x, blk_y, block, c);
			free_mat(blk_x, block);
		}
	}
}