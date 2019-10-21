#include "linear.h"

fmat_t* init_fmat(uint32_t m, uint32_t n) {
	fmat_t* fmat = (fmat_t*) malloc (sizeof(fmat_t));

	fmat->lines = m;
	fmat->cols = n;
	fmat->transposed = false;
	fmat->mat = zero_mat_f(m, n);

	return fmat;
}

float* zero_arr_f(uint32_t n) {
	float* arr;
	arr = (float *) calloc(n, sizeof(float));
	return arr;	
}

float* rand_arr_f(uint32_t n) {
	float* arr;
	
	arr = (float *) malloc(n * sizeof(float));
	
	for(uint32_t i = 0; i < n; ++i) {
		arr[i] = rand();
	}

	return arr;
}

float** zero_mat_f(uint32_t m, uint32_t n) {
	float** mat;

	mat = (float**) malloc(m * sizeof(float*));

	for(uint32_t i = 0; i < m; ++i) {
		mat[i] = zero_arr_f(n);
	}

	return mat;
}

float** rand_mat_f(uint32_t m, uint32_t n) {
	float** mat;

	mat = (float**) malloc(m * sizeof(float*));

	for(uint32_t i = 0; i < m; ++i) {
		mat[i] = rand_arr_f(n);
	}

	return mat;
}

void free_fmat(fmat_t* mat) {
	free_mat(mat->lines, mat->mat);
	free(mat);
}

// Private
void free_mat(uint32_t n, float** mat) {
	for(uint32_t i = 0; i < n; i++) {
		free(mat[i]);
	}
	free(mat);
}

void print_fmat(fmat_t* mat) {
	uint32_t m = mat->lines;
	uint32_t n = mat->cols;

	for(uint32_t i = 0; i < m; ++i) {
		for(uint32_t j = 0; j < n; ++j) {
			printf("%f ", mat->mat[i][j]);
		}
		printf("\n");
	}	
}

fmat_t* mul_fmat(fmat_t* a, fmat_t* b) {
	bool is_b_transposed = b->transposed;
	fmat_t* res = init_fmat(a->cols, b->lines);

	if(is_b_transposed == false) {
		transpose_fmat(b);
	}

	for(uint32_t i = 0; i < a->lines; ++i) {
		for(uint32_t j = 0; j < b->cols; ++j) {
			res->mat[i][j] = mul_arr(res->cols, a->mat[i], b->mat[j]);
		}
	}

	if(is_b_transposed == false) {
		transpose_fmat(b);
	}

	return res;
}

void rand_fmat(fmat_t* mat) {
	free_mat(mat->lines, mat->mat);
	mat->mat = rand_mat_f(mat->lines, mat->cols);
}

// Deprecated
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


void transpose_fmat(fmat_t* fmat) {
	float aux;
	uint32_t temp;

	for(uint32_t i = 0; i < fmat->lines; ++i) {
		for(uint32_t j = 0; j < i; ++j) {
			aux = fmat->mat[i][j];
			fmat->mat[i][j] = fmat->mat[j][i];
			fmat->mat[j][i] = aux;
		}
	}

	temp = fmat->lines;
	fmat->lines = fmat->cols;
	fmat->cols = temp;
	fmat->transposed = !fmat->transposed;
}

bool cmp_fmat(fmat_t* a, fmat_t* b) {
	if(a->lines != b->lines || a->cols != b->cols) {
		return false;
	}

	for(uint32_t i = 0; i < a->lines; ++i) {
		for(uint32_t j = 0; j < b->cols; ++j) {
			if(a->mat[i][j] != b->mat[i][j]) {
				return true;
			}
		}
	}

	return false;
}


// Deprecated
//float** blk_mul(uint32_t n, uint32_t blk_x, uint32_t blk_y, float **a, float **b) {
//	float **block;
//
//	block = (float**) malloc(blk_x * sizeof(float*));
//	for(uint32_t i = 0; i < blk_x; ++i) {
//		block[i] = (float*) calloc(blk_y, sizeof(float*));
//	}
//
//	for(uint32_t i = 0; i < blk_x; ++i) {
//		for(uint32_t j = 0; j < blk_y; ++j) {
//			block[i][j] = mul_arr(n, a[i], b[j]);
//		}
//	}
//
//	return block;
//}

fmat_t* blk_mul(float **a, float **b, uint32_t blk_x, uint32_t blk_y, uint32_t arr_size) {
	fmat_t* block = init_fmat(blk_x, blk_y);

	for(uint32_t i = 0; i < blk_x; ++i) {
		for(uint32_t j = 0; j < blk_y; ++j) {
			block->mat[i][j] = mul_arr(arr_size, a[i], b[j]);
		}
	}

	return block;
}

// Deprecated
void set_blk(uint32_t i, uint32_t j, uint32_t blk_x, uint32_t blk_y, float **blk, float **c) {
	for(uint32_t k = 0; k < blk_x; k++) {
		for(uint32_t l = 0; l < blk_y; l++) {
			c[i+k][j+l] = blk[k][l];
		}
	}
}

// Deprecated
/*
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
*/