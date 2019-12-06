#include "linear.h"

fmat_t* init_fmat(uint32_t m, uint32_t n) {
	fmat_t* fmat = (fmat_t*) malloc (sizeof(fmat_t));

	fmat->lines = m;
	fmat->cols = n;
	fmat->transposed = false;
	fmat->mat = zero_mat_f(m, n);

	return fmat;
}

void free_fmat(fmat_t* mat) {
	free_mat(mat->lines, mat->mat);
	free(mat);
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

void transpose_fmat(fmat_t* fmat) {
	uint32_t temp;
	float** t_mat = zero_mat_f(fmat->cols, fmat->lines);

	for(uint32_t i = 0; i < fmat->lines; ++i) {
		for(uint32_t j = 0; j < fmat->cols; ++j) {
			t_mat[j][i] = fmat->mat[i][j];
		}
	}

	free_mat(fmat->lines, fmat->mat);
	fmat->mat = t_mat;

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
				return false;
			}
		}
	}

	return true;
}

fmat_t* slice_fmat_lines(fmat_t* fmat, uint32_t lines) {
	fmat_t* slice = init_fmat(lines, fmat->cols);

	return slice;
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

// Float matrix utilities

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


void free_mat(uint32_t n, float** mat) {
	for(uint32_t i = 0; i < n; i++) {
		free(mat[i]);
	}
	free(mat);
	mat = NULL;
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

fmat_t* blk_mul(float **a, float **b, uint32_t blk_x, uint32_t blk_y, uint32_t arr_size) {
	fmat_t* block = init_fmat(blk_x, blk_y);

#pragma omp parallel for schedule(dynamic) num_threads(OPENMP_THREADS)
	for(uint32_t i = 0; i < blk_x; ++i) {
		for(uint32_t j = 0; j < blk_y; ++j) {
			//printf("i = %d j = %d\n", i, j);
			block->mat[i][j] = mul_arr(arr_size, a[i], b[j]);
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


// MPI functions
void MPI_Send_fmat(fmat_t* fmat, int dest) {
	MPI_Send(&fmat->lines, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(&fmat->cols, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

	for(uint32_t i = 0; i < fmat->lines; ++i) {
		MPI_Send(&fmat->mat[i][0], fmat->cols, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
	}
}

fmat_t* MPI_Recv_fmat(int src) {
	MPI_Status stat;

	fmat_t* fmat = init_fmat(1,1);
	free_mat(fmat->lines, fmat->mat);

	MPI_Recv(&fmat->lines, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &stat);
	MPI_Recv(&fmat->cols, 1, MPI_INT, src, 0, MPI_COMM_WORLD, &stat);
	
	float** new_mat = zero_mat_f(fmat->lines, fmat->cols);

	for(uint32_t i = 0; i < fmat->lines; ++i) {
		MPI_Recv(&new_mat[i][0], fmat->cols, MPI_FLOAT, src, 0, MPI_COMM_WORLD, &stat);
	}

	fmat->mat = new_mat;

	return fmat;
}
