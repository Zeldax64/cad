#include "blockmul.h"

fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t b_transposed = b->transposed;
	
	mm_t* tasks = init_mm(a, b, blk_height, blk_width);
	
	if(!b_transposed) {
		transpose_fmat(b);
	}

	for(uint32_t i = 0; i < tasks->res->lines; i+=tasks->blk_height) {
		for(uint32_t j = 0; j < tasks->res->cols; j+=tasks->blk_width) {
			create_task();
			send_task(tasks, i, j);
		}
	}

	finalize();

	if(!b_transposed) {
		transpose_fmat(b);
	}

	free_mm(tasks);
	return tasks->res;
}

mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t res_lines, res_cols;

	res_lines = a->lines;
	res_cols = b->cols;

	mm_t* mm = (mm_t*) malloc(sizeof(mm_t));
	mm->a = a;
	mm->b = b;
	mm->res = init_fmat(res_lines, res_cols);
	
	// Assert block height and width.
	if(blk_height > a->lines || blk_height > b->lines) {
		mm->blk_height = a->lines < b->lines? a->lines : b->lines;
	}
	else {
		mm->blk_height = blk_height;
	}
	if(blk_width > a->cols || blk_width > b->cols) {
		mm->blk_width = a->cols < b->cols? a->cols : b->lines;
	}
	else {
		mm->blk_width = blk_width;	
	}

	return mm;
}

void free_mm(mm_t* mm) {
	free(mm);
}


void create_task() {
	int dest = 1;
	int task = 1;
	MPI_Send(&task, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
}


void finalize() {
	int dest = 1;
	int task = 0;
	MPI_Send(&task, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
}

void send_task(mm_t* mm, uint32_t i, uint32_t j) {
	uint32_t blk_height, blk_width;
	fmat_t* block;

	float **a_slc, **b_slc;

	a_slc = &mm->a->mat[i];
	b_slc = &mm->b->mat[j];

	// Assert block height and width.
	blk_height = (mm->res->lines-i) > mm->blk_height ? mm->blk_height : mm->res->lines-i;
	blk_width = (mm->res->cols-j) > mm->blk_width ? mm->blk_width : mm->res->cols-j;

	// Create slices of matrices A and B.
	fmat_t *a, *b;

	a = init_fmat(blk_height, mm->a->cols);
	b = init_fmat(blk_width, mm->b->cols);

	for(uint32_t line = 0; line < blk_height; ++line) {
		for(uint32_t col = 0; col < a->cols; ++col) {
			a->mat[line][col] = a_slc[line][col];
		}
	}

	for(uint32_t line = 0; line < blk_width; ++line) {
		for(uint32_t col = 0; col < b->cols; ++col) {
			b->mat[line][col] = b_slc[line][col];
		}
	}

	// Send data and wait response.
	MPI_Send(&blk_height, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	MPI_Send(&blk_width, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

	MPI_Send_fmat(a, 1);
	MPI_Send_fmat(b, 1);

	block = MPI_Recv_fmat(1);

	set_blk(i, j, blk_height, blk_width, block->mat, mm->res->mat);
	
	free_fmat(block);	
}

int has_task() {
	int has_task;
	MPI_Status stat;

	MPI_Recv(&has_task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	
	return has_task;	
}

void recv_task() {
	uint32_t blk_height, blk_width;
	MPI_Status stat;

	fmat_t *a, *b, *block;
	MPI_Recv(&blk_height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	MPI_Recv(&blk_width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);

	a = MPI_Recv_fmat(0);
	b = MPI_Recv_fmat(0);

	block = blk_mul(a->mat, b->mat, blk_height, blk_width, a->cols);
	MPI_Send_fmat(block, 0);

	free_fmat(block);
}
