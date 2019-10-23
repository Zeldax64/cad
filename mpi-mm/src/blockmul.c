#include "blockmul.h"

fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t b_transposed = b->transposed;
	mm_t* tasks = init_mm(a, b, blk_height, blk_width);
	
	if(!b_transposed) {
		transpose_fmat(b);
	}

	//print_fmat(tasks->res);
	for(uint32_t i = 0; i < tasks->res->lines; i+=tasks->blk_height) {
		for(uint32_t j = 0; j < tasks->res->cols; j+=tasks->blk_width) {
			send_task(tasks, i, j);
		}
	}

	if(!b_transposed) {
		transpose_fmat(b);
	}

	// free tasks;
	return tasks->res;
}

mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t res_lines, res_cols;
	//uint32_t blk_lines, blk_cols;
	res_lines = a->lines;
	res_cols = b->cols;

	mm_t* mm = (mm_t*) malloc(sizeof(mm_t));
	mm->a = a;
	mm->b = b;
	mm->res = init_fmat(res_lines, res_cols);
	
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

void send_task(mm_t* mm, uint32_t i, uint32_t j) {
	uint32_t blk_height, blk_width;
	fmat_t* block;

	// Assert block height and width.
	blk_height = (mm->res->lines-i) > mm->blk_height ? mm->blk_height : mm->res->lines-i;
	blk_width = (mm->res->cols-j) > mm->blk_width ? mm->blk_width : mm->res->cols-j;
	
	// Send stub for matrix A...
	//for(uint32_t height = 0; height < blk_height; ++height) {
		//MPI_Send(mm->a->mat[height+i]);
	//}

	// Send stub for matrix B...
	//for(uint32_t width = 0; width < blk_width; ++width) {
		//MPI_Send(mm->b->mat[width+j]);
	//}

	//TODO: receive calculated matrix.

	//TODO: copy received matrix to mm->res.
	// Here is a temporary workaround while not using MPI.

	block = blk_mul(&mm->a->mat[i], &mm->b->mat[j], blk_height, blk_width, mm->a->cols);
	set_blk(i, j, blk_height, blk_width, block->mat, mm->res->mat);
	
	free_fmat(block);	

}

// Not used
//block_t* init_block(uint32_t a_size, uint32_t b_size) {}

void free_mm(mm_t* mm) {

}

void free_block(block_t* block) {

}