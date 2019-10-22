#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "linear.h"

typedef struct block_t {
	uint32_t a_size;
	uint32_t b_size;
	float** a_lines;
	float** b_lines;
	fmat_t* res;
} block_t;


typedef struct mm_t {
	block_t** blocks;
	uint32_t blk_height, blk_width;
	fmat_t *a, *b, *res;	
} mm_t;

fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);

mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);
block_t* init_block(uint32_t a_size, uint32_t b_size);

void send_task(mm_t* mm, uint32_t i, uint32_t j); 

void free_mm(mm_t* mm);
void free_block(block_t* block);