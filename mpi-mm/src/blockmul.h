#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "linear.h"

typedef struct mm_t {
	uint32_t blk_height, blk_width;
	fmat_t *a, *b, *res;	
} mm_t;

fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);
mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);

// Master
void create_task();
void finalize();
void send_task(mm_t* mm, uint32_t i, uint32_t j); 

// Slave
int has_task();
void recv_task();

void free_mm(mm_t* mm);
