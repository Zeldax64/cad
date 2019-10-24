#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "linear.h"

typedef struct block_t {
	uint32_t i, j;
	int done;
} block_t;

typedef struct mm_t {
	uint32_t blk_height, blk_width;
	fmat_t *a, *b, *res;
	block_t *blocks;
} mm_t;

// Master
fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);
mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width);
void init_workers();
int  has_worker();
void receive_block(mm_t* mm, int src);
void wait_any(mm_t* mm);
void wait_all(mm_t* mm);
void update_status();
void create_task();
void send_task(mm_t* mm, uint32_t i, uint32_t j, int dest); 
void finalize(mm_t* mm);

// Slave
void slave_loop();
int has_task();
void recv_task();

void free_mm(mm_t* mm);
