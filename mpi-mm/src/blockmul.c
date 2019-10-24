#include "blockmul.h"

static const int MASTER_RANK = 0;
static int WORLD_SIZE = 0;
static MPI_Request *workers = NULL;
static int workers_free = 0;


// Master
fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t b_transposed = b->transposed;
	
	init_workers();
	mm_t* tasks = init_mm(a, b, blk_height, blk_width);
	
	if(!b_transposed) {
		transpose_fmat(b);
	}

	for(uint32_t i = 0; i < tasks->res->lines; i+=tasks->blk_height) {
		for(uint32_t j = 0; j < tasks->res->cols; j+=tasks->blk_width) {
			if(has_worker()) {
				create_task();
				int dest = 1;
				send_task(tasks, i, j, dest);
			}
			else {
				wait_any(tasks);
				j -= tasks->blk_width;
			}
		}
	}

	if(!b_transposed) {
		transpose_fmat(b);
	}

	finalize(tasks);
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

	mm->blocks = (block_t*) malloc(WORLD_SIZE*sizeof(block_t));

	return mm;
}

void free_mm(mm_t* mm) {
	free(mm->blocks);
	free(mm);
}

void init_workers() {
	//MPI_Comm_size(MPI_COMM_WORLD, &WORLD_SIZE);
	WORLD_SIZE = 2;
	workers = (MPI_Request *) malloc(WORLD_SIZE * sizeof(MPI_Request));
	workers_free = WORLD_SIZE-1;
}

int has_worker() {
	return workers_free > 0? 1:0;
}

void receive_block(mm_t* mm, int src) {
	fmat_t *block;

	block = MPI_Recv_fmat(src);
	set_blk(mm->blocks[src].i, mm->blocks[src].j, block->lines, block->cols, block->mat, mm->res->mat);
	mm->blocks[src].done = 1;

	free_fmat(block);	
}

void wait_any(mm_t* mm) {
	int worker = 0;
	MPI_Status stat;

	MPI_Waitany(WORLD_SIZE-1, &workers[1], &worker, &stat);
	receive_block(mm, worker+1);

	workers_free++;
}

void wait_all(mm_t* mm) {
	int worker = 0;
	MPI_Status *stats = (MPI_Status*) malloc((WORLD_SIZE-1) * sizeof(MPI_Status));

	MPI_Waitall(WORLD_SIZE-1, &workers[1], stats);
	
	// Loop here
	for(uint32_t i = 1; i < WORLD_SIZE; ++i) {
		if(mm->blocks[i].done == 0) {
			receive_block(mm, i);
		}
	}

	workers_free = WORLD_SIZE-1;
}

void update_status() {

}

void create_task() {
	int dest = 1;
	int task = 1;
	MPI_Send(&task, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

	workers_free -= 1;
}

void send_task(mm_t* mm, uint32_t i, uint32_t j, int dest) {
	uint32_t blk_height, blk_width;
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
	MPI_Send(&blk_height, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(&blk_width, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

	MPI_Send_fmat(a, dest);
	MPI_Send_fmat(b, dest);

	mm->blocks[dest].i = i;
	mm->blocks[dest].j = j;

	MPI_Irecv(&mm->blocks[dest].done, 1, MPI_INT, dest, 0, MPI_COMM_WORLD, &workers[dest]);

	free_fmat(a);
	free_fmat(b);
}

void finalize(mm_t* mm) {
	int dest = 1;
	int task = 0;

	wait_all(mm);

	MPI_Send(&task, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	free(workers);
}

// Slave
void slave_loop() {
	while(1) {
		if(has_task()) {
			recv_task();
		}
		else {
			break;
		}
	}
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
	int done = 0;

	MPI_Recv(&blk_height, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD, &stat);
	MPI_Recv(&blk_width, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD, &stat);
	a = MPI_Recv_fmat(MASTER_RANK);
	b = MPI_Recv_fmat(MASTER_RANK);

	block = blk_mul(a->mat, b->mat, blk_height, blk_width, a->cols);

	// Tell master that the block is done.
	MPI_Send(&done, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD);

	// Send block.
	MPI_Send_fmat(block, MASTER_RANK);

	free_fmat(block);
}
