#include "blockmul.h"

// Private variables
static const int MASTER_RANK = 0;
static int WORLD_SIZE = 0;
static MPI_Request *workers = NULL;
static int workers_free = 0;

// Private Master Functions
static void inc_workers_free() {
	workers_free++;
}

static void dec_workers_free() {
	workers_free--;
}

static int get_workers_free() {
	return workers_free;
}

// Master functions.
fmat_t* mpi_blk_mul(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	int worker;
	uint32_t b_transposed = b->transposed;
	int master_threads;
	
	init_workers(); // Init MPI workers.
	mm_t* tasks = init_mm(a, b, blk_height, blk_width);
	
	// Transpose matrix B to achieve better performance.
	if(!b_transposed) {
		transpose_fmat(b);
	}

	master_threads = (OPENMP_THREADS == 1 && WORLD_SIZE == 1) ? 2 : OPENMP_THREADS;
#pragma omp parallel num_threads(master_threads)
{
	#pragma omp single
	for(uint32_t i = 0; i < tasks->res->lines; i+=tasks->blk_height) {
		for(uint32_t j = 0; j < tasks->res->cols; j+=tasks->blk_width) {			
			if(has_worker()) { // Is any worker free?				
				// Get a free worker and give him a task.
				worker = get_free_worker(tasks);
				allocate_task(tasks, worker, i, j);
			}
			else {
				// Wait until a free worker appears.				
				wait_any(tasks);
				j -= tasks->blk_width;
			}
		}
	}
}

	// Transpose matrix B back if necessary.
	if(!b_transposed) {
		transpose_fmat(b);
	}

	// Finalize other MPI processes and free objects.
	finalize(tasks);
	free_mm(tasks);

	return tasks->res;
}

mm_t* init_mm(fmat_t* a, fmat_t* b, uint32_t blk_height, uint32_t blk_width) {
	uint32_t res_lines, res_cols;

	mm_t* mm = (mm_t*) malloc(sizeof(mm_t));

	// Set matrices in MM object
	mm->a = a;
	mm->b = b;
	res_lines = a->lines;
	res_cols = b->cols;
	mm->res = init_fmat(res_lines, res_cols);
	
	// Assert block height and width to avoid errors.
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

	// Init list of working blocks.
	mm->blocks = (block_t*) malloc(WORLD_SIZE*sizeof(block_t));
	for(int i = 0; i < WORLD_SIZE; ++i) {
		mm->blocks[i].status = DONE;
	}

	return mm;
}

void free_mm(mm_t* mm) {
	free(mm->blocks);
	free(mm);
}

void init_workers() {
	MPI_Comm_size(MPI_COMM_WORLD, &WORLD_SIZE);
	workers = (MPI_Request *) malloc(WORLD_SIZE * sizeof(MPI_Request));
	for(uint32_t i = 0; i < WORLD_SIZE; ++i) { 
		workers[i] = MPI_REQUEST_NULL;
	}
	workers_free = WORLD_SIZE;
}

int has_worker() {
	return workers_free > 0? 1:0;
}

// Iterate trhough the MM object and return a free worker.
int get_free_worker(mm_t *mm) {
	int worker = 0;

	while(mm->blocks[worker].status != DONE) {
		worker=worker+1;
	}

	return worker;
}

void receive_block(mm_t* mm, int src) {
	fmat_t *block;

	block = MPI_Recv_fmat(src);
	set_blk(mm->blocks[src].i, mm->blocks[src].j, block->lines, block->cols, block->mat, mm->res->mat);
	mm->blocks[src].status = DONE;

	free_fmat(block);	
}

void wait_any(mm_t* mm) {
	int worker = 0;
	MPI_Status stat;

	MPI_Waitany(WORLD_SIZE, &workers[0], &worker, &stat);

	if(worker == MASTER_RANK) {
	
	}
	else {
		receive_block(mm, worker);
	}

	inc_workers_free();
}

void wait_all(mm_t* mm) {
	MPI_Status *stats = (MPI_Status*) malloc((WORLD_SIZE-1) * sizeof(MPI_Status));

	MPI_Waitall(WORLD_SIZE-1, &workers[1], stats);
	
	for(uint32_t i = 1; i < WORLD_SIZE; ++i) {
		if(mm->blocks[i].status == WAITING) {
			receive_block(mm, i);
		}
	}

	workers_free = WORLD_SIZE-1;
}

// This function takes the destiny process and allocate
// a task to it.
void allocate_task(mm_t* mm, int dest, uint32_t i, uint32_t j) {
	
	dec_workers_free();
	// Update block status
	mm->blocks[dest].status = PROCESSING;
	mm->blocks[dest].i = i;
	mm->blocks[dest].j = j;

	// Assert block height and width.
	mm->blocks[dest].blk_height = (mm->res->lines-i) > mm->blk_height ? mm->blk_height : mm->res->lines-i;
	mm->blocks[dest].blk_width = (mm->res->cols-j) > mm->blk_width ? mm->blk_width : mm->res->cols-j;

	if(dest == 0) {
		#pragma omp task 
		{ master_worker(mm); }
		// Set worker as pending in workers' array.
		MPI_Irecv(&mm->blocks[dest].status, 1, MPI_INT, dest, 0, MPI_COMM_WORLD, &workers[dest]);
	}
	else {
		// Send task to another process.
		send_task(mm, dest);		
	}

}

void send_task(mm_t* mm, int dest) {
	uint32_t i = mm->blocks[dest].i;
	uint32_t j = mm->blocks[dest].j;
	uint32_t blk_height = mm->blocks[dest].blk_height; 
	uint32_t blk_width = mm->blocks[dest].blk_width;
	
	float **a_slc, **b_slc;

	a_slc = &mm->a->mat[i];
	b_slc = &mm->b->mat[j];

	// TODO: A improvement could be made here if 
	// the slicing operation didn't need to create
	// a new fmat object. 
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

	// Send data to worker.

	// Tell another process that he will receive a task.
	int task = 1;
	MPI_Send(&task, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

	// Send block dimensions.
	MPI_Send(&blk_height, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	MPI_Send(&blk_width, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

	// Send matrices.
	MPI_Send_fmat(a, dest);
	MPI_Send_fmat(b, dest);

	// Set worker as pending in workers' array.
	MPI_Irecv(&mm->blocks[dest].status, 1, MPI_INT, dest, 0, MPI_COMM_WORLD, &workers[dest]);

	free_fmat(a);
	free_fmat(b);
}

// Use master process to also multiply the matrix.
void master_worker(mm_t* mm) {
	fmat_t *res;
	block_t *task = &(mm->blocks[MASTER_RANK]);
	float **a = &mm->a->mat[task->i];
	float **b = &mm->b->mat[task->j];
	int done = DONE;

	res = blk_mul(a, b, task->blk_height, task->blk_width, mm->a->cols);

	set_blk(task->i, task->j, res->lines, res->cols, res->mat, mm->res->mat);
	task->status = DONE;

	// Tell master that the block is done.
	MPI_Send(&done, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD);

	free_fmat(res);	
}

void finalize(mm_t* mm) {
	int task = 0;

	wait_all(mm);

	for(int i = 1; i < WORLD_SIZE; ++i) {
		MPI_Send(&task, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	}
	free(workers);
}

// Slave
void slave_loop() {
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

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
	int task;
	MPI_Status stat;

	MPI_Recv(&task, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
	return task;	
}

void recv_task() {
	uint32_t blk_height, blk_width;
	MPI_Status stat;
	fmat_t *a, *b, *block;
	int done = WAITING;

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
	free_fmat(a);
	free_fmat(b);
}
