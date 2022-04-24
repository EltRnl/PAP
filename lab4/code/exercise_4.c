/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"
#include <math.h>
#include <assert.h>

int find_factor(int n){
	int sqr = sqrt(n);
	for(int d=sqr; d>1; d--){
		if(n%d==0)
			return d;
	}
	return 1;
}

/****************************************************/
void lbm_comm_init_ex4(lbm_comm_t * comm, int total_width, int total_height)
{
	//
	// DONE: calculate the splitting parameters for the current task.
	//

	//get infos
	int rank;
	int comm_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	// DONE: calculate the number of tasks along X axis and Y axis.
	if(total_width<total_height){
		comm->nb_x = find_factor(comm_size);
		comm->nb_y = comm_size/comm->nb_x;
	}else{
		comm->nb_y = find_factor(comm_size);
		comm->nb_x = comm_size/comm->nb_y;
	}
	assert(total_width%comm->nb_x == 0 && total_height%comm->nb_y == 0);

	// DONE: calculate the current task position in the splitting
	comm->rank_x = rank%comm->nb_x;
	comm->rank_y = rank/comm->nb_x;

	// DONE : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = (total_width/comm->nb_x) + 2;
	comm->height = (total_height/comm->nb_y) + 2;

	// DONE : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = comm->rank_x*(comm->width-2);
	comm->y = comm->rank_y*(comm->height-2);

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	comm->buffer_recv_down=malloc(sizeof(double)*comm->width*DIRECTIONS); 
	comm->buffer_recv_up=malloc(sizeof(double)*comm->width*DIRECTIONS); 
	comm->buffer_send_down=malloc(sizeof(double)*comm->width*DIRECTIONS); 
	comm->buffer_send_up=malloc(sizeof(double)*comm->width*DIRECTIONS);
	
	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4(lbm_comm_t * comm)
{
	free(comm->buffer_recv_down);
	free(comm->buffer_recv_up);
	free(comm->buffer_send_down);
	free(comm->buffer_send_up);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate. 
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	// Send left and right
	double* send_left = lbm_mesh_get_cell(mesh,1,0);
	double* recv_left = lbm_mesh_get_cell(mesh,0,0);

	double* send_right = lbm_mesh_get_cell(mesh,comm->width-2,0);
	double* recv_right = lbm_mesh_get_cell(mesh,comm->width-1,0);

	MPI_Status status;
	
	if(comm->rank_x!=0) MPI_Recv(recv_left,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x-1,0,MPI_COMM_WORLD,&status);
	if(comm->rank_x!=comm->nb_x-1) MPI_Send(send_right,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x+1,0,MPI_COMM_WORLD);

	if(comm->rank_x!=comm->nb_x-1) MPI_Recv(recv_right,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x+1,0,MPI_COMM_WORLD,&status);
	if(comm->rank_x!=0) MPI_Send(send_left,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x-1,0,MPI_COMM_WORLD);

	// Send up and down
	for(int i=0; i<comm->width*DIRECTIONS; i+=DIRECTIONS){
		comm->buffer_send_up[i] = *lbm_mesh_get_cell(mesh,i/DIRECTIONS,0);
		comm->buffer_send_down[i] =  *lbm_mesh_get_cell(mesh,i/DIRECTIONS,comm->height-1);
	}


}
