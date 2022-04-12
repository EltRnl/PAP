/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// DONE: Implement the 1D communication with non-blocking MPI functions.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	double* send_left = lbm_mesh_get_cell(mesh,1,0);
	double* recv_left = lbm_mesh_get_cell(mesh,0,0);

	double* send_right = lbm_mesh_get_cell(mesh,comm->width-2,0);
	double* recv_right = lbm_mesh_get_cell(mesh,comm->width-1,0);

	MPI_Request request[4];
	int count = 0;
	
	if(comm->rank_x!=0) MPI_Irecv(recv_left,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x-1,0,MPI_COMM_WORLD,request+(count++));
	if(comm->rank_x!=comm->nb_x-1) MPI_Isend(send_right,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x+1,0,MPI_COMM_WORLD,request+(count++));

	if(comm->rank_x!=comm->nb_x-1) MPI_Irecv(recv_right,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x+1,0,MPI_COMM_WORLD,request+(count++));
	if(comm->rank_x!=0) MPI_Isend(send_left,DIRECTIONS*comm->height,MPI_DOUBLE,comm->rank_x-1,0,MPI_COMM_WORLD,request+(count++));

	MPI_Waitall(count,request,MPI_STATUS_IGNORE);

}
