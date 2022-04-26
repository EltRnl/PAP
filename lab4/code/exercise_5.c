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
//      8 neighbors using MPI types for non contiguous
//      side.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
// NEW:
//     - >>> MPI type for non contiguous cells <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"
#include <stdbool.h>

/****************************************************/
void lbm_comm_init_ex5(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5 execpt for type creation
	lbm_comm_init_ex4(comm, total_width, total_height);

	//DONE: create MPI type for non contiguous side in comm->type.
	
	MPI_Type_vector(comm->width, DIRECTIONS, DIRECTIONS*comm->height, MPI_DOUBLE, &comm->type);
	MPI_Type_commit(&comm->type);
}

/****************************************************/
void lbm_comm_release_ex5(lbm_comm_t * comm)
{
	//we use the same implementation than ex5 except for type destroy
	lbm_comm_release_ex4(comm);

	//DONE: release MPI type created in init.
	MPI_Type_free(&comm->type);
}

/****************************************************/
int rank_from_coords(int x, int y, int w, int h)
{
	if(x<0 || x>=w || y<0 || y>=h)
		return -1;
	return y * w + x;
}

/****************************************************/
void lbm_comm_ghost_exchange_ex5(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// DONE: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - use MPI type for non contiguous side 
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

	//DONE:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	bool do_corners = false;

	/************* Calculating neighboors rank *************/

	// Direct neighboors
	int n_l, n_r, n_u, n_d; 
	MPI_Cart_shift(comm->communicator, 0, 1, &n_l, &n_r);
	MPI_Cart_shift(comm->communicator, 1, -1, &n_d, &n_u);

	// Diagonal neighboors
	int n_ul = MPI_PROC_NULL; int n_ur = MPI_PROC_NULL; int n_dl = MPI_PROC_NULL; int n_dr = MPI_PROC_NULL;
	if(do_corners){
		int ncoords[2] = {comm->rank_x-1, comm->rank_y-1};
		if(!(ncoords[0]<0 || ncoords[1]<0 || ncoords[0]>=comm->nb_x || ncoords[1]>=comm->nb_y)) MPI_Cart_rank(comm->communicator,ncoords,&n_ul);
		ncoords[0] += 2;
		if(!(ncoords[0]<0 || ncoords[1]<0 || ncoords[0]>=comm->nb_x || ncoords[1]>=comm->nb_y)) MPI_Cart_rank(comm->communicator,ncoords,&n_ur);
		ncoords[1] += 2; 
		if(!(ncoords[0]<0 || ncoords[1]<0 || ncoords[0]>=comm->nb_x || ncoords[1]>=comm->nb_y)) MPI_Cart_rank(comm->communicator,ncoords,&n_dr);
		ncoords[0] -= 2; 
		if(!(ncoords[0]<0 || ncoords[1]<0 || ncoords[0]>=comm->nb_x || ncoords[1]>=comm->nb_y)) MPI_Cart_rank(comm->communicator,ncoords,&n_dl);
	}

	/************* Sending left and right *************/

	// Getting data pointers
	double* send_left = lbm_mesh_get_cell(mesh,1,0);
	double* recv_left = lbm_mesh_get_cell(mesh,0,0);

	double* send_right = lbm_mesh_get_cell(mesh,comm->width-2,0);
	double* recv_right = lbm_mesh_get_cell(mesh,comm->width-1,0);

	double* send_up = lbm_mesh_get_cell(mesh,0,1);
	double* recv_up = lbm_mesh_get_cell(mesh,0,0);

	double* send_down = lbm_mesh_get_cell(mesh,0,comm->height-2);
	double* recv_down = lbm_mesh_get_cell(mesh,0,comm->height-1);

	MPI_Status status;
	
	
	MPI_Recv(recv_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,comm->communicator,&status);
	MPI_Send(send_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,comm->communicator);

	MPI_Recv(recv_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,comm->communicator,&status);
	MPI_Send(send_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,comm->communicator);

	/************* Sending up and down *************/

	MPI_Recv(recv_up,1,comm->type,n_u,0,comm->communicator,&status);
	MPI_Send(send_down,1,comm->type,n_d,0,comm->communicator);
	

	MPI_Recv(recv_down,1,comm->type,n_d,0,comm->communicator,&status);
	MPI_Send(send_up,1,comm->type,n_u,0,comm->communicator);

	/************* Sending diagonaly *************/

	// Getting data pointers
	if(do_corners){
		double* send_ul = lbm_mesh_get_cell(mesh,1,1);
		double* recv_ul = lbm_mesh_get_cell(mesh,0,0);

		double* send_ur = lbm_mesh_get_cell(mesh,comm->width-2,1);
		double* recv_ur = lbm_mesh_get_cell(mesh,comm->width-1,0);

		double* send_dl = lbm_mesh_get_cell(mesh,1,comm->height-2);
		double* recv_dl = lbm_mesh_get_cell(mesh,0,comm->height-1);

		double* send_dr = lbm_mesh_get_cell(mesh,comm->width-2,comm->height-2);
		double* recv_dr = lbm_mesh_get_cell(mesh,comm->width-1,comm->height-1);	

		MPI_Recv(recv_ur,DIRECTIONS,MPI_DOUBLE,n_ur,0,comm->communicator,&status);
		MPI_Send(send_dl,DIRECTIONS,MPI_DOUBLE,n_dl,0,comm->communicator);

		MPI_Recv(recv_dl,DIRECTIONS,MPI_DOUBLE,n_dl,0,comm->communicator,&status);
		MPI_Send(send_ur,DIRECTIONS,MPI_DOUBLE,n_ur,0,comm->communicator);
		
		MPI_Recv(recv_ul,DIRECTIONS,MPI_DOUBLE,n_ul,0,comm->communicator,&status);
		MPI_Send(send_dr,DIRECTIONS,MPI_DOUBLE,n_dr,0,comm->communicator);

		MPI_Recv(recv_dr,DIRECTIONS,MPI_DOUBLE,n_dr,0,comm->communicator,&status);
		MPI_Send(send_ul,DIRECTIONS,MPI_DOUBLE,n_ul,0,comm->communicator);
	}
}
