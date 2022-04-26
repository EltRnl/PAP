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
#include <stdbool.h>
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
	int comm_size;
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

	// DONE: calculate the number of tasks along X axis and Y axis.
	int dims[2] = {0, 0};
    MPI_Dims_create(comm_size, 2, dims);
    int periods[2] = {false, false};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &comm->communicator);
 
	comm->nb_x = dims[0];
	comm->nb_y = dims[1];

	assert(total_width%comm->nb_x == 0 && total_height%comm->nb_y == 0);

	int rank;
	MPI_Comm_rank(comm->communicator, &rank);
	// DONE: calculate the current task position in the splitting
    int my_coords[2];
    MPI_Cart_coords(comm->communicator, rank, 2, my_coords);

	comm->rank_x = my_coords[0];
	comm->rank_y = my_coords[1];

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
	#ifndef NDEBUG
	lbm_comm_print(comm);
	#endif
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
void copy_cell(double * source, double * target){
	for(int i=0; i<DIRECTIONS; i++){
		target[i] = source[i];
	}
}

/****************************************************/
void lbm_comm_ghost_exchange_ex4(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// DONE: Implement the 2D communication with :
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

	MPI_Status status;
	
	if(n_l!=-1) MPI_Recv(recv_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,comm->communicator,&status);
	if(n_r!=-1) MPI_Send(send_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,comm->communicator);

	if(n_r!=-1) MPI_Recv(recv_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,comm->communicator,&status);
	if(n_l!=-1) MPI_Send(send_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,comm->communicator);


	/************* Sending up and down *************/
	
	// Filling buffers with our data 

	for(int i=0; i<comm->width; i++){
		if(n_u!=-1) copy_cell(lbm_mesh_get_cell(mesh,i,1),comm->buffer_send_up+(i*DIRECTIONS));
		if(n_d!=-1) copy_cell(lbm_mesh_get_cell(mesh,i,comm->height-2),comm->buffer_send_down+(i*DIRECTIONS));
	}

	MPI_Recv(comm->buffer_recv_up,DIRECTIONS*comm->width,MPI_DOUBLE,n_u,0,comm->communicator,&status);
	MPI_Send(comm->buffer_send_down,DIRECTIONS*comm->width,MPI_DOUBLE,n_d,0,comm->communicator);


	MPI_Recv(comm->buffer_recv_down,DIRECTIONS*comm->width,MPI_DOUBLE,n_d,0,comm->communicator,&status);
	MPI_Send(comm->buffer_send_up,DIRECTIONS*comm->width,MPI_DOUBLE,n_u,0,comm->communicator);

	for(int i=0; i<comm->width; i++){
		if(n_u!=-1) copy_cell(comm->buffer_recv_up+(i*DIRECTIONS),lbm_mesh_get_cell(mesh,i,0));
		if(n_d!=-1) copy_cell(comm->buffer_recv_down+(i*DIRECTIONS),lbm_mesh_get_cell(mesh,i,comm->height-1));
	}

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
