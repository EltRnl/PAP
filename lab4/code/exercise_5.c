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

/****************************************************/
void lbm_comm_init_ex5(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5 execpt for type creation
	lbm_comm_init_ex4(comm, total_width, total_height);

	//TODO: create MPI type for non contiguous side in comm->type.

	MPI_Type_vector(comm->width, 1, comm->height, MPI_DOUBLE, &comm->type);
	MPI_Type_commit(&comm->type);
}

/****************************************************/
void lbm_comm_release_ex5(lbm_comm_t * comm)
{
	//we use the same implementation than ex5 except for type destroy
	lbm_comm_release_ex4(comm);

	//TODO: release MPI type created in init.
	MPI_Type_free(&comm->type);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex5(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
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

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	/************* Calculating neighboors rank *************/

	int n_l = rank_from_xy(comm->rank_x-1,comm->rank_y,comm->nb_x,comm->nb_y);
	int n_r = rank_from_xy(comm->rank_x+1,comm->rank_y,comm->nb_x,comm->nb_y);
	int n_u = rank_from_xy(comm->rank_x,comm->rank_y-1,comm->nb_x,comm->nb_y);
	int n_d = rank_from_xy(comm->rank_x,comm->rank_y+1,comm->nb_x,comm->nb_y);
	//printf("Process (%d,%d)\n",n_u,n_d);

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
	
	
	if(n_l!=-1) MPI_Recv(recv_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,MPI_COMM_WORLD,&status);
	if(n_r!=-1) MPI_Send(send_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,MPI_COMM_WORLD);

	if(n_r!=-1) MPI_Recv(recv_right,DIRECTIONS*comm->height,MPI_DOUBLE,n_r,0,MPI_COMM_WORLD,&status);
	if(n_l!=-1) MPI_Send(send_left,DIRECTIONS*comm->height,MPI_DOUBLE,n_l,0,MPI_COMM_WORLD);

	/************* Sending up and down *************/

	//printf("Process (%d,%d) before u/d\n", n_u,n_d);
	if(n_u!=-1) MPI_Recv(recv_up,comm->width,comm->type,n_u,0,MPI_COMM_WORLD,&status);
	//printf("Process (%d,%d) after recv u\n", n_u,n_d);
	if(n_d!=-1) MPI_Send(send_down,comm->width,comm->type,n_d,0,MPI_COMM_WORLD);
	//printf("Process (%d,%d) after send d", n_u,n_d);
	

	if(n_d!=-1) MPI_Recv(recv_down,comm->width,comm->type,n_d,0,MPI_COMM_WORLD,&status);
	//printf("Process (%d,%d) after recv d", n_u,n_d);
	if(n_u!=-1) MPI_Send(send_up,comm->width,comm->type,n_u,0,MPI_COMM_WORLD);
	//printf("Process (%d,%d) after send u", n_u,n_d);
	
}
