#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include <cmath>
using namespace std;
#define	INDEX(i,j,n) i*n+j // Array Indexing (i,j)


void initialize(float* x, long nbound, long mbound){
	// Initializing Internal Random Array
	for(int i = 0; i < nbound; i++) {
		for(int j = 0; j < mbound + 2; j++) {
		  x[INDEX(i,j,mbound)] = 1.0 * rand()/RAND_MAX; 
		}
	}
}

void smooth(float* x, float* y, long n, long m, float a, float b, float c){ 
  for(int i = 1; i < n-1; i++){ 
	 	for(int j = 1; j < m-1; j++){
			y[INDEX(i,j,m)] = a*(x[INDEX((i-1),(j-1),m)] + x[INDEX((i-1),(j+1),m)] + x[INDEX((i+1),(j-1),m)] + x[INDEX((i+1),(j+1),m)]) + b*(x[INDEX((i-1),j,m)] + x[INDEX((i+1),j,m)] + x[INDEX(i,(j-1),m)] + x[INDEX(i,(j+1),m)]) + c*x[INDEX(i,j,m)];
		}
	}	
}

void count(float* array, long n, long m, float t, long long &count) {
  int i,j;
  for(i = 1; i < n + 1 ; ++i){
		for(j = 1; j < m + 1; ++j){
			if(array[INDEX(i,j,m)] < t){ // Checking Threshold Values 
				count++;   
			}
		}
	}
}

double get_wall_time() {
  struct timeval time; 
  gettimeofday(&time,NULL);
  return (double)time.tv_sec + (double)time.tv_usec * .000001; 
}

int main(int argc, char* argv[]){
	
	// Intializing MPI
	int ierr, nranks, rank; 
	ierr = MPI_Init(&argc, &argv);
	MPI_Comm comm = MPI_COMM_WORLD;        // comminicator, comm & comm_cart are used in this program, descriptions on comm_cart are located below
	ierr = MPI_Comm_size(comm, &nranks);   // Totle Number of Processes 
	ierr = MPI_Comm_rank(comm, &rank);     // Unique "process ID" within a communicator (0 ~np -1)

	MPI_Status status;
	MPI_Request request; 

	// =================================== Tile Dimensions ===================================
	// Finding Task Dimensions -Rectangle: (N x M) 
	// * Preference given to a Square-Arrangment
	int n, m;  
	float taskTemp = sqrt(nranks);
	//printf("Sqrt value %g \n",taskTemp);
	if (fmod(taskTemp,1.0) == 0 && nranks != 2) {
		n = (int)taskTemp;
		m = n;
	} else {
		//printf("Numer of Tasks is not a square %g\n",ceil(nranks / 2.0));
		for (float i = 2.0; i < ceil(nranks / 2.0) ; i = i + 1.0) {
			taskTemp = nranks / i;
			if (fmod(taskTemp,1.0) == 0) {
				n = (int)taskTemp;
				m = i;
			}
		}
	}
	if (rank == 0 ){ printf("Task Tile Arrangement:  %d x %d\n", n,m);}
	ierr = MPI_Barrier(comm);

	// Creating Size of Each Tile
	long total_size = 10000;
	long n_size = 15810; 
	long m_size = 15810;

	// Individual Tile Dimensions 
	long nbound = n_size + 2; 
	long mbound = m_size + 2; 

	if (rank == 0) {
		printf("Tile Dimensions:  %d x %d\n", nbound, mbound);
	}
 
	// ================================== Declaring Variables ===================================
	float *x, *y;
	float a = 0.0, b = 0.2, c = 0.2, t = 0.5;    
	double start,time_initialize,time_smooth,time_count_x,time_count_y,time_alloc_x,time_alloc_y;
	long long *count_x, *count_y, *elem_count_x, *elem_count_y;

	// ================================= Global Data Array (MPI) =================================
	// MPI_Datatype init_array; 
	// MPI_Type_contiguous(nbound * mbound, MPI_FLOAT, &init_array); 
	// MPI_Type_commit(&init_array);

	// ==================================== X & Y Allocation ====================================
	start = get_wall_time();
	x = (float *) malloc( nbound * mbound * sizeof(float));
	time_alloc_x = get_wall_time() - start; 

	start = get_wall_time();
	y = (float *) malloc( n_size * m_size * sizeof(float));
	time_alloc_y = get_wall_time() - start; 

	// ================================== Array Intialization ===================================
	// Rank = Sends a Random Array of Numbers 
	 if (nranks != 1){
	 	if (rank == 0) {
	 		printf("Sending x array to tasks\n");
	 		start = get_wall_time();
	 		for (int i = 1; i < nranks; i++) { 
	 			initialize(x,nbound,mbound); 																
	 			ierr = MPI_Isend(&x[0],1000,MPI_FLOAT,i, 1, comm, &request);
	 			printf("sent to %d \n",i);
	 		}
	time_initialize = get_wall_time() - start; 
		}
	} else {
		printf("THe number of tasks is below 2\n");
	}
	ierr = MPI_Barrier(comm);
	ierr = MPI_Irecv(&x[0], 1000, MPI_FLOAT, 0, 1, comm, &request);
	printf("%d received \n",rank);
	
	// ================================= Moving Rows & Columns =================================
	// Create Cartesian Topology 
	// tdim: Tile dimension, pdim: periodicity 
	ierr = MPI_Barrier(comm);
	if (rank == 0 ){ printf("Setting Cartesian Topology\n");}
	int ndim =2;
	int dim[2] = {n,m};
	int pdim[2] = {0,1};
	int coordinates[2];
	int rank2d;
	MPI_Comm comm_cart;
	ierr = MPI_Cart_create(MPI_COMM_WORLD,ndim,dim,pdim,1,&comm_cart);
	if (ierr < 0) {
		perror("unable to create");
	}
	ierr = MPI_Cart_coords(comm_cart,rank,ndim,coordinates);
	ierr = MPI_Cart_rank(comm_cart,coordinates,&rank2d);
	printf("Rank: %d (original: %d) : (%d,%d) \n",rank,rank2d,coordinates[0],coordinates[1]);

	// MPI Type for Column Transfers (Vector)
	float col_recv[nbound], col_send[nbound];
	// MPI_Datatype n_col;
	// MPI_Type_vector(nbound,1,mbound,MPI_FLOAT,&n_col);  	 // Column Vectors for column movement
	// MPI_Type_commit(&n_col);

	// ================================== Cartesian Shift ==================================
	// Setting Cartesian source & receive of ranks (src_{action}, dest_{action}
	int src_up, dest_up, src_down, dest_down, src_left, dest_left, src_right, dest_right;
	printf("Setting Cart Shift: %d \n",rank2d);
	ierr = MPI_Cart_shift(comm_cart,0,-1,&src_up, &dest_up); 		 	// Row Up
	ierr = MPI_Cart_shift(comm_cart,0,1,&src_down, &dest_down);	 	// Row Down 
	ierr = MPI_Cart_shift(comm_cart,1,-1,&src_left, &dest_left);  // Column Left 
	ierr = MPI_Cart_shift(comm_cart,1,1,&src_right, &dest_right); // Column Right
	// direction 0 goes down on y, direction 1 goes right on x 
	// displacements are relative to directions 

	ierr = MPI_Barrier(comm);
	ierr = MPI_Barrier(comm_cart);
	printf("cart shift finished: %d \n", rank2d);

	// ================================= Transferring Data =================================
	// Row Up ===========================================================================
	printf("Row Up: %d\n",rank);
	if (dest_up != -1) {
		ierr = MPI_Isend(&x[0], mbound, MPI_FLOAT, dest_up, 1, comm, &request);
		MPI_Wait(&request, &status);
	}

	if (src_up != -1) {
		ierr = MPI_Irecv(&x[INDEX(nbound-1,0,mbound)], mbound, MPI_FLOAT, src_up, 1, comm, &request);
		MPI_Wait(&request,&status);
	}

	// Row Down ===========================================================================
	printf("Row down: %d\n",rank);
	if (dest_down != -1) {
		ierr = MPI_Isend(&x[INDEX(nbound-1,0,mbound)],mbound, MPI_FLOAT, dest_down, 1, comm, &request);
		MPI_Wait(&request, &status);
	}

	if (src_down != -1) {
		ierr = MPI_Irecv(&col_recv[0], mbound, MPI_FLOAT, src_down, 1, comm, &request);
		MPI_Wait(&request, &status);
	}

	// Column Left ========================================================================
	printf("Sending Column left: %d\n",rank);
	if (dest_left != -1) {
		for (int i = 0; i < nbound; i++) {
			col_send[i] = x[INDEX(i,0,mbound)];
		}
		ierr = MPI_Isend(&col_send[0], nbound, MPI_FLOAT, dest_left, 1, comm, &request);
		ierr = MPI_Wait(&request,&status);
	}

		printf("Receiving col left %d\n",rank);	
	if (src_left != -1) {
		ierr = MPI_Irecv(&col_recv[0], nbound, MPI_FLOAT, src_left, 1, comm, &request);
		MPI_Wait(&request, &status);
		// Inserting received in Position 
		for (int i = 0; i < nbound; i++) {
			x[INDEX(i,mbound-1,mbound)] = col_recv[i];
		}
	}
	
	// Column Right =======================================================================
	printf("Column right: %d\n",rank);
	if (dest_right != -1) {
		for (int i = 0; i < nbound; i++) {
			col_send[i] = x[INDEX(i,mbound-1,mbound)];
		}
		ierr = MPI_Isend(&x[0], nbound, MPI_FLOAT, dest_right, 1, comm, &request);
		ierr = MPI_Wait(&request,&status);
	}

	if (src_right != -1) {
		ierr = MPI_Irecv(&col_recv[0], nbound, MPI_FLOAT, src_right, 1, comm, &request);
		ierr = MPI_Wait(&request, &status);
		for (int i = 0; i < nbound; i++) {
			x[INDEX(i,0,mbound)] = col_recv[i];
		}
	}
	printf("transfer finished %d \n",rank);

	// ====================================== Smoothing ======================================
	// Smoothing Function (Serial)
	printf("Smoothing: %d\n",rank);
	start = get_wall_time();
	smooth(x, y, nbound, mbound, a, b, c);
	time_smooth = get_wall_time() - start; 
	ierr = MPI_Barrier(comm);

	// ====================================== Counting =======================================
	printf("Counting X: %d \n",rank);
	start = get_wall_time();
	count(x,nbound, mbound, t, count_x);
	time_count_x = get_wall_time() - start; 

	printf("Counting Y: %d \n",rank);
	start = get_wall_time();
	count(y,nbound, mbound, t, count_y);
	time_count_y = get_wall_time() - start; 

	// Adding up all local counts -> global counts 
	MPI_Reduce(&elem_count_x,&count_x,nranks,MPI_INT, MPI_SUM, 0, comm);
	MPI_Reduce(&elem_count_y,&count_y,nranks,MPI_INT, MPI_SUM, 0, comm);
	printf("Counting finished: %d",rank);
	ierr = MPI_Barrier(comm);
	

	// ======================================= Output ==========================================
	// Output by thread 0
	if (rank == 0){
		printf("\n");
		printf("Summary\n");
		printf("-------\n");
		printf("Number of task                           :: %d\n", nranks);
		printf("Number of elements in a row/column       :: %d x %d \n", nbound, mbound);
		printf("Number of inner elements in a row/column :: %d x %d \n", n_size, m_size);
		printf("Total number of elements                 :: %d\n", nbound * mbound * nranks);
		printf("Total number of inner elements           :: %d\n", n_size * m_size);
		printf("Memory (GB) used per array               :: %g\n", (nbound)*(mbound)*sizeof(float) / (float)(1024*1024*1024));
		printf("Threshold                                :: %g\n", t);
		printf("Smoothing constants (a, b, c)            :: %g, %g, %g\n", a,b,c);
		printf("Number of elements below threshold (X)   :: %d\n", count_x);
		printf("Fraction of elements below threshold     :: %g\n", count_x / (double)(nbound * mbound * nranks));
		printf("Number of elements below threshold (Y)   :: %d\n", count_y);
		printf("Fraction of elements below threshold     :: %g\n", count_y / (double)nbound * mbound * nranks); 
		
		printf("CPU: Alloc-X :: %1.6f\n", time_alloc_x);
		printf("CPU: Alloc-Y :: %1.6f\n", time_alloc_y);
		printf("CPU: Init-X :: %1.6f\n", time_initialize);
		printf("CPU: Smooth :: %1.6f\n", time_smooth);
		printf("CPU: Count-X :: %1.6f\n", time_count_x);
		printf("CPU: Count-Y :: %1.6f\n", time_count_y);
	}
	
	ierr = MPI_Finalize();
}