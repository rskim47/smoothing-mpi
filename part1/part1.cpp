#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include <cmath>
using namespace std;
#define	INDEX(i,j,n) i*n+j // Array Indexing (i,j)


void initialize(float* x, long nbound, long mbound){
	// Initializing Internal Random Array
	for(int i = 0; i < nbound; i++) {
		for(int j = 0; j < mbound; j++) {
		  x[INDEX(i,j,mbound)] = 1.0 * rand()/RAND_MAX; 
		}
	}
}

void smooth(float* x, float* y, long nbound, long mbound, float a, float b, float c){ 
  for(int i = 1; i < nbound - 1; i++){ 
	 	for(int j = 1; j < mbound - 1; j++){
				y[INDEX(i,j,mbound)] = a*(x[INDEX((i-1),(j-1),mbound)] + x[INDEX((i-1),(j+1),mbound)] + x[INDEX((i+1),(j-1),mbound)] + x[INDEX((i+1),(j+1),mbound)]) + b*(x[INDEX((i-1),j,mbound)] + x[INDEX((i+1),j,mbound)] + x[INDEX(i,(j-1),mbound)] + x[INDEX(i,(j+1),mbound)]) + c*x[INDEX(i,j,mbound)];
		}
	}	
}

void count(float* array, long nbound, long mbound, float t, uint64_t * count = 0) {
	uint64_t below = 0; 
  for(int i = 1; i < nbound - 1; i++){
		for(int j = 1; j < mbound - 1; j++){
			if(array[INDEX(i,j,mbound)] < t){ // Checking Threshold Values 
				below++;
			}
		}
	}
	*count = below;
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
	if (nranks > 2) {
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
	} else if (nranks == 2){
		n = 1; m = 2;
	} else if (nranks == 1) {
		n = 1; m = 1;
	}
	if (rank == 0 ){ printf("Task Tile Arrangement:  %d x %d\n", n,m);}
	ierr = MPI_Barrier(comm);

	// Creating Size of Each Tile
	long total_size = 16000;

	long n_size = total_size / n; 
	long m_size = total_size / m;

	// Individual Tile Dimensions 
	long nbound = n_size + 2; 
	long mbound = m_size + 2; 

	if (rank == 0) {
		printf("Tile Dimensions:  %d x %d\n", nbound, mbound);
	}
 
	// ================================== Declaring Variables ===================================
	float *x, *y;
	float a = 0.0, b = 0.2, c = 0.2, t = 0.5;    
	double start,start2,time_initialize,time_intialize2,time_smooth,time_count_x,time_count_y,time_alloc_x,time_alloc_y;
	uint64_t count_x, count_y, elem_count_x, elem_count_y;

	// ==================================== X & Y Allocation ====================================
	start = get_wall_time();
	x = (float *) malloc( nbound * mbound * sizeof(float)); 
	ierr = MPI_Barrier(comm);
	time_alloc_x = get_wall_time() - start;

	start = get_wall_time();
	y = (float *) malloc( nbound * mbound * sizeof(float));
	ierr = MPI_Barrier(comm);
	time_alloc_y = get_wall_time() - start; 

	// ================================== Array Intialization ===================================
	// Rank = Sends a Random Array of Numbers 
	start = get_wall_time();
	if (nranks > 1){
		if (rank == 0) {
			printf("Sending x array to tasks\n");
			start2 = get_wall_time();
			for (int i = 1; i < nranks; i++) { 
				start2 = get_wall_time();
				initialize(x,nbound,mbound);
				time_intialize2 = get_wall_time() - start2;																
				ierr = MPI_Isend(&x[0],1000,MPI_FLOAT,i, 1, comm, &request);
				printf("sent to %d \n",i);
			}
		}
	} 
	ierr = MPI_Irecv(&x[0], 1000, MPI_FLOAT, 0, 1, comm, &request);
	ierr = MPI_Barrier(comm);
	time_initialize = get_wall_time() - start; 
	//printf("%d received \n",rank);
	
	// ================================= Moving Rows & Columns =================================
	// Create Cartesian Topology 
	// tdim: Tile dimension, pdim: periodicity 
	ierr = MPI_Barrier(comm);
	if (nranks > 1) {
		//if (rank == 0 ){ printf("Setting Cartesian Topology\n");}
		int ndim =2;
		int dim[2] = {n,m};
		int pdim[2] = {0,1};
		int coordinates[2];
		int rank2d;
		MPI_Comm comm_cart;
		ierr = MPI_Cart_create(MPI_COMM_WORLD,ndim,dim,pdim,1,&comm_cart);
		ierr = MPI_Cart_coords(comm_cart,rank,ndim,coordinates);
		ierr = MPI_Cart_rank(comm_cart,coordinates,&rank2d);
		//printf("Rank: %d (original: %d) : (%d,%d) \n",rank,rank2d,coordinates[0],coordinates[1]);

		// MPI Type for Column Transfers (Vector)
		float col_recv[nbound], col_send[nbound];

		// ================================== Cartesian Shift ==================================
		// Setting Cartesian source & receive of ranks (src_{action}, dest_{action}
		int src_up, dest_up, src_down, dest_down, src_left, dest_left, src_right, dest_right;
		//printf("Setting Cart Shift: %d \n",rank2d);
		ierr = MPI_Cart_shift(comm_cart,0,-1,&src_up, &dest_up); 		 	// Row Up
		ierr = MPI_Cart_shift(comm_cart,0,1,&src_down, &dest_down);	 	// Row Down 
		ierr = MPI_Cart_shift(comm_cart,1,-1,&src_left, &dest_left);  // Column Left 
		ierr = MPI_Cart_shift(comm_cart,1,1,&src_right, &dest_right); // Column Right
		// direction 0 goes down on y, direction 1 goes right on x 
		// displacements are relative to directions 

		ierr = MPI_Barrier(comm);
		ierr = MPI_Barrier(comm_cart);
		//printf("cart shift finished: %d \n", rank2d);

		// ================================= Transferring Data =================================
		// Row Up ===========================================================================
		//printf("Row Up: %d\n",rank);
		if (dest_up != -1) {
			ierr = MPI_Isend(&x[0], mbound, MPI_FLOAT, dest_up, 1, comm, &request);
			MPI_Wait(&request, &status);
		}

		if (src_up != -1) {
			ierr = MPI_Irecv(&x[INDEX(nbound-1,0,mbound)], mbound, MPI_FLOAT, src_up, 1, comm, &request);
			MPI_Wait(&request,&status);
		}

		// Row Down ===========================================================================
		//printf("Row down: %d\n",rank);
		if (dest_down != -1) {
			ierr = MPI_Isend(&x[INDEX(nbound-1,0,mbound)],mbound, MPI_FLOAT, dest_down, 1, comm, &request);
			MPI_Wait(&request, &status);
		}

		if (src_down != -1) {
			ierr = MPI_Irecv(&col_recv[0], mbound, MPI_FLOAT, src_down, 1, comm, &request);
			MPI_Wait(&request, &status);
		}

		// Column Left ========================================================================
		//printf("Sending Column left: %d\n",rank);
		if (dest_left != -1) {
			for (int i = 0; i < nbound; i++) {
				col_send[i] = x[INDEX(i,0,mbound)];
			}
			ierr = MPI_Isend(&col_send[0], nbound, MPI_FLOAT, dest_left, 1, comm, &request);
			ierr = MPI_Wait(&request,&status);
		}

		//printf("Receiving col left %d\n",rank);	
		if (src_left != -1) {
			ierr = MPI_Irecv(&col_recv[0], nbound, MPI_FLOAT, src_left, 1, comm, &request);
			MPI_Wait(&request, &status);
			// Inserting received in Position 
			for (int i = 0; i < nbound; i++) {
				x[INDEX(i,mbound-1,mbound)] = col_recv[i];
			}
		}
		
		// Column Right =======================================================================
		//printf("Column right: %d\n",rank);
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
		//printf("transfer finished %d \n",rank);
	}

	ierr = MPI_Barrier(comm);
	// ====================================== Smoothing ======================================
	// Smoothing Function (Serial)
	printf("Smoothing: %d\n",rank);
	start = get_wall_time();
	smooth(x, y,nbound,mbound, a, b, c);
	printf("Smoothing finished: %d\n",rank);
	ierr = MPI_Barrier(comm);
	time_smooth = get_wall_time() - start; 
	

	// ====================================== Counting =======================================
	printf("Counting X: %d \n",rank);
	start = get_wall_time();
	count(x,nbound, mbound, t, &count_x);
	ierr = MPI_Barrier(comm);
	time_count_x = get_wall_time() - start; 

	printf("Counting Y: %d \n",rank);
	start = get_wall_time();
	count(y,nbound, mbound, t, &count_y);
	ierr = MPI_Barrier(comm);
	time_count_y = get_wall_time() - start; 

	// Adding up all local counts -> global counts 
	MPI_Reduce(&elem_count_x,&count_x,nranks,MPI_INT, MPI_SUM, 0, comm);
	MPI_Reduce(&elem_count_y,&count_y,nranks,MPI_INT, MPI_SUM, 0, comm);
	printf("Counting finished: %d \n",rank);
	ierr = MPI_Barrier(comm);
	ierr = MPI_Finalize();

	// ======================================= Output ==========================================
	if (rank == 0){
		printf("\n");
		printf("Summary\n");
		printf("-------\n");
		printf("Number of task                           :: %d\n", nranks);
		printf("Number of elements in a row/column       :: %d x %ld \n", nbound, mbound);
		printf("Number of inner elements in a row/column :: %d x %d \n", n_size, m_size);
		printf("Total number of elements                 :: %d\n", nbound * mbound * nranks);
		printf("Total number of inner elements           :: %d\n", n_size * m_size * nranks);
		printf("Memory (GB) used per array               :: %f\n", (float)(nbound)*(mbound)*sizeof(float) / (1024*1024*1024));
		printf("Threshold                                :: %g\n", t);
		printf("Smoothing constants (a, b, c)            :: %g, %g, %g\n", a,b,c);
		printf("Number of elements below threshold (X)   :: %d\n", count_x);
		printf("Fraction of elements below threshold     :: %g\n", count_x / (double)(nbound * mbound * nranks));
		printf("Number of elements below threshold (Y)   :: %d\n", count_y);
		printf("Fraction of elements below threshold     :: %g\n", count_y / (double)(nbound * mbound * nranks)); 
		
		printf("CPU: Alloc-X :: %1.6f\n", time_alloc_x);
		printf("CPU: Alloc-Y :: %1.6f\n", time_alloc_y);
		printf("CPI: new     :: %1.6f\n", time_intialize2);
		printf("CPU: Init-X  :: %1.6f\n", time_initialize);
		printf("CPU: Smooth  :: %1.6f\n", time_smooth);
		printf("CPU: Count-X :: %1.6f\n", time_count_x);
		printf("CPU: Count-Y :: %1.6f\n", time_count_y);
}
		
		
	
	

}