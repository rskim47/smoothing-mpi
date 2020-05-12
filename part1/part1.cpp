#include <iostream>
#include <sys/time.h>
#include <mpi.h>
#include <cmath>
using namespace std;
#define	INDEX(i,j,n) i*n+j // Array Indexing (i,j)


void initialize(float* x, long n, long m){
	// Initializing Internal Random Array
	for(int i = 0; i < n + 2; i++) {
		for(int j = 0; j < m + 2; j++) {
		  x[INDEX(i,j,m)] = 1.0 * rand()/RAND_MAX; 
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

void count(float* array, long n, long m, float t, double &count) {
  int i,j;
  for(i = 1; i < n + 1 ; ++i){
		for(j = 1; j < m + 1; ++j){
			if(array[INDEX(i,j,m)] < t){ // Checking Threshold Values 
				count++;   
			}
		}
	}
}

void parllelRegionCheck(int rank){
	// Outputs Thread # 
	cout << "Thread " << rank << "| Staus: Active" << endl;	
	MPI_Barrier(MPI_COMM_WORLD);
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
	MPI_Comm comm_old = MPI_COMM_WORLD;        // comminicator, comm_old & comm_cart are used in this program, descriptions on comm_cart are located below
	ierr = MPI_Comm_size(comm_old, &nranks);   // Totle Number of Processes 
	ierr = MPI_Comm_rank(comm_old, &rank);     // Unique "process ID" within a communicator (0 ~np -1)
	
	MPI_Comm comm_cart;
	MPI_Status status;
	MPI_Request request; 

	
	// =================================== Tile Dimensions ===================================
	// Finding Task Dimensions -Rectangle: (N x M) 
	// * Preference given to a Square-Arrangment
	int n, m;  
	float taskTemp = sqrt(nranks);
	if (rank == 0) {
		printf("Sqrt value %g \n",taskTemp);
		if (fmod(taskTemp,1.0) == 0 && nranks != 2) {
			n = (int)taskTemp;
			m = n;
		} else {
			printf("Numer of Tasks is not a square %g\n",ceil(nranks / 2.0));
			for (float i = 2.0; i < ceil(nranks / 2.0) ; i = i + 1.0) {
				taskTemp = nranks / i;
				if (fmod(taskTemp,1.0) == 0) {
					n = (int)taskTemp;
					m = i;
				}
			}
		}
		printf("Task Tile Arrangement:  %d x %d\n", n,m);
	}
	
	ierr = MPI_Barrier(comm_old);

	// Creating Size of Each Tile
	long total_size = 98306;
	long n_size = total_size / n; 
	long m_size = total_size / m;
	long nbound = n_size + 2; 
	long mbound = m_size + 2; 
	if (rank == 0) {
		printf("Tile Dimensions:  %d x %d\n", nbound, mbound);
	}
 
	// ================================ Intializing Variables =================================
	float *x, *y;
	float a = 0.0, b = 0.2, c = 0.2, t = 0.5;    
	double start,time_initialize,time_smooth,time_count_x,time_count_y,time_alloc_x,time_alloc_y;
	double count_x, count_y, elem_count_x, elem_count_y;

	// ================================= Global Data Array (MPI) =================================
	// MPI_Datatype init_array; 
	// MPI_Type_contiguous(nbound * mbound, MPI_FLOAT, &init_array); 
	// MPI_Type_commit(&init_array);

	// ==================================== X & Y Allocation ====================================
	// Allocation made for all tasks.
	start = get_wall_time();
	x = (float *) malloc( nbound * mbound * sizeof(float));
	time_alloc_x = get_wall_time() - start; 

	start = get_wall_time();
	y = (float *) malloc( nbound * mbound * sizeof(float));
	time_alloc_y = get_wall_time() - start; 

	// ================================= Array Intialization =================================
	// Scheme 1 : Every MPI tasks initializes its own array
	if (nrank != 1){
		if (rank == 0) {
			start = get_wall_time();
			for (int i = 1; i < nranks; i++) { 
				initialize(x,nbound,mbound); 																// Creating Array
				ierr = MPI_Isend(&x[0], nbound*mbound ,MPI_FLOAT,i, 1, comm_old, &request);
				ierr = MPI_recv(&x[0], nbound*mbound, MPI_FLOAT, 0, 1, comm_old, &request);
				ierr = MPI_Wait(&request,&statusus);
			}
			time_initialize = get_wall_time() - start; 
		}
	} else {
		printf("THe number of tasks is below 2");
	}


	// ================================= Moving Rows & Columns =================================
	// Create Cartesian Topology 
	// tdim: Tile dimension, pdim: periodicity 
	int dim[2] = {n,m};
	int pdim[2] = {0,0};
	ierr = MPI_Cart_create(comm_old,2,dim,pdim,0,&comm_cart);
	int MPI_Cart_get(comm_cart)

	// Rank Locations 
	 // Coordinates of each rank in Cartesian Topology
	int coordinates[2];
	ierr = MPI_Cart_coords(comm_cart,rank,2,coordinates);
	ierr = MPI_Cart_rank(comm_cart,coordinates,&rank);
	//printf("Rank %d : %d x %d",rank,coordinates[0],coordinates[1]);

	// MPI Type for Column Transfers (Vector)
	float col_recv[nbound]; 
	MPI_Datatype n_col;
	MPI_Type_vector(nbound,1,mbound,MPI_FLOAT,&n_col);  	 // Column Vectors for column movement
	MPI_Type_commit(&n_col);

	// ================================== Cartesian Shift ==================================
	// Setting Cartesian source & receive of ranks (src_{action}, dest_{action}
	int src_up, dest_up, src_down, dest_down, src_left, dest_left, src_right, dest_right;
	
	ierr = MPI_Cart_shift(comm_cart,0,-1,&src_up, &dest_up); 		 	// Row Up
	ierr = MPI_Cart_shift(comm_cart,0,1,&src_down, &dest_down);	 	// Row Down 
	ierr = MPI_Cart_shift(comm_cart,1,-1,&src_left, &dest_left);  // Column Left 
	ierr = MPI_Cart_shift(comm_cart,1,1,&src_right, &dest_right); // Column Right
	// direction 0 goes down on y, direction 1 goes right on x 
	// displacements are relative to directions 

	// ================================= Transferring Data =================================
	// Row Up ===========================================================================
	if (dest_up != -1) {
		ierr = MPI_Isend(&x[0], mbound, MPI_FLOAT, dest_up, 1, comm_old, &request);
		MPI_Wait(&request, &status);
	}

	if (src_up != -1) {
		ierr = MPI_Irecv(&x[INDEX(nbound,0,mbound)], mbound, MPI_FLOAT, src_up, 1, comm_old, &request);
		MPI_Wait(&request,&status);
	}

	// Row Down ===========================================================================
	if (dest_down != -1) {
		ierr = MPI_Isend(&x[INDEX(nbound,0,mbound)],mbound, MPI_FLOAT, dest_down, 1, comm_old, &request);
		MPI_Wait(&request, &status);
	}

	if (src_down != -1) {
		ierr = MPI_Irecv(&x[0], mbound, MPI_FLOAT, src_down, 1, comm_old, &request);
		MPI_Wait(&request, &status);
	}

	// Column Left ========================================================================
	if (dest_left != -1) {
		ierr = MPI_Isend(&x[nbound+1], 1, n_col, dest_left, 1, comm_old, &request);
		ierr = MPI_Wait(&request,&status);
	}

	if (src_left != -1) {
		ierr = MPI_Irecv(&col_recv[0], nbound, MPI_FLOAT, src_left, 1, comm_old, &request);
		MPI_Wait(&request, &status);
		// Inserting received in Position 
		for (int i = (2*mbound)-1, j = 0; i < mbound*(mbound-1); i += mbound, j++) {
			x[i] = col_recv[j];
		}
	}

	// Column Right =======================================================================
	if (dest_right != -1) {
		ierr = MPI_Isend(&x[(2*mbound-2)], 1, n_col, dest_right, 1, comm_old, &request);
	}

	if (src_right != -1) {
		ierr = MPI_Irecv(&col_recv[0], nbound, MPI_FLOAT, src_right, 1, comm_old, &request);
		MPI_Wait(&request, &status);
		for (int i = mbound, j = 0; i < mbound*(nbound-1); i += mbound, j++) {
			x[i] = col_recv[j];
		}
	}
	ierr = MPI_Barrier(comm_old);

	// ====================================== Smoothing ======================================
	// Smoothing Function (Serial)
	start = get_wall_time();
	smooth(x, y, nbound, mbound, a, b, c);
	time_smooth = get_wall_time() - start; 
	ierr = MPI_Barrier(comm_old);

	// ====================================== Counting =======================================
	start = get_wall_time();
	count(x,nbound, mbound, t, count_x);
	time_count_x = get_wall_time() - start; 
	ierr = MPI_Barrier(comm_old);

	start = get_wall_time();
	count(x,nbound, mbound, t, count_y);
	time_count_y = get_wall_time() - start; 

	// Adding up all local counts -> global counts 
	MPI_Reduce(&elem_count_x,&count_x,nranks ,MPI_INT, MPI_SUM, 0, comm_old);
	MPI_Reduce(&elem_count_y,&count_y,nranks ,MPI_INT, MPI_SUM, 0, comm_old);
	ierr = MPI_Barrier(comm_old);

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