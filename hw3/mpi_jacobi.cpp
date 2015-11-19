// There's an old C++ bug in MPI that requires it to be included before stdlib.h
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <my_timer.h>

#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
T *Allocate (const int n)
{
   T *ptr = new T[n];
   return ptr;
}

template <typename T>
T *Deallocate (const T *ptr)
{
   if (ptr)
      delete [] ptr;

   return NULL;
}

template <typename T>
T ** AllocateMesh (const int M, const int N)
{
   T **data = new T* [M];
   data[0] = new T [M*N];
   for (int i = 1; i < M; ++i)
      data[i] = data[0] + i*N;

   return data;
}

template <typename T>
T ** DeallocateMesh (T ** ptr)
{
   if (ptr)
   {
      if (ptr[0])
         delete [] ptr[0];

      delete [] ptr;
   }

   return NULL;
}

// template <typename T>
// T * AllocateArray (const int N) {
//    T *data = new T* [N];

// }

// template <typename T>
// T * DeallocateArray (T * ptr) {

// }

int numProcs = -1;
int myRank = -1;
int **ProcMap = NULL;
int iProc = -1;
int jProc = -1;
int numProcs_i = -1;
int numProcs_j = -1;
int **iStart = NULL;
int **iEnd   = NULL;
int **jStart = NULL;
int **jEnd   = NULL;
const int MPI_TAG_EXCHANGE = 0;
MPI_Comm Comm = MPI_COMM_NULL;

// Send/recv the edges of each grid to the north/south/east/west
// neighbors.
void exchange_boundaries (double **x, int Ni, int Nj)
{
   /* MPI_PROC_NULL - This rank may be used to send or receive from no-one.
    *  http://linux.die.net/man/3/mpi_proc_null
    */
   int leftRank = MPI_PROC_NULL;
   int rightRank = MPI_PROC_NULL;
   int upperRank = MPI_PROC_NULL;
   int downRank = MPI_PROC_NULL;

   int column = 0;
   int row = 0;

   for (int i=0; i<numProcs_i; i++)
      for (int j=0; j<numProcs_j; j++)
         if (ProcMap[i][j] == myRank) {
            column = i;
            row = j;
            break;
         }

   int size = iEnd[column][row] - iStart[column][row] + 1;

   if (column>0) { //send left
      leftRank = myRank - numProcs_j;
      printf("left> myRank is %d and I'm exchanging boudary with rank %d\n", myRank, leftRank);
      // MPI_Sendrecv(&x[1], size, MPI_DOUBLE, leftRank, MPI_TAG_EXCHANGE, &x[Ni-1], size, MPI_DOUBLE, leftRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }

   if (column<numProcs_i-1) { //send right
      rightRank = myRank + numProcs_j;
      printf("right> myRank is %d and I'm exchanging boudary with rank %d\n", myRank, rightRank);
      // MPI_Sendrecv(&x[Ni-2], size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, &x[0], size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }

   if (row>0) { //send up
      upperRank = myRank - 1;
      printf("upper> myRank is %d and I'm exchanging boudary with rank %d\n", myRank, upperRank);
   }

   if (row<numProcs_j-1) { //send down
      downRank = myRank + 1;
      printf("down> myRank is %d and I'm exchanging boudary with rank %d\n", myRank, downRank);
   }



            // if (i>0) { //send left
            //    leftRank = j - numProcs_j + i * numProcs_j;
            //    MPI_Sendrecv(&x[1], size, MPI_DOUBLE, leftRank, MPI_TAG_EXCHANGE, &x[Ni-1], size, MPI_DOUBLE, leftRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // }

            // if (i<numProcs_i-1) { //send right
            //    rightRank = j + (i + 1) * numProcs_j;
            //    MPI_Sendrecv(&x[Ni-2], size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, &x[0], size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // }

            // size = iEnd[i][j] - iStart[i][j] + 1;

            // if (j>0) { //send up
            //    upperRank = myRank-1;
            //    double *arraySend = new double* [size];
            //    double *arrayReceive = new double* [size];
            //    for (int k=1; k<Ni-1; k++) {
            //       arraySend[k-1] = x[k][1];
            //       arrayReceive[k-1] = x[k][Nj-1];
            //    }
            //    MPI_Sendrecv(arraySend, size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, arrayReceive, size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               
            //    if (arraySend) delete [] arraySend;
            //    if (arrayReceive) delete [] arrayReceive;

            // }

            // if (j<numProcs_j-1) { //send down
            //    downRank = myRank+1;
            //    double *arraySend = new double* [size];
            //    double *arrayReceive = new double* [size];
            //    for (int k=1; k<Ni-1; k++) {
            //       arraySend[k-1] = x[k][Nj-1];
            //       arrayReceive[k-1] = x[k][1];
            //    }
            //    MPI_Sendrecv(arraySend, size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, arrayReceive, size, MPI_DOUBLE, rightRank, MPI_TAG_EXCHANGE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               
            //    if (arraySend) delete [] arraySend;
            //    if (arrayReceive) delete [] arrayReceive;
            // }
}

int main (int argc, char* argv[])
{
   int numThreads = 1;
#ifdef _OPENMP
   numThreads = omp_get_max_threads();
#endif

   int mpi_error;
   mpi_error = MPI_Init (&argc, &argv);

   MPI_Comm_dup (MPI_COMM_WORLD, &Comm);

   //int myRank, numProcs;
   mpi_error = MPI_Comm_rank (Comm, &myRank);
   mpi_error = MPI_Comm_size (Comm, &numProcs);

   printf("I am %d of %d %d\n", myRank, numProcs, numThreads);

   int N = 10; // 10 x 10 global mesh.
   if (argc > 1)
      if (isdigit(*argv[1]))
         N = atoi(argv[1]);

   int maxIterations = 100; // Maximum # of iterations.
   if (argc > 2)
      if (isdigit(*argv[2]))
         maxIterations = atoi(argv[2]);

   double maxResidual = 1e-4; // Maximum residual before terminating.
   if (argc > 3)
      if (isdigit(*argv[3]))
         maxResidual = atof(argv[3]);

   // Partition the mesh across the processes.

   // Number of partitions in the x and y directions.
   numProcs_i = 1;
   numProcs_j = numProcs;

   // Try to find a nice partition if even or square.
   int imin = 1, jmin = 1, min = 2*numProcs;
   for (int i = 1; i < numProcs; i *= 2)
   {
      if (numProcs % i == 0)
      {
         int ni = i;
         int nj = numProcs / ni;
         if (ni+nj < min)
         {
            imin = ni;
            jmin = nj;
            min = ni+nj;
         }
      }
   }
   numProcs_i = imin;
   numProcs_j = jmin;

   if (myRank == 0)
      printf("numProcs, i, j = %d, %d, %d\n", numProcs, numProcs_i, numProcs_j);

   // Create a mapping of processes onto the 2d mesh.
   //int **ProcMap = AllocateMesh<int>(numProcs_i, numProcs_j);
   ProcMap = AllocateMesh<int>(numProcs_i, numProcs_j);

   // Where am I in the process grid?
   //int iProc, jProc;

   for (int i = 0; i < numProcs_i; ++i)
      for (int j = 0; j < numProcs_j; ++j)
      {
         int rank = j + i * numProcs_j;
         ProcMap[i][j] = rank;

         if (rank == myRank)
         {
            iProc = i;
            jProc = j;
         }
      }

   // Translate the process coordinates into mesh coordinates.
   iStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   iEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);
   jStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   jEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);

   for (int i = 0; i < numProcs_i; ++i)
   {
      for (int j = 0; j < numProcs_j; ++j)
      {
         int Npts_i = (N-2) / numProcs_i;
         if (i < (N-2) % numProcs_i)
            Npts_i++;

         if (i == 0)
            iStart[i][j] = 1;

         iEnd[i][j] = iStart[i][j] + Npts_i - 1;
         if (i != numProcs_i-1)
            iStart[i+1][j] = iEnd[i][j] + 1;

         int Npts_j = (N-2) / numProcs_j;
         if (j < (N-2) % numProcs_j)
            Npts_j++;

         if (j == 0)
            jStart[i][j] = 1;

         jEnd[i][j] = jStart[i][j] + Npts_j - 1;
         if (j != numProcs_j-1)
            jStart[i][j+1] = jEnd[i][j] + 1;
      }
   }

   int Ni = iEnd[iProc][jProc] - iStart[iProc][jProc] + 3;
   int Nj = jEnd[iProc][jProc] - jStart[iProc][jProc] + 3;

   printf("rank=%d, iProc,jProc=%d,%d, iStart,iEnd,Ni=%d,%d,%d, jStart,jEnd,Nj=%d,%d,%d\n", myRank, iProc, jProc, iStart[iProc][jProc], iEnd[iProc][jProc], Ni, jStart[iProc][jProc], jEnd[iProc][jProc], Nj);

   double **x = AllocateMesh<double>(Ni, Nj);
   double **xtemp = AllocateMesh<double>(Ni, Nj);

   // x[i][j] is initially zero everywhere expect ...
   for (int i = 0; i < Ni; ++i)
      for (int j = 0; j < Nj; ++j)
         x[i][j] = 0;

   // x[i][0] is the lower boundary = 1
   if (jProc == 0)
      for (int i = 0; i < Ni; ++i)
         x[i][0] = 1;

   // x[0][j] is the left boundary = 1
   if (iProc == 0)
      for (int j = 0; j < Nj; ++j)
         x[0][j] = 1;

   // Set xtemp = x so the boundaries are consistent.
   for (int i = 0; i < Ni; ++i)
      for (int j = 0; j < Nj; ++j)
         xtemp[i][j] = x[i][j];

   double mpi_p2p_time = 0;
   double mpi_coll_time = 0;

   myTimer_t total_timer = getTimeStamp();

   // Iterate for some number of steps or until we converge.
   int iteration = 0;
   double residual = 1;

   while (residual > maxResidual and iteration < maxIterations)
   {
      myTimer_t mpi_p2p_timer = getTimeStamp();

      if (numProcs > 1)
      {
         exchange_boundaries (x, Ni, Nj);
      }

      mpi_p2p_time += getElapsedTime( mpi_p2p_timer, getTimeStamp() );

      residual = 0;
      #pragma omp parallel for reduction(+:residual)
      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
         {
            xtemp[i][j] = (x[i+1][j] + x[i-1][j] +
                           x[i][j+1] + x[i][j-1]) / 4.0;
            double delta = xtemp[i][j] - x[i][j];
            residual += delta*delta;
         }

      #pragma omp parallel for
      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
            x[i][j] = xtemp[i][j];

      myTimer_t mpi_coll_timer = getTimeStamp();
      
      double local_residual = residual;

      // Add a method to reduce the residual values on each process to a single value
      // and make sure that all processes have the same value.
      MPI_Allreduce(&local_residual, &residual, 1, MPI_DOUBLE, MPI_SUM, Comm);

      mpi_coll_time += getElapsedTime( mpi_coll_timer, getTimeStamp() );

      residual = sqrt(residual);
      iteration++;
      if (myRank == 0 and iteration % 1 == 0)
         printf("%d %4d: %e\n", myRank, iteration, residual);
   }

   double total_time = getElapsedTime( total_timer, getTimeStamp() );
   printf("rank=%d, timers = %f, %f, %f\n", myRank, total_time, mpi_p2p_time, mpi_coll_time);

   double local_times[3] = {1000*total_time, 1000*mpi_p2p_time, 1000*mpi_coll_time};
   double global_times[3];

   MPI_Reduce (local_times, global_times, 3, MPI_DOUBLE, MPI_MAX, 0, Comm);

   /*
    * The last thing is to find the maximum run time for each process. I have added timers to 
    *  measure the total time, the point-to-point communication time, and the collection 
    *  communication time. Reduce these three values to find the maximum across all of the 
    *  processes.
    */
   double local_total_time = 1000*total_time;
   double global_total_time;

   MPI_Reduce(&local_total_time, &global_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, Comm);

   if (myRank == 0)
      printf("N = %d Iterations = %d  residual = %e time = %f %f %f Max time = %f Procs = %d %d %d %d\n", N, iteration, residual, global_times[0], global_times[1], global_times[2], global_total_time, numProcs, numProcs_i, numProcs_j, numThreads);

   if (N < 100)
   {
      char flname[12];
      sprintf(flname, "jacobi%d.out", myRank);
      FILE *f = fopen(flname,"w");
      for (int i = 0; i < Ni; ++i)
         for (int j = 0; j < Nj; ++j)
            fprintf(f,"%e %e %e\n", (iStart[iProc][jProc]+i-1)/double(N-1), (jStart[iProc][jProc]+j-1)/double(N-1), x[i][j]);
            //fprintf(f,"%d %d %e\n", (iStart[iProc][jProc]+i-1), (jStart[iProc][jProc]+j-1), x[i][j]);
      fclose(f);
   }

   DeallocateMesh(xtemp);
   DeallocateMesh(x);

   MPI_Finalize();

   return 0;
}
