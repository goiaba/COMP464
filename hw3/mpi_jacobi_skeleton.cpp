// There's an old C++ bug in MPI that requires it to be included before stdlib.h
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <my_timer.h>

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
MPI_Comm Comm = MPI_COMM_NULL;

void exchange_boundaries (double **x, int Ni, int Nj)
{
   return;
}

int main (int argc, char* argv[])
{
   int mpi_error;
   mpi_error = MPI_Init (&argc, &argv);

   MPI_Comm_dup (MPI_COMM_WORLD, &Comm);

   //int myRank, numProcs;
   mpi_error = MPI_Comm_rank (Comm, &myRank);
   mpi_error = MPI_Comm_size (Comm, &numProcs);

   printf("I am %d of %d\n", myRank, numProcs);

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

   // Partition the mesh across the processes in the x and y directions.

   // Number of partitions in the x and y directions.
   numProcs_i = 1;
   numProcs_j = numProcs;

   // Try to find a nice partition if even or square.
   for (int i = 1; i < numProcs; i *= 2)
   {
      if (numProcs % i == 0)
      {
         numProcs_i = i;
         numProcs_j = numProcs / numProcs_i;
      }
   }
   if (myRank == 0)
      printf("numProcs i,j = %d, %d, %d\n", numProcs, numProcs_i, numProcs_j);

   // Create a mapping of processes onto a 2d mesh.
   ProcMap = AllocateMesh<int>(numProcs_i, numProcs_j);

   // Where am I in the process grid?

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
   // Each partition has a start and end index. These indices
   // are the points computed and do not include the halo or
   // boundary points.
   iStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   iEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);
   jStart = AllocateMesh<int>(numProcs_i, numProcs_j);
   jEnd   = AllocateMesh<int>(numProcs_i, numProcs_j);

   {
      for (int i = 0; i < numProcs_i; ++i)
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

   // x[][] is initially zero everywhere expect ...
   // x[][0] is the lower boundary = 1
   // x[0][] is the left boundary = 1

   for (int i = 0; i < Ni; ++i)
      for (int j = 0; j < Nj; ++j)
         x[i][j] = 0;

   if (iProc == 0)
      for (int j = 0; j < Nj; ++j)
         x[0][j] = 1;

   if (jProc == 0)
      for (int i = 0; i < Ni; ++i)
         x[i][0] = 1;

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
         // Exchange boundary information between the neighboring
         // partitions.
         exchange_boundaries (x, Ni, Nj);
      }

      mpi_p2p_time += getElapsedTime( mpi_p2p_timer, getTimeStamp() );

      residual = 0;
      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
         {
            xtemp[i][j] = (x[i+1][j] + x[i-1][j] +
                           x[i][j+1] + x[i][j-1]) / 4.0;
            double delta = xtemp[i][j] - x[i][j];
            residual += delta*delta;
         }

      for (int i = 1; i < Ni-1; ++i)
         for (int j = 1; j < Nj-1; ++j)
            x[i][j] = xtemp[i][j];

      myTimer_t mpi_coll_timer = getTimeStamp();

      // Each process has a residual based on it's local solution.
      // Combine all of the residuals into a single (global) residual.

      mpi_coll_time += getElapsedTime( mpi_coll_timer, getTimeStamp() );

      residual = sqrt(residual);
      iteration++;
      //if (myRank == 0 and iteration % 100 == 0)
      //   printf("%d %4d: %e\n", myRank, iteration, residual);
   }

   double total_time = getElapsedTime( total_timer, getTimeStamp() );
   printf("rank=%d, timers = %f, %f, %f\n", myRank, total_time, mpi_p2p_time, mpi_coll_time);

   // Each process has local timers. Combine all three timers to get the
   // maximum across all of the processes.

   if (myRank == 0)
      printf("N = %d, Iterations = %d,  residual = %e, time = %f %f %f Procs = %d\n", N, iteration, residual, total_time, mpi_p2p_time, mpi_coll_time, numProcs);

   if (N < 100)
   {
      // Write out a solutions if the mesh isn't too large.
      char flname[12];
      sprintf(flname, "jacobi%d.out", myRank);
      FILE *f = fopen(flname,"w");
      for (int i = 0; i < Ni; ++i)
         for (int j = 0; j < Nj; ++j)
            fprintf(f,"%e %e %e\n", (iStart[iProc][jProc]+i-1)/double(N-1), (jStart[iProc][jProc]+j-1)/double(N-1), x[i][j]);
      fclose(f);
   }

   DeallocateMesh(xtemp);
   DeallocateMesh(x);

   MPI_Finalize();

   return 0;
}
