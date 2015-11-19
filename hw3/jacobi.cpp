// There's an old C++ bug in MPI that requires it to be included before stdlib.h
//#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <my_timer.h>

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

int main (int argc, char* argv[])
{
   int N = 10; // 10 x 10 mesh.
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

   double **x = AllocateMesh<double>(N, N);
   double **xtemp = AllocateMesh<double>(N, N);

   // x[][] is initially zero everywhere expect ...
   // x[][0] is the lower boundary = 1
   // x[0][] is the left boundary = 1

   for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
         x[i][j] = 0;

   for (int j = 0; j < N; ++j)
      x[0][j] = 1;

   for (int i = 0; i < N; ++i)
      x[i][0] = 1;

   // Set xtemp = x so the boundaries are consistent.
   for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
         xtemp[i][j] = x[i][j];

   myTimer_t timer_start = getTimeStamp();

   // Iterate for some number of steps or until we converge.
   int iteration = 0;
   double residual = 1;

   while (residual > maxResidual and iteration < maxIterations)
   {
      residual = 0;
      for (int i = 1; i < N-1; ++i)
         for (int j = 1; j < N-1; ++j)
         {
            xtemp[i][j] = (x[i+1][j] + x[i-1][j] +
                           x[i][j+1] + x[i][j-1]) / 4.0;
            double delta = xtemp[i][j] - x[i][j];
            residual += delta*delta;
         }

      for (int i = 1; i < N-1; ++i)
         for (int j = 1; j < N-1; ++j)
            x[i][j] = xtemp[i][j];

      residual = sqrt(residual);
      iteration++;
      //printf("%4d: %e\n", iteration, residual);
   }

   myTimer_t timer_stop = getTimeStamp();
   printf("N = %d, Iterations = %d, residual = %e, time = %f\n", N, iteration, residual, 1000*getElapsedTime( timer_start, timer_stop));

   FILE *f = fopen("jacobi.out","w");
   for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
         fprintf(f,"%e %e %e\n", i/double(N-1), j/double(N-1), x[i][j]);
   fclose(f);

   DeallocateMesh(xtemp);
   DeallocateMesh(x);

   return 0;
}
