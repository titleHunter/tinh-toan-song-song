#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#define M 500
#define N 500

int main()
{
  double epsilon = 0.001;
  double diff = epsilon; // difference between new and old value of a node
  int i, j, iterations, iterations_print; // loop variables
  double mean; // mean value of the boundary 
  double diffInThread; // diff in each thread
  double u[M][N]; // for storing the old values of the grid
  double w[M][N]; // for storing the current values of the grid
  double wtime; // for estimating computed time
  double nThreads;

  printf("\n");
  printf("HEATED_PLATE_OPENMP\n");
  printf("  C/OpenMP version\n");
  printf("  A program to solve for the steady state temperature distribution\n");
  printf("  over a rectangular plate.\n");
  printf("\n");
  printf("  Spatial grid of %d by %d points.\n", M, N);
  printf("  The iteration will be repeated until the change is <= %e\n", epsilon);
  printf("  Number of processors available = %d\n", omp_get_num_procs());
  printf("  Number of threads =              %d\n", omp_get_max_threads());

  
  mean = 0.0;

#pragma omp parallel shared(w) private(i, j)
  {
  // fix boundary condition
#pragma omp for
    for (i = 1; i < M - 1; i++)
    {
      w[i][0] = 100.0;
    }
#pragma omp for
    for (i = 1; i < M - 1; i++)
    {
      w[i][N - 1] = 100.0;
    }
#pragma omp for
    for (j = 0; j < N; j++)
    {
      w[M - 1][j] = 100.0;
    }
#pragma omp for
    for (j = 0; j < N; j++)
    {
      w[0][j] = 0.0;
    }
    
	//	average boundary values
#pragma omp for reduction(+ : mean)
    for (i = 1; i < M - 1; i++)
    {
      mean = mean + w[i][0] + w[i][N - 1];
    }
#pragma omp for reduction(+ : mean)
    for (j = 0; j < N; j++)
    {
      mean = mean + w[M - 1][j] + w[0][j];
    }
  }

  mean = mean / (double)(2 * M + 2 * N - 4);
  printf("\n");
  printf("  MEAN = %f\n", mean);

	// initialize temperature array
#pragma omp parallel shared(mean, w) private(i, j)
  {
#pragma omp for
    for (i = 1; i < M - 1; i++)
    {
      for (j = 1; j < N - 1; j++)
      {
        w[i][j] = mean;
      }
    }
  }
  
  // main loop
  iterations = 0;
  iterations_print = 1;
  printf("\n");
  printf(" Iteration Change\n");
  printf("\n");
  wtime = omp_get_wtime();

  diff = epsilon;

  while (epsilon <= diff)
  {
#pragma omp parallel shared(u, w) private(i, j)
    {

  // Save the old solution in U.

#pragma omp for
      for (i = 0; i < M; i++)
      {
        for (j = 0; j < N; j++)
        {
          u[i][j] = w[i][j];
        }
      }

  // Determine the new estimate of the solution at the interior points.

#pragma omp for
      for (i = 1; i < M - 1; i++)
      {
        for (j = 1; j < N - 1; j++)
        {
          w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
        }
      }
    }
 
    // Compute the maximun diff
    diff = 0.0;
#pragma omp parallel shared(diff, u, w) private(i, j, diffInThread)
    {
      diffInThread = 0.0;
#pragma omp for
      for (i = 1; i < M - 1; i++)
      {
        for (j = 1; j < N - 1; j++)
        {
          if (diffInThread < fabs(w[i][j] - u[i][j]))
          {
            diffInThread = fabs(w[i][j] - u[i][j]);
          }
        }
      }
#pragma omp critical
      {
        if (diff < diffInThread)
        {
          diff = diffInThread;
        }
      }
    }

    iterations++;
    if (iterations == iterations_print)
    {
      printf("  %8d  %f\n", iterations, diff);
      iterations_print = 2 * iterations_print;
    }
  }
  wtime = omp_get_wtime() - wtime;

  printf("\n");
  printf("  %8d  %f\n", iterations, diff);
  printf("\n");
  printf("  Error tolerance achieved.\n");
  printf("  Wallclock time = %f\n", wtime);
  /*
  Terminate.
*/
  printf("\n");
  printf("HEATED_PLATE_OPENMP:\n");
  printf("  Normal end of execution.\n");

  return 0;
}
