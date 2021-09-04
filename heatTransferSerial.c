#include <omp.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define M 500
#define N 500

int main()
{
  double epsilon = 0.001;
  double mean = 0.0;
  double diff; // difference between new and old value of a node
  double u[M][N]; // for storing the old value of the grid
  double w[M][N]; // for storing the current value of the grid
  double wtime;   // for estimate computed time
  int i, j, iterations, iterations_print; // loop variables

  printf("\n");
  printf("HEATED_PLATE_OPENMP\n");
  printf("  C/Serial version\n");
  printf("  A program to solve for the steady state temperature distribution\n");
  printf("  over a rectangular plate.\n");
  printf("\n");
  printf("  Spatial grid of %d by %d points.\n", M, N);
  printf("  The iteration will be repeated until the change is <= %e\n", epsilon);

  // fix boundary conditions
  for (i = 1; i < M - 1; i++)
  {
    w[i][0] = 100.0;
    w[i][N - 1] = 100.0;
  }
  for (j = 0; j < N; j++)
  {
    w[0][j] = 100.0;
    w[M - 1][j] = 0.0;
  }


  // average the boundary values
  for (i = 1; i < M - 1; i++)
  {
    mean = mean + w[i][0] + w[i][N - 1];
  }
  for (j = 0; j < N; j++)
  {
    mean = mean + w[M - 1][j] + w[0][j];
  }

  mean = mean / (double)(2 * M + 2 * N - 4);
  printf("\n");
  printf("  MEAN = %f\n", mean);

  // initialize temperature array
  for (i = 1; i < M - 1; i++)
  {
    for (j = 1; j < N - 1; j++)
    {
      w[i][j] = mean;
    }
  }

  // main loop
  iterations = 0;
  iterations_print = 1;
  printf("\n");
  printf(" Iteration  Change\n");
  printf("\n");
  wtime = omp_get_wtime();
  diff = epsilon;

  while (epsilon <= diff)
  {
    // Save the old solution in U
    for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
      {
        u[i][j] = w[i][j];
      }
    }

  	// Determine the new estimate of the solution at the interior points.
    for (i = 1; i < M - 1; i++)
    {
      for (j = 1; j < N - 1; j++)
      {
        w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
      }
    }

    // Compute the maximun diff
    diff = 0.0;
    for (i = 1; i < M - 1; i++)
    {
      for (j = 1; j < N - 1; j++)
      {
        if (diff < fabs(w[i][j] - u[i][j]))
        {
          diff = fabs(w[i][j] - u[i][j]);
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
  printf("HEATED_PLATE_SERIAL:\n");
  printf("  Normal end of execution.\n");
}
