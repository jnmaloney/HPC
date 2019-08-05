#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(int argc, char** argv)
{
  #pragma omp parallel
  {
      
    #pragma omp for
    for (int i = 0; i < 10; ++i)
    {

      struct drand48_data rand_data;
      int seed = omp_get_thread_num();
      srand48_r(seed, &rand_data);
      double x;
      drand48_r(&rand_data, &x);
      printf("Random number: %f\tThread number: %i\n", x, omp_get_thread_num());

    }
  }

  return 0;
}
