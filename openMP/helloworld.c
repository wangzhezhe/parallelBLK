#include "stdio.h"
#include "stdlib.h"

int main()
{
#pragma omp parallel
    {
        // Code inside this region runs in parallel.
        int threadlocal = 1;
        printf("Parallel Hello private %d!\n", threadlocal);
    }
    printf("non parallel part\n");

    extern int parallelism_enabled;
    int n = 10;
    int c = 0;
    int num_threads = omp_get_num_threads();
    printf("num_threads before parallel %d\n", num_threads);
//if the extern parameter is 1, the parallel part will be used
//and the output will be the chaos order
#pragma omp parallel for if (parallelism_enabled)
    for (c = 0; c < n; ++c)
    {
        int num_threads = omp_get_num_threads();
        printf("num_threads during parallel %d\n", num_threads);
        int this_thread = omp_get_thread_num();
        printf("thread number %d\n", this_thread);
        printf("handle c %d\n", c);
    }
}
