#include "stdio.h"
#include "stdlib.h"

int main()
{
    int n = 10;
    int c = 0;
    int nt=5;
#pragma omp parallel num_threads(nt)
    {
        {
            int this_thread = omp_get_thread_num();
            printf("thread id %d\n", this_thread);

#pragma omp for
            for (c = 0; c < n; ++c)
            {

                printf("handle c %d\n", c);
            }
        }
    }
}