#include <omp.h>
#include "stdlib.h"
#include "stdio.h"
#include "time.h"

#define COUNT 100000

int main()
{
    int sum = 100; // Assign an initial value.
    long double tCPU=clock()/(double)CLOCKS_PER_SEC;
    long double tmr = omp_get_wtime();

    printf("CLOCKS_PER_SEC %g\n",CLOCKS_PER_SEC);
    
#pragma omp parallel for reduction(+: sum)
    for (int i = 0; i < COUNT; i++)
    {
        sum += i;
    }
    printf("Sum: %d\n", sum);
    long double tmrange = omp_get_wtime() - tmr;
    long double tCPUEnd=clock()/(double)CLOCKS_PER_SEC-tCPU;
    printf("omp range omp %g\n", tmrange);
    printf("omp range cpu clock %g\n",tCPUEnd);
    return 0;
}