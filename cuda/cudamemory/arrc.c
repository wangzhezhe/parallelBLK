#include <stdio.h>
#define MAX 20
typedef int Values[MAX];

int changeArr(Values vals2) 
{
    vals2[0] = 200;
    vals2[1] = 100;
    printf("---in function---\n");
    printf("%d and %d\n", vals2[0],vals2[1]);
    
    return 0;
}   

int main (int argc, char *argv[]) 
{
    Values vals;
    changeArr(&vals);
    printf("%d and %d\n", vals[0],vals[1]);
    return 0;
}