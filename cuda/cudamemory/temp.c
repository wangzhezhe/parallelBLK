#include "stdio.h"
#include "stdlib.h"

void getGuardRecord(int** gRecord)
{
    int i;
    int gindex = 0;
    //(*gRecord)[0]=0;
    //(*gRecord)[1]=1;
    //(*gRecord)[2]=2;
    int index=0;
    for (i = 0; i < 5; i++)
    {
           (*gRecord)[index]=i;
           printf("guard %d i %d index %d\n",(*gRecord)[index],i,index);
           index++;
        
    }
  
}

void getrecord(int* gRecord)
{
    int i;
    int gindex = 0;
    gRecord[0]=1;
    gRecord[1]=2;
}

int main(){
    int* g=malloc(sizeof(int)*5);
    getrecord(g);
    int i;
    for (i = 0; i < 5; i++)
    {
        printf("guard %d i %d \n",*(g+i),i);     
    }

    printf("%d\n",1%1);
    return 0;
}