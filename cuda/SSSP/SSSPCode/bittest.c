#include "stdio.h"
#include "stdlib.h"


int main(){
    int mask = 5;
    int count=2;
    int i,j;
    int num=0;
    for(i=0;i<32;i++){
        if((1<<i & mask)!=0){
            num++;
            printf("position %d num %d\n",i,num);
           
        }
    }



    return 0;
}