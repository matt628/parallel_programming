#include<omp.h>
#include<stdio.h>
#include <time.h>
#include <stdlib.h>

int main () {
    const int MAX_SIZE = 1000;
    srand(time(NULL));   // Initialization, should only be called once.
    // int r = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.
    int arr[MAX_SIZE+10];
    int start =  omp_get_wtime();
    #pragma omp parallel 
    {
        printf("Liczba wątków: %d\n", omp_get_num_threads());
        printf("Jestem wątek nr: %d\n", omp_get_thread_num());
        int i;
        #pragma omp for 
        for(i = 1; i < MAX_SIZE/4; i++) {
            arr[(i*omp_get_thread_num())] = 1;
        }
    }
    int end =  omp_get_wtime();

    printf("time %d\n", end-start);
    // int i;
    // for(i = 0; i < MAX_SIZE; ++i) {
    //     printf("%d\n", arr[i]);
    // }
    return 0;
}