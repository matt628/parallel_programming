#include<omp.h>
#include<stdio.h>
#include <time.h>
#include <stdlib.h>

int main () {
    const int MAX_SIZE = 1000000;
    srand(time(NULL));   // Initialization, should only be called once.
    // int r = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.
    int arr[MAX_SIZE+10];
    omp_set_num_threads(1);
    double start =  omp_get_wtime();
    #pragma omp parallel
    {
        int i;
        #pragma omp for private (i)
        for(i = 1; i < MAX_SIZE/4; i++) {
            arr[(i*omp_get_thread_num())] = 1;
        }
    }
    double end =  omp_get_wtime();

    printf("time %f\n", end-start);
    // int i;
    // for(i = 0; i < MAX_SIZE; ++i) {
    //     printf("%d\n", arr[i]);
    // }
    return 0;
}