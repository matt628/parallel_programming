#include<omp.h>
#include<stdio.h>

int main () {
    int v1, v2, v3;

    printf("Liczba wątków: %d\n", omp_get_num_threads());
    printf("Jestem wątek nr: %d\n", omp_get_thread_num())

    return 0;
}