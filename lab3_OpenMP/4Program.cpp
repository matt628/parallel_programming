#include<omp.h>
#include<stdio.h>

int main () {
    int v1, v2, v3;

    printf("Liczba wątków: %d", omp_get_num_thread());

    return 0;
}