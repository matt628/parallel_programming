#include <vector>
#include <random>
#include <stdio.h>
#include <omp.h>
#include <algorithm>
#include "argh/argh.h"
#include <iostream>

#ifndef SCHEDULE
#define SCHEDULE schedule(static)
#define SCHEDULE_STR "schedule(static)"
#endif

int threads = 1, size = 1e6, repeat = 1;

template<int min=0, int max=1>
void uniform_fill(std::vector<double>& array) {
  int size = array.size();
  #pragma omp parallel num_threads(threads)
  {
    std::uniform_real_distribution<double> distribution(min, max);
    std::default_random_engine generator; 
    int t_thread = omp_get_thread_num();
    generator.seed(t_thread * time(NULL) + 17);
    
    #pragma omp for SCHEDULE
    for (int i = 0; i < size; i++) {
      array[i] = distribution(generator);
    }
  }
}

void bucket_sort(std::vector<double>& array, int no_buckets) {

}

int main(int argc, char* argv[]) { 
    std::vector<double> data(size);  



    // cmdl({ "-t", "--threads"}) >> threads;
    // cmdl({ "-s", "--size" }) >> size;
    // cmdl({ "-r", "--repeat" }) >> repeat;


    double fill_time_0 = omp_get_wtime();
    uniform_fill(data);
    double fill_time = omp_get_wtime() - fill_time_0;

    argh::parser cmdl(argv);

    std::vector<double> original = data;

    double bucket_sort_1 = omp_get_wtime();
    bucket_sort(data, 8);
    double bucket_sort_time = omp_get_wtime() - bucket_sort_1;


    printf("fill_time, bucket_time\n");
    printf("%lf, %lf\n", fill_time, bucket_sort_time);

}