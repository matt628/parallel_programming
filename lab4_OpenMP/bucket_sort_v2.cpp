#include <vector>
#include <random>
#include <stdio.h>
#include <omp.h>
#include <algorithm>
#include <argh/argh.h>
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

int main(int argc, char* argv[]) { 
    std::vector<double> data(size);  


    double fill_time_0 = omp_get_wtime();
    uniform_fill(data);
    double fill_time = omp_get_wtime() - fill_time_0;

    cout << 'd';

}