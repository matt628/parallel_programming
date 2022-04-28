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

// GLOBAL 
int threads = 1, size = 1e6, repeat = 1, 	bucket_size = 50;

template<int min = 0, int max = 1>
void uniform_fill(std::vector<double>& array) {
  printf("THREADS %d", threads);
#pragma omp parallel num_threads(threads)
  {
	std::uniform_real_distribution<double> distribution(min, max);
	std::default_random_engine generator;
	int t_thread = omp_get_thread_num();
	generator.seed(t_thread * time(NULL) + 17);

#pragma omp for SCHEDULE
	for (size_t i = 0; i < array.size(); i++) {
	  array[i] = distribution(generator);
	}
  }
}

void bucket_sort(std::vector<double>& array, int no_buckets) {
    
}


bool verify(std::vector<double>& supposedly_sorted, std::vector<double>& original) {
  std::sort(original.begin(), original.end());

  bool are_equal = supposedly_sorted == original;

  if(are_equal) {
    // printf("Sorted Successful\n");
  } else {
    // printf("!!! NOT Sorted !!!\n");
  }
  return are_equal;
} 

int main(int argc, char* argv[]) { 
    argh::parser cmdl(argv);

    std::vector<double> data(size);  

    threads = cmdl["threads"];
    size = cmdl["size"];
    repeat = cmdl["repeat"];
    bucket_size = cmdl["bucket"];
    for(int i = 0; i<repeat; i++){
      double fill_time_0 = omp_get_wtime();
      uniform_fill(data);
      double fill_time = omp_get_wtime() - fill_time_0;

      std::vector<double> original = data;

      double bucket_sort_1 = omp_get_wtime();
      bucket_sort(data, 8);
      double bucket_sort_time = omp_get_wtime() - bucket_sort_1;


      bool isSorted = verify(data, original);
      // printf("thread_number, task_array_size, bucket_size, repeat, fill_time, bucket_time, total_time, is_sorted\n");
      printf("%d, %d, %d, %d, %lf, %lf, %lf, %d\n", threads, size, bucket_size, repeat, fill_time, bucket_sort_time, fill_time+bucket_sort_time, isSorted);
    }
}