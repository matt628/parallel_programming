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

template<typename Function>
double inline timeit(Function&& timed_function) {
  double time_0 = omp_get_wtime();
  timed_function();
  return omp_get_wtime() - time_0;
}

// GLOBAL 
int threads, size, repeat, bucket_size;

template<int min = 0, int max = 1>
void uniform_fill(std::vector<double>& array) {
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

void generate_perfect_array(std::vector<double>& array) {
  int seed = 21567;
  std::shuffle(array.begin(), array.end(), std::default_random_engine(seed));
}

void bucket_sort(std::vector<double>& array, int no_buckets) {
    
}


template<int max = 1>
void parallel_bucket_sort_1(std::vector<double>& array) {
  // allocate memory for buckets.
  int no_buckets = size / bucket_size;
  int buckets_per_thread = no_buckets / threads;
  int estimated_bucket_size = std::max((int)array.size() / no_buckets, 1);
  std::vector<std::vector<double>> buckets(no_buckets);
  for (auto bucket : buckets) {
	  bucket.reserve(estimated_bucket_size);
  }

#pragma omp parallel shared(buckets) firstprivate(no_buckets) num_threads(threads)
  {
	int tid = omp_get_thread_num();

	// we start by populating buckets
	// each thread fills its own buckets.
	// double split_to_buckets_time = timeit([&] {
	  
	// });
  for (size_t i = tid; i < tid + array.size(); i++) {
		int bucket_index = std::min((int)(no_buckets * array[i % array.size()] / max), no_buckets - 1);
		if (tid * buckets_per_thread <= bucket_index && (bucket_index < (tid + 1) * buckets_per_thread || tid == threads - 1)) {
		  buckets[bucket_index].push_back(array[i % array.size()]);
		}
  }
	// now each thread sorts its share of buckets.
	double sort_buckets_time = timeit([&] {
#pragma omp for schedule(static)
	  for (int bucket_index = 0; bucket_index < no_buckets; bucket_index++) {
		std::sort(buckets[bucket_index].begin(), buckets[bucket_index].end());
	  }
	});

	// after the buckets have been sorted
	double write_sorted_buckets_time = timeit([&] {

	  // we compute indices where to start writing in the original array.
	  std::vector<int> bucket_idx_to_array_idx_table(no_buckets);
	  for (int bucket_index = 1; bucket_index < no_buckets; bucket_index++) {
		bucket_idx_to_array_idx_table[bucket_index] =
			buckets[bucket_index - 1].size() + bucket_idx_to_array_idx_table[bucket_index - 1];
	  }

	  // finally, we can write the result.
#pragma omp for schedule(static)
	  for (int bucket_index = 0; bucket_index < no_buckets; bucket_index++) {
		int start_idx = bucket_idx_to_array_idx_table[bucket_index];
		for (size_t i = 0; i < buckets[bucket_index].size(); i++) {
		  array[start_idx + i] = buckets[bucket_index][i];
		}
	  }
	});

	// // update measurements at the end.
	// if (tid == 0) {
	//   measurement.split_to_buckets_time = split_to_buckets_time;
	//   measurement.sort_buckets_time = sort_buckets_time;
	//   measurement.write_sorted_buckets_time = write_sorted_buckets_time;
	// }
  // }
}
}

//a) czy potrzebna jest jakaś ochrona danych wspólnych (
  // tablica początkowa: przy odczycie i przy zapisie; 
  //kubełki: przy zapisie, sortowaniu  kubełka, - nie kady wątek ma swoj kubelek
  // odczycie 
//b) jaki jest rząd złożoności obliczeniowej algorytmu, a jaka jest praca algorytmu równoległego, czy algorytm jest sekwencyjnie-efektywny?

void perfect_bucket_sort(std::vector<double>& array){
  bucket_size =1; 
  int no_buckets = size / bucket_size;
  int buckets_per_thread = no_buckets / threads;
  std::vector<int> buckets(no_buckets);

  // each thread reads the whole array
  #pragma omp parallel shared(buckets) firstprivate(no_buckets) num_threads(threads)
  {
	int tid = omp_get_thread_num();
    for (size_t i = tid; i < tid + array.size(); i++) {
      size_t index = i % array.size();
      int bucket_index = std::min((int)(no_buckets * array[i % array.size()] / max), no_buckets - 1);
      if (tid * buckets_per_thread <= bucket_index && (bucket_index < (tid + 1) * buckets_per_thread || tid == threads - 1)) {
        buckets[bucket_index].push_back(array[index]);
      }
  }
  
  }
  // put numbers into own buckets

  // sort own bucket

  // put results together
}

bool verify(std::vector<double>& supposedly_sorted, std::vector<double>& original) {
  std::sort(original.begin(), original.end());
  bool are_equal = supposedly_sorted == original;
  return are_equal;
} 

int main(int argc, char* argv[]) { 
    // argh::parser cmdl(argv);

    std::vector<double> data(size);  

    threads = atoi(argv[1]);
    size = atoi(argv[2]);
    repeat = atoi(argv[3]);
    bucket_size = atoi(argv[4]);
    char use_perfect_data_string = argv[5][0];

    bool use_perfect_data = use_perfect_data_string == 't' ? true : false;
    for(int i = 0; i<repeat; i++){
      double fill_time_0 = omp_get_wtime();
      if(use_perfect_data) {
        generate_perfect_array(data);
      } else {
        uniform_fill(data);
      }
      double fill_time = omp_get_wtime() - fill_time_0;

      std::vector<double> original = data;

      double bucket_sort_1 = omp_get_wtime();
      if(use_perfect_data){
        perfect_bucket_sort(data);
      } else {
        parallel_bucket_sort_1(data);
      }
      // parallel_bucket_sort_1(data);
      double bucket_sort_time = omp_get_wtime() - bucket_sort_1;


      bool isSorted = verify(data, original);
      // printf("thread_number, task_array_size, bucket_size, repeat, fill_time, bucket_time, total_time, is_sorted\n");
      printf("%d, %d, %d, %d, %lf, %lf, %lf, %s,%d\n", threads, size, bucket_size, repeat, fill_time, bucket_sort_time, fill_time+bucket_sort_time, use_perfect_data, isSorted);
    }
}