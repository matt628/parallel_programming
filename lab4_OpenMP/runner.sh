#!/bin/bash

make clean
make build

size=100000000

mkdir -p temp

function measure_alg() {
  results="temp/results.txt"
  echo "thread_number, task_array_size, bucket_size, repeat, fill_time, bucket_time, total_time, is_sorted" > "$results"
  for i in {1..20}; do
    ./bucketSort.x ${i} ${size} 4 50 "false" | tee -a "$results"
  done
}

measure_alg