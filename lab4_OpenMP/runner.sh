#!/bin/bash

make clean
make build

size=100000000

mkdir -p temp

function measure_alg() {
  results="temp/results.txt"
  echo "thread_number, task_array_size, bucket_size, repeat, fill_time, bucket_time, total_time, is_sorted" > "$results"
  for i in {1..8}; do
    ./bucketSort.x --threads="${i}" --size="${size}" --repeat=4 --bucket=50 | tee -a "$results"
  done
}

measure_alg 3