#!/bin/bash

make clean
make build

size=10000000

mkdir -p temp/algorithm

function measure_alg() {
  version=$1
  results="temp/res_${version}.txt"
  echo "fill_time, bucket_time, total_time, is_sorted" > "$results"
  for i in {1..8}; do
    ./build/measure --threads="${i}" --size="${size}" --repeat=3 --version="$version" --bucket-size=50 | tee -a "$results"
  done
}