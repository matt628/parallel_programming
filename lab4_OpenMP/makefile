all: build run

build:
	g++ -std=c++11 -Wall -o bucketSort.x bucket_sort_v2.cpp -fopenmp
run:
	./bucketSort.x --threads=1 --size=1000000 --repeat=5 --bucket=50

clean:
	rm -rf *.x
	