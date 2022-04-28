all: build run

build:
	g++ -std=c++11 -Wall -o bucketSort.x bucket_sort_v2.cpp -fopenmp
run:
	./bucketSort.x 1 1000000 5 50

clean:
	rm -rf *.x
	