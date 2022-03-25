#include "parameters.h"

double sender(int size) {
    char *buff = malloc(size);
    int i;
    double start = MPI_Wtime();
    for (i = 0; i < N; i++) {
        // int MPI Ssend(void *buf, int count, MPI Datatype datatype, int dest, int tag, MPIP Comm comm)
        MPI_SSend(buff, size, MPI_BYTE, RECEIVER, 0, MPI_COMM_WORLD);
    }
    return MPI_Wtime() - start;
}

void receiver(int size) {
    char *buff = malloc(size);
    int i;
    for (i = 0; i < N; i++) {
        MPI_Recv(buff, size, MPI_BYTE, SENDER, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void test(int rank) {
    int size;
    printf("Bandwidth[Mbit/s], size[B]")
    for (size = 1; size <= MAX_SIZE; size *= 2) {
        if (rank == SENDER) {
            double time = sender(size);
            printf("%.5f, %d\n", N*size/time/1000000*2, size);
        } else if (rank == RECEIVER) {
            receiver(size);
        }
    }
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv); // manual initialisation

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //init of default communicator

    test(rank);

    MPI_Finalize(); // manual finalisation
    return 0;
}