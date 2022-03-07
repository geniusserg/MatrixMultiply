

#include <stdlib.h>
#include <iostream>
#include "mpi.h"
#include "math.h"

int size = 0, rank;

int ndim;
int ndims[2] = { 0,0 };
int period[2] = {0, 0};
int coords[2];
MPI_Comm com;

int main(int argc, char** argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!((size == 1) || (size == 2) || (size == 4) || (size == 9) || (size == 16))) {
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	ndim = (int)sqrt(size);
	ndims[0] = ndim;
	ndims[1] = ndim;
	MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, period, 1, &com);

	MPI_Cart_coords(com, rank, ndim, coords);

	std::cout << coords[0] << coords[1] << std::endl;

	MPI_Finalize();
	return 0;
}
