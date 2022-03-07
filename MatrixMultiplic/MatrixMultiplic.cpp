#define N 4000

#include <stdlib.h>
#include <iostream>
#include "mpi.h"
#include "math.h"

int size, rank;

int ndim;
int ndims[2] = { 0,0 };
int period[2] = {0, 0};
int coords[2];
int submatrix_size = 7;
MPI_Comm com;

int* matrixA;
int* matrixB;
int* matrixC;

int* matrixAoriginal;
int* matrixBoriginal;
int* matrixCoriginal;

void multiply_matrix() {
	for (int i = 0; i < submatrix_size; i++) {
		for (int j = 0; j < submatrix_size; j++) {
			matrixC[i * submatrix_size + j] = matrixA[i * submatrix_size + j] * matrixB[i + j * submatrix_size];
		}
	}
}

void print_matrix() {
	for (int i = 0; i < submatrix_size; i++) {
		for (int j = 0; j < submatrix_size; j++) {
			std::cout << matrixC[i * submatrix_size + j] << " ";
		}
		std::cout << std::endl;
	}
}

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
	MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, period, 0, &com);

	MPI_Cart_coords(com, rank, 2, coords);	

	submatrix_size = N / ndim;

	matrixA = new int[submatrix_size * submatrix_size];
	matrixB = new int[submatrix_size * submatrix_size];
	matrixC = new int[submatrix_size * submatrix_size];

	matrixAoriginal = new int[N * N];
	matrixBoriginal = new int[N * N];
	matrixCoriginal = new int[N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matrixCoriginal[i * N + j] = 0;
			matrixBoriginal[i * N + j] = 10;
			matrixAoriginal[i * N + j] = 10;
		}
	}

	long long index = 0;
	for (int i = submatrix_size*coords[0]; i < submatrix_size * (coords[0]+1); i++) {
		for (int j = submatrix_size * coords[1]; j < submatrix_size * (coords[1]+1); j++) {
			matrixC[index] = 0;
			matrixB[index] = matrixBoriginal[i * submatrix_size + j];
			matrixA[index] = matrixAoriginal[i * submatrix_size + j];
			index++;
		}
	}

	
	multiply_matrix();

	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixC;
	delete[] matrixAoriginal;
	delete[] matrixBoriginal;
	delete[] matrixCoriginal;

	MPI_Finalize();
	return 0;
}
