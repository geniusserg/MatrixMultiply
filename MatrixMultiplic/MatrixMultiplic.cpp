#define N 16

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
MPI_Comm row_comm;
MPI_Comm col_comm;

int* matrixA;
int* matrixB;
int* matrixC;

int* matrixAoriginal;
int* matrixBoriginal;
int* matrixCoriginal;

void multiply_matrix() {
	for (int i = 0; i < submatrix_size; i++) {
		for (int j = 0; j < submatrix_size; j++) {
			for (int k = 0; k < submatrix_size; k++) {
				matrixC[i * submatrix_size + j] += matrixA[i * submatrix_size + k] * matrixB[j + k * submatrix_size];
			}
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

void prepare_block(int a, int b, int* matrix, int* matrixOriginal) {
	long long index = 0;
	for (int i = submatrix_size * a; i < submatrix_size * (a + 1); i++) {
		for (int j = submatrix_size * b; j < submatrix_size * (b + 1); j++) {
			matrix[index] = matrixOriginal[i * N + j];
			index++;
		}
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
	ndims[0] = ndims[1] = ndim;
	submatrix_size = N / ndim;
	MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, period, 0, &com);

	int subdims[2] = {0, 1};
	MPI_Cart_sub(com, subdims, &row_comm);
	subdims[0] = 1;
	subdims[1] = 0;
	MPI_Cart_sub(com, subdims, &col_comm);

	MPI_Cart_coords(com, rank, 2, coords);	

	//prepare data
	if (rank == 0) {
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
	}

	matrixA = new int[submatrix_size * submatrix_size];
	matrixB = new int[submatrix_size * submatrix_size];
	matrixC = new int[submatrix_size * submatrix_size];

	// A B distribution for all processes
	if (rank == 0) {
		for (int i = 0; i < ndim; i++) {
			for (int j = 0; j < ndim; j++) {
				if (i == 0 && j == 0) {
					continue;
				}
				std::cout << "RANK ONE SENDING BLOCK" << std::endl;
				prepare_block(i, j, matrixA, matrixAoriginal);
				prepare_block(i, j, matrixB, matrixBoriginal);
				int trank = 0;
				int tcoords[2] = { i, j };
				MPI_Cart_rank(com, tcoords, &trank);
				MPI_Send(matrixA, submatrix_size * submatrix_size, MPI_DOUBLE, trank, 0, com);
				std::cout << "sent block" << std::endl;
				MPI_Send(matrixB, submatrix_size * submatrix_size, MPI_DOUBLE, trank, 1, com);
				std::cout << "RANK ONE SENDING BLOCK ..." << std::endl;
			}
		}
		prepare_block(0, 0, matrixA, matrixAoriginal);
		prepare_block(0, 0, matrixB, matrixBoriginal);
	}
	else {
		MPI_Status mStatus;
		MPI_Recv(matrixA, submatrix_size * submatrix_size, MPI_DOUBLE, rank, 0, com, &mStatus);
		MPI_Recv(matrixB, submatrix_size * submatrix_size, MPI_DOUBLE, rank, 1, com, &mStatus);
		std::cout << "received" << std::endl;
	}

	

	// main loop 
	/*
	for (int iter = 0; iter < ndim; iter++) {
		// A` distribution
		if (coords[1] == (coords[0] + iter) % ndim) {
			MPI_Bcast(matrixA, submatrix_size * submatrix_size, MPI_DOUBLE, coords[1], row_comm);
		}

		// multiply matrix
		multiply_matrix();

		// B` distribution
		int next = (coords[0] + 1)%(ndim - 1);
		int prev = (coords[0] - 1);
		if (coords[0] == 0) {
			prev = ndim - 1;
		}
		MPI_Sendrecv_replace(matrixB, submatrix_size * submatrix_size, MPI_DOUBLE, next, 0, prev, 0, col_comm, &mStatus);
	}
	*/

	if (rank == 0){
		std::cout << "finalize" << std::endl;
	}

	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixC;
	delete[] matrixAoriginal;
	delete[] matrixBoriginal;
	delete[] matrixCoriginal;

	MPI_Finalize();
	return 0;
}
