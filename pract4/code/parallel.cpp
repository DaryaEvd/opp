#include <cstdio>
#include <iostream>
#include <mpi.h>

struct MyMatrix {
  double *data = nullptr;
  size_t colmns;
  size_t rows;

  MyMatrix(size_t rows, size_t colmns) {
    data = new double[colmns * rows];
    this->rows = rows;
    this->colmns = colmns;
  }
};

void initMatrix(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      matrix.data[i * matrix.colmns + j] = rand() % 100 + 15;
    }
  }
}

void zerofyMatrix(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      matrix.data[i * matrix.colmns + j] = 0;
    }
  }
}

void freeMatrix(MyMatrix matrix) { delete[] matrix.data; }

void printMatrix(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.rows; ++i) {
    for (size_t j = 0; j < matrix.colmns; ++j) {
      std::cout << matrix.data[i * matrix.colmns + j] << "\t";
    }
    std::cout << std::endl;
  }
}

void multimplyMtrices(MyMatrix m1, MyMatrix m2, MyMatrix mRes) {
  for (size_t i = 0; i < m1.rows; i++) {
    for (size_t j = 0; j < m2.colmns; j++) {
      for (size_t k = 0; k < m1.colmns; k++) {
        mRes.data[i * m2.colmns + j] +=
            m1.data[i * m1.colmns + k] * m2.data[k * m2.colmns + j];
      }
    }
  }
}

void testWithRandomValues(MyMatrix m1, MyMatrix m2) {
  initMatrix(m1);
  initMatrix(m2);
}

void testWithConstantValues(MyMatrix m1, MyMatrix m2) {
  for (size_t i = 0; i < m1.colmns * m1.rows; ++i) {
    m1.data[i] = i + 10;
  }

  for (size_t j = 0; j < m2.colmns * m2.rows; ++j) {
    m2.data[j] = j + 20;
  }
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cout
        << "Error! Enter rows and columns amount for 2 matrixes!"
        << std::endl;
    return 0;
  }

  const size_t dim1 = atoi(argv[1]); // m1.rows
  const size_t dim2 = atoi(argv[2]); // m1.columns = m2.rows
  const size_t dim3 = atoi(argv[3]); // m2.columns

  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();

  int amountOfProcs, rankOfCurrProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  std::cout << "amountOfProcs: " << amountOfProcs << std::endl;

  if (rankOfCurrProc == 0) {
    // testWithRandomValues(m1, m2);
    testWithConstantValues(m1, m2);

    std::cout << "Matrix1 : " << std::endl;
    printMatrix(m1);
    std::cout << "Matrix2 : " << std::endl;
    printMatrix(m2);
  }

  int reorder = 0;
  int dims[2] = {3, 8};
  int periods[2] = {0, 0};

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm2d);

  int coords[2];

  MPI_Cart_coords(comm2d, rankOfCurrProc, 2, coords[]);

  int currNewRank;

  MPI_Comm comm2d;
  MPI_Dims_create(amountOfProcs, 2, dims);

  int sizeX = dims[0];
  int sizeY = dims[1];

  MPI_Comm_rank(comm2d, &currNewRank);

  MPI_Cart_get(comm2d, 2, dims, periods, coords);
  int rankX = coords[0], rankY = coords[1];

  MyMatrix m1 = MyMatrix(dim1, dim2);
  MyMatrix m2 = MyMatrix(dim2, dim3);

  MyMatrix mRes = MyMatrix(dim1, dim3);

  if (m1.colmns != m2.rows || m1.rows != mRes.rows ||
      m2.colmns != mRes.colmns) {
    std::cerr << "Error! You've entered a wrong dimention to matrices"
              << std::endl;
    MPI_Finalize();
    return 0;
  }

  zerofyMatrix(mRes);
  multimplyMtrices(m1, m2, mRes);
  std::cout << "Matrix mRes : " << std::endl;
  printMatrix(mRes);

  double endTime = MPI_Wtime();
  if (rankOfCurrProc == 0) {
    std::cout << "Time taken: " << endTime - startTime << " sec"
              << std::endl;
  }

  freeMatrix(m1);
  freeMatrix(m2);
  freeMatrix(mRes);
  MPI_Finalize();
}
