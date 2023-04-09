#include <cstdio>
#include <iostream>
#include <mpi.h>

struct MyMatrix {
  double *data = nullptr;
  size_t colmns;
  size_t rows;

  MyMatrix(size_t rows, size_t colmns) {
    data = new double[colmns * rows]();
    this->rows = rows;
    this->colmns = colmns;
  }
};

struct MyGrid {
  MPI_Comm entireGridComm;
  MPI_Comm rowComm;
  MPI_Comm colmnsComm;

  int totalProcsNUmber;
  int reorderOfGrid;
  int currRowNumber;
  int currColmnNumber;
  int currRankNumber;
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

void multimplyMtrices(MyMatrix matrixA, MyMatrix matrixB,
                      MyMatrix mRes) {
  for (size_t i = 0; i < matrixA.rows; i++) {
    for (size_t j = 0; j < matrixB.colmns; j++) {
      for (size_t k = 0; k < matrixA.colmns; k++) {
        mRes.data[i * matrixB.colmns + j] +=
            matrixA.data[i * matrixA.colmns + k] *
            matrixB.data[k * matrixB.colmns + j];
      }
    }
  }
}

void testWithRandomValues(MyMatrix matrixA, MyMatrix matrixB) {
  initMatrix(matrixA);
  initMatrix(matrixB);
}

void testWithConstantValues(MyMatrix matrixA, MyMatrix matrixB) {
  for (size_t i = 0; i < matrixA.colmns * matrixA.rows; ++i) {
    matrixA.data[i] = i + 10;
  }

  for (size_t j = 0; j < matrixB.colmns * matrixB.rows; ++j) {
    matrixB.data[j] = j + 20;
  }
}

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cout
        << "Error! Enter rows and columns amount for 2 matrixes!"
        << std::endl;
    MPI_Finalize();

    return 0;
  }

  const size_t dim1 = atoi(argv[1]); // matrixA.rows
  const size_t dim2 = atoi(argv[2]); // matrixA.columns = matrixB.rows
  const size_t dim3 = atoi(argv[3]); // matrixB.columns

  const int gridHeight = atoi(argv[4]); // height of 2d grid
  const int gridWeight = atoi(argv[5]); // weight of 2d grid

  MyMatrix matrixA = MyMatrix(dim1, dim2);
  MyMatrix matrixB = MyMatrix(dim2, dim3);

  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();

  int amountOfProcs, rankOfCurrProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  if (amountOfProcs != gridHeight * gridWeight) {
    if (rankOfCurrProc == 0) {
      std::cout << "Bad input! Amount of processes should be equal "
                   "to gridHeight * gridWeight"
                << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  MPI_Comm gridEntireComm;
  int dimOfGrid = 2;
  int dims[2] = {gridHeight,
                 gridWeight}; // массив, содержащий количество
                              // процессов в каждом измерении
  int periods[2] = {0, 0};
  int reorder = 0;

  MPI_Cart_create(MPI_COMM_WORLD, dimOfGrid, dims, periods, reorder,
                  &gridEntireComm);

  int ndims = 2;
  int coords[ndims];
  // определяет координаты процесса по его номеру
  MPI_Cart_coords(gridEntireComm, rankOfCurrProc, ndims, coords);
  int rankX_rows = coords[0];
  int rankY_columns = coords[1];

  if (rankX_rows == 0 && rankY_columns == 0) {
    initMatrix(matrixA);
    printMatrix(matrixA);
    std::cout << "============" << std::endl;
    initMatrix(matrixB);
  }

  MPI_Comm rowsComm;
  int coordsRowRemain[2] = {
      0, 1}; // выбрасываем первую координату, оставляя вторую
  MPI_Cart_sub(gridEntireComm, coordsRowRemain, &rowsComm);

  MPI_Comm colmnsComm;
  int coordsColmnsRemain[2] = {
      1, 0}; // выбрасываем вторую координау, оставляя первую
  MPI_Cart_sub(gridEntireComm, coordsColmnsRemain, &colmnsComm);

  // double partA = new double[matrixA.rows / amountOfProcs *
  // matrixA.colmns];

  MyMatrix partA =
      MyMatrix(matrixA.rows / amountOfProcs, matrixA.colmns);
  if (coords[1] == 0) {
    MPI_Scatter(matrixA.data,
                matrixA.rows / amountOfProcs * matrixA.colmns,
                MPI_DOUBLE, partA.data, partA.rows * partA.colmns,
                // matrixA.rows / amountOfProcs * matrixA.colmns,
                MPI_DOUBLE, 0, colmnsComm);
  }
  MPI_Bcast(partA.data, partA.rows * partA.colmns, MPI_DOUBLE, 0,
            rowsComm);

  MPI_Datatype bColmns;

  int numberOfBlocks = matrixB.rows;
  int numberOfElemsInEachBlock = matrixB.colmns / amountOfProcs;
  int stride =
      matrixB.colmns; // nmber of elems between start of each block
  MPI_Type_vector(numberOfBlocks, numberOfElemsInEachBlock, stride,
                  MPI_DOUBLE, &bColmns);
  MPI_Type_commit(&bColmns);

  MPI_Type_free(bColmns);
  // MPI_Comm_free(&entireGridComm);
  MPI_Finalize();

  return 0;
}
