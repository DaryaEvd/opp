#include <iostream>
#include <mpi.h>

/* Axises are
|-----Y
|
|
X
*/

static int dimOfAnyGrid = 2;
static int X_AXIS = 0;
static int Y_AXIS = 1;

struct MyMatrix {
  size_t row;
  size_t column;
  double *data = nullptr;

  MyMatrix(size_t row, size_t column) {
    this->row = row;
    this->column = column;
    data = new double[row * column]();
  }
};

void freeMatrix(MyMatrix matrix) { delete[] matrix.data; }

void initMatrix(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.row; ++i) {
    for (size_t j = 0; j < matrix.column; ++j) {
      matrix.data[i * matrix.column + j] = rand() % 100 + 15;
    }
  }
}

void multiplyMtrices(MyMatrix m1, MyMatrix m2, MyMatrix mRes) {
  for (size_t i = 0; i < m1.row; i++) {
    for (size_t j = 0; j < m2.column; j++) {
      for (size_t k = 0; k < m1.column; k++) {
        mRes.data[i * m2.column + j] +=
            m1.data[i * m1.column + k] * m2.data[k * m2.column + j];
      }
    }
  }
}

void printMat(MyMatrix matrix) {
  for (size_t i = 0; i < matrix.row; i++) {
    for (size_t j = 0; j < matrix.column; j++) {
      std::cout << matrix.data[i * matrix.column + j] << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cout
        << "Error! Enter rowsInC and columns amount for 2 matrixes!"
        << std::endl;
    return 0;
  }

  srand(time(nullptr)); // just for random vales each time
  // srand(0); // for testint the same values in one task

  const size_t dim1 = atoi(argv[1]);
  const size_t dim2 = atoi(argv[2]);
  const size_t dim3 = atoi(argv[3]);

  const int rowsGrid = atoi(argv[4]);    // height of 2d grid
  const int columnsGrid = atoi(argv[5]); // weight of 2d grid

  MPI_Init(&argc, &argv);
  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  if (amountOfProcs != rowsGrid * columnsGrid) {
    if (rankOfCurrProc == 0) {
      std::cout << "Bad input! Amount of processes should be equal "
                   "to rowsGrid * columnsGrid"
                << std::endl;
    }
    MPI_Finalize();
    return 0;
  }

  int dimSize[] = {rowsGrid, columnsGrid};
  int periods[] = {0};

  MPI_Comm commGrid;

  double startt = MPI_Wtime();

  MPI_Cart_create(MPI_COMM_WORLD, dimOfAnyGrid, dimSize, periods, 0,
                  &commGrid); // reorder = 0

  int coordsOfCurrProc[dimOfAnyGrid];
  MPI_Cart_coords(commGrid, rankOfCurrProc, dimOfAnyGrid,
                  coordsOfCurrProc);

  MPI_Comm commRow;
  int remainX[] = {X_AXIS, Y_AXIS};
  MPI_Cart_sub(commGrid, remainX, &commRow);

  MPI_Comm commColumn;
  int remainY[] = {Y_AXIS, X_AXIS};
  MPI_Cart_sub(commGrid, remainY, &commColumn);

  MyMatrix A(dim1, dim2);
  MyMatrix B(dim2, dim3);
  MyMatrix C(dim1, dim3);

  if (rankOfCurrProc == 0) {
    initMatrix(A);
    initMatrix(B);

    // std::cout << "A[" << A.row << " x " << A.column << "]"
    //           << std::endl;
    // printMat(A);
    // std::cout << "B[" << B.row << " x " << B.column << "]"
    //           << std::endl;
    // printMat(B);
  }

  MyMatrix partA(dim1 / dimSize[X_AXIS], dim2);

  if (coordsOfCurrProc[Y_AXIS] == 0) {
    MPI_Scatter(A.data, partA.row * partA.column, MPI_DOUBLE,
                partA.data, partA.row * partA.column, MPI_DOUBLE, 0,
                commColumn);
  }
  MPI_Bcast(partA.data, partA.row * partA.column, MPI_DOUBLE, 0,
            commRow);

  MyMatrix partB(dim2, dim3 / dimSize[Y_AXIS]);

  MPI_Datatype bSendType; // type of columns

  // number of blocks, numbers of elems in each block, number of
  // elems between the start of each block, old type, new type
  MPI_Type_vector(dim2, partB.column, dim3, MPI_DOUBLE, &bSendType);
  MPI_Type_commit(&bSendType);

  if (coordsOfCurrProc[X_AXIS] == 0 &&
      coordsOfCurrProc[Y_AXIS] == 0) {
    for (int i = 1; i < columnsGrid; i++) {
      MPI_Send(B.data + partB.column * i, 1, bSendType, i, 11,
               commRow);
    }

    for (int i = 0; i < B.row; i++) {
      for (int j = 0; j < partB.column; j++) {
        partB.data[i * partB.column + j] = B.data[i * B.column + j];
      }
    }
  }

  else if (coordsOfCurrProc[X_AXIS] == 0) {
    MPI_Recv(partB.data, partB.row * partB.column, MPI_DOUBLE, 0, 11,
             commRow, MPI_STATUS_IGNORE);
  }

  MPI_Bcast(partB.data, dim2 * partB.column, MPI_DOUBLE, 0,
            commColumn);
  MyMatrix partC(dim1 / dimSize[X_AXIS], dim3 / dimSize[Y_AXIS]);
  multiplyMtrices(partA, partB, partC);

  MPI_Datatype cRecvType;
  MPI_Type_vector(partC.row, partC.column, dim3, MPI_DOUBLE,
                  &cRecvType);
  MPI_Type_commit(&cRecvType);

  int offset[amountOfProcs];
  for (int procRank = 0; procRank < amountOfProcs; procRank++) {
    MPI_Cart_coords(commGrid, procRank, dimOfAnyGrid,
                    coordsOfCurrProc);

    // define the location of sqaure (partC)
    // relatively from the start of full matrix C
    offset[procRank] =
        coordsOfCurrProc[X_AXIS] * partC.row * C.column +
        coordsOfCurrProc[Y_AXIS] * partC.column;
  }

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < partC.row; i++) {
      for (int j = 0; j < partC.column; j++) {
        C.data[i * C.column + j] = partC.data[i * partC.column + j];
      }
    }
    for (int i = 1; i < amountOfProcs; i++) {
      MPI_Recv(C.data + offset[i], 1, cRecvType, i, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
  } else {
    MPI_Send(partC.data, partC.row * partC.column, MPI_DOUBLE, 0, 0,
             MPI_COMM_WORLD);
  }

  MPI_Type_free(&bSendType);
  MPI_Type_free(&cRecvType);

  MPI_Comm_free(&commGrid);
  MPI_Comm_free(&commColumn);
  MPI_Comm_free(&commRow);

  double endt = MPI_Wtime();

  if (rankOfCurrProc == 0) {
    // std::cout << "C[" << C.row << " x " << C.column
    // << "]"
    //           << std::endl;
    // printMat(C);

    std::cout << "Time taken: " << endt - startt << " sec"
              << std::endl;
  }

  freeMatrix(A);
  freeMatrix(B);
  freeMatrix(C);

  freeMatrix(partA);
  freeMatrix(partB);
  freeMatrix(partC);

  MPI_Finalize();
  return 0;
}
