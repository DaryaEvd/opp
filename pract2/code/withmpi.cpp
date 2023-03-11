#include <cmath>
#include <iostream>
#include <mpi.h>

void printMatrix(double *matrix, const size_t sizeInput) {
  std::cout << "Full matrix A is: " << std::endl;
  for (size_t i = 0; i < sizeInput; ++i) {
    for (size_t j = 0; j < sizeInput; ++j) {
      std::cout << matrix[j + i * sizeInput] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void printVector(const double *vector, const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; ++i) {
    std::cout << std::fixed << vector[i] << " ";
  }
  std::cout << std::endl;
}

void fillRandomMatrix(double *fullMatrA, size_t sizeInput) {
  srand(0);

  for (size_t i = 0; i < sizeInput; i++) {
    for (size_t j = 0; j < sizeInput; j++) {
      fullMatrA[i * sizeInput + j] =
          (i == j) ? 9999 : rand() % 500 + 150.15;
    }
  }
}

void fillConstantMatrix(double *fullMatrA, size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    for (size_t j = 0; j < sizeInput; j++) {
      fullMatrA[i * sizeInput + j] = (i == j) ? 2.0 : 1.0;
    }
  }
}

void fillRandomVector(double *vector, size_t sizeInput) {
  srand(0);

  for (size_t i = 0; i < sizeInput; i++) {
    vector[i] = rand() % 500;
  }
}

void fillConstantVector(double *vector, size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    vector[i] = sizeInput + 1;
  }
}

void multimplyMatrixOnVector(const double *matrix,
                             const double *vector, double *res,
                             size_t matrixRowsNumber,
                             size_t matrixColumnsNumber) {
  for (size_t i = 0; i < matrixRowsNumber; i++) {
    res[i] = 0;
    for (size_t j = 0; j < matrixColumnsNumber; j++) {
      res[i] += matrix[i * matrixColumnsNumber + j] * vector[j];
    }
  }
}

double *fillVectorU(const int N) {
  double *u = new double[N];
  for (int i = 0; i < N; i++) {
    u[i] = sin(2 * M_PI * (i) / N);
  }
  return u;
}

void fillVectorB(double *b, const double *fullMatrixA,
                 size_t sizeInput) {
  double *u = fillVectorU((int)sizeInput);
  std::cout << "vector u is: " << std::endl;
  printVector(u, sizeInput);

  multimplyMatrixOnVector(fullMatrixA, u, b, sizeInput, sizeInput);
  std::cout << "vector b is: " << std::endl;
  printVector(b, sizeInput);
  delete[] u;
}

double countScalarMult(const double *vector1, const double *vector2,
                       const size_t size) {

  double res = 0;
  for (size_t i = 0; i < size; i++) {
    res += vector1[i] * vector2[i];
  }

  return res;
}

double countVectorLength(const double *vector, const size_t size) {
  return sqrt(countScalarMult(vector, vector, size));
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    res[i] = scalar * vector[i];
  }
}

void subtractVectors(const double *vector1, const double *vector2,
                     double *res, const size_t firstVectorSize,
                     const int firstVectorElementsOffset,
                     const int secondVectorElementsOffset) {
  for (size_t i = 0; i < firstVectorSize; i++) {
    res[i] = vector1[firstVectorElementsOffset + i] -
             vector2[secondVectorElementsOffset + i];
  }
}

void copyVectors(const double *src, double *dst, const size_t size) {
  for (size_t i = 0; i < size; i++) {
    dst[i] = src[i];
  }
}

int *countElemsNumInEachProc(size_t amountOfProcs, size_t sizeInput) {
  int *elemsNum = new int[amountOfProcs];

  int basicRowsCount = sizeInput / amountOfProcs;
  int restRowsCount = sizeInput % amountOfProcs;

  for (size_t i = 0; i < amountOfProcs; i++) {
    elemsNum[i] = basicRowsCount * sizeInput;
    if (restRowsCount > 0) {
      elemsNum[i] += sizeInput;
      --restRowsCount;
    }
  }
  return elemsNum;
}

int *countRowsInEachProcess(const int *elementsNumberArr,
                            int amountOfProcs, size_t sizeInput) {
  int *rowsNum = new int[amountOfProcs];
  for (size_t i = 0; i < amountOfProcs; i++) {
    rowsNum[i] = elementsNumberArr[i] / sizeInput;
  }
  return rowsNum;
}

int *createElemsOffsetArr(const int *elemsNum, int amountOfProcs) {
  int *elementsOffsetArray = new int[amountOfProcs];
  int elementsOffset = 0;
  for (size_t i = 0; i < amountOfProcs; i++) {
    elementsOffsetArray[i] = elementsOffset;
    elementsOffset += elemsNum[i];
  }
  return elementsOffsetArray;
}

int *createRowsOffsetArr(const int *elementsOffsetArray,
                         int amountOfProcs, size_t sizeInput) {
  int *rowsOffsetArray = new int[amountOfProcs];
  for (size_t i = 0; i < amountOfProcs; i++) {
    rowsOffsetArray[i] = elementsOffsetArray[i] / sizeInput;
  }

  return rowsOffsetArray;
}

int main(int argc, char **argv) {

  const size_t sizeInput = atoi(argv[1]);
  const double epsilon = 10e-7;

  double critCurrentEnd = 1;
  double prevCritEnd = 1;
  double prevPrevCritEnd = 1;

  double tau;
  bool isEndOfAlgo = false;
  size_t iterationCounts = 0;
  const size_t maxIterationCounts = 50000;
  int convergenceCount = 0;
  const size_t maxConvergenceCount = 5;

  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();

  int amountOfProcs, rankOfCurrProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  double *fullMatrA = new double[sizeInput * sizeInput];
  double *b = new double[sizeInput];
  double *xCurr = new double[sizeInput];

  // if (rankOfCurrProc == 0) {
  //   fillConstantMatrix(fullMatrA, sizeInput);
  //   // printMatrix(fullMatrA, sizeInput);

  //   fillVectorB(b, fullMatrA, sizeInput);
  //   // std::cout << "vector b is: " << std::endl;
  //   // printVector(b, sizeInput);

  //   std::fill(xCurr, xCurr + sizeInput, 0);
  //   // std::cout << "vector X at start is: " << std::endl;
  //   // printVector(xCurr, sizeInput);
  // }

  if (rankOfCurrProc == 0) {
    fillRandomMatrix(fullMatrA, sizeInput);
    // printMatrix(fullMatrA, sizeInput);

    fillRandomVector(b, sizeInput);
    // std::cout << "vector b is: " << std::endl;
    // printVector(b, sizeInput);

    std::fill(xCurr, xCurr + sizeInput, 0);
    // std::cout << "vector X at start is: " << std::endl;
    // printVector(xCurr, sizeInput);
  }

  double bNorm;
  if (rankOfCurrProc == 0) {
    bNorm = sqrt(countScalarMult(b, b, sizeInput));
  }

  MPI_Bcast(b, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // amount of elems, handling each process
  int *elemsNum = countElemsNumInEachProc(amountOfProcs, sizeInput);

  // amount of rows, which each process handles
  int *rowsNum =
      countRowsInEachProcess(elemsNum, amountOfProcs, sizeInput);

  // count, from which element data sends to each process
  int *elementsOffsetArray =
      createElemsOffsetArr(elemsNum, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray = createRowsOffsetArr(
      elementsOffsetArray, amountOfProcs, sizeInput);

  double *partMatrA = new double[elemsNum[rankOfCurrProc]];
  MPI_Scatterv(fullMatrA, elemsNum, elementsOffsetArray, MPI_DOUBLE,
               partMatrA, elemsNum[rankOfCurrProc], MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

  double *partMultResultVector_Ax =
      new double[rowsNum[rankOfCurrProc]];
  double *partMultResultVector_Ay =
      new double[rowsNum[rankOfCurrProc]];
  double *partVectorY = new double[rowsNum[rankOfCurrProc]];
  double *partNextX = new double[rowsNum[rankOfCurrProc]];
  double *partMultVectorByScalar_TauY =
      new double[rowsNum[rankOfCurrProc]];
  double *vectorY = new double[sizeInput];
  double *xNext = new double[sizeInput];

  while (1) {
    MPI_Bcast(xCurr, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    multimplyMatrixOnVector(partMatrA, xCurr, partMultResultVector_Ax,
                            rowsNum[rankOfCurrProc], sizeInput);

    subtractVectors(partMultResultVector_Ax, b, partVectorY,
                    rowsNum[rankOfCurrProc], 0,
                    rowsOffsetArray[rankOfCurrProc]);

    // y_n = A * x_n - b
    MPI_Allgatherv(partVectorY, rowsNum[rankOfCurrProc], MPI_DOUBLE,
                   vectorY, rowsNum, rowsOffsetArray, MPI_DOUBLE,
                   MPI_COMM_WORLD);

    // A * y_n
    multimplyMatrixOnVector(partMatrA, vectorY,
                            partMultResultVector_Ay,
                            rowsNum[rankOfCurrProc], sizeInput);

    // || A * x_n - b ||
    double yNorm = 0;
    if (rankOfCurrProc == 0) {
      double squareMetricVectorY =
          countScalarMult(vectorY, vectorY, sizeInput);
      yNorm = sqrt(squareMetricVectorY);
    }

    // start counting tau
    double numeratorTau, denominatorTau;
    double partScalarProduct_YAy =
        countScalarMult(partVectorY, partMultResultVector_Ay,
                        rowsNum[rankOfCurrProc]);
    // (y_n, A * y_n)
    MPI_Reduce(&partScalarProduct_YAy, &numeratorTau, 1, MPI_DOUBLE,
               MPI_SUM, 0, MPI_COMM_WORLD);

    double partScalarProduct_AyAy = countScalarMult(
        partMultResultVector_Ay, partMultResultVector_Ay,
        rowsNum[rankOfCurrProc]);
    // (A * y_n, A * y_n)
    MPI_Reduce(&partScalarProduct_AyAy, &denominatorTau, 1,
               MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rankOfCurrProc == 0) {
      tau = numeratorTau / denominatorTau;
    }
    MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    countVectorMultNumber(partVectorY, tau,
                          partMultVectorByScalar_TauY,
                          rowsNum[rankOfCurrProc]);

    subtractVectors(xCurr, partMultVectorByScalar_TauY, partNextX,
                    rowsNum[rankOfCurrProc],
                    rowsOffsetArray[rankOfCurrProc], 0);

    // x_n+1 = x_n - tau * y
    MPI_Allgatherv(partNextX, rowsNum[rankOfCurrProc], MPI_DOUBLE,
                   xNext, rowsNum, rowsOffsetArray, MPI_DOUBLE,
                   MPI_COMM_WORLD);

    critCurrentEnd = yNorm / bNorm;
    if (rankOfCurrProc == 0) {
      if (critCurrentEnd < epsilon && critCurrentEnd < prevCritEnd &&
          critCurrentEnd < prevPrevCritEnd) {
        isEndOfAlgo = true;
      }
      if (iterationCounts > maxIterationCounts) {
        isEndOfAlgo = true;
      }
    }
    MPI_Bcast(&isEndOfAlgo, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);

    if (isEndOfAlgo) {
      break;
    }

    copyVectors(xNext, xCurr, sizeInput);

    prevPrevCritEnd = prevCritEnd;
    prevCritEnd = critCurrentEnd;

    if (rankOfCurrProc == 0) {
      iterationCounts++;
    }

    MPI_Bcast(&iterationCounts, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rankOfCurrProc == 0) {
      std::cout << "iteration " << iterationCounts
                << " ended ====" << std::endl;
    }
  }

  double endTime = MPI_Wtime();
  if (rankOfCurrProc == 0) {
    if (iterationCounts < maxIterationCounts) {
      std::cout << "Iteration amount in total: " << iterationCounts
                << std::endl;
      std::cout << "Time taken: " << endTime - startTime << " sec"
                << std::endl;
      // std::cout << "Solution (vector X is): " << std::endl;
      // printVector(xCurr, sizeInput);
    }

    else {
      std::cout << "There are no solutions =( Change input \n";
    }

    delete[] fullMatrA;
  }

  delete[] partMatrA;
  delete[] partMultVectorByScalar_TauY;
  delete[] partVectorY;
  delete[] partNextX;
  delete[] partMultResultVector_Ax;
  delete[] partMultResultVector_Ay;

  delete[] b;
  delete[] xCurr;
  delete[] xNext;
  delete[] vectorY;
  delete[] rowsOffsetArray;
  delete[] elementsOffsetArray;
  delete[] rowsNum;
  delete[] elemsNum;

  MPI_Finalize();

  return 0;
}
