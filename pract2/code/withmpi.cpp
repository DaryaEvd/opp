#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mpi.h>

double *fillRandomMatrix(const size_t sizeInput) {
  double *fullMatrA = new double[sizeInput * sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    // srand(i);
    for (size_t j = 0; j < sizeInput; ++j) {
      if (i == j) {
        fullMatrA[i * sizeInput + j] = 9999;
      } else {
        fullMatrA[i * sizeInput + j] = rand() % 500 + 150.15;
        // fullMatrA[i * sizeInput + j] = 100;
      }
    }
  }
  return fullMatrA;
}

void fillConstantMatrix(double *fullMatrA, const size_t sizeInput) {
  // double *fullMatrA = new double[sizeInput * sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    for (size_t j = 0; j < sizeInput; ++j) {
      if (i == j) {
        fullMatrA[i * sizeInput + j] = 2.0;
      } else {
        fullMatrA[i * sizeInput + j] = 1.0;
      }
    }
  }
  // return fullMatrA;
}

double *fillRandomVector(const size_t sizeInput) {
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    // srand(i);
    // vector[i] = rand() % 500;
    vector[i] = 300;
  }
  return vector;
}

double *fillConstantVector(const size_t sizeInput) {
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = sizeInput + 1.0;
  }
  return vector;
}

void fillVectorU(double *vector, const size_t sizeInput) {
  // double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = sin(2 * 3.1415 * i / sizeInput);
  }
  // return vector;
}

void printMatrix(double *matrix, const size_t sizeInput) {
  std::cout << "Matrix A is: " << std::endl;
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

void multimplyMatrixOnVector(const double *matrix,
                             const double *vector, double *res,
                             const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; ++i) {
    res[i] = 0;
    for (size_t j = 0; j < sizeInput; ++j) {
      res[i] += matrix[i * sizeInput + j] * vector[j];
    }
  }
}

void parallelMultMatrixOnVector(const double *matrix,
                                const double *vector, double *res,
                                size_t matrixRowsNumber,
                                size_t matrixColumnsNumber) {
  for (size_t i = 0; i < matrixRowsNumber; i++) {
    res[i] = 0;
    for (size_t j = 0; j < matrixColumnsNumber; j++) {
      int sum_part = matrix[i * matrixColumnsNumber + j] * vector[j];
      res[i] += sum_part;
    }
  }
}

double countScalarMult(const double *vector1, const double *vector2,
                       const size_t sizeInput) {
  double res = 0;
  for (size_t i = 0; i < sizeInput; ++i) {
    res += vector1[i] * vector2[i];
  }
  return res;
}

double countVectorLength(const size_t sizeInput,
                         const double *vector) {
  return sqrt(countScalarMult(vector, vector, sizeInput));
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; ++i) {
    res[i] = scalar * vector[i];
  }
}

void substructVectors(const double *vector1, const double *vector2,
                      double *res, const size_t sizeInput) {
  for (size_t j = 0; j < sizeInput; ++j) {
    res[j] = vector1[j] - vector2[j];
  }
}

void parallelSubstrucVectors(const double *vector1,
                             const double *vector2, double *res,
                             const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    res[i] = vector1[i] - vector2[i];
  }
}

void zerofyVectors(double *vector, const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    vector[i] = 0;
  }
}

void copyVectors(const double *src, double *dst,
                 const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    dst[i] = src[i];
  }
}

int main(int argc, char **argv) {
  const int sizeInput = atoi(argv[1]);
  const double epsilon = 1e-10;

  const size_t maxIterationCounts = 100;
  size_t iterationCounts = 0;
  double tau = 0;
  // double *fullMatrA = fillRandomMatrix(sizeInput);
  // double *fullMatrA = fillConstantMatrix(sizeInput);
  // printMatrix(fullMatrA, sizeInput);
  // double *b = fillConstantVector(sizeInput);
  // double *b = fillRandomVector(sizeInput);
  // std::cout << "vector b is: " << std::endl;
  // printVector(b, sizeInput);
  // zerofyVectors(xCurr, sizeInput);
  // std::cout << "vector X is: " << std::endl;
  // printVector(xCurr, sizeInput);
  // std::cout << "matrix A is" << std::endl;
  // printMatrix(fullMatrA, sizeInput);
  // std::cout << "vector b is" << std::endl;
  // printVector(b, sizeInput);

  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();

  int amountOfProcs, rankOfCurrProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);

  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);
  std::cout << "rankOfCurrProc: " << rankOfCurrProc << std::endl;

  if (rankOfCurrProc == 0) {
    std::cout << "amountOfProcs: " << amountOfProcs << std::endl;
  }
  double *fullMatrA = new double[sizeInput * sizeInput];
  double *u = new double[sizeInput];
  double *b = new double[sizeInput];
  double *xCurr = new double[sizeInput];
  zerofyVectors(xCurr, sizeInput);

  if (rankOfCurrProc == 0) {

    fillConstantMatrix(fullMatrA, sizeInput);
    printMatrix(fullMatrA, sizeInput);

    fillVectorU(u, sizeInput);
    std::cout << "vector u is" << std::endl;
    printVector(u, sizeInput);

    multimplyMatrixOnVector(fullMatrA, u, b, sizeInput);
    std::cout << "vector b is" << std::endl;
    printVector(b, sizeInput);
  }

  double bNorm = countVectorLength(sizeInput, b);
  // std::cout << "b norm is: " << bNorm << ",  ";

  /* data for counting rows amount*/
  int basicRowsCount = sizeInput / amountOfProcs;
  int restRowsCount = sizeInput % amountOfProcs;

  int *rowsNum = new int[amountOfProcs];
  int *sendCounts = new int[amountOfProcs];
  int *recvCounts = new int[amountOfProcs];

  for (int i = 0; i < amountOfProcs; i++) {
    rowsNum[i] = basicRowsCount +
                 (i < restRowsCount); // how many rows in each proc
    sendCounts[i] =
        rowsNum[i] * sizeInput; // how many elems in each proc

    recvCounts[i] = basicRowsCount + (i < restRowsCount);
  }

  int *displs = new int[amountOfProcs];
  displs[0] = 0;
  for (int i = 1; i < amountOfProcs; i++) {
    displs[i] = displs[i - 1] + sendCounts[i - 1];
  }

  int *offset = new int[amountOfProcs];
  offset[0] = 0;
  for (int i = 1; i < amountOfProcs; i++) {
    recvCounts[i] = basicRowsCount + (i < restRowsCount);
    offset[i] = offset[i - 1] + recvCounts[i - 1];
  }

  double *partMatrA = new double[sendCounts[rankOfCurrProc]];
  MPI_Scatterv(fullMatrA, sendCounts, displs, MPI_DOUBLE, partMatrA,
               sendCounts[rankOfCurrProc], MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
  // ура, разослали каждому процессу элементы (их разное количество в
  // процессах)

  MPI_Bcast(b, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // MPI_Bcast(xCurr, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  int widthPartA = sizeInput;
  int heightPartA = rowsNum[rankOfCurrProc];

  double *partMatrA_X = new double[rowsNum[rankOfCurrProc]];
  double *yPart = new double[rowsNum[rankOfCurrProc]];
  zerofyVectors(yPart, rowsNum[rankOfCurrProc]);

  double *yFull = new double[sizeInput];
  zerofyVectors(yFull, sizeInput);

  double *partAY = new double[rowsNum[rankOfCurrProc]];
  double *tauY = new double[rowsNum[rankOfCurrProc]];

  double *partXNext = new double[rowsNum[rankOfCurrProc]];

  double *xNext = new double[sizeInput];
  zerofyVectors(xNext, sizeInput);

  double *tmpAMulX = new double[sizeInput];

  bool isContinue = true;

  while (1) {
    MPI_Bcast(xCurr, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // std::cout << "After B cast " << std::endl;

    parallelMultMatrixOnVector(partMatrA, xCurr, partMatrA_X,
                               rowsNum[rankOfCurrProc], sizeInput);
    // std::cout << "After mult matrix on vector" << std::endl;

    parallelSubstrucVectors(partMatrA_X, b, yPart,
                            rowsNum[rankOfCurrProc]);
    // std::cout << "After substru matrix on vector" << std::endl;

    MPI_Allgatherv(yPart, rowsNum[rankOfCurrProc], MPI_DOUBLE, yFull,
                   recvCounts, offset, MPI_DOUBLE, MPI_COMM_WORLD);
    std::cout << "y full is: ";
    printVector(yFull, sizeInput);

    double yNorm = countVectorLength(sizeInput, yFull);
    if (rankOfCurrProc == 0) {
      std::cout << "Y norm is: " << yNorm << ",  ";
    }

    parallelMultMatrixOnVector(partMatrA, yFull, partAY,
                               rowsNum[rankOfCurrProc], sizeInput);

    double partScalarY_Ay =
        countScalarMult(yPart, partAY, rowsNum[rankOfCurrProc]);
    double numerator;
    MPI_Reduce(&partScalarY_Ay, &numerator, 1, MPI_DOUBLE, MPI_SUM, 0,
               MPI_COMM_WORLD);
    // std::cout << "numerator tau is: " << numerator << std::endl;

    double partScalarAy_Ay =
        countScalarMult(partAY, partAY, rowsNum[rankOfCurrProc]);
    double denumerator;
    MPI_Reduce(&partScalarAy_Ay, &denumerator, 1, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    // std::cout << "denumerator tau is: " << denumerator <<
    // std::endl;

    if (rankOfCurrProc == 0) {
      tau = numerator / denumerator;
    }

    MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    countVectorMultNumber(yPart, tau, tauY, rowsNum[rankOfCurrProc]);

    parallelSubstrucVectors(xCurr, tauY, partXNext,
                            rowsNum[rankOfCurrProc]);

    MPI_Allgatherv(partXNext, rowsNum[rankOfCurrProc], MPI_DOUBLE,
                   xNext, recvCounts, offset, MPI_DOUBLE,
                   MPI_COMM_WORLD);

    double endCriteria = yNorm / bNorm;
    if (rankOfCurrProc == 0) {
      std::cout << "endcriteria is: " << std::fixed << endCriteria
                << " ";
    }


    if (iterationCounts > maxIterationCounts) {
      std::cout << "Too many iterations. Change init values"
                << std::endl;
      delete[] partXNext;
      delete[] partAY;
      delete[] tmpAMulX;
      delete[] yPart;
      delete[] yFull;
      delete[] partMatrA;
      delete[] displs;
      delete[] rowsNum;
      delete[] sendCounts;

      delete[] xCurr;
      delete[] xNext;
      delete[] partMatrA_X;
      delete[] tauY;
      delete[] fullMatrA;
      delete[] b;

      MPI_Finalize();
      return 0;
    }

    if (rankOfCurrProc == 0) {
      if (endCriteria < epsilon) {
        std::cout << "endCriteria < epsilon !!!!!!!!!!!!!!!!!!!!"
                  << std::endl;
        isContinue = false;
        // MPI_Bcast(&isContinue, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        // break;
      }
    }
    MPI_Bcast(&isContinue, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    
    if (isContinue == false) {
      break;
    }

    copyVectors(xNext, xCurr, sizeInput);

    iterationCounts++;

    if (rankOfCurrProc == 0) {
      std::cout << "iteration " << iterationCounts
                << " ended ====" << std::endl;
    }
  }

  double endTime = MPI_Wtime();

  if (rankOfCurrProc == 0) {
    std::cout << "Iteration amount in total: " << iterationCounts
              << "\n";

    std::cout << "Time taken: " << endTime - startTime << " sec";
    std::cout << std::endl;

    std::cout << "Solution (vector x is): " << std::endl;
    printVector(xNext, sizeInput);
  }

  delete[] partXNext;
  delete[] partAY;
  delete[] tmpAMulX;
  delete[] yPart;
  delete[] yFull;
  delete[] partMatrA;
  delete[] displs;
  delete[] rowsNum;
  delete[] sendCounts;

  delete[] xCurr;
  delete[] xNext;
  delete[] partMatrA_X;
  delete[] tauY;
  delete[] fullMatrA;
  delete[] b;

  MPI_Finalize();
  return 0;
}
