#include <cmath>
#include <cstdlib>
#include <iostream>

void printMatrix(double *matrix, const size_t sizeInput) {
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

double *fillRandomMatrix(const size_t sizeInput) {
  double *matrixA = new double[sizeInput * sizeInput];
  srand(0);
  for (size_t i = 0; i < sizeInput; ++i) {
    for (size_t j = 0; j < sizeInput; ++j) {
      if (i == j) {
        matrixA[i * sizeInput + j] = 9999;
      } else {
        matrixA[i * sizeInput + j] = rand() % 500 + 150.15;
      }
    }
  }
  return matrixA;
}

double *fillConstantMatrix(const size_t sizeInput) {
  double *matrixA = new double[sizeInput * sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    for (size_t j = 0; j < sizeInput; ++j) {
      if (i == j) {
        matrixA[i * sizeInput + j] = 2.0;
      } else {
        matrixA[i * sizeInput + j] = 1.0;
      }
    }
  }
  return matrixA;
}

double *fillRandomVector(const size_t sizeInput) {
  srand(0);
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = rand() % 500;
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

double *fillVectorU(const size_t sizeInput) {
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = sin(2 * 3.1415 * i / sizeInput);
  }
  return vector;
}

void fillVectorB(double *b, const double *fullMatrixA,
                 size_t sizeInput) {
  double *vectorU = fillVectorU((int)sizeInput);
  // std::cout << "vector U is: " << std::endl;
  // printVector(vectorU, sizeInput);
  multimplyMatrixOnVector(fullMatrixA, vectorU, b, sizeInput);

  delete[] vectorU;
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

void copyVector(const double *src, double *dst,
                const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    dst[i] = src[i];
  }
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Bad input! Enter matrix size" << std::endl;
  }
  const size_t sizeInput = atoi(argv[1]);
  const double accuracy = 1e-10;
  const double epsilon = 1e-10;

  const size_t maxIterationCounts = 50000;
  const size_t maxConvergenceCount = 5;
  size_t iterationCounts = 0;
  int convergenceCount = 0;

  double critCurrentEnd = 1;
  double prevCritEnd = 1;
  double prevPrevCritEnd = 1;

  bool isEndOfAlgo = false;

  struct timespec endt, startt;
  clock_gettime(CLOCK_MONOTONIC_RAW, &startt);

  double *Atmp = new double[sizeInput];
  double *y = new double[sizeInput];
  double *tauY = new double[sizeInput];

  double *xCurr = new double[sizeInput];
  double *xNext = new double[sizeInput];

  ///* for testing data with RANDOM values ========= */
  double *matrixA = fillRandomMatrix(sizeInput);
  // printMatrix(matrixA, sizeInput);
  double *b = fillRandomVector(sizeInput);
  // std::cout << "vector b is: " << std::endl;
  // printVector(b, sizeInput);
  std::fill(xCurr, xCurr + sizeInput, 0);
  // std::cout << "vector X at start is: " << std::endl;
  // printVector(xCurr, sizeInput);
  ///* for testing data with RANDOM values ========= */

  ///* for testing when vector b uses sin ========= */
  // double *matrixA = fillConstantMatrix(sizeInput);
  // // printMatrix(matrixA, sizeInput);

  // double *b = new double[sizeInput];
  // fillVectorB(b, matrixA, sizeInput);
  // // std::cout << "vector b is" << std::endl;
  // // printVector(b, sizeInput);

  // std::fill(xCurr, xCurr + sizeInput, 0);
  // // std::cout << "vector X at start is: " << std::endl;
  // // printVector(xCurr, sizeInput);
  ///* for testing when vector b uses sin ========= */

  double bNorm = countVectorLength(sizeInput, b);

  while (1) {
    multimplyMatrixOnVector(matrixA, xCurr, Atmp,
                            sizeInput);      // A * x_n
    substructVectors(Atmp, b, y, sizeInput); // y_n = A * x_n - b
    double yNorm =
        countVectorLength(sizeInput, y); // || A * x_n - b ||

    std::fill(Atmp, Atmp + sizeInput, 0);
    multimplyMatrixOnVector(matrixA, y, Atmp, sizeInput); // A * y_n
    double numeratorTau =
        countScalarMult(y, Atmp, sizeInput); // (y_n, A * y_n)
    double denominatorTau =
        countScalarMult(Atmp, Atmp, sizeInput); // (A * y_n, A * y_n)

    double tau = numeratorTau / denominatorTau;

    countVectorMultNumber(y, tau, tauY, sizeInput); // tau * y
    substructVectors(xCurr, tauY, xNext,
                     sizeInput); // x_n+1 = x_n - tau * y

    double critCurrentEnd = yNorm / bNorm;

    if (iterationCounts > maxIterationCounts) {
      std::cout << "Too many iterations. Change init values"
                << std::endl;
      delete[] xCurr;
      delete[] xNext;
      delete[] Atmp;
      delete[] y;
      delete[] tauY;
      delete[] matrixA;
      delete[] b;

      return 0;
    }
    if (critCurrentEnd < epsilon && critCurrentEnd < prevCritEnd &&
        critCurrentEnd < prevPrevCritEnd) {

      isEndOfAlgo = true;
      break;
    }

    copyVector(xNext, xCurr, sizeInput);

    prevPrevCritEnd = prevCritEnd;
    prevCritEnd = critCurrentEnd;

    iterationCounts++;

    std::cout << "iteration " << iterationCounts
              << " ended ===" << std::endl;
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &endt);

  if (isEndOfAlgo) {
    std::cout << "Iteration amount in total: " << iterationCounts
              << std::endl;
    std::cout << "Time taken: "
              << endt.tv_sec - startt.tv_sec +
                     accuracy * (endt.tv_nsec - startt.tv_nsec)
              << " sec" << std::endl;

    // std::cout << "Solution (vector X is): " << std::endl;
    // printVector(xNext, sizeInput);

  } else {
    std::cout << "There are no solutions =( Change input \n";
  }

  delete[] xCurr;
  delete[] xNext;
  delete[] Atmp;
  delete[] y;
  delete[] tauY;
  delete[] matrixA;
  delete[] b;

  return 0;
}
