#include <cmath>
#include <cstdlib>
#include <iostream>

double *fillRandomMatrix(const size_t size) {
  double *matrixA = new double[size * size];
  for (size_t i = 0; i < size; ++i) {
    // srand(i);
    for (size_t j = 0; j < size; ++j) {
      if (i == j) {
        matrixA[i * size + j] = 9999;
      } else {
        matrixA[i * size + j] = rand() % 500 + 150.15;
        // matrixA[i * size + j] = 100;
      }
    }
  }
  return matrixA;
}

double *fillConstantMatrix(const size_t size) {
  double *matrixA = new double[size * size];
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      if (i == j) {
        matrixA[i * size + j] = 2.0;
      } else {
        matrixA[i * size + j] = 1.0;
      }
    }
  }
  return matrixA;
}

double *fillRandomVector(const size_t size) {
  double *vector = new double[size];
  for (size_t i = 0; i < size; ++i) {
    // srand(i);
    // vector[i] = rand() % 500;
    vector[i] = 300;
  }
  return vector;
}

double *fillConstantVector(const size_t size) {
  double *vector = new double[size];
  for (size_t i = 0; i < size; ++i) {
    vector[i] = size + 1.0;
  }
  return vector;
}

double *fillVectorU(const size_t size) {
  double *vector = new double[size];
  for (size_t i = 0; i < size; ++i) {
    vector[i] = sin(2 * 3.1415 * i / size);
  }
  return vector;
}

void printMatrix(double *matrix, const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    for (size_t j = 0; j < size; ++j) {
      std::cout << matrix[j + i * size] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void printVector(const double *vector, const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    std::cout << std::fixed << vector[i] << " ";
  }
  std::cout << std::endl;
}

void multimplyMatrixOnVector(const double *matrix,
                             const double *vector, double *res,
                             const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    res[i] = 0;
    for (size_t j = 0; j < size; ++j) {
      res[i] += matrix[i * size + j] * vector[j];
    }
  }
}

double countScalarMult(const double *vector1, const double *vector2,
                       const size_t size) {
  double res = 0;
  for (size_t i = 0; i < size; ++i) {
    res += vector1[i] * vector2[i];
  }
  return res;
}

double countVectorLength(const size_t size, const double *vector) {
  return sqrt(countScalarMult(vector, vector, size));
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    res[i] = scalar * vector[i];
  }
}

void substructVectors(const double *vector1, const double *vector2,
                      double *res, const size_t size) {
  for (size_t j = 0; j < size; ++j) {
    res[j] = vector1[j] - vector2[j];
  }
}

void zerofyVectors(double *vector, const size_t size) {
  for (size_t i = 0; i < size; i++) {
    vector[i] = 0;
  }
}

void copyVectors(const double *src, double *dst, const size_t size) {
  for (size_t i = 0; i < size; i++) {
    dst[i] = src[i];
  }
}

int main(int argc, char *argv[]) {
  const size_t size = atoi(argv[1]);
  const double accuracy = 1e-10;
  const double epsilon = 1e-10;

  const size_t maxIterationCounts = 50000;
  const size_t maxConvergenceCount = 5;
  size_t iterationCounts = 0;
  int convergenceCount = 0;

  double *Atmp = new double[size];
  double *y = new double[size];
  double *tauY = new double[size];

  double *xCurr = new double[size];
  double *xNext = new double[size];
  zerofyVectors(xNext, size);

  double prevEndCrit = 0, endCrit = 0;
  struct timespec endt, startt;

  double *matrixA = fillRandomMatrix(size);
  // double *matrixA = fillConstantMatrix(size);
  // printMatrix(matrixA, size);

  // double *b = fillConstantVector(size);
  double *b = fillRandomVector(size);
  // std::cout << "vector b is: " << std::endl;
  // printVector(b, size);

  zerofyVectors(xCurr, size);
  // std::cout << "vector X is: " << std::endl;
  // printVector(xCurr, size);

  // std::cout << "matrix A is" << std::endl;
  // printMatrix(matrixA, size);
  // std::cout << "vector b is" << std::endl;
  // printVector(b, size);




  // double *matrixA = fillConstantMatrix(size);
  // printMatrix(matrixA, size);
  // double *u = fillVectorU(size);
  // std::cout << "vector u is" << std::endl;
  // printVector(u, size);
  // double *b = new double[size];
  // multimplyMatrixOnVector(matrixA, u, b, size);
  // std::cout << "vector b is" << std::endl;
  // printVector(b, size);

  clock_gettime(CLOCK_MONOTONIC_RAW, &startt);
  double bNorm = countVectorLength(size, b);
  // std::cout << "b norm is: " << bNorm << ",  ";

  while (1) {
    multimplyMatrixOnVector(matrixA, xCurr, Atmp, size); // A * x_n
    substructVectors(Atmp, b, y, size);        // y_n = A * x_n - b
    double yNorm = countVectorLength(size, y); // || A * x_n - b ||
    std::cout << "Y norm is: " << yNorm << ",  ";
    zerofyVectors(Atmp, size);
    multimplyMatrixOnVector(matrixA, y, Atmp, size); // A * y_n
    double numerator =
        countScalarMult(y, Atmp, size); // (y_n, A * y_n)
    double denumerator =
        countScalarMult(Atmp, Atmp, size); // (A * y_n, A * y_n)

    double tau = numerator / denumerator;

    countVectorMultNumber(y, tau, tauY, size); // tau * y
    substructVectors(xCurr, tauY, xNext,
                     size); // x_n+1 = x_n - tau * y

    double endCriteria = yNorm / bNorm;
    std::cout << "endcrit is: " << std::fixed << endCriteria << " ";
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
    if (endCriteria < epsilon) {
      std::cout << "endCriteria < epsilon !!!!!!!!!!!!!!!!!!!!"
                << std::endl;
      break;
    }

    // if (convergenceCount > maxConvergenceCount) {
    //   std::cout << "convergenceCount > maxConvergenceCount "
    //                "!!!!!!!!!!!!!!!!!!!!"
    //             << std::endl;
    //   break;
    // } else if (endCrit < prevEndCrit) {
    //   std::cout << "I'm in endCrit < prevEndCrit" << std::endl;
    //   convergenceCount++;
    // } else if (endCrit > prevEndCrit) {
    //   std::cout << "I'm in endCrit > prevEndCrit" << std::endl;
    //   convergenceCount = 0;
    // }

    copyVectors(xNext, xCurr, size);
    prevEndCrit = endCrit;
    iterationCounts++;

    std::cout << "iteration " << iterationCounts
              << " ended ===" << std::endl;
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &endt);

  if (convergenceCount <= maxConvergenceCount) {
    std::cout << "Iteration amount in total: " << iterationCounts
              << "\n";

  } else {
    std::cout << "There are no solutions!\n";
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &endt);

  std::cout << "Time taken: "
            << endt.tv_sec - startt.tv_sec +
                   accuracy * (endt.tv_nsec - startt.tv_nsec)
            << " sec";
  std::cout << std::endl;

  // std::cout << "Solution (vector x is): " << std::endl;
  // printVector(xNext, size);

  delete[] xCurr;
  delete[] xNext;
  delete[] Atmp;
  delete[] y;
  delete[] tauY;
  delete[] matrixA;
  delete[] b;

  return 0;
}
