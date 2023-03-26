#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>

/* for printinf type of scheduling*/
#define xstr(x) str(x)
#define str(x) #x

// #define TYPE auto
#define TYPE static
// #define TYPE dynamic
// #define TYPE guided
// #define TYPE runtime

#define CHUNK_SIZE 200

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

double *fillRandomVector(const size_t sizeInput) {
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = rand() % 500;
  }
  return vector;
}

void multimplyMatrixOnVector(const double *matrix,
                             const double *vector, double *res,
                             const size_t sizeInput) {
#pragma omp for schedule(TYPE, CHUNK_SIZE)
  for (int i = 0; i < sizeInput; i++) {
    res[i] = 0;
  }
#pragma omp for collapse(2) schedule(TYPE, CHUNK_SIZE)
  for (size_t i = 0; i < sizeInput; ++i) {
    for (size_t j = 0; j < sizeInput; ++j) {
#pragma omp atomic

      res[i] += matrix[i * sizeInput + j] * vector[j];
    }
  }
}

void countScalarMult(const double *vector1, const double *vector2,
                     const size_t sizeInput, double *res) {
#pragma omp single
  { *res = 0; }
#pragma omp barrier

#pragma omp for reduction(+ : res[0]) schedule(TYPE, CHUNK_SIZE)
  for (size_t i = 0; i < sizeInput; ++i) {
    *res += vector1[i] * vector2[i];
  }
}

void countVectorLength(const size_t sizeInput, const double *vector,
                       double *res) {
#pragma omp single
  { *res = 0; }
#pragma omp barrier

#pragma omp for reduction(+ : res[0]) schedule(TYPE, CHUNK_SIZE)
  for (size_t i = 0; i < sizeInput; ++i) {
    *res += vector[i] * vector[i];
  }
#pragma omp single
  { *res = sqrt(*res); }
#pragma omp barrier
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t sizeInput) {
#pragma omp for
  for (size_t i = 0; i < sizeInput; ++i) {
    res[i] = scalar * vector[i];
  }
}

void substructVectors(const double *vector1, const double *vector2,
                      double *res, const size_t sizeInput) {
#pragma omp for schedule(TYPE, CHUNK_SIZE)
  for (size_t j = 0; j < sizeInput; ++j) {
    res[j] = vector1[j] - vector2[j];
  }
}

void copyVector(const double *src, double *dst,
                const size_t sizeInput) {
#pragma omp for schedule(TYPE, CHUNK_SIZE)
  for (size_t i = 0; i < sizeInput; i++) {
    dst[i] = src[i];
  }
}

int main(int argc, char *argv[]) {
  // In this code the 1st arg is sizeMatrix
  // the 2nd arg is amount of threads

  if (argc != 3) {
    std::cout << "Bad input! Enter matrix size" << std::endl;
  }
  srand(0);

  const size_t sizeInput = atoi(argv[1]);
  const int amountOfThreads = atoi(argv[2]);

  const double precision = 1e-10;
  const double epsilon = 1e-10;

  const size_t maxIterationCounts = 50000;
  size_t iterationCounts = 0;

  double critCurrentEnd = 1;
  double prevCritEnd = 1;
  double prevPrevCritEnd = 1;

  bool isEndOfAlgo = false;

  double startt = omp_get_wtime();

  double *Atmp = new double[sizeInput];
  double *y = new double[sizeInput];
  double *tauY = new double[sizeInput];

  double *xCurr = new double[sizeInput];
  double *xNext = new double[sizeInput];

  ///* for testing data with RANDOM values, starts ========= */
  double *matrixA = fillRandomMatrix(sizeInput);
  // printMatrix(matrixA, sizeInput);
  double *b = fillRandomVector(sizeInput);
  // std::cout << "vector b is: " << std::endl;
  // printVector(b, sizeInput);
  std::fill(xCurr, xCurr + sizeInput, 0);
  // std::cout << "vector X at start is: " << std::endl;
  // printVector(xCurr, sizeInput);
  ///* for testing data with RANDOM values, ended ========= */

  double bNorm, yNorm;
  countVectorLength(sizeInput, b, &bNorm);

  double numeratorTau, denominatorTau;

#pragma omp parallel num_threads(amountOfThreads)
  {
    while (1) {
      multimplyMatrixOnVector(matrixA, xCurr, Atmp,
                              sizeInput); // A * x_n

      substructVectors(Atmp, b, y, sizeInput); // y_n = A * x_n - b

      countVectorLength(sizeInput, y,
                        &yNorm); // || A * x_n - b ||

      multimplyMatrixOnVector(matrixA, y, Atmp, sizeInput); // A * y_n

      countScalarMult(y, Atmp, sizeInput,
                      &numeratorTau); // (y_n, A * y_n)
      countScalarMult(Atmp, Atmp, sizeInput,
                      &denominatorTau); // (A * y_n, A * y_n)

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

        abort();
      }

      if (critCurrentEnd < epsilon && critCurrentEnd < prevCritEnd &&
          critCurrentEnd < prevPrevCritEnd) {

        isEndOfAlgo = true;
        break;
      }

      copyVector(xNext, xCurr, sizeInput);

#pragma omp barrier

      prevPrevCritEnd = prevCritEnd;
      prevCritEnd = critCurrentEnd;

#pragma omp single
      iterationCounts++;

#pragma omp barrier
    }
  }

  double endt = omp_get_wtime();

  if (isEndOfAlgo) {
    std::cout << "Iteration amount in total: " << iterationCounts
              << std::endl;
    std::cout << "Time taken: " << endt - startt << " sec"
              << std::endl;

    std::cout << "Type of scheduling is: " << xstr(TYPE) << std::endl;
    std::cout << "Chunk size is: " << xstr(CHUNK_SIZE) << std::endl;

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
