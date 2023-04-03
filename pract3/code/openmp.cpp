#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>

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
#pragma omp for
  for (size_t i = 0; i < sizeInput; ++i) {
    res[i] = 0;
    for (size_t j = 0; j < sizeInput; ++j) {

      /*
      Почему здесь не нужен atomic?
      Потому что res считается отдельно для каждого потока
      и здесь не происходит гонки данных, ведь
      данные не зависят друг от друга
      (то есть res[i] принадлежит только одному потоку).


      Как работает atomic?
      -Он не создает копии локальных переменных
      (в отличие от reduction)
      -Гарантирует, что во время работы одной нити у других
      нитей не будет доступа к res[i]
      */
      res[i] += matrix[i * sizeInput + j] * vector[j];
    }
  }
}

/*============ only for test, starts ============*/
/* for testing when vector b uses sin*/
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
  std::cout << "vector U is: " << std::endl;
  printVector(vectorU, sizeInput);
  multimplyMatrixOnVector(fullMatrixA, vectorU, b, sizeInput);

  delete[] vectorU;
}
/*============ only for test, ended ============*/

void countScalarMult(const double *vector1, const double *vector2,
                     const size_t sizeInput, double *res) {
#pragma omp single
  { *res = 0; }

#pragma omp for reduction(+ : res[0])
  for (size_t i = 0; i < sizeInput; ++i) {
    *res += vector1[i] * vector2[i];
  }
}

void countVectorLength(const size_t sizeInput, const double *vector,
                       double *res) {
#pragma omp single
  { *res = 0; }
#pragma omp for reduction(+ : res[0])
  for (size_t i = 0; i < sizeInput; ++i) {
    *res += vector[i] * vector[i];
  }
#pragma omp single
  { *res = sqrt(*res); }
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t sizeInput) {
#pragma omp for
  for (size_t i = 0; i < sizeInput; ++i) {
    res[i] = scalar * vector[i];
  }
}

/* Что происходит в методах countScalarMult(...) и
countVectorLength(...) ? В мэйне создается переменная, и бы работаем с
ней. То есть в методы передаем указатель на эту переменную, зануляем
её, потом производим с ней расчеты. И записываем результат в *res.

Почему нельзя написать '...reduction(+ : *res)' ?
Ответ: потому что компилятор (как минимум на сервере) не позволяет так
сделать. Поэтому можно считать, что res - это массив из одной
переменной, поэтому и записываем в нулевую ячейку.

А можно ли было обойтись без введения переменной res извне? Да.

--------
double countScalarMult(const double *vector1, const double *vector2,
                       const size_t sizeInput, double *res) {
   double res = 0;
  #pragma omp parallel for reduction(+ : res)
    for (size_t i = 0; i < sizeInput; ++i) {
      res += vector1[i] * vector2[i];
    }
    return res;
}

double countVectorLength(const size_t sizeInput,
                         const double *vector) {
  return sqrt(countScalarMult(vector, vector, sizeInput));
}
----------
Вот здесь '#pragma omp parallel for reduction(+ : res)' мы
создаем параллельную секцию, в которой res не приватная, и с
переменной res можно спокойно работать.
Если убрать слово 'parallel', то программа не скомпилируется,
как раз-таки потому что тогда res будет приватной в цикле)
*/

void substructVectors(const double *vector1, const double *vector2,
                      double *res, const size_t sizeInput) {
#pragma omp for
  for (size_t j = 0; j < sizeInput; ++j) {
    res[j] = vector1[j] - vector2[j];
  }
}

void copyVector(const double *src, double *dst,
                const size_t sizeInput) {
#pragma omp for
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

  ///* for testing when vector b uses sin starts ========= */
  // double *matrixA = fillConstantMatrix(sizeInput);
  // // printMatrix(matrixA, sizeInput);

  // double *b = new double[sizeInput];
  // fillVectorB(b, matrixA, sizeInput);
  // // std::cout << "vector b is" << std::endl;
  // // printVector(b, sizeInput);

  // std::fill(xCurr, xCurr + sizeInput, 0);
  // // std::cout << "vector X at start is: " << std::endl;
  // // printVector(xCurr, sizeInput);
  ///* for testing when vector b uses sin ended ========= */

  double bNorm, yNorm;
  countVectorLength(sizeInput, b, &bNorm);

  double numeratorTau, denominatorTau;

#pragma omp parallel num_threads(amountOfThreads)                    \
    shared(iterationCounts, critCurrentEnd, prevPrevCritEnd,         \
           isEndOfAlgo, Atmp, y, tauY, xCurr, xNext, matrixA, b,     \
           bNorm, yNorm, numeratorTau, denominatorTau)
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

#pragma omp single
      {
        prevPrevCritEnd = prevCritEnd;
        prevCritEnd = critCurrentEnd;
        iterationCounts++;
      }
    }
  }

  double endt = omp_get_wtime();

  if (isEndOfAlgo) {
    std::cout << "Iteration amount in total: " << iterationCounts
              << std::endl;
    std::cout << "Time taken: " << endt - startt << " sec"
              << std::endl;

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