#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mpi.h>

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

double *fillConstantVector(const size_t sizeInput) {
  double *vector = new double[sizeInput];
  for (size_t i = 0; i < sizeInput; ++i) {
    vector[i] = sizeInput + 1.0;
  }
  return vector;
}

void zerofyDoubleVectors(double *vector, const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    vector[i] = 0;
  }
}
void zerofyIntVectors(int *vector, const size_t sizeInput) {
  for (size_t i = 0; i < sizeInput; i++) {
    vector[i] = 0;
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

int countRows(int amountOfAvailableProc, int rankOfCurrentProc,
              int sizeInput) {
  int basicCount = sizeInput / amountOfAvailableProc;
  int restCount = sizeInput % amountOfAvailableProc;
  return basicCount + (rankOfCurrentProc < restCount ? 1 : 0);
}

void multimplyMatrixOnVector(const double *matrix,
                             const double *vector, double *res,
                             const size_t height,
                             const size_t width) {
  for (size_t i = 0; i < height; ++i) {
    // res[i] = 0;
    for (size_t j = 0; j < width; ++j) {
      res[i] += matrix[i * width + j] * vector[j];
    }
  }
}

void parallelMultMatrixOnVector(const double *matrixA_part,
                                const double *vector, double *res,
                                const size_t widthPartMatrA,
                                const size_t heightPartMatrA,
                                int amountOfAvailableProc,
                                int rankOfCurrentProc) {
  // matrixA_part = width * rowsNum[rankOfCurrentProc];
  // res - vector which size dimension is heightPartMatrA
  double *tmpRes = new double[heightPartMatrA];
  zerofyDoubleVectors(tmpRes, heightPartMatrA);
  multimplyMatrixOnVector(matrixA_part, vector, tmpRes, heightPartMatrA,
                          widthPartMatrA);

  int *recvCounts = new int[amountOfAvailableProc];
  zerofyIntVectors(recvCounts, widthPartMatrA);
  for (size_t i = 0; i < amountOfAvailableProc; i++) {
    recvCounts[i] = widthPartMatrA / amountOfAvailableProc;
    recvCounts[i] += (i < widthPartMatrA % amountOfAvailableProc);
  }

  int *displs = new int[amountOfAvailableProc];
  zerofyIntVectors(displs, widthPartMatrA);
  displs[0] = 0;
  for (size_t i = 1; i < amountOfAvailableProc; i++) {
    displs[i] = displs[i - 1] + recvCounts[i - 1];
  }

  /*** MPI_Allgatherv - собирает блоки с разным числом элементов от каждого процесса 
   * @param sendbuf [in] - starting address of send buffer
   * @param sendcounts [in] - amount of elems in send buffer
   * @param sendtype - data type of send buffer elements
   * @param recvbuf [out] - address of receive buffer
   * @param recvcounts - array (=group size) containing the number of
   elements that are to be received from each process
   * @param displs - array (=group size). Entry i
   specifies the displacement (relative to recvbuf) at which to place
   the incoming data from process i
   * @param recvtype - data type of recieve buffer elements
   * @param Comm - communicator
   */
  MPI_Allgatherv(tmpRes, heightPartMatrA, MPI_DOUBLE, res,
                 recvCounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);

  delete[] tmpRes;
  delete[] recvCounts;
  delete[] displs;
}

void substructVectors(const double *vector1, const double *vector2,
                      double *res, const size_t size) {
  for (size_t j = 0; j < size; ++j) {
    res[j] = vector1[j] - vector2[j];
  }
}

void countVectorMultNumber(const double *vector, double scalar,
                           double *res, const size_t size) {
  for (size_t i = 0; i < size; ++i) {
    res[i] = scalar * vector[i];
  }
}

void copyVectors(const double *src, double *dst, const size_t size) {
  for (size_t i = 0; i < size; i++) {
    dst[i] = src[i];
  }
}

int main(int argc, char *argv[]) {
  const size_t sizeInput = atoi(argv[1]);

  /*    constants declararion   */
  const double accuracy = 1e-10;
  const double epsilon = 1e-10;

  const size_t maxIterationCounts = 50000;
  const size_t maxConvergenceCount = 5;
  size_t iterationCounts = 0;
  int convergenceCount = 0;

  /* turn on mpi.h*/
  MPI_Init(&argc, &argv);
  double startTime = MPI_Wtime();

  /* count amount of proccesses and ranks in proc*/
  int amountOfAvailableProc;
  MPI_Comm_size(MPI_COMM_WORLD, &amountOfAvailableProc);

  int rankOfCurrentProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrentProc);

  double *matrixA_full = nullptr;
  if (rankOfCurrentProc == 0) {
    matrixA_full = fillConstantMatrix(sizeInput);
  }
  // const int rowsAmount = countRows(amountOfAvailableProc,
  // rankOfCurrentProc, sizeInput);

  int basicCount = sizeInput / amountOfAvailableProc;
  int restCount = sizeInput % amountOfAvailableProc;

  int *sendCounts = new int[amountOfAvailableProc];
  int *rowsNum = new int[amountOfAvailableProc];
  zerofyIntVectors(sendCounts, amountOfAvailableProc);
  zerofyIntVectors(rowsNum, amountOfAvailableProc);

  for (size_t i = 0; i < amountOfAvailableProc; i++) {
    rowsNum[i] = basicCount +
                 (i < restCount); // amount of rows in each process
    sendCounts[i] =
        rowsNum[i] * sizeInput; // amount of elements in each process
  }

  int *displs = new int[amountOfAvailableProc];
  displs[0] = 0;
  for (int i = 1; i < sizeInput; i++) {
    displs[i] = displs[i - 1] + sendCounts[i - 1];
  }
  double *matrixA_part = new double[sendCounts[rankOfCurrentProc]];

  /*** MPI_Scatterv - рапределяет блоки данных по всем процессам
   * @param sendbuf [in] - address of send buffer (significant only at
   * root)
   * @param sendcounts [in] - array (=group size)
   * specifying the number of elements to send to each processor
   * @param displs [in] - array (=group size).
   * Entry i specifies the displacement (relative to sendbuf from
   * which to take the outgoing data to process i)
   * i-oe значение посылает данные в i-й блок
   * @param sendtype - type of sending data
   * @param recvbuf [out] - address of reciever buffer's starting
   * @param recvcounts - amount of elems that we recieve
   * @param recvtype - type of receiving data
   * @param root - number of provess-reciever
   * @param Comm - communicator
   */
  MPI_Scatterv(matrixA_full, sendCounts, displs, MPI_DOUBLE,
               matrixA_part, sendCounts[rankOfCurrentProc],
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  size_t widthPartA = sizeInput; // ширина матрицы
  size_t heightPartA =
      rowsNum[rankOfCurrentProc]; // количество строк
                                  // в данном процессе
  
  double *xCurr = new double[sizeInput];
  if (rankOfCurrentProc == 0) {
    zerofyDoubleVectors(xCurr, sizeInput);
  }

  double *b = nullptr;
  if (rankOfCurrentProc == 0) {
    b = fillConstantVector(sizeInput);
  }

  //посылает сообщения всем процессам группы, включая себя
  MPI_Bcast(b, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // size_t b_width = 1;
  size_t b_height = sizeInput;
  double bNorm = countVectorLength(b_height, b);

  double *Atmp = new double[sizeInput];
  zerofyDoubleVectors(Atmp, sizeInput);
  double *y = new double[sizeInput];
  zerofyDoubleVectors(y, sizeInput);
  double *tauY = new double[sizeInput];
  zerofyDoubleVectors(tauY, sizeInput);
  double *xNext = new double[sizeInput];
  zerofyDoubleVectors(xNext, sizeInput);

  double prevEndCrit = 0, endCrit = 0;

  while (1) {
    MPI_Bcast(xCurr, sizeInput, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // we need to count A * x_n (= Atmp)
    // it's 'fat operation', so need to parallel it :)
    parallelMultMatrixOnVector(matrixA_part, b, Atmp, widthPartA,
                               heightPartA, amountOfAvailableProc,
                               rankOfCurrentProc);

    substructVectors(Atmp, b, y, sizeInput); // y_n = A * x_n - b
    double yNorm =
        countVectorLength(sizeInput, y); // || A * x_n - b ||

    zerofyDoubleVectors(Atmp, sizeInput);
    parallelMultMatrixOnVector(matrixA_part, y, Atmp, widthPartA,
                               heightPartA, amountOfAvailableProc,
                               rankOfCurrentProc); // A * y_n
    double numerator =
        countScalarMult(y, Atmp, sizeInput); // (y_n, A * y_n)
    double denumerator =
        countScalarMult(Atmp, Atmp, sizeInput); // (A * y_n, A * y_n)
    
    double tau = numerator / denumerator;

    countVectorMultNumber(y, tau, tauY, sizeInput); // tau * y
    substructVectors(xCurr, tauY, xNext,
                     sizeInput); // x_n+1 = x_n - tau * y

    double endCriteria = yNorm / bNorm;

    iterationCounts++;

    if (iterationCounts > maxIterationCounts) {
      std::cout << "Too many iterations. Change init values"
                << std::endl;
      delete[] xCurr;
      delete[] xNext;
      delete[] Atmp;
      delete[] y;
      delete[] tauY;
      delete[] matrixA_part;
      delete[] matrixA_full;
      delete[] b;

      return 0;
    }

    if (endCriteria < epsilon) {
      break;
    }

    if (convergenceCount > maxConvergenceCount) {
      break;
    } else if (endCrit < prevEndCrit) {
      convergenceCount++;
    } else if (endCrit > prevEndCrit) {
      convergenceCount = 0;
    }

    copyVectors(xNext, xCurr, sizeInput);
    prevEndCrit = endCrit;

    std::cout << "iteration " << iterationCounts
              << " ended ===" << std::endl;
  }  

  if (rankOfCurrentProc == 0) {
    double endTime = MPI_Wtime();
    if (convergenceCount <= maxConvergenceCount) {
      std::cout << "Iteration amount in total: " << iterationCounts
                << "\n";
      std::cout << "Time taken: " << endTime - startTime << " sec";
      std::cout << std::endl;
    } else {
      std::cout << "There are no solutions!\n";
    }
  }

  delete[] xCurr;
  delete[] xNext;
  delete[] Atmp;
  delete[] y;
  delete[] tauY;
  delete[] matrixA_part;
  delete[] matrixA_full;
  delete[] b;

  // MPI_Barrier(MPI_COMM_WORLD);

  /*don't forget to finalize*/
  MPI_Finalize();
  return 0;
}

