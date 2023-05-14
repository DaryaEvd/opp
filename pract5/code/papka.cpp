#include <cmath>
#include <fstream>
#include <iostream>
#include <mpi.h>

int *countElemsNumInEachProc(int amountOfProcs, int rows,
                             int columns) {
  int *elemsNumArray = new int[amountOfProcs];

  int basicRowsCount = rows / amountOfProcs;
  int restRowsCount = rows % amountOfProcs;

  for (int i = 0; i < amountOfProcs; i++) {
    elemsNumArray[i] = basicRowsCount * columns;
    if (restRowsCount > 0) {
      elemsNumArray[i] += columns;
      --restRowsCount;
    }
  }
  return elemsNumArray;
}

int *countRowsInEachProcess(const int *elementsNumberArr,
                            int amountOfProcs, int columns) {
  int *rowsNumArray = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsNumArray[i] = elementsNumberArr[i] / columns;
  }
  return rowsNumArray;
}

int *createElemsOffsetArr(const int *elemsNumArray,
                          int amountOfProcs) {
  int *elementsOffsetArray = new int[amountOfProcs];
  int elementsOffset = 0;
  for (int i = 0; i < amountOfProcs; i++) {
    elementsOffsetArray[i] = elementsOffset;
    elementsOffset += elemsNumArray[i];
  }
  return elementsOffsetArray;
}

int *createRowsOffsetArr(const int *elementsOffsetArray,
                         int amountOfProcs, int rows) {
  int *rowsOffsetArray = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsOffsetArray[i] = elementsOffsetArray[i] / rows;
  }

  return rowsOffsetArray;
}

void generateGlider(int *data, int rows, int columns) {
  data[0 * columns + 1] = 1;
  data[1 * columns + 2] = 1;
  data[2 * columns + 0] = 1;
  data[2 * columns + 1] = 1;
  data[2 * columns + 2] = 1;
}

void printMatrixToFile(int *data, int rows, int columns,
                       std::fstream &file) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      file << data[i * columns + j] << " ";
    }
    file << "\n";
  }
}

void copyMatrix(int *dest, int *src, int rows, int columns) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {
      dest[i * columns + j] = src[i * columns + j];
    }
  }
}

// x - rows, y - columns
int countNeighbors(int *data, int rows, int columns, int xMatr,
                   int yMatr) {
  int sum = 0;

  for (int i = -1; i < 2; ++i) {   // rows
    for (int j = -1; j < 2; ++j) { // columns

      int currRow = (xMatr + i + rows) % rows;
      int currColumn = (yMatr + j + columns) % columns;

      sum += data[currRow * columns + currColumn];
    }
  }

  sum -= data[xMatr * columns + yMatr]; // cell itself
  // std::cout << "SUMMA is: " << sum << std::endl << std::endl;

  return sum;
}

void computeNextGeneration(int *oldData, int *nextData, int rows,
                           int columns) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < columns; ++j) {

      int state = oldData[i * columns + j];

      // std::cout << i << " " << j << "   ";
      int neighborsAmount =
          countNeighbors(oldData, rows, columns, i, j);

      if (state == 0 && neighborsAmount == 3) {
        nextData[i * columns + j] = 1;
      } else if (state == 1 &&
                 (neighborsAmount < 2 || neighborsAmount > 3)) {
        nextData[i * columns + j] = 0;
      } else {
        nextData[i * columns + j] = oldData[i * columns + j] = state;
      }
    }
  }
}

bool equalsToPrevEvolution(int *prevMatr, int *currMatrix, int size) {
  // for (int i = 0; i < rows; i++) {
  //   for (int j = 0; j < columns; j++) {
  //     if (currMatrix[i * columns + j] != prevMatr[i * columns + j])
  //     {
  //       return false;
  //     }
  //   }
  // }
  // return true;

  for (int i = 0; i < size; i++) {
    if (prevMatr[i] != currMatrix[i]) {
      return false;
    }
  }
  return true;
}

bool *calcVector(int **historyOfEvolution, int *currentGen, int iter,
                 int size) {
  bool *vectorStop = new bool[iter];
  for (int j = 0; j < iter; j++) {
    vectorStop[j] = equalsToPrevEvolution(historyOfEvolution[j],
                                          currentGen, size);
  }
  return vectorStop;
}

// iter - x
// amountOfProcs - y
bool isEnd(bool *stopMatrix, int iter, int amountOfProcs) {
  for (int x = 0; x < iter; x++) {
    int count = 0;

    for (int y = 0; y < amountOfProcs; y++) {
      if (stopMatrix[y * iter + x]) {
        count++;
      }
    }

    if (count == amountOfProcs) {
      return true;
    }
  }
  return false;
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Bad amount of arguments!\n"
              << "Enter rows amount, then enter columns amount.\n"
              << std::endl;
    return 0;
  }

  const int rowsAmount = atoi(argv[1]);
  const int columnsAmount = atoi(argv[2]);
  const long maxIterations = 1000;

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  double startt, endt;

  int **historyOfEvolution =
      new int *[maxIterations](); // TODO: don't forget to init mem

  // amount of elems, handling each process
  int *elemsNumArray = countElemsNumInEachProc(
      amountOfProcs, rowsAmount, columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArray = countRowsInEachProcess(
      elemsNumArray, amountOfProcs, columnsAmount);

  historyOfEvolution[0] =
      new int[(rowsNumArray[rankOfCurrProc] + 2) * columnsAmount]();

  std::fstream inputFile;
  if (rankOfCurrProc == 0) {

    inputFile.open("begin.txt",
                   std::ios::in | std::ios::out | std::ios::trunc);
    if (!inputFile) {
      std::cout << "Can't open input file\n";
      return 0;
    }

    generateGlider(historyOfEvolution[0], rowsAmount, columnsAmount);

    printMatrixToFile(historyOfEvolution[0], rowsAmount,
                      columnsAmount, inputFile);
  }

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < amountOfProcs; i++) {
      std::cout << i << ": " << rowsNumArray[i] << std::endl;
    }
  }

  int *currentGen; // = new int[rowsAmount * columnsAmount]();
  int *nextGen;    // = new int[rowsAmount * columnsAmount]();

  int rankPrev = (amountOfProcs + rankOfCurrProc - 1) % amountOfProcs;
  int rankNext = (amountOfProcs + rankOfCurrProc + 1) % amountOfProcs;

  int tagFirstLine = 0;
  int tagLastLine = 1;

  MPI_Request requestSendFirstLine;
  MPI_Request requestSendLastLine;
  MPI_Request requestGetLastLine;
  MPI_Request requestGetFirstLine;
  MPI_Request requestVector;

  MPI_Status status;

  bool *vectorStopFlagPerIter = new bool[maxIterations];
  bool *allStopVectors = new bool[maxIterations * amountOfProcs];

  int iterCurr = 0;
  bool repeated = false;

  if (rankOfCurrProc == 0) {
    startt = MPI_Wtime();
  }
  while (!repeated) {

    historyOfEvolution[iterCurr + 1] =
        new int[((rowsNumArray[rankOfCurrProc] + 2) * columnsAmount)];

    // currentGen = historyOfEvolution[iterCurr];
    // nextGen = historyOfEvolution[iterCurr + 1];

    historyOfEvolution[iterCurr] = currentGen;
    currentGen = nextGen;

    // // 1 - initiation of sending first line to the prev core
    MPI_Isend(&currentGen[columnsAmount], columnsAmount, MPI_INT,
              rankPrev, tagFirstLine, MPI_COMM_WORLD,
              &requestSendFirstLine);

    // // 2 - initiation of sending last line the to next core
    MPI_Isend(
        &currentGen[(rowsNumArray[rankOfCurrProc]) * columnsAmount],
        columnsAmount, MPI_INT, rankNext, tagLastLine, MPI_COMM_WORLD,
        &requestSendLastLine);

    // // 3 - initiation of receiving last line from the previous core
    MPI_Irecv(currentGen, columnsAmount, MPI_INT, rankPrev,
              tagLastLine, MPI_COMM_WORLD, &requestGetLastLine);

    // // 4 - initiation of receiving first line from the next core
    MPI_Irecv(&currentGen[(rowsNumArray[rankOfCurrProc] + 1) *
                          columnsAmount],
              columnsAmount, MPI_INT, rankNext, tagFirstLine,
              MPI_COMM_WORLD, &requestGetFirstLine);

    // // 5 - count vector of stop flags

    bool *stopMatrix;
    bool *stopVector;
    if (iterCurr > 1) {
      stopMatrix = new bool[(iterCurr - 1) * amountOfProcs];

      stopVector =
          calcVector(historyOfEvolution, currentGen, iterCurr - 1,
                     elemsNumArray[rankOfCurrProc]);

      // // 6 - init changing of stop vectors with all cores
      MPI_Iallgather(stopVector, iterCurr - 1, MPI_C_BOOL, stopMatrix,
                     iterCurr - 1, MPI_C_BOOL, MPI_COMM_WORLD,
                     &requestVector);
    }

    // // 7 - count stages of rows, except first and last line
    computeNextGeneration(
        &currentGen[columnsAmount], &nextGen[columnsAmount],
        rowsNumArray[rankOfCurrProc] - 2, columnsAmount);

    // // 8 - wait end sending 1st line to prev core
    MPI_Wait(&requestSendFirstLine, &status);

    // 9 - wait end of receiving from the 3rd step
    MPI_Wait(&requestGetLastLine, &status);

    // // 10 - count stages of the first line
    // // computeNextGenerationInFirstLine()
    computeNextGeneration(currentGen, nextGen, 3, columnsAmount);

    // // 11 - wait end sending last line to the next core
    MPI_Wait(&requestSendLastLine, &status);

    // // 12 - wait end receiving
    MPI_Wait(&requestGetFirstLine, &status);

    // // 13 - count stages of the last line
    computeNextGeneration(
        &currentGen[(rowsNumArray[rankOfCurrProc] - 3) *
                    columnsAmount],
        &nextGen[(rowsNumArray[rankOfCurrProc] - 3) * columnsAmount],
        3, columnsAmount);

    // std::copy(nextGen, &nextGen[elemsNumArray[rankOfCurrProc]],
    //           currentGen);

    if (iterCurr > 1) {
      // 14 - wait end of exchanginf stop vectors with each other
      // process
      MPI_Wait(&requestVector, MPI_STATUS_IGNORE);

      // 15 - compare vectors of stop
      if (isEnd(stopMatrix, iterCurr - 1, amountOfProcs)) {
        repeated = true;
        // break;
      }

      delete[] stopVector;
      delete[] stopMatrix;
    }

    // printMatrixToFile(currentGen, rowsNumArray[rankOfCurrProc],
    //                   columnsAmount, inputFile);

    iterCurr++;
  }

  if (repeated) {
    std::cout << "First repeat of cell automat stage after: "
              << iterCurr << " iterCurr" << std::endl;
  } else {
    std::cout << "Finished after iteration: " << iterCurr
              << std::endl;
  }

  inputFile.close();

  delete[] elemsNumArray;
  delete[] rowsNumArray;

  for (int i = 0; i < maxIterations; i++) {
    delete[] historyOfEvolution[i];
  }
  delete[] historyOfEvolution;

  MPI_Finalize();
  return 0;
}
