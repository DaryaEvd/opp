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

      std::cout << i << " " << j << "   ";
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

bool equalsToPrevEvolution(int *prevMatr, int *currMatrix, int rows,
                           int columns) {
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < columns; j++) {
      if (currMatrix[i * columns + j] != prevMatr[i * columns + j]) {
        return false;
      }
    }
  }
  return true;
}

void countStopVector(int **history, int *matrix, int *stopFlags,
                     int iteration, int rows, int columns) {
  int componentVector;

  for (int i = 0; i < iteration; i++) {
    componentVector = 1;

    if (!equalsToPrevEvolution(history[i], matrix, rows, columns)) {
      componentVector = 0;
      break;
    }

    if (componentVector == 0) {
      break;
    }

    stopFlags[i] = componentVector;
  }
}

void countStopFlags(int **historyOfEvolution, bool *stopFlag,
                    int *rowsNumArray, int iter, int rankOfCurrProc,
                    int amountOfProcs, int rows, int columns) {
  if (iter == 0) {
    for (int proc = 0; proc < amountOfProcs; proc++) {
      stopFlag[proc] = 0;
      return;
    }
  }

  int *now = historyOfEvolution[iter];

  for (int i = 0; i < iter; i++) {
    stopFlag[iter * rankOfCurrProc + i] = 1;
    if (!equalsToPrevEvolution(historyOfEvolution[i], now, rows,
                               columns)) {
      stopFlag[iter * rankOfCurrProc + i] = 0;
    }
  }
}

bool checkComparison(bool *allStopVectors, int iteration,
                     int amountOfProcs, int rankOfCurrProc) {
  for (int i = 0; i < iteration; i++) {
    if (allStopVectors[i] == false) {
      return false;
    }
  }
  return true;
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

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  const long maxIterations = 1000;
  int **historyOfEvolution =
      new int *[maxIterations](); // TODO: don't forget to init mem

  // amount of elems, handling each process
  int *elemsNumArray = countElemsNumInEachProc(
      amountOfProcs, rowsAmount, columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArray = countRowsInEachProcess(
      elemsNumArray, amountOfProcs, columnsAmount);

  historyOfEvolution[0] = new int[rowsNumArray[rankOfCurrProc] + 2]();

  if (rankOfCurrProc == 0) {
    std::fstream inputFile;
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

  MPI_Request requestFirstLineSend;
  MPI_Request requestLastLineSend;
  MPI_Request requestGetLastLine;
  MPI_Request requestGetFirstLine;
  MPI_Request requestVector;

  MPI_Status status;

  bool *vectorStopFlagPerIter = new bool[maxIterations];
  bool *allStopVectors = new bool[maxIterations * amountOfProcs];

  int iterCurr = 0;
  bool repeated = false;

  while (iterCurr < maxIterations && !repeated) {

    historyOfEvolution[iterCurr + 1] =
        new int[((rowsNumArray[rankOfCurrProc] + 2) * columnsAmount)];

    currentGen = historyOfEvolution[iterCurr];
    nextGen = historyOfEvolution[iterCurr + 1];

    // 1 - initiation of sending first line to the prev core
    MPI_Isend(&currentGen[columnsAmount], columnsAmount, MPI_INT,
              rankPrev, tagFirstLine, MPI_COMM_WORLD,
              &requestFirstLineSend);

    // 2 - initiation of sending last line the to next core
    MPI_Isend(
        &currentGen[(rowsNumArray[rankOfCurrProc]) * columnsAmount],
        columnsAmount, MPI_INT, rankNext, tagLastLine, MPI_COMM_WORLD,
        &requestLastLineSend);

    // 3 - initiation of receiving last line from the previous core
    MPI_Irecv(currentGen, columnsAmount, MPI_INT, rankPrev,
              tagLastLine, MPI_COMM_WORLD, &requestGetLastLine);

    // 4 - initiation of receiving first line from the next core
    MPI_Irecv(&currentGen[(rowsNumArray[rankOfCurrProc] + 1) *
                          columnsAmount],
              columnsAmount, MPI_INT, rankNext, tagFirstLine,
              MPI_COMM_WORLD, &requestGetFirstLine);

    // 5 - count vector of stop flags
    countStopFlags(historyOfEvolution, vectorStopFlagPerIter,
                   rowsNumArray, iterCurr, rankOfCurrProc,
                   amountOfProcs, rowsNumArray[rankOfCurrProc],
                   columnsAmount);

    // 6 - init changing of stop vectors with all cores
    MPI_Ialltoall(vectorStopFlagPerIter, iterCurr, MPI_C_BOOL,
                  allStopVectors, iterCurr, MPI_C_BOOL,
                  MPI_COMM_WORLD, &requestVector);

    // 7 - count stages of rows, except first and last line
    computeNextGeneration(&currentGen[columnsAmount], nextGen,
                          rowsNumArray[rankOfCurrProc] - 1,
                          columnsAmount);

    // 8 - wait end sending 1st line to prev core
    MPI_Wait(&requestFirstLineSend, &status);

    // 9 - wait end of receiving from the 3rd step
    MPI_Wait(&requestGetLastLine, &status);

    // 10 - count stages of the first line
    // computeNextGenerationInFirstLine()
    computeNextGeneration(currGen, nextGen, 3, columnsAmount);

    // 11 - wait end sending last line to the next core
    MPI_Wait(&requestLastLineSend, &status);

    // 12 - wait end receiving
    MPI_Wait(&requestGetFirstLine, &status);

    // 13 - count stages of the last line
    computeNextGeneration(
        &currentGen[(rowsNumArray[rankOfCurrProc] - 3) *
                    columnsAmount],
        &nextGen[(rowsNumArray[rankOfCurrProc] - 3) * columnsAmount],
        3, columnsAmount);
    // 14 - wait end of exchanginf stop vectors with each other
    // process
    MPI_Wait(&requestVector, &status);

    // 15 - compare vectors of stop
    if (checkComparison((allStopVectors, iterCurr, rankOfCurrProc,
                         amountOfProcs))) {
      break;
    }

    iterCurr++;
  }

  MPI_Finalize();
  return 0;
}
