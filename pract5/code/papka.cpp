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

// bool *countStopVector(int lengthVector) {
//   bool *vector = new bool[lengthVector]();
// }

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

void countStopVector() {}

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

  int *currentGen = new int[rowsAmount * columnsAmount]();

  if (rankOfCurrProc == 0) {
    std::fstream inputFile;
    inputFile.open("begin.txt",
                   std::ios::in | std::ios::out | std::ios::trunc);
    if (!inputFile) {
      std::cout << "Can't open input file\n";
      return 0;
    }

    generateGlider(currentGen, rowsAmount, columnsAmount);

    printMatrixToFile(currentGen, rowsAmount, columnsAmount,
                      inputFile);
  }

  // amount of elems, handling each process
  int *elemsNumArray = countElemsNumInEachProc(
      amountOfProcs, rowsAmount, columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArray = countRowsInEachProcess(
      elemsNumArray, amountOfProcs, columnsAmount);

  // count, from which element data sends to each process
  int *elementsOffsetArray =
      createElemsOffsetArr(elemsNumArray, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray = createRowsOffsetArr(
      elementsOffsetArray, amountOfProcs, rowsAmount);

  if (rankOfCurrProc == 0) {
    for (int i = 0; i < amountOfProcs; i++) {
      std::cout << i << ": " << rowsNumArray[i] << std::endl;
    }
  }

  int *partCurrGen = new int[elemsNumArray[rankOfCurrProc]]();

  MPI_Scatterv(currentGen, elemsNumArray, elementsOffsetArray,
               MPI_INT, partCurrGen, elemsNumArray[rankOfCurrProc],
               MPI_INT, 0, MPI_COMM_WORLD);

  int *upperLine = new int[columnsAmount]();
  int *lowerLine = new int[columnsAmount]();

  const long maxIterations = 1000;

  int **historyOfEvolution =
      new int *[maxIterations](); // TODO: don't forget to init mem

  int iterCurr = 0;
  bool repeated = false;

  // determine lower and upper nrighbour

  int rankPrev = (amountOfProcs + rankOfCurrProc - 1) % amountOfProcs;
  int rankNext = (amountOfProcs + rankOfCurrProc + 1) % amountOfProcs;

  int *nextGen = new int[rowsAmount * columnsAmount]();
  int *partNextGen = new int[elemsNumArray[rankOfCurrProc]]();

  int tagFirstLine = 0;
  int tagLastLine = 1;

  MPI_Request requestFirstLineSend;
  MPI_Request requestLastLineSend;
  MPI_Request requestGetLastLine;
  MPI_Request requestGetFirstLine;

  MPI_Status status;

  bool *matrixStopFlag = new bool[maxIterations * amountOfProcs];
  bool *vectorStopFlag = new bool[maxIterations];

  for (int i = 0; i < maxIterations; i++) {
    vectorStopFlag[i] = false;
  }

  while (iterCurr < maxIterations && !repeated) {

    // 1 - initiation of sending first line to the prev core
    MPI_Isend(partCurrGen, columnsAmount, MPI_INT, rankPrev,
              tagFirstLine, MPI_COMM_WORLD, &requestFirstLineSend);

    // 2 - initiation of sending last line the to next core
    MPI_Isend(
        // partCurrGen[elemsNumArray[rankOfCurrProc] - columnsAmount],
        &partCurrGen[rowsOffsetArray[rankOfCurrProc] -
                     1], // mb not -1, but -2, idk yet ...
        columnsAmount, MPI_INT, rankNext, tagLastLine, MPI_COMM_WORLD,
        &requestLastLineSend);

    // 3 - initiation of receiving last line from the previous core
    MPI_Irecv(upperLine, columnsAmount, MPI_INT, rankPrev,
              tagLastLine, MPI_COMM_WORLD, &requestGetLastLine);

    // 4 - initiation of receiving first line from the next core
    MPI_Irecv(lowerLine, columnsAmount, MPI_INT, rankNext,
              tagFirstLine, MPI_COMM_WORLD, &requestGetFirstLine);

    historyOfEvolution[iterCurr] =
        new int[elemsNumArray[rankOfCurrProc]]();

    copyMatrix(historyOfEvolution, partCurrGen,
               rowsNumArray[iterCurr], columnsAmount);

    // countStopVector();
    for (int i = iterCurr; i > -1; i--) {
      if (equalsToPrevEvolution(
              historyOfEvolution[iterCurr], partCurrGen,
              rowsNumArray[rankOfCurrProc], columnsAmount)) {
          vectorStopFlag[i] = 1;
      }
      else {
        vectorStopFlag[i] = 0;
      }
    }

    // 7 - count stages of rows, except first and last line
    computeNextGeneration(partCurrGen, partNextGen, rowsAmount,
                          columnsAmount);

    // 8 - wait end receiving from the 2nd step
    MPI_Wait(&requestLastLineSend, &status);

    // 9 - wait end of receiving from the 3rd step
    MPI_Wait(&requestGetLastLine, &status);

    // 10 - count stages of the first line
  }

  MPI_Finalize();
  return 0;
}
