#include "mpi.h"
#include <iostream>
#include <iterator>
#include <vector>

void generateGlider(bool *extendedPartMatr, int columnsAmount) {
  extendedPartMatr[2 * columnsAmount + 0] = true;
  extendedPartMatr[2 * columnsAmount + 1] = true;
  extendedPartMatr[2 * columnsAmount + 2] = true;
  extendedPartMatr[1 * columnsAmount + 2] = true;
  extendedPartMatr[0 * columnsAmount + 1] = true;
}

bool equalsMatrices(const bool *first, const bool *second,
                    size_t from, size_t to) {
  for (size_t i = from; i < to; ++i) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

void calcStopVectors(std::vector<bool *> historyOfEvolution,
                     bool *stopVector, bool *extendedPartMatr,
                     int rowsAmount, int columnsAmount) {
  size_t vector_size = historyOfEvolution.size() - 1;
  auto it = historyOfEvolution.begin();
  for (int i = 0; i < vector_size; ++i) {
    stopVector[i] =
        equalsMatrices(*it, extendedPartMatr, columnsAmount,
                       columnsAmount * (rowsAmount + 1));
    it++;
  }
}

bool isStop(int rowsAmount, int columnsAmount,
            const bool *stopMatrix) {
  for (int i = 0; i < columnsAmount; ++i) {
    bool stop = true;
    for (int j = 0; j < rowsAmount; ++j) {
      stop &= stopMatrix[j * columnsAmount + i];
    }
    if (stop)
      return true;
  }
  return false;
}

// x - rows, y - columns
int countNeighbors(bool *oldData, int columnsAmount, int i, int j) {
  int cnt = 0;
  if (oldData[i * columnsAmount + (j + 1) % columnsAmount])
    cnt++;
  if (oldData[i * columnsAmount +
              (j + columnsAmount - 1) % columnsAmount])
    cnt++;
  if (oldData[(i + 1) * columnsAmount + (j + 1) % columnsAmount])
    cnt++;
  if (oldData[(i + 1) * columnsAmount +
              (j + columnsAmount - 1) % columnsAmount])
    cnt++;
  if (oldData[(i - 1) * columnsAmount + (j + 1) % columnsAmount])
    cnt++;
  if (oldData[(i - 1) * columnsAmount +
              (j + columnsAmount - 1) % columnsAmount])
    cnt++;
  if (oldData[(i + 1) * columnsAmount + j])
    cnt++;
  if (oldData[(i - 1) * columnsAmount + j])
    cnt++;
  return cnt;
}

void computeNextGeneration(bool *oldData, bool *nextData,
                           int rowsAmount, int columnsAmount) {
  for (int i = 1; i < rowsAmount - 1; ++i) {
    for (int j = 0; j < columnsAmount; ++j) {

      int state = oldData[i * columnsAmount + j];

      int neighborsAmount =
          countNeighbors(oldData, columnsAmount, i, j);

      if (state == 0 && neighborsAmount == 3) {
        nextData[i * columnsAmount + j] = 1;
      } else if (state == 1 &&
                 (neighborsAmount < 2 || neighborsAmount > 3)) {
        nextData[i * columnsAmount + j] = 0;
      } else {
        nextData[i * columnsAmount + j] =
            oldData[i * columnsAmount + j] = state;
      }
    }
  }
}

int *countElemsNumInEachProc(int amountOfProcs, int rows,
                             int columns) {
  int *elemsNum = new int[amountOfProcs];

  int basicRowsCount = rows / amountOfProcs;
  int restRowsCount = rows % amountOfProcs;

  for (int i = 0; i < amountOfProcs; i++) {
    elemsNum[i] = basicRowsCount * columns;
    if (restRowsCount > 0) {
      elemsNum[i] += columns;
      --restRowsCount;
    }
  }
  return elemsNum;
}

int *countRowsInEachProcess(const int *elementsNumberArr,
                            int amountOfProcs, int columns) {
  int *rowsNumArr = new int[amountOfProcs];
  for (int i = 0; i < amountOfProcs; i++) {
    rowsNumArr[i] = elementsNumberArr[i] / columns;
  }
  return rowsNumArr;
}

int *createElemsOffsetArr(const int *elemsNum, int amountOfProcs) {
  int *elementsOffsetArray = new int[amountOfProcs];
  int elementsOffset = 0;
  for (int i = 0; i < amountOfProcs; i++) {
    elementsOffsetArray[i] = elementsOffset;
    elementsOffset += elemsNum[i];
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

void startLife(int amountOfProcs, int rankOfCurrProc,
               bool *startMatrix, int rowsAmount, int columnsAmount) {
  // amount of elems, handling each process
  int *counts = countElemsNumInEachProc(amountOfProcs, rowsAmount,
                                        columnsAmount);

  // amount of rows, which each process handles
  int *rowsNumArr =
      countRowsInEachProcess(counts, amountOfProcs, columnsAmount);

  // count, from which element data sends to each process
  int *offsets = createElemsOffsetArr(counts, amountOfProcs);

  // count, from which row data sends to each process
  int *rowsOffsetArray =
      createRowsOffsetArr(offsets, amountOfProcs, rowsAmount);

  // if (rankOfCurrProc == 0) {
  //   for (int i = 0; i < amountOfProcs; i++) {
  //     std::cout << i << ": " << rowsNumArr[i] << std::endl;
  //   }
  // }

  bool *extendedPartMatr =
      new bool[counts[rankOfCurrProc] + columnsAmount * 2];
  bool *basePartMatr = extendedPartMatr + columnsAmount;

  MPI_Scatterv(startMatrix, counts, offsets, MPI_C_BOOL, basePartMatr,
               counts[rankOfCurrProc], MPI_C_BOOL, 0, MPI_COMM_WORLD);

  int prevRank = (rankOfCurrProc + amountOfProcs - 1) % amountOfProcs;
  int nextRank = (rankOfCurrProc + 1) % amountOfProcs;

  std::vector<bool *> historyOfEvolution;

  int iteration = 0;
  bool stop = false;
  while (!stop) {

    bool *extendedNextPartMatr =
        new bool[counts[rankOfCurrProc] + 2 * columnsAmount];
    bool *baseNextPartMatr = extendedNextPartMatr + columnsAmount;
    historyOfEvolution.push_back(extendedPartMatr);

    iteration++;

    MPI_Request requestSendFirstLine;
    // 1 - initiation of sending first line to the prev core
    MPI_Isend(basePartMatr, columnsAmount, MPI_C_BOOL, prevRank, 1,
              MPI_COMM_WORLD, &requestSendFirstLine);

    MPI_Request requestSendLastLine;
    // 2 - initiation of sending last line the to next core
    MPI_Isend(basePartMatr + counts[rankOfCurrProc] - columnsAmount,
              columnsAmount, MPI_C_BOOL, nextRank, 0, MPI_COMM_WORLD,
              &requestSendLastLine);

    MPI_Request requestGetLastLine;
    // 3 - initiation of receiving last line from the previous
    // core
    MPI_Irecv(extendedPartMatr, columnsAmount, MPI_C_BOOL, prevRank,
              0, MPI_COMM_WORLD, &requestGetLastLine);

    MPI_Request requestGetFirstLine;
    // 4 - initiation of receiving first line from the next core
    MPI_Irecv(basePartMatr + counts[rankOfCurrProc], columnsAmount,
              MPI_C_BOOL, nextRank, 1, MPI_COMM_WORLD,
              &requestGetFirstLine);

    // 5 - count vector of stop flags
    MPI_Request flagsReq;
    bool *stopVector;
    bool *stopMatrix;
    int vector_size = historyOfEvolution.size() - 1;
    if (vector_size > 1) {
      stopVector = new bool[vector_size];
      calcStopVectors(historyOfEvolution, stopVector,
                      extendedPartMatr, rowsNumArr[rankOfCurrProc],
                      columnsAmount);
      stopMatrix = new bool[vector_size * amountOfProcs];

      // 6 - init changing of stop vectors with all cores

      MPI_Iallgather(stopVector, (int)vector_size, MPI_C_BOOL,
                     stopMatrix, (int)vector_size, MPI_C_BOOL,
                     MPI_COMM_WORLD, &flagsReq);
    }

    // 7 - count stages of rows, except first and last line
    computeNextGeneration(basePartMatr, baseNextPartMatr,
                          rowsNumArr[rankOfCurrProc], columnsAmount);

    // 8 - wait end sending 1st line to prev core
    MPI_Wait(&requestSendFirstLine, MPI_STATUS_IGNORE);

    // 9 - wait end of receiving from the 3rd step
    MPI_Wait(&requestGetLastLine, MPI_STATUS_IGNORE);

    // 10 - count stages of the first line
    computeNextGeneration(extendedPartMatr, extendedNextPartMatr, 3,
                          columnsAmount);

    11 - wait end sending last line to the next core MPI_Wait(
             &requestSendLastLine, MPI_STATUS_IGNORE);

    12 - wait end receiving MPI_Wait(&requestGetFirstLine,
                                     MPI_STATUS_IGNORE);

    13 - count stages of the last line computeNextGeneration(
             basePartMatr +
                 (rowsNumArr[rankOfCurrProc] - 2) * columnsAmount,
             baseNextPartMatr +
                 (rowsNumArr[rankOfCurrProc] - 2) * columnsAmount,
             3, columnsAmount);

    if (vector_size > 1) {
      // 14 - wait end of exchanging stop vectors with each other
      // process

      MPI_Wait(&flagsReq, MPI_STATUS_IGNORE);

      // 15 - compare vectors of stop
      stop = isStop(amountOfProcs, (int)vector_size, stopMatrix);

      delete[] stopVector;
      delete[] stopMatrix;
    }

    if (stop)
      break;

    extendedPartMatr = extendedNextPartMatr;
    basePartMatr = baseNextPartMatr;
  }
  if (rankOfCurrProc == 0) {
    std::cout << "Repeat after: " << iteration - 1 << " iterations"
              << std::endl;
  }

  for (bool *matrixDump : historyOfEvolution) {
    delete[] (matrixDump);
  }

  delete[] offsets;
  delete[] counts;
}

int main(int argc, char *argv[]) {

  if (argc != 3) {
    std::cout
        << "Bad amount of arguments!\n"
        << "Enter rowsAmount amount, then enter columns amount.\n"
        << std::endl;
    return 0;
  }

  const int rowsAmount = std::atoi(argv[1]);
  const int columnsAmount = std::atoi(argv[2]);

  MPI_Init(&argc, &argv);

  int amountOfProcs, rankOfCurrProc;

  MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rankOfCurrProc);

  bool *startMatrix = nullptr;
  if (rankOfCurrProc == 0) {
    startMatrix = new bool[rowsAmount * columnsAmount];
    std::fill(startMatrix, startMatrix + columnsAmount * rowsAmount,
              false);
    generateGlider(startMatrix, columnsAmount);
  }

  double startt = MPI_Wtime();
  startLife(amountOfProcs, rankOfCurrProc, startMatrix, rowsAmount,
            columnsAmount);
  double endt = MPI_Wtime();

  if (rankOfCurrProc == 0) {
    std::cout << "Time taken: " << endt - startt << " sec."
              << std::endl;
  }

  if (rankOfCurrProc == 0)
    delete[] startMatrix;
  MPI_Finalize();
  return 0;
}